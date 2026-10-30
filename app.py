from flask import Flask, request, jsonify, send_from_directory
from scipy import stats
import numpy as np
import math
import os

app = Flask(__name__, static_folder=None)

MAX_SHAPIRO_N = 5000
MAX_TOTAL_N = 100000


def fail(message, code=400):
    return jsonify({"error": message}), code


def valid_alpha(value):
    try:
        alpha = float(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError("Nível de significância inválido.")
    if not math.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("O nível de significância deve estar entre 0 e 1.")
    return alpha


def numeric_values(values, label):
    if not isinstance(values, list) or not values:
        raise ValueError(f"{label}: informe pelo menos uma observação.")
    try:
        result = [float(v) for v in values]
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f"{label}: há valores inválidos.")
    if not all(math.isfinite(v) for v in result):
        raise ValueError(f"{label}: todos os valores devem ser números finitos.")
    return result


@app.get("/")
def index():
    return send_from_directory(".", "index.html")


@app.get("/api/health")
def health():
    import scipy
    return jsonify({
        "status": "ok",
        "engine": "SciPy",
        "scipy_version": scipy.__version__
    })


@app.post("/api/calculate")
def calculate():
    data = request.get_json(silent=True)
    if not isinstance(data, dict):
        return fail("Envie os dados em formato JSON válido.")

    test = data.get("test")

    try:
        alpha = valid_alpha(data.get("alpha", 0.05))

        # SHAPIRO–WILK
        if test == "shapiro":
            values = data.get("values")
            if not isinstance(values, list) or len(values) < 3:
                return fail("O Shapiro–Wilk exige pelo menos 3 observações.")
            if len(values) > MAX_SHAPIRO_N:
                return fail(f"Para este aplicativo, o limite é de {MAX_SHAPIRO_N:,} observações por amostra.".replace(",", "."))
            x = numeric_values(values, "A amostra")
            if len(set(x)) < 2:
                return fail("A amostra não pode ter todos os valores iguais.")

            result = stats.shapiro(x)
            sorted_x = np.sort(np.asarray(x, dtype=float))
            n = len(sorted_x)

            # Quantis teóricos de Filliben, usados no gráfico Q–Q
            probabilities = (np.arange(1, n + 1) - 0.375) / (n + 0.25)
            theoretical = stats.norm.ppf(probabilities)

            # Reta de referência baseada nos quartis
            q25, q75 = np.percentile(sorted_x, [25, 75])
            norm_q25, norm_q75 = stats.norm.ppf([0.25, 0.75])
            slope = (q75 - q25) / (norm_q75 - norm_q25)
            intercept = q25 - slope * norm_q25
            trendline = slope * theoretical + intercept

            return jsonify({
                "test": "Shapiro–Wilk",
                "statistic": float(result.statistic),
                "p_value": float(result.pvalue),
                "n": n,
                "alpha": alpha,
                "summary": {
                    "mean": float(np.mean(x)),
                    "median": float(np.median(x)),
                    "std": float(np.std(x, ddof=1)),
                    "iqr": float(stats.iqr(x))
                },
                "chart": {
                    "type": "qqplot",
                    "observed": sorted_x.tolist(),
                    "theoretical": theoretical.tolist(),
                    "trendline": trendline.tolist()
                },
                "method": "scipy.stats.shapiro"
            })

        # KRUSKAL–WALLIS
        if test == "kruskal":
            groups = data.get("groups")
            if not isinstance(groups, list) or len(groups) < 2:
                return fail("Informe pelo menos dois grupos independentes.")

            clean_groups, labels = [], []
            for i, group in enumerate(groups):
                if not isinstance(group, dict):
                    return fail(f"O formato do grupo {i + 1} é inválido.")
                try:
                    vals = numeric_values(group.get("values"), f"Grupo {i + 1}")
                except ValueError as exc:
                    return fail(str(exc))
                clean_groups.append(vals)
                labels.append(str(group.get("name") or f"Grupo {i + 1}")[:80])

            total_n = sum(map(len, clean_groups))
            if total_n < 3:
                return fail("Informe pelo menos três observações no total.")
            if total_n > MAX_TOTAL_N:
                return fail(f"O limite total é de {MAX_TOTAL_N:,} observações.".replace(",", "."))

            pooled = np.asarray([v for group in clean_groups for v in group], dtype=float)
            if np.all(pooled == pooled[0]):
                return fail("O teste não pode ser calculado quando todas as observações são iguais.")

            result = stats.kruskal(*clean_groups, nan_policy="raise")
            _, counts = np.unique(pooled, return_counts=True)
            tie_sum = sum(int(t) ** 3 - int(t) for t in counts if t > 1)
            correction = 1.0 - tie_sum / (total_n ** 3 - total_n)
            ranks = stats.rankdata(pooled, method="average")

            details, chart_groups, offset = [], [], 0
            for label, group in zip(labels, clean_groups):
                group_ranks = ranks[offset:offset + len(group)]
                mean_rank = float(np.mean(group_ranks))
                median = float(np.median(group))
                q25, q75 = np.percentile(group, [25, 75])

                details.append({
                    "name": label,
                    "n": len(group),
                    "mean_rank": mean_rank,
                    "median": median,
                    "min": float(np.min(group)),
                    "max": float(np.max(group)),
                    "q25": float(q25),
                    "q75": float(q75)
                })
                chart_groups.append({
                    "name": label,
                    "mean_rank": mean_rank,
                    "median": median,
                    "q25": float(q25),
                    "q75": float(q75)
                })
                offset += len(group)

            return jsonify({
                "test": "Kruskal–Wallis",
                "statistic": float(result.statistic),
                "p_value": float(result.pvalue),
                "df": len(clean_groups) - 1,
                "n": total_n,
                "tie_correction": float(correction),
                "groups": details,
                "alpha": alpha,
                "chart": {
                    "type": "rank_means",
                    "data": chart_groups
                },
                "method": "scipy.stats.kruskal"
            })

        # TESTE EXATO DE FISHER
        if test == "fisher":
            table = data.get("table")
            if (
                not isinstance(table, list)
                or len(table) != 2
                or any(not isinstance(row, list) or len(row) != 2 for row in table)
            ):
                return fail("A tabela deve ter exatamente 2 linhas e 2 colunas.")

            cells = []
            for row in table:
                clean_row = []
                for value in row:
                    if isinstance(value, bool):
                        return fail("As quatro frequências devem ser números inteiros.")
                    try:
                        number = float(value)
                    except (TypeError, ValueError, OverflowError):
                        return fail("As quatro frequências devem ser números inteiros.")
                    if not math.isfinite(number) or not number.is_integer():
                        return fail("As quatro frequências devem ser números inteiros.")
                    if number < 0:
                        return fail("As frequências não podem ser negativas.")
                    clean_row.append(int(number))
                cells.append(clean_row)

            total = sum(sum(row) for row in cells)
            if total == 0:
                return fail("A tabela não pode estar vazia.")

            alternative = {
                "two": "two-sided",
                "greater": "greater",
                "less": "less"
            }.get(data.get("alternative", "two"))
            if alternative is None:
                return fail("Hipótese alternativa inválida.")

            result = stats.fisher_exact(cells, alternative=alternative)
            odds = float(result.statistic)

            total_g1 = cells[0][0] + cells[0][1]
            total_g2 = cells[1][0] + cells[1][1]
            pct_g1 = cells[0][0] / total_g1 * 100 if total_g1 else 0.0
            pct_g2 = cells[1][0] / total_g2 * 100 if total_g2 else 0.0

            return jsonify({
                "test": "Exato de Fisher",
                "odds_ratio": odds if math.isfinite(odds) else None,
                "p_value": float(result.pvalue),
                "n": total,
                "alpha": alpha,
                "alternative": alternative,
                "chart": {
                    "type": "proportions",
                    "g1_event": cells[0][0],
                    "g1_no_event": cells[0][1],
                    "g2_event": cells[1][0],
                    "g2_no_event": cells[1][1],
                    "pct_g1": float(pct_g1),
                    "pct_g2": float(pct_g2)
                },
                "method": "scipy.stats.fisher_exact"
            })

        return fail("Teste não reconhecido.")

    except ValueError as exc:
        return fail(str(exc))
    except Exception:
        app.logger.exception("Falha inesperada durante o cálculo estatístico.")
        return fail("Não foi possível concluir o cálculo. Verifique os dados e tente novamente.", 422)


if __name__ == "__main__":
    app.run(
        host="0.0.0.0",
        port=int(os.environ.get("PORT", 5000)),
        debug=False
    )
