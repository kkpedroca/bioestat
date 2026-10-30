from flask import Flask, request, jsonify, send_from_directory
from scipy import stats
import math

app = Flask(__name__, static_folder=None)

def fail(message, code=400):
    return jsonify({"error": message}), code

def valid_alpha(value):
    try:
        alpha = float(value)
    except (TypeError, ValueError):
        raise ValueError("Nível de significância inválido.")
    if not 0 < alpha < 1:
        raise ValueError("O nível de significância deve estar entre 0 e 1.")
    return alpha

@app.get("/")
def index():
    return send_from_directory(".", "index.html")

@app.get("/api/health")
def health():
    import scipy
    return jsonify({"status": "ok", "engine": "SciPy", "scipy_version": scipy.__version__})

@app.post("/api/calculate")
def calculate():
    data = request.get_json(silent=True) or {}
    test = data.get("test")
    try:
        alpha = valid_alpha(data.get("alpha", 0.05))
        if test == "shapiro":
            values = data.get("values")
            if not isinstance(values, list) or len(values) < 3:
                return fail("O Shapiro–Wilk exige pelo menos 3 observações.")
            if len(values) > 5000:
                return fail("Para este aplicativo, o limite é de 5.000 observações por amostra.")
            try: x = [float(v) for v in values]
            except (TypeError, ValueError): return fail("A amostra contém valores inválidos.")
            if any(not math.isfinite(v) for v in x): return fail("Todos os valores devem ser números finitos.")
            if len(set(x)) < 2: return fail("A amostra não pode ter todos os valores iguais.")
            result = stats.shapiro(x)
            return jsonify({"test":"Shapiro-Wilk","statistic":float(result.statistic),"p_value":float(result.pvalue),"n":len(x),"alpha":alpha,"method":"scipy.stats.shapiro"})

        if test == "kruskal":
            groups = data.get("groups")
            if not isinstance(groups, list) or len(groups) < 2: return fail("Informe pelo menos dois grupos.")
            clean_groups, labels = [], []
            for i, group in enumerate(groups):
                if not isinstance(group, dict): return fail("Formato de grupo inválido.")
                values = group.get("values")
                if not isinstance(values, list) or not values: return fail(f"O grupo {i+1} está vazio.")
                try: vals = [float(v) for v in values]
                except (TypeError, ValueError): return fail(f"O grupo {i+1} contém valores inválidos.")
                if any(not math.isfinite(v) for v in vals): return fail(f"O grupo {i+1} contém valores não finitos.")
                clean_groups.append(vals); labels.append(str(group.get("name") or f"Grupo {i+1}")[:80])
            if sum(map(len, clean_groups)) < 3: return fail("Informe pelo menos três observações no total.")
            result = stats.kruskal(*clean_groups, nan_policy="raise")
            pooled = [v for group in clean_groups for v in group]
            N = len(pooled)
            _, counts = __import__('numpy').unique(pooled, return_counts=True)
            tie_sum = sum(int(t)**3-int(t) for t in counts if t > 1)
            correction = 1.0 - tie_sum/(N**3-N) if N > 1 else 1.0
            ranks = stats.rankdata(pooled, method="average")
            details, offset = [], 0
            for label, group in zip(labels, clean_groups):
                gr = ranks[offset:offset+len(group)]
                details.append({"name":label,"n":len(group),"mean_rank":float(sum(gr)/len(group))})
                offset += len(group)
            return jsonify({"test":"Kruskal-Wallis","statistic":float(result.statistic),"p_value":float(result.pvalue),"df":len(clean_groups)-1,"n":N,"tie_correction":float(correction),"groups":details,"alpha":alpha,"method":"scipy.stats.kruskal"})

        if test == "fisher":
            table = data.get("table")
            if not isinstance(table,list) or len(table)!=2 or any(not isinstance(row,list) or len(row)!=2 for row in table):
                return fail("A tabela deve ter exatamente 2 linhas e 2 colunas.")
            try: cells = [[int(v) for v in row] for row in table]
            except (TypeError,ValueError,OverflowError): return fail("As quatro frequências devem ser números inteiros.")
            if any(v<0 for row in cells for v in row): return fail("As frequências não podem ser negativas.")
            if any(float(table[i][j])!=cells[i][j] for i in range(2) for j in range(2)): return fail("As frequências devem ser números inteiros.")
            if sum(map(sum,cells))==0: return fail("A tabela não pode estar vazia.")
            alternative = {"two":"two-sided","greater":"greater","less":"less"}.get(data.get("alternative","two"))
            if alternative is None: return fail("Hipótese alternativa inválida.")
            result = stats.fisher_exact(cells, alternative=alternative)
            odds = float(result.statistic)
            return jsonify({"test":"Exato de Fisher","odds_ratio":odds if math.isfinite(odds) else None,"p_value":float(result.pvalue),"n":sum(map(sum,cells)),"alpha":alpha,"alternative":alternative,"method":"scipy.stats.fisher_exact"})
        return fail("Teste não reconhecido.")
    except ValueError as exc:
        return fail(str(exc))
    except Exception as exc:
        return fail(f"Não foi possível calcular: {str(exc)}",422)

if __name__ == "__main__":
    app.run(host="127.0.0.1",port=5000,debug=False)
