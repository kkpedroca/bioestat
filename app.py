from flask import Flask, request, jsonify, send_from_directory
from scipy import stats
from scipy.stats import norm
import numpy as np
import math, os

app = Flask(__name__, static_folder=None)
MAX_SHAPIRO_N = 5000
MAX_TOTAL_N = 100000

def fail(message, code=400):
    return jsonify({"error": message}), code

def valid_alpha(value):
    try: alpha = float(value)
    except (TypeError, ValueError, OverflowError): raise ValueError("Nível de significância inválido.")
    if not math.isfinite(alpha) or not 0 < alpha < 1: raise ValueError("O nível de significância deve estar entre 0 e 1.")
    return alpha

def numeric_values(values, label):
    if not isinstance(values, list) or not values: raise ValueError(f"{label}: informe pelo menos uma observação.")
    try: result = [float(v) for v in values]
    except (TypeError, ValueError, OverflowError): raise ValueError(f"{label}: há valores inválidos.")
    if not all(math.isfinite(v) for v in result): raise ValueError(f"{label}: todos os valores devem ser números finitos.")
    return result

def calculation_modes(simple, complete):
    return {"calculation": {"simplified": simple, "complete": complete}}

@app.get("/")
def index(): return send_from_directory(".", "index.html")

@app.get("/api/health")
def health():
    import scipy
    return jsonify({"status":"ok", "engine":"SciPy", "scipy_version":scipy.__version__})

@app.post("/api/calculate")
def calculate():
    data = request.get_json(silent=True)
    if not isinstance(data, dict): return fail("Envie os dados em formato JSON válido.")
    test = data.get("test")
    try:
        alpha = valid_alpha(data.get("alpha", .05))

        if test == "shapiro":
            values = data.get("values")
            if not isinstance(values, list) or len(values) < 3: return fail("O Shapiro–Wilk exige pelo menos 3 observações.")
            if len(values) > MAX_SHAPIRO_N: return fail(f"O limite é de {MAX_SHAPIRO_N} observações por amostra.")
            x = numeric_values(values, "A amostra")
            if len(set(x)) < 2: return fail("A amostra não pode ter todos os valores iguais.")
            result = stats.shapiro(x)
            ordered = np.sort(np.asarray(x, dtype=float)); n = len(ordered)
            probs = (np.arange(1,n+1)-.375)/(n+.25)
            theoretical = stats.norm.ppf(probs)
            q25,q75 = np.percentile(ordered,[25,75]); nq25,nq75=stats.norm.ppf([.25,.75])
            slope=(q75-q25)/(nq75-nq25); intercept=q25-slope*nq25
            mean=float(np.mean(ordered)); deviations=(ordered-mean); ss=float(np.sum(deviations**2))
            # Numerador equivalente reconstruído da identidade W = numerator² / SS.
            numerator=math.sqrt(max(0.0,float(result.statistic)*ss))
            simple=[
                f"1. Ordene os dados: {', '.join(f'{v:g}' for v in ordered)}.",
                f"2. Média amostral: x̄ = {mean:.10g}.",
                f"3. Soma dos quadrados: Σ(xᵢ − x̄)² = {ss:.10g}.",
                f"4. Estatística: W = (Σ aᵢx₍ᵢ₎)² / Σ(xᵢ − x̄)² = {float(result.statistic):.10g}.",
                f"5. p-valor calculado pelo algoritmo de Shapiro–Wilk (SciPy): {float(result.pvalue):.10g}; compare com α = {alpha:g}."
            ]
            complete=[
                "Hipóteses: H₀ — a amostra é compatível com uma distribuição normal; H₁ — a amostra não é compatível com uma distribuição normal.",
                f"n = {n}; dados ordenados: {', '.join(f'{v:.8g}' for v in ordered)}.",
                f"Média x̄ = Σxᵢ/n = {mean:.12g}.",
                f"Desvios (xᵢ − x̄): {', '.join(f'{v:.6g}' for v in deviations)}.",
                f"Denominador Σ(xᵢ − x̄)² = {ss:.12g}.",
                "A estatística W utiliza coeficientes aᵢ derivados dos valores esperados e da matriz de covariâncias dos estatísticos de ordem normais. O SciPy aplica o algoritmo AS R94; os coeficientes internos não são expostos pela API pública.",
                f"Numerador equivalente (Σaᵢx₍ᵢ₎)² = W × denominador = {float(result.statistic)*ss:.12g}; raiz positiva equivalente = {numerator:.12g}.",
                f"W = {float(result.statistic):.12g}; p-valor = {float(result.pvalue):.12g}; α = {alpha:g}.",
                "Decisão: " + ("rejeitar H₀." if result.pvalue < alpha else "não rejeitar H₀. Não rejeitar H₀ não comprova normalidade; indica apenas que não houve evidência suficiente contra H₀ neste teste.")
            ]
            return jsonify({"test":"Shapiro–Wilk","statistic":float(result.statistic),"p_value":float(result.pvalue),"n":n,"alpha":alpha,
                "summary":{"mean":mean,"median":float(np.median(x)),"std":float(np.std(x,ddof=1)),"iqr":float(stats.iqr(x))},
                "chart":{"type":"qqplot","observed":ordered.tolist(),"theoretical":theoretical.tolist(),"trendline":(slope*theoretical+intercept).tolist()},
                "normality_caution":True,"method":"scipy.stats.shapiro",**calculation_modes(simple,complete)})

        if test == "kruskal":
            groups=data.get("groups")
            if not isinstance(groups,list) or len(groups)<2: return fail("Informe pelo menos dois grupos independentes.")
            clean=[]; labels=[]
            for i,g in enumerate(groups):
                if not isinstance(g,dict): return fail(f"Formato do grupo {i+1} inválido.")
                try: vals=numeric_values(g.get("values"),f"Grupo {i+1}")
                except ValueError as exc: return fail(str(exc))
                clean.append(vals); labels.append(str(g.get("name") or f"Grupo {i+1}")[:80])
            N=sum(map(len,clean))
            if N<3: return fail("Informe pelo menos três observações no total.")
            if N>MAX_TOTAL_N: return fail(f"O limite total é de {MAX_TOTAL_N} observações.")
            pooled=np.asarray([v for g in clean for v in g],dtype=float)
            if np.all(pooled==pooled[0]): return fail("O teste não pode ser calculado quando todas as observações são iguais.")
            ranks=stats.rankdata(pooled,method="average"); _,counts=np.unique(pooled,return_counts=True)
            tie_sum=sum(int(t)**3-int(t) for t in counts if t>1); C=1-tie_sum/(N**3-N)
            result=stats.kruskal(*clean,nan_policy="raise")
            details=[]; chart=[]; offset=0; rank_sums=[]
            for label,g in zip(labels,clean):
                gr=ranks[offset:offset+len(g)]; mr=float(np.mean(gr)); rank_sums.append(float(np.sum(gr)))
                med=float(np.median(g)); q1,q3=np.percentile(g,[25,75])
                details.append({"name":label,"n":len(g),"mean_rank":mr,"rank_sum":float(np.sum(gr)),"median":med,"min":float(np.min(g)),"max":float(np.max(g)),"q25":float(q1),"q75":float(q3)})
                chart.append({"name":label,"mean_rank":mr,"median":med,"q25":float(q1),"q75":float(q3)})
                offset+=len(g)
            # Dunn pairwise comparisons, tie-adjusted variance and Holm familywise correction.
            dunn=[]
            if result.pvalue < alpha:
                variance=(N*(N+1)/12)*C
                raw=[]
                for i in range(len(clean)):
                    for j in range(i+1,len(clean)):
                        se=math.sqrt(variance*(1/len(clean[i])+1/len(clean[j])))
                        z=(details[i]["mean_rank"]-details[j]["mean_rank"])/se if se else 0.0
                        p=float(2*norm.sf(abs(z)))
                        raw.append({"group1":labels[i],"group2":labels[j],"mean_rank_difference":details[i]["mean_rank"]-details[j]["mean_rank"],"z":z,"p_raw":p})
                order=sorted(range(len(raw)),key=lambda k:raw[k]["p_raw"]); m=len(raw); running=0.0; adjusted=[0.0]*m
                for rank_idx,k in enumerate(order):
                    running=max(running,min(1.0,(m-rank_idx)*raw[k]["p_raw"]))
                    adjusted[k]=running
                for k,item in enumerate(raw):
                    dunn.append({**item,"p_holm":adjusted[k],"significant":adjusted[k]<alpha})
            simple=[
                f"1. Reúna os N = {N} valores e atribua postos conjuntos; empates recebem postos médios.",
                f"2. Soma dos postos por grupo: {', '.join(f'{labels[i]} = {rank_sums[i]:.8g}' for i in range(len(labels)))}.",
                f"3. H = [12/(N(N+1)) × Σ(Rᵢ²/nᵢ)] − 3(N+1), com correção de empates C = {C:.10g}.",
                f"4. H corrigido = {float(result.statistic):.10g}; gl = {len(clean)-1}; p = {float(result.pvalue):.10g}; α = {alpha:g}."
            ]
            complete=["Hipóteses: H₀ — as distribuições dos grupos são iguais; H₁ — pelo menos uma distribuição difere.",
                f"N = {N}; fator de correção de empates C = 1 − Σ(t³−t)/(N³−N) = 1 − {tie_sum}/{N**3-N} = {C:.12g}.",
                "Postos conjuntos por observação: " + "; ".join(f"{labels[i]}: " + ", ".join(f"{v:g}" for v in ranks[sum(len(g) for g in clean[:i]):sum(len(g) for g in clean[:i+1])]) for i in range(len(clean))),
                "Somas dos postos e tamanhos: " + "; ".join(f"{d['name']}: R={d['rank_sum']:.10g}, n={d['n']}, posto médio={d['mean_rank']:.10g}" for d in details),
                f"H sem correção = {float(result.statistic)*C:.12g}; H corrigido = H/C = {float(result.statistic):.12g}.",
                f"Aproximação qui-quadrado: gl = k−1 = {len(clean)-1}; p-valor = {float(result.pvalue):.12g}; α = {alpha:g}.",
                "Decisão: " + ("rejeitar H₀; executar comparações múltiplas." if result.pvalue<alpha else "não rejeitar H₀; o pós-teste não é acionado.")]
            if dunn:
                complete.append("Pós-teste de Dunn: erro-padrão = √[(N(N+1)C/12)(1/nᵢ+1/nⱼ)]; z = diferença dos postos médios/erro-padrão. P-valores bilaterais ajustados pelo método de Holm.")
                complete += [f"{d['group1']} × {d['group2']}: diferença={d['mean_rank_difference']:.8g}; z={d['z']:.8g}; p bruto={d['p_raw']:.8g}; p-Holm={d['p_holm']:.8g}; {'significativo' if d['significant'] else 'não significativo'}." for d in dunn]
            return jsonify({"test":"Kruskal–Wallis","statistic":float(result.statistic),"p_value":float(result.pvalue),"df":len(clean)-1,"n":N,"tie_correction":float(C),"groups":details,"alpha":alpha,"dunn":dunn,"posthoc_performed":bool(dunn),"chart":{"type":"rank_means","data":chart},"method":"scipy.stats.kruskal + Dunn com ajuste de Holm",**calculation_modes(simple,complete)})

        if test == "fisher":
            table=data.get("table")
            if not isinstance(table,list) or len(table)!=2 or any(not isinstance(row,list) or len(row)!=2 for row in table): return fail("A tabela deve ter exatamente 2 linhas e 2 colunas.")
            cells=[]
            for row in table:
                cr=[]
                for v in row:
                    if isinstance(v,bool): return fail("As frequências devem ser números inteiros não negativos.")
                    try: z=float(v)
                    except (TypeError,ValueError,OverflowError): return fail("As frequências devem ser números inteiros não negativos.")
                    if not math.isfinite(z) or not z.is_integer() or z<0: return fail("As frequências devem ser números inteiros não negativos.")
                    cr.append(int(z))
                cells.append(cr)
            total=sum(map(sum,cells))
            if total==0: return fail("A tabela não pode estar vazia.")
            alt={"two":"two-sided","greater":"greater","less":"less"}.get(data.get("alternative","two"))
            if alt is None: return fail("Hipótese alternativa inválida.")
            result=stats.fisher_exact(cells,alternative=alt); odds=float(result.statistic)
            # Intervalo condicional exato de 95% (não-Wald), quando estimável.
            ci=None; conditional_or=None; ci_note=None
            try:
                cond=stats.contingency.odds_ratio(cells,kind="conditional")
                conditional_or=float(cond.statistic); interval=cond.confidence_interval(confidence_level=.95)
                lo,hi=float(interval.low),float(interval.high)
                if math.isfinite(lo) and math.isfinite(hi): ci={"low":lo,"high":hi,"level":.95,"estimate":conditional_or,"method":"condicional exato"}
                else: ci_note="O intervalo condicional inclui limite infinito; não foi possível exibir dois limites finitos."
            except Exception: ci_note="O intervalo de confiança não pôde ser estimado para esta configuração da tabela."
            r1=cells[0][0]+cells[0][1]; r2=cells[1][0]+cells[1][1]; c1=cells[0][0]+cells[1][0]; c2=cells[0][1]+cells[1][1]
            pct1=cells[0][0]/r1*100 if r1 else 0; pct2=cells[1][0]/r2*100 if r2 else 0
            simple=[f"1. OR amostral = (a×d)/(b×c) = ({cells[0][0]}×{cells[1][1]})/({cells[0][1]}×{cells[1][0]}) = {odds:.10g}.",
                f"2. Teste exato de Fisher: p-valor = {float(result.pvalue):.10g}; α = {alpha:g}.",
                (f"3. IC 95% condicional exato do OR: [{ci['low']:.8g}; {ci['high']:.8g}]." if ci else f"3. IC 95%: {ci_note or 'não disponível'}")]
            complete=["Hipóteses bilaterais: H₀ — odds ratio = 1 (sem associação); H₁ — odds ratio ≠ 1. Para alternativas unilaterais, a hipótese alternativa segue a opção selecionada.",
                f"Tabela observada: [[{cells[0][0]}, {cells[0][1]}], [{cells[1][0]}, {cells[1][1]}]]; total N={total}.",
                f"Margens fixas: linhas ({r1}, {r2}); colunas ({c1}, {c2}).",
                f"OR amostral = ad/bc = ({cells[0][0]}×{cells[1][1]})/({cells[0][1]}×{cells[1][0]}) = {odds:.12g}.",
                f"Fisher calcula a probabilidade hipergeométrica da tabela observada e soma as probabilidades das tabelas tão ou mais extremas, conforme alternativa '{alt}'. p = {float(result.pvalue):.12g}.",
                (f"Estimativa condicional do OR = {conditional_or:.12g}; IC 95% condicional exato = [{ci['low']:.12g}; {ci['high']:.12g}]." if ci else f"IC 95%: {ci_note or 'não disponível'}"),
                f"Taxas de evento: Grupo 1 = {pct1:.8g}%; Grupo 2 = {pct2:.8g}%; decisão ao nível α={alpha:g}: {'rejeitar H₀' if result.pvalue<alpha else 'não rejeitar H₀'}."]
            return jsonify({"test":"Exato de Fisher","odds_ratio":odds if math.isfinite(odds) else None,"p_value":float(result.pvalue),"n":total,"alpha":alpha,"alternative":alt,"confidence_interval":ci,"ci_note":ci_note,"conditional_odds_ratio":conditional_or,"chart":{"type":"proportions","g1_event":cells[0][0],"g1_no_event":cells[0][1],"g2_event":cells[1][0],"g2_no_event":cells[1][1],"pct_g1":float(pct1),"pct_g2":float(pct2)},"method":"scipy.stats.fisher_exact + odds_ratio condicional",**calculation_modes(simple,complete)})
        return fail("Teste não reconhecido.")
    except ValueError as exc: return fail(str(exc))
    except Exception:
        app.logger.exception("Falha inesperada durante cálculo estatístico.")
        return fail("Não foi possível concluir o cálculo. Verifique os dados e tente novamente.",422)

if __name__ == "__main__":
    app.run(host="0.0.0.0",port=int(os.environ.get("PORT",5000)),debug=False)
