# BioEstat — versão aprimorada

Aplicação educacional de testes estatísticos com backend Python/SciPy. Os cálculos são determinísticos; não se usa IA generativa para produzir resultados estatísticos.

## Recursos
- Shapiro–Wilk com gráfico Q–Q, resumo descritivo e aviso interpretativo sobre a hipótese nula.
- Kruskal–Wallis com postos, correção para empates e pós-teste de Dunn com ajuste de Holm quando o teste global é significativo.
- Exato de Fisher com OR amostral e intervalo de confiança condicional exato de 95% quando estimável.
- Para cada teste, o usuário escolhe entre demonstração simplificada e completa do procedimento.
- Validação de entradas no servidor; limite de 5.000 observações para Shapiro–Wilk.

## Executar localmente (Windows)
```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
py -m pip install -r requirements.txt
py app.py
```
Abra http://127.0.0.1:5000.

## Publicar no Render
- Build Command: `pip install -r requirements.txt`
- Start Command: `gunicorn app:app --bind 0.0.0.0:$PORT`
- Runtime: Python 3
- Root Directory: deixe vazio se os arquivos estiverem na raiz do repositório.

## Observações metodológicas
- Shapiro–Wilk: não rejeitar H₀ não prova normalidade. A API pública do SciPy não expõe os coeficientes internos do algoritmo AS R94; a demonstração completa explicita essa limitação e mostra a identidade matemática de W e seus componentes calculáveis.
- Kruskal–Wallis: o pós-teste de Dunn é executado apenas quando p global < α. Os p-valores bilaterais são ajustados pelo método sequencial de Holm, controlando o erro familiar.
- Fisher: o IC exibido é o intervalo condicional exato de 95%, calculado a partir da estimativa condicional do odds ratio. Ele não é o mesmo estimador pontual do OR amostral ad/bc exibido no resultado.
- Em dados com células zero, o OR amostral ou os limites do IC podem ser infinitos/indefinidos; o aplicativo informa quando não pode apresentar limites finitos.
- Para pesquisa/publicação, registre versões do software, defina previamente hipóteses e alternativas e avalie os pressupostos do delineamento.
