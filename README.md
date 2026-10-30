# BioEstat

Aplicação web de testes estatísticos com backend Python, Flask e SciPy.
A interface inclui gráficos interativos com Chart.js.

## Testes disponíveis
- Shapiro–Wilk: estatística W, p-valor, resumo descritivo e gráfico Q–Q.
- Kruskal–Wallis: estatística H, p-valor, graus de liberdade, correção para empates, postos médios, medianas e gráfico.
- Exato de Fisher: odds ratio, p-valor bilateral ou unilateral e gráfico de proporções.

## Executar localmente (Windows)
1. Extraia o ZIP.
2. Abra o terminal dentro desta pasta.
3. Crie o ambiente: `py -m venv .venv`
4. Ative: `.venv\Scripts\Activate.ps1`
5. Instale: `py -m pip install -r requirements.txt`
6. Inicie: `py app.py`
7. Acesse `http://127.0.0.1:5000`.

## Publicar no Render
Crie um Web Service conectado ao repositório GitHub que contenha estes arquivos na raiz.

- Language: Python 3
- Root Directory: deixe vazio
- Build Command: `pip install -r requirements.txt`
- Start Command: `gunicorn app:app --bind 0.0.0.0:$PORT`
- Instance Type: Free ou outro plano desejado

O endpoint `/api/health` informa se o backend está ativo e mostra a versão do SciPy.

## Observações estatísticas
- Shapiro–Wilk: o aplicativo limita a amostra a 5.000 observações; para amostras acima desse limite, a precisão do p-valor do SciPy também merece atenção.
- Kruskal–Wallis: a significância indica diferença entre pelo menos dois grupos; não identifica quais pares diferem. Um pós-teste apropriado pode ser necessário.
- Exato de Fisher: a interpretação do odds ratio depende da organização das linhas e colunas da tabela. A opção bilateral é a configuração padrão.
- A ferramenta auxilia a análise, mas não substitui a avaliação dos pressupostos, do desenho do estudo e do plano estatístico.

## Arquivos
- `app.py`: API Flask e cálculos SciPy.
- `index.html`: interface e gráficos.
- `requirements.txt`: dependências.
