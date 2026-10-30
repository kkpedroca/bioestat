# BioEstat — versão com backend estatístico SciPy

Esta versão mantém a interface HTML e transfere os cálculos para um backend Python.
Não utiliza IA generativa: testes estatísticos devem ser calculados por algoritmos
determinísticos e validados, não por um modelo de linguagem. O backend utiliza SciPy.

## Requisitos
- Python 3.10 ou superior
- Internet na primeira instalação para baixar Flask, SciPy e NumPy

## Executar no Windows
1. Extraia o ZIP e abra o terminal dentro da pasta `bioestat_robusto`.
2. Crie o ambiente: `py -m venv .venv`
3. Ative no PowerShell: `.venv\Scripts\Activate.ps1`
   ou no CMD: `.venv\Scripts\activate.bat`
4. Instale: `py -m pip install -r requirements.txt`
5. Inicie: `py app.py`
6. Abra: http://127.0.0.1:5000

## macOS/Linux
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python app.py
```
Acesse http://127.0.0.1:5000.

## Implementações
- Shapiro–Wilk: `scipy.stats.shapiro` (até 5.000 observações).
- Kruskal–Wallis: `scipy.stats.kruskal`, com correção para empates.
- Exato de Fisher: `scipy.stats.fisher_exact`, com opções bilateral e unilateral.
- Validação de entradas no servidor e mensagens de erro.

Importante: não abra `index.html` diretamente com duplo clique. Inicie o servidor
e acesse pelo endereço local acima. Para pesquisa/publicação, registre as versões
do software e avalie os pressupostos e o plano estatístico do estudo.
