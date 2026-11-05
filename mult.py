
from itertools import product
from pathlib import Path

# Lê o arquivo input.txt que está na mesma pasta deste programa
caminho = Path(__file__).with_name("mult.txt")

with open(caminho, "r", encoding="utf-8") as arquivo:
    linhas = arquivo.readlines()

# Extrai as sequências do formato FASTA
sequencias = []
atual = ""

for linha in linhas:
    linha = linha.strip()

    if not linha:
        continue

    if linha.startswith(">"):
        if atual:
            sequencias.append(atual)
            atual = ""
    else:
        atual += linha

if atual:
    sequencias.append(atual)

if len(sequencias) != 4:
    raise ValueError(
        f"Esperadas 4 sequências, mas foram encontradas {len(sequencias)}."
    )

# Comprimentos das quatro sequências
tamanhos = tuple(len(seq) for seq in sequencias)

# Calcula a pontuação de uma coluna
def pontuar_coluna(caracteres):
    pontos = 0

    for i in range(4):
        for j in range(i + 1, 4):
            if caracteres[i] != caracteres[j]:
                pontos -= 1

    return pontos


# Programação dinâmica em quatro dimensões
dp = {(0, 0, 0, 0): 0}
caminho = {}

# Percorre todos os estados possíveis
for estado in product(
    range(tamanhos[0] + 1),
    range(tamanhos[1] + 1),
    range(tamanhos[2] + 1),
    range(tamanhos[3] + 1)
):

    if estado == (0, 0, 0, 0):
        continue

    melhor_pontuacao = float("-inf")
    melhor_anterior = None

    # Cada máscara representa quais sequências avançam
    # de posição nesta coluna do alinhamento.
    for mascara in range(1, 16):

        anterior = list(estado)
        coluna = []
        valido = True

        for i in range(4):

            if mascara & (1 << i):
                if anterior[i] == 0:
                    valido = False
                    break

                coluna.append(
                    sequencias[i][anterior[i] - 1]
                )
                anterior[i] -= 1

            else:
                coluna.append("-")

        if not valido:
            continue

        anterior = tuple(anterior)

        if anterior not in dp:
            continue

        pontos_coluna = pontuar_coluna(coluna)
        pontuacao = dp[anterior] + pontos_coluna

        if pontuacao > melhor_pontuacao:
            melhor_pontuacao = pontuacao
            melhor_anterior = (anterior, coluna)

    dp[estado] = melhor_pontuacao
    caminho[estado] = melhor_anterior


# Recupera o alinhamento encontrado
estado = tamanhos
alinhamentos = [[], [], [], []]

while estado != (0, 0, 0, 0):

    anterior, coluna = caminho[estado]

    for i in range(4):
        alinhamentos[i].append(coluna[i])

    estado = anterior

# Inverte as sequências, pois a recuperação foi feita de trás para frente
for i in range(4):
    alinhamentos[i].reverse()

# Exibe o resultado
print(dp[tamanhos])

for alinhamento in alinhamentos:
    print("".join(alinhamento))
