
from Bio.Align import substitution_matrices

# Carregar matriz BLOSUM62
blosum = substitution_matrices.load("BLOSUM62")

GAP = -5


def alinhamento_global(seq1, seq2):

    m = len(seq1)
    n = len(seq2)

    # Criar matriz
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Inicializar primeira coluna
    for i in range(1, m + 1):
        dp[i][0] = dp[i - 1][0] + GAP

    # Inicializar primeira linha
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j - 1] + GAP

    # Preencher matriz
    for i in range(1, m + 1):
        for j in range(1, n + 1):

            diagonal = (
                dp[i - 1][j - 1]
                + blosum[seq1[i - 1], seq2[j - 1]]
            )

            acima = dp[i - 1][j] + GAP
            esquerda = dp[i][j - 1] + GAP

            dp[i][j] = max(diagonal, acima, esquerda)

    return dp[m][n]


# Ler arquivo TXT
with open("input.txt", "r") as arquivo:
    linhas = arquivo.readlines()

# Extrair sequências do FASTA
sequencias = []

for linha in linhas:
    linha = linha.strip()

    if linha and not linha.startswith(">"):
        sequencias.append(linha)

print("Sequência 1:", sequencias[0])
print("Sequência 2:", sequencias[1])
print("Comprimento 1:", len(sequencias[0]))
print("Comprimento 2:", len(sequencias[1]))

# Executar alinhamento
resultado = alinhamento_global(
    sequencias[0],
    sequencias[1]
)

print("Pontuação máxima:", resultado)
