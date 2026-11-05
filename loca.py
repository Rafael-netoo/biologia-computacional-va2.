
from Bio.Align import substitution_matrices

# ==========================================
# 1. CONFIGURAÇÕES
# ==========================================

# Carregar matriz PAM250
pam = substitution_matrices.load("PAM250")

# Penalidade por lacuna
GAP = -5


# ==========================================
# 2. ALGORITMO SMITH-WATERMAN
# ==========================================

def alinhamento_local(seq1, seq2):

    m = len(seq1)
    n = len(seq2)

    # Criar matriz de pontuação
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Matriz para armazenar os caminhos
    caminho = [[None] * (n + 1) for _ in range(m + 1)]

    maior_pontuacao = 0
    melhor_i = 0
    melhor_j = 0

    # Preencher a matriz
    for i in range(1, m + 1):
        for j in range(1, n + 1):

            # Pontuação diagonal
            diagonal = (
                dp[i - 1][j - 1]
                + pam[seq1[i - 1], seq2[j - 1]]
            )

            # Pontuação vertical
            acima = dp[i - 1][j] + GAP

            # Pontuação horizontal
            esquerda = dp[i][j - 1] + GAP

            # Escolher a melhor pontuação
            melhor = max(0, diagonal, acima, esquerda)

            dp[i][j] = melhor

            # Registrar o caminho
            if melhor == 0:
                caminho[i][j] = None
            elif melhor == diagonal:
                caminho[i][j] = "D"
            elif melhor == acima:
                caminho[i][j] = "C"
            else:
                caminho[i][j] = "E"

            # Identificar a maior pontuação
            if melhor > maior_pontuacao:
                maior_pontuacao = melhor
                melhor_i = i
                melhor_j = j

    # ==========================================
    # 3. RECUPERAR OS TRECHOS ALINHADOS
    # ==========================================

    i = melhor_i
    j = melhor_j

    alinhamento1 = []
    alinhamento2 = []

    while i > 0 and j > 0 and dp[i][j] > 0:

        movimento = caminho[i][j]

        if movimento == "D":

            alinhamento1.append(seq1[i - 1])
            alinhamento2.append(seq2[j - 1])

            i -= 1
            j -= 1

        elif movimento == "C":

            alinhamento1.append(seq1[i - 1])
            alinhamento2.append("-")

            i -= 1

        elif movimento == "E":

            alinhamento1.append("-")
            alinhamento2.append(seq2[j - 1])

            j -= 1

    # Inverter os trechos recuperados
    alinhamento1.reverse()
    alinhamento2.reverse()

    trecho1 = "".join(alinhamento1).replace("-", "")
    trecho2 = "".join(alinhamento2).replace("-", "")

    return maior_pontuacao, trecho1, trecho2


# ==========================================
# 4. LER ENTRADA DO ARQUIVO TXT
# ==========================================

with open("loca.txt", "r", encoding="utf-8") as arquivo:

    linhas = arquivo.readlines()

sequencias = []

for linha in linhas:

    linha = linha.strip()

    # Ignorar identificadores FASTA
    if linha and not linha.startswith(">"):
        sequencias.append(linha)

# Verificar quantidade de sequências
if len(sequencias) != 2:
    raise ValueError(
        "O arquivo input.txt deve conter exatamente duas sequências."
    )


# ==========================================
# 5. EXECUTAR E MOSTRAR RESULTADO
# ==========================================

pontuacao, trecho1, trecho2 = alinhamento_local(
    sequencias[0],
    sequencias[1]
)

print(pontuacao)
print(trecho1)
print(trecho2)
