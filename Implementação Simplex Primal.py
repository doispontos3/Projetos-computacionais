# Antes de implementar nossa função do Simplex propriamente dita, vamos implementar funções básicas como uma função de Eliminação de Gauss pra resolver
# sistemas lineares e outras operações basicas.
def prod_interno(x, y):
    if len(x) != len(y):
        return False
    z = 0
    for i in range(len(x)):
        z += x[i]*y[i]
    return z
# Função de transpor matrizes:
def mat_tranposta(A):
    A_transposta = []
    for j in range(len(A[0])):
        At_linha = []
        for i in range(len(A)):
            At_linha.append(A[i][j])
        A_transposta.append(At_linha)
    return A_transposta
# Função de Eliminação de Gauss:
def Gauss_piv(Ax, B):
    Ax_copia = [linha[:] for linha in Ax]
    B_copia = B[:]
    n = len(Ax_copia)
    e1 = 1e-9 #Vamos definir um limite de significancia para os valores de B. (importante pois muitas vezes um valor que era pra ser nulo é
    # não é nulo por causa de erros de arredondamento da máquina da base decimal pra binária)
    for k in range(n):
        p = 0
        linha = 0
        trocou = False
        for i in range(n-k):
            if abs(Ax_copia[i+k][k]) > p:
                p = abs(Ax_copia[i+k][k])
                linha = i + k
                trocou = True
        #Botando o maior pivô de uma linha em sua devida posição:
        if trocou == True:
            Pivo = Ax_copia[linha]
            pivob = B_copia[linha]
            Ax_copia[linha] = Ax_copia[k]
            B_copia[linha] = B_copia[k]
            Ax_copia[k] = Pivo
            B_copia[k] = pivob
        #Zerando os elementos da coluna abaixo do pivô de sua respectiva linha:
        for i in range(n-k-1):
            if abs(Ax_copia[k][k]) < e1: continue
            m_fator = Ax_copia[i+k+1][k]/Ax_copia[k][k]
            B_copia[i+k+1] = B_copia[i+k+1] - m_fator * B_copia[k]
            if abs(B_copia[i+k+1]) < e1:
                B_copia[i+k+1] = 0
            for j in range(len(Ax_copia[i])):
                Ax_copia[i+k+1][j] = Ax_copia[i+k+1][j] - m_fator * Ax_copia[k][j]
                if abs(Ax_copia[i+k+1][j]) < e1:
                    Ax_copia[i+k+1][j] = 0
    x=[]
    for i in range(n):
        x.append(0)
    #Obtendo o valor das variáveis:
    for k in range(n):
        somatorio = 0
        for i in range(n):
            somatorio += Ax_copia[n-k-1][n-i-1] * x[n-i-1]
        if abs(Ax_copia[n-k-1][n-k-1]) > e1:
            x[n-k-1] = (B_copia[n-k-1] - somatorio)/Ax_copia[n-k-1][n-k-1]
    return x

# Vamos computar nossas entradas:
m = int(input()) # Número de restrições (linhas);
n = int(input()) # Número de variáveis (colunas);
c = [float(custo) for custo in input().split()] # Vetor dos custos na funcao objetivo a ser minimizada no formato "c_1 c_2 ... c_n";
b = [float(recurso) for recurso in input().split()] # Vetor dos recursos no formato "b_1 b_2 ... b_m";
A = [] # Matriz com os coeficientes das nossas restrições.
for i in range(m):
    A_linha = [float(a_ij) for a_ij in input().split()] # Linha i da matriz com os coeficientes das restrições no formato "a_i1 a_i2 ... a_in"
    A.append(A_linha)

# Agora vamos para a implementação do método Simplex. Primeiro passo é determinar uma partição básica e calcular o vetor x_B em B_B*x_B = b_B
# com x_N = 0 em primeiro momento. Para a escolha da nossa base, assumiremos nossas ultimas m colunas da matriz A como sendo nossas variáveis de folga e
# portanto nossas variáveis básicas:

def Simplex(m, n, c, b, A):
    PARE = False
# Vamos montar nossa partição básica e não-básica inicial. Salvaremos as matrizes B e N e os índices de suas colunas em relação a matriz A:
    B = []
    B_indice = []
    b_B = []
    c_B = []
    for i in range(m):
        B_linha = []
        for j in range(m):
            B_linha.append(A[i][n-m+j])
        B.append(B_linha)
        B_indice.append(n-m+i)
        b_B.append(b[i])
        c_B.append(c[n-m+i])
    N = []
    N_indice = []
    c_N = []
    for i in range(m):
        N_linha = []
        for j in range(n-m):
            N_linha.append(A[i][j])
        N.append(N_linha)
    for j in range(n-m):
        c_N.append(c[j])
        N_indice.append(j)
    while not PARE == True:
        c_relativo = float('inf')
# Passo 1): Calculamos então nosso x_B:
        x_B = Gauss_piv(B, b_B)
# Passo 2.1): Em seguida calculamos o multiplicador simplex:
        B_transposto = mat_tranposta(B)
        vetor_lambda = Gauss_piv(B_transposto, c_B)
# Passos 2.2) e 2.3) Custos relativos, vamos calcular os custos pra cada variável não básica. Calculado todas, nós salvamos a menor delas
# e seu índice se ela for menor que 0, caso contrário, todas são maiores que 0 e portanto chegamos na solução ótima. Do contrário, nossa variável
# com menor custo relativo entra na base.
        menor_N = 0
        for i in range(n-m):
            N_coluna_i = [N[k][i] for k in range(m)]
            lambda_a_Nj = prod_interno(vetor_lambda, N_coluna_i)
            c_Nrelativo = c_N[i] - lambda_a_Nj
            if c_Nrelativo <= c_relativo:
                c_relativo = c_Nrelativo
                menor_N = i
# Passo 3) Teste de otimalidade (Nossa condição de parada principal que diz se encontramos nossa solução ótima):
        if c_relativo >= -1e-9: # Por questões vistas em cálculo numérico de precisão de máquina, usaremos um valor infimo invés de 0 para
            # garantir que imprecisões não afetem muito o algoritmo. Faremos isso em todas comparações pra 0.
            fx_otimo = prod_interno(c_B, x_B)
            print(B_indice)
            print(x_B)
            print(fx_otimo)
# As saídas sairão de maneira que a primeira linha nos diz qual a varíavel utilizada na base e em seguida seu valor. A terceira linha será o valor
# da função objetivo nesse ponto.
            return x_B
# Passos 4) e 5) Faremos o cálculo da direção Simplex e verificação de direção e determinaremos quem sai da base:
        N_coluna_menor = [N[k][menor_N] for k in range(m)]
        y = Gauss_piv(B, N_coluna_menor)
        limitado = False
        for yi in y:
            if yi > 1e-9:
                limitado = True
        if limitado == False:
            print("A solução factivel é ilimitada, portanto a solução é infinita")
            return
        menor_epsilon = float('inf')
        menor_B = 0
        for j in range(m):
            if y[j] > 1e-9:
                epsilon = x_B[j]/y[j]
                if epsilon < menor_epsilon:
                    menor_epsilon = epsilon
                    menor_B = j
# Passo 6) Vamos trocar a j-esima coluna de B pela i-esima coluna de N (Fazemos isso também nos nossos vetores de índices):
        for i in range(m):
            B_copia = B[i][menor_B]
            B[i][menor_B] = N[i][menor_N]
            N[i][menor_N] = B_copia
        c_copia = c_B[menor_B]
        B_indicecopia = B_indice[menor_B]
        B_indice[menor_B] = N_indice[menor_N]
        N_indice[menor_N] = B_indicecopia
        c_B[menor_B] = c_N[menor_N]
        c_N[menor_N] = c_copia
# e o algoritmo se repete até determinar que a região é ilimitada ou o menor custo relativo seja não-negativo, mostrando que chegamos na solução ótima.
