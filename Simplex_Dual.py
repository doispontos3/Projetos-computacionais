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
def mat_transposta(A):
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

# Vamos computar nossas entradas na forma primal e depois converte-las para a forma dual:
m = int(input()) # Número de restrições (linhas);
n = int(input()) # Número de variáveis (colunas);
c = [float(custo) for custo in input().split()] # Vetor dos custos na funcao objetivo a ser minimizada no formato "c_1 c_2 ... c_n";
b = [float(recurso) for recurso in input().split()] # Vetor dos recursos no formato "b_1 b_2 ... b_m";
A = [] # Matriz com os coeficientes das nossas restrições.
for i in range(m):
    A_linha = [float(a_ij) for a_ij in input().split()] # Linha i da matriz com os coeficientes das restrições no formato "a_i1 a_i2 ... a_in"
    A.append(A_linha)
# Nesse código vamos implementar a versão dual do problema. Primeiro passo é converter o problema pra sua forma dual e, logo em seguida, executar o simplex dual.
def primal_para_dual(m, n, c, b, A):
    A = mat_transposta(A)
    m_copia = m
    m = n
    n = m_copia
    b_copia = b
    b = c
    c = b_copia
    return
Dual = primal_para_dual(m, n, c, b, A)
# Vamos implementar a função Simplex_dual. Ela será uma alteração da de simplex_dual entregue anteriormente.
def Simplex_Dual(m, n, c, b, A):
    PARE = False
# Vamos montar nossa partição básica e não-básica inicial (assumindo que a base dual inicial dada seja factivel):
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
#Passo 1) callculo da solução básica dual (vetor_lambda) e lista de custos relativos:
        B_transposto = mat_transposta(B)
        vetor_lambda = Gauss_piv(B_transposto, c_B)
        c_Nrelativo_lista = []
        for j in range(n-m):
            N_coluna_j = [N[k][j] for k in range(m)]
            lambda_a_Nj = prod_interno(vetor_lambda, N_coluna_j)
            c_Nrelativo = c_N[j] - lambda_a_Nj
            c_Nrelativo_lista.append(c_Nrelativo)
#Passo 2.1) calculo da solução básica primal (x_B):
        x_B = Gauss_piv(B, b_B)
#Passo 2.2) determinaremos qual variável a sair da base e o teste de otimalidade:
        menor_x_B = float('inf')
        menor_B = 0
        for i in range(m):
            if x_B[i] < menor_x_B:
                menor_x_B = x_B[i]
                menor_B = i
        if menor_x_B >= -1e-9:
            fx_otimo = prod_interno(c_B, x_B)
            print(B_indice)
            print(x_B)
            print(fx_otimo)
            return x_B
#Passo 3) calculo da direção dual simplex:
        menos_e_l = [0.0] * m
        menos_e_l[menor_B] = -1.0
        eta_l = Gauss_piv(B_transposto, menos_e_l)
#Passo 4) vamos determinar o passo e a variável a entrar na base:
        factivel = False
        menor_delta = float('inf')
        menor_N = 0
        for j in range(n-m):
            N_coluna_j = [N[k][j] for k in range(m)]
            eta_a_Nj = prod_interno(eta_l, N_coluna_j)
            
            if eta_a_Nj > 1e-9:
                factivel = True
                delta = c_Nrelativo_lista[j] / eta_a_Nj
                if delta < menor_delta:
                    menor_delta = delta
                    menor_N = j
        if factivel == False:
            print("O problema primal é infactível, logo o dual é ilimitado")
            return False
#Passo 5) iremos trocar a menor_B-ésima coluna de B pela menor_N-ésima coluna de N:
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
Resposta = Simplex_Dual(m, n, c, b, A)