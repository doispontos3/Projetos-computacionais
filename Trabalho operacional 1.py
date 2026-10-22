'''
Nome: Felipe de Freitas Percinoto RA: 275170
Docente: Eduardo Cardoso de Abreu
Trabalho Computacional 1 de MS211 - Calculo Numérico
'''
# Questão 1
x_1 = 0.5
x_2 = 0.11
n = 30000
# Vamos primeiro obter o valor exato do somatório para x_1 e x_2 e então comparar com o valor obtido por meio de sua iteração
Valor_exato_x1 = x_1 * n
Valor_exato_x2 = x_2 * n

# Agora montando o somatório por meio de n iterações:
S_1 = 0
S_2 = 0
for i in range(n):
    S_1 += x_1
    S_2 += x_2
# Agora vamos calcular os erros obtidos das iterações
Erro_absoluto1 = abs(Valor_exato_x1 - S_1)
Erro_absoluto2 = abs(Valor_exato_x2 - S_2)
Erro_relativo1 = Erro_absoluto1/Valor_exato_x1
Erro_relativo2 = Erro_absoluto2/Valor_exato_x2

print("O Valor exato de S_1 é {0}, a iteração retorna {1}, com erro absoluto de {2} e erro relativo de {3}".format(Valor_exato_x1, S_1, Erro_absoluto1, Erro_relativo1))
print("O Valor exato de S_2 é {0}, a iteração retorna {1}, com erro absoluto de {2} e erro relativo de {3}".format(Valor_exato_x2, S_2, Erro_absoluto2, Erro_relativo2))

# A diferença dos erros é oriundo da forma como 0.5 e 0.11 são escritos na base binária, a base utilizada pelo computador. 
# 0.5 = 0.1 na base binária, enquanto 0.11 possui uma representação infinita periódica, e a abreviação dos digitos leva a esse erro acumulado nas somas.
# Quanto mais digitos abreviados da representação periódica, maior o erro acumulado (e menor o poder computacional utilizado).

# Questão 2
# a)
x0 = 0.5
e1 = 1*(10**(-4))
e2 = e1
xk = x0
xk_1= 0
def fxk(xk):
    fx = xk**3 - 9*xk + 3
    return fx
dfxk = 3*xk**2 - 9
kn = 0
while abs(fxk(xk)) > e1 and abs(xk - xk_1)/abs(xk) > e2:
    kn+=1
    xk_1= xk
    xk = xk - fxk(xk)/dfxk
    dfxk = 3*xk**2 - 9
print("O valor aproximado da raiz de fxk é {0}, com um total de {1} iterações realizadas.".format(xk, kn))

# b)
x0 = 0
x1 = 1
e1 = 5*(10**(-4))
e2 = e1
def secante(x0, x1):
    xk = x1
    xk_1 = x0
# Aqui no método da secante, a principal diferença se encontra no fato de que utilizamos uma aproximação local da derivada invés da derivada exata,
# para a função do problema é desnecessário essa aproximação uma vez que ela é contínua, fácil de calcular e fica longe de se anular na vizinhança
# da raiz de fxk. Para determinar x0 e x1, usamos o teorema do valor intermediário, que irá nos garantir que se f(x0)f(x1) < 0 => há uma raiz em [x0,x1].
    dfxk = (fxk(xk) - fxk(xk_1))/(xk - xk_1)
    k = 0 
    while abs(fxk(xk)) > e1 and abs(xk - xk_1)/abs(xk) > e2:
        k+=1
        xk_1a = xk
        xk = xk - fxk(xk)*((xk - xk_1)/(fxk(xk)-fxk(xk_1)))
        xk_1 = xk_1a
        dfxk = (fxk(xk) - fxk(xk_1))/(xk - xk_1)
    return xk, k
xk, k= secante(x0,x1)
print("Pelo método da secante, obtemos a raiz aproximada {0} com um total de {1} iterações.".format(xk, k))
# Para outros valores iniciais de x0 e x1, precisamos apenas garantir que há uma raiz entre x0 e x1, isto é, f(x0)f(x1) < 0.
# Para f(x), podemos verificar com um gráfico que teremos raízes em [0,1], [-4,-3] e [2,3]. Tomando então x0 = -4 e x1= -3 e então 
# x0 = 2 e x1 = 3. Para averiguar a existencia de raízes, podemos montar diferentes intervalos onde xk+1= xk + 1 e analisar f(xk)f(xk+1) 
# e ver se há alguma troca de sinal. Observada a troca de sinal, tomamos os dois valores que limitam o intervalo como nossos x0 e x1 para 
# obter nossa raiz.
x0 = -4
x1 = -3
xk, k= secante(x0,x1)
print("Pelo método da secante, obtemos a raiz aproximada {0} com um total de {1} iterações.".format(xk, k))
x0 = 2
x1 = 3
xk, k= secante(x0,x1)
print("Pelo método da secante, obtemos a raiz aproximada {0} com um total de {1} iterações.".format(xk, k))
# c)
# O algoritmo da bissecção trabalha com base no teorema do valor intermediário, isto é, se f(x) muda de sinal em um intervalo contínuo [a,b],
# então existe ao menos uma raiz no intervalo. Observamos pelos metódos anteriores que a raiz de f(x) pertence ao intervalo [0,1]. A convergência da 
# bissecção é certa se a raiz estiver num intervalo compacto [a,b], porém o metódo realiza muitas iterações e pode ser muito custoso em relação aos 
# anteriores.
a = 0
b = 1 
xk = 1
kb = 0
e1 = 1*(10**(-4))
e2 = e1
xk_1 = 0
while abs(fxk(xk)) > e1 and abs(xk - xk_1)/abs(xk) > e2:
    kb+=1
    if (fxk(a)*fxk(b)) < 0:
        xk = (a+b)/2
    if (fxk(a)*fxk(xk)) < 0:
        b = xk
    else:
        a = xk
g = kb - kn
print("Pelo método da bissecção, obtemos a raiz aproximada {0} com um total de {1} iterações, um total de {2} a mais em relação ao método de Newton-Raphson.".format(xk, kb, g))

# d) O uso dos dois critérios de parada permite uma aproximação mais confiável das iterações, o critério e1 garante a convergência correta para a raiz 
# de f(x), já e2 garante e estabilidade da sequência de xk, evitando oscilações em torno da raiz e uma parada prematura na iteração.

# Questão 3
# a) O objetivo do algoritmo de Gauss com Pivoteamento Parcial é transformar uma matriz normal em uma matriz triangular superior, que é mais fácil de se
# se trabalhar e manipular.

Ax = [[2, 2, 1, 1], [1, -1, 2, -1], [3, 2, -3, -2], [4, 3, 2 ,1]]
B = [7, 1, 4, 12]

# Iremos, primeiramente, verificar o pivô de cada linha e então averiguar que o pivô com maior absoluto esteja na primeira linha:
def Gauss_piv(Ax, B):
    n = len(Ax)
    e1 = 1*(10**(-15)) #Vamos definir um limite de significancia para os valores de B.
    for k in range(n):
        p = 0
        linha = 0
        trocou = False
        for i in range(n-k):
            if abs(Ax[i+k][k]) > p:
                p = abs(Ax[i+k][k])
                linha = i + k
                trocou = True
        #Botando o maior pivô de uma linha em sua devida posição:
        if trocou == True:
            Pivo = Ax[linha]
            pivob = B[linha]
            Ax[linha] = Ax[k]
            B[linha] = B[k]
            Ax[k] = Pivo
            B[k] = pivob
        #Zerando os elementos da coluna abaixo do pivô de sua respectiva linha:
        for i in range(n-k-1):
            m = Ax[i+k+1][k]/Ax[k][k]
            B[i+k+1] = B[i+k+1] - m * B[k]
            if abs(B[i+k+1]) < e1:
                B[i+k+1] = 0 
            for j in range(len(Ax[i])):
                Ax[i+k+1][j] = Ax[i+k+1][j] - m * Ax[k][j]
                if abs(Ax[i+k+1][j]) < e1:
                    Ax[i+k+1][j] = 0 
    x=[]
    for i in range(n):
        x.append(0)
    #Obtendo o valor das variáveis:
    for k in range(n):
        somatorio = 0
        for i in range(n):
            somatorio += Ax[n-k-1][n-i-1] * x[n-i-1]
        x[n-k-1] = (B[n-k-1] - somatorio)/Ax[n-k-1][n-k-1]
    return Ax, B, x
Ax, B, x= Gauss_piv(Ax, B)      
print(Ax)
print(B)
print(x)
# b) Vamos primeiramente programar a matriz de Hilbert
def Hilbert(n):
    Hn = []
    Bn = []
    for i in range(n):
        linha = []
        Bi = 0 
        for j in range(n):
            linha.append(1/(i+j+1))
            Bi += 1/(i+j+1)
        Hn.append(linha)
        Bn.append(Bi)
    return Hn, Bn
Hn, Bn = Hilbert(5)
Hn, Bn, xHn = Gauss_piv(Hn, Bn)
print(Hn, Bn, xHn)
Hn, Bn = Hilbert(10)
Hn, Bn, xHn = Gauss_piv(Hn, Bn)
print(Hn, Bn, xHn)
Hn, Bn = Hilbert(15)
Hn, Bn, xHn = Gauss_piv(Hn, Bn)
print(xHn)
Ax = [[2, 2, 1, 1], [1, -1, 2, -1], [3, 2, -3, -2], [4, 3, 2 ,1]]
B = [7, 1, 4, 12]
# Observe que como a matriz é extremamente mal condicionada e cheia de números com dizimas periódicas, diversos arredondamentos devem ser feitos 
# para que a máquina consiga fazer as iterações numéricas. Esses arredondamentos geram perdas de precisão que, conforme n aumenta, se tornam cada
# vez maiores, ao ponto que o algoritmo falhe devido os sucessivos arredondamentos. Fazendo alguns testes, é possível observar que n = 11 é o último 
# valor para n que as incógnitas mantem um valor próximo de 1 com uma imprecisão de no máximo +- 0.01. Para n = 12, os valores de x já dispersam 
# e começam a divergir da solução analítica = (1, 1, ..., 1)^T por +-0.3 devido aos erros de arredondamento acumulados. Apartir de n = 42 minha máquina
# desistiu de dar prosseguimento no algoritmo e começou a dividir por 0. Os erros se tornam tão significativos que para n = 20 é obtida uma divisão 
# por 0, porém para n = 25 esse erro já não ocorre pois o erro extra em relação a n = 20 cancela a divisão por 0.

# c) Primeiro passo do algoritmo é verificar se a matriz segue as condições de Gauss-Jacobi: |aᵢᵢ| ≥ ∑|aᵢⱼ| com i =/= j e aᵢᵢ =/= 0 
# para todo i = 1, ..., n. Esses critérios garantem a convergência do método, se o primeiro critério falhar, não há garantia de convergência ou a 
# convergência pode ser muito lenta.
def diag_teste(Ax):
    n = len(Ax)
    for k in range(n):
        p = 0
        linha = 0
        trocou = False
        for i in range(n-k):
            if abs(Ax[i+k][k]) > p:
                p = abs(Ax[i+k][k])
                linha = i + k
                trocou = True
        #Botando o maior pivô de uma linha em sua devida posição:
        if trocou == True:
            Pivo = Ax[linha]
            Ax[linha] = Ax[k]
            Ax[k] = Pivo
    diagonal_dom = 0
    diag_quase_dom = 0
    for i in range(n):
        soma = 0
        for j in range(n):        
            if j != i:
                soma += abs(Ax[i][j])  
        if abs(Ax[i][i]) > soma:
            diagonal_dom +=  1
        elif abs(Ax[i][i]) == soma:
            diag_quase_dom += 1
    return diagonal_dom, diag_quase_dom
diag_dom, diag_q_dom = diag_teste(Ax)
Hn, Bn = Hilbert(5)
diag_dom_h, diag_q_dom_h = diag_teste(Hn)
print (diag_dom, diag_q_dom)
print (diag_dom_h, diag_q_dom_h)
# Feito o teste, quanto mais linhas com diagonais dominantes, mais provável é que o algoritmo convirja. Note, porém, que a matriz inicial e a matriz de Hilbert possuem 0 
# diagonais dominantes. Além do mais, conforme n cresce, os coeficientes ann -> 0 na matriz de Hilbert, falhando em ambos critérios de Gauss-Jacobi.
def Gauss_Jacobi(Ax, B, max):
    e1 = 1*(10**(-4))
    n = len(Ax)
    x = [0] * n
    xk = [0] * n
    erro_max = 0
    for k in range(max):
        for i in range(n):
            soma = 0
            for j in range(n):
                if j != i:
                    soma += Ax[i][j]*x[j]
            xk[i] = (B[i]-soma)/Ax[i][i]
            erro = abs(x[i] - xk[i])
        if erro > erro_max:
            erro_max = erro
        if erro_max < e1:
           return x, erro_max, True
        x = xk
    return x, erro_max, False
# O algoritmo possui dois critérios de parada, um e1 que para o algoritmo quando atinge o critério de tolerância, e um de número maximo de iterações. Esse segundo é importante
# quando estamos lidando com matrizes que falham no critério de diagonais dominantes, pois não há garantia de convergência para essas matrizes.
Ax = [[2, 2, 1, 1], [1, -1, 2, -1], [3, 2, -3, -2], [4, 3, 2 ,1]]
B = [7, 1, 4, 12]
max = 15
x, erro_max, convergiu = Gauss_Jacobi(Ax, B, max)
x_hn, erro_max_hn, convergiu_hn = Gauss_Jacobi(Hn, Bn, max)
print(x, erro_max, convergiu)
print(x_hn, erro_max_hn, convergiu_hn)

# d) A principal diferença de Gauss-Seidel para Gauss-Jacobi está na formula iterativa. Enquanto Gauss-Jacobi utiliza a iteração anterior das variáveis 
# do sistema, Gauss-Seidel utiliza a iteração atual, ou seja, suponhamos que estejamos calculando a iteração 2 de x_2, Gauss-Jacobi utilizará os valores
# da iteração 1 pra todas variáveis, já Gauss-Seidel utizará o valor da iteração 2 da variável x_1, os valores da iteração 2 de x_1 e x_2 para calcular
# x_3 até x_n. Como em Gauss-Seidel nós imediatamente utilizamos valores mais recentes em suas iterações, o algoritmo é um pouco mais dificil de ser
# implementado, porém é mais rápido que Gauss-Jacobi e converge em menos iterações. Gauss-Seidel também não é capaz de guardar dois vetores de uma 
# iteração k e k+1, então o calculo de erro exige a criação de uma cópia das iterações anteriores para ser efetuada, enquanto que em Gauss-Jacobi os dois
# vetores já estão sempre salvos.

# e) Sim, matrizes mal condicionadas, como a matriz de Hilbert, geram erros cada vez maiores que se acumulam a cada iteração e divergem. Métodos mais 
# robustos como a própria eliminação de Gauss divergem devidos esses erros, o mesmo ocorre em uma taxa ainda maior para metódos mais restritos e sensíveis
# Gauss-Jacobi e Gauss-Seidel.

# Questão 4)
# a) Note que o método de Newton é a versão do algoritmo de Newton-Raphson aplicado em um sistema não linear com n equações e n icógnitas
# Primeiramente, vamos definir a matriz contendo as funções $F(x)=0$ onde $x = [x_1, ..., x_n]$. Em seguida, no lugar de f'(x) calculado em 
# Newton-Raphson, utilizaremos a matriz jacobiana $J(x)$ que carrega as derivadas parciais de cada $x_i$ em relação a cada equação j com $i, j = 1, ..., n$ 
def F(x_1, x_2):
    f1 = 4*x_1 - x_1**3 + x_2
    f2 = (-x_1**2)/9 + (4*x_2 - x_2**2)/4 + 1
    FX = [f1, f2] 
    return FX
def Jacob(x_1, x_2):
    df1x_1 = 4 - 3*x_1**2
    df1x_2 = 1
    df2x_1 = (-2*x_1)/9
    df2x_2 = 2 - x_2/2 
    J = [[df1x_1, df1x_2], [df2x_1, df2x_2]]
    return J
# Definidas ambas funções, vamos adaptar nosso algoritmo de Newton-Raphson feito em 2) a) para receber matrizes e vetores.
def norma(Vetor):
    n = len(Vetor)
    Norma = 0
    for i in range(n):
        Norma += (Vetor[i])**2
    Norma = Norma**(1/2)
    return Norma
xk = [0] * len(F(x_1, x_2))
e1 = 10**(-4) 
max_iter = 500
def Newton(max_iter, e1, xk):
    n = len(xk)
    xk_1 = xk
    u = 1
    v = 1
    k = 0
    while u > e1 and v > e1:
        k += 1
        xk_1 = xk
        Fxk_1 = F(xk_1[0],xk_1[1])
        Norma_Fx = norma(Fxk_1)
        Jxk = Jacob(xk_1[0],xk_1[1])
        Fxk_1 = [(-i) for i in Fxk_1]
    # Vamos utilizar nossa função de Eliminação de Gauss para resolver o sistema &\delta xJ(x_k) = -F(x_k)
        AX, B, delta_x = Gauss_piv(Jxk, Fxk_1)
        xk = []
        for j in range(len(delta_x)):
            xk.append(xk_1[j] + delta_x[j])
        if k >= max_iter:
            break
        u = abs(norma(xk) - norma(xk_1)) 
        v = norma(F(xk[0], xk[1]))
    return xk, k
xk, k = Newton(max_iter, e1, xk)
print(xk, k)
# Para $e_1$ obtemos um total de 11 iterações para obter nossas raízes.
# b) 
xk = [0,2]
xk, k = Newton(max_iter, e1, xk)
print(xk, k)

xk = [-5,5]
xk, k = Newton(max_iter, e1, xk)
print(xk, k)

xk = [1000000,1000000]
xk, k = Newton(max_iter, e1, xk)
print(xk, k)

# c)
xk = [0,0]
e1 = 10**(-12)
xk, k = Newton(max_iter, e1, xk)
print(xk, k)

# d)
xk = [0,0]
e1 = 10**(-25)
xk, k = Newton(max_iter, e1, xk)
print(xk, k)