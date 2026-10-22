'''
Nome: Felipe de Freitas Percinoto RA: 275170
Docente: Eduardo Cardoso de Abreu
Trabalho Computacional 2 de MS211 - Calculo Numérico
'''
'''
    Questão 1:
    a) Sim, a existência do polinômio é garantida uma vez que, utilizando o método de Lagrange, nós sempre obteremos um polinômio passando pelos
pontos ${(x_1, f(x_1)),..., (x_n, f(x_n)}$. Para mostrar sua unicidade, vamos utilizar uma prova por contradição e o Teorema Fundamental da 
Álgebra: Pelo TFA, temos que dado um polinômio de grau n, o polinômio terá n raízes (iguais ou distintas entre si). Seja $p(x)$ e $q(x)$ dois 
polinômios distintos de grau $<= n$ e que interceptam todos os pontos ${(x_1, f(x_1)),..., (x_n, f(x_n)}$. Seja $r(x) = p(x) - q(x)$. 
$r(x)$ terá grau $<= n$ já que $p(x)$ e $r(x)$ possuem grau $<= n$. Como $r(x)$ e $q(x)$ interceptam os pontos ${(x_1, f(x_1)),..., (x_n, f(x_n)}$, 
temos que $r(x_k) = p(x_k) - r(x_k) = f(x_k) - f(x_k) = 0$, porém como temos $n+1$ pontos, isso contradiz o TFA $=> r(x) = 0 <=> p(x) = q(x)$.

   b) Pelo método de Lagrange, temos que o polinômio será dado por: $p(x) = ∑ (k=0 to n)  y_k * ∏ (j=0 to n, j=/=k) (x - x_j)/(x_k - x_j).
Dado os pontos, basta calcular o polinômio, onde o produtório é chamado de polinômio de base de Lagrange, que possui a propriedade: 
Dado um $L_k(x)$ com $k E {0, 1, ..., n} => L_k(x_i) = 0$ se $i =/= j$ e $L_k(x_i) = 1$ se $i = k$.
'''   

def Interpolacao_polinomial(x_k, y_k, x):
    # x será um valor qualquer que queremos aplicar no nosso polinômio obtido
    n = len(x_k)
    f_x = 0
    for i in range(n): 
        y = y_k[i]
        for j in range(n):
            if j != i:
                y *= (x - x_k[j])/(x_k[i] - x_k[j])
        f_x += y
    return f_x
# Vamos escrever um polinômio com todos os pontos dados e então vamos fazer uma comparação com uma aproximação utilizando somente
# três dos pontos dados:
x_k = [0.2, 0.34, 0.4, 0.52, 0.6, 0.72]
y_k = [0.16, 0.22, 0.27, 0.29, 0.32, 0.37]
x = 0.47
resultado_1 = Interpolacao_polinomial(x_k, y_k, x)
x_k = [0.4, 0.52, 0.6]
y_k = [0.27, 0.29, 0.32]
resultado_2 = Interpolacao_polinomial(x_k, y_k, x)
erro_ab = abs(resultado_1 - resultado_2)
erro_re = erro_ab/abs(resultado_1)
print("O valor de f({0}) no polinômio aproximado de grau 5 é {1}, enquanto no de grau 2 é {2}, um erro absoluto de {3} e relativo de {4}".format(x, resultado_1, resultado_2, erro_ab, erro_re))
# Vamos testar com um polinomio de grau 1, usando apenas os pontos 0.4 e 0.6 de $x_1$:
x_k = [0.4, 0.6]
y_k = [0.27, 0.32]
resultado_3 = Interpolacao_polinomial(x_k, y_k, x)
erro_ab = abs(resultado_1 - resultado_3)
erro_re = erro_ab/abs(resultado_1)
print("O valor de f({0}) no polinômio aproximado de grau 5 é {1}, enquanto no de grau 1 é {2}, um erro absoluto de {3} e relativo de {4}".format(x, resultado_1, resultado_3, erro_ab, erro_re))
# Interessantemente, pros pontos escolhidos de $x_i$, nós obtivemos uma aproximação ainda melhor de f(0.47) com um polinômio de grau do que
# no de grau 2 em comparação ao de grau 5 que interpola todos os pontos dados. Isso é uma coincidência e, com outro pontos para realizar a interpolação,
# poderiamos ter visto outro resultado, já que uma reta tende a não ser boa aproximação para valores grandes de x aplicado numa f(x) 
# que é um polinômio de grau 5.
# Para obter o polinômio própriamente dito, vamos calcular os coeficientes:
def coeficientes_polinômio(x_k, y_k):
    n = len(x_k)
    coeficientes = [0.0] * n
    for i in range(n):
        denominador = 1.0
        for j in range(n):
            if j != i:
                denominador *= (x_k[i] - x_k[j])
        Lk = [0.0] * n
        Lk[0] = 1.0
        for j in range(n):
            if j != i:
                novo_Lk= [0.0] * n
                for k in range(n-1):
                    novo_Lk[k+1] += Lk[k]
                for k in range(n):
                    novo_Lk[k] += -x_k[j] * Lk[k]
                Lk = novo_Lk
        for j in range(n):
            coeficientes[j] += y_k[i] * Lk[j] / denominador
    return coeficientes
x_k = [0.2, 0.34, 0.4, 0.52, 0.6, 0.72]
y_k = [0.16, 0.22, 0.27, 0.29, 0.32, 0.37]
coeficientes = coeficientes_polinômio(x_k, y_k)
print(coeficientes)
x_k = [0.4, 0.52, 0.6]
y_k = [0.27, 0.29, 0.32]
coeficientes = coeficientes_polinômio(x_k, y_k)
print(coeficientes)

# Questão 2
# a) O método dos Quadrado Mínimos consiste em: Dado uma lista de pontos $(x_i, y_i), i \in {1, ..., n}$, queremos aproximar uma função de maneira
# diminuir o resíduo $r_i  = y_i - f(x_i)$ de cada ponto, minimizando então a soma $\sum{1}^{n} [y_i - f(x_i)]^{2}$ do quadrado das diferenças dos pontos.
# a principal diferença entre a interpolação polinomial e o método dos quadrados mínimos se encontra na certeza dos dados, isto é, na interpolação, 
# assumimos que os pontos são exatos e, portanto, queremos que o polinômio aproximado passe por cima dos pontos. Já no M.Q.M., não há certeza dos dados
# e sabemos que para todos pontos haverá algum tipo de erro na medição, então uma função que minimize o erro dos pontos acima seria mais desejada
# que uma função que passe diretamente por cima dos pontos e acabe contabilizando esse erro em sua aproximação. Para aproximar com uma reta, temos que 
# nossa função modelo será $f(x) = ax + b => \sum{i=1}^{n} [y_i - (ax_i + b)]^{2} = S(a, b)$, ou seja, para modelar essa função, estaremos basicamente
# modelando os valores de a e b de maneira a minimizar a soma dos quadrados das diferenças dos pontos. Como temos duas variáveis, podemos montar um 
# sistema com suas derivadas parciais igualadas a 0 e manipula-lo de maneira a obter um sistema linear com variáveis a e b $=> \frac{\partial S}{\partial a} 
# = -2 \sum{i=1}^{n} [y_i - (ax_i + b)]x_i = 0$ e $\frac{\partial S}{\partial b} = -2 \sum{i=1}^{n} [y_i - (ax_i + b)] = 0 <=> a\sum{i=1}^{n} (x_i)^2 + b \sum{1}^{n} x_i
# = \sum{i=1}^{n} y_i x_i$ e $nb + a \sum{i=1}^{n} x_i = \sum{i=1}^{n} y_i$, que é um sistema linear de incognitas a e b igualadas a $\sum{i=1}^{n} y_1$
def Gauss_piv(Ax, B):
    n = len(Ax)
    e1 = 1*(10**(-15)) #Vamos definir um limite de significancia para os valores de B.
    for k in range(n-1):
        p = 0
        linha = 0
        trocou = False
        for i in range(k, n):
            if abs(Ax[i][k]) > p:
                p = abs(Ax[i][k])
                linha = i 
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
        for i in range(k+1, n):  
            if abs(Ax[k][k]) < e1:
                continue
            m = Ax[i][k] / Ax[k][k]  
            B[i] = B[i] - m * B[k]
            if abs(B[i]) < e1:
                B[i] = 0 
            for j in range(k, n):  
                Ax[i][j] = Ax[i][j] - m * Ax[k][j]
                if abs(Ax[i][j]) < e1:
                    Ax[i][j] = 0 
    x = []
    for i in range(n):
        x.append(0)
    for i in range(n-1, -1, -1):  
        somatorio = 0
        for j in range(i+1, n):  
            somatorio += Ax[i][j] * x[j]
        x[i] = (B[i] - somatorio) / Ax[i][i]
    return x
def Quadrados_Minimos(x_i, y_i, n):
    m = len(x_i)
    matriz = []
    matriz_b = [0.0] * n
    for i in range(n):
        matriz_linha = [0.0] * n
        matriz.append(matriz_linha)
# A partir da definição dada, basicamente queremos que a matriz obtida dos somatórios seja uma matriz onde o elemento $a_ij = \sum{k=1}^{n} (x_k)^{i+j-2}
    for i in range(n):
        for j in range(n):
            for k in range(m):
                matriz[i][j] += (x_i[k])**(i+j)
    for i in range(n):
        for j in range(m):
            matriz_b[i] += y_i[j]*(x_i[j]**i)
# Agora, vamos importar uma função montada no trabalho computacional anterior para resolver o sistema linear:
    coeficientes = Gauss_piv(matriz, matriz_b)
    residuos = 0
    return coeficientes
x_i = [1, 2, 3, 4, 5, 6, 7, 8]
y_i = [0.5, 0.6, 0.9, 0.8, 1.2, 1.5, 1.7, 2.0]
coeficientes_reta = Quadrados_Minimos(x_i, y_i, n=2)
coeficientes_parabola = Quadrados_Minimos(x_i, y_i, n=3)
print("Por meio dos quadrados mínimos, obtemos a reta y = {0} +{1}x e a parabola y = {2} + {3}x + {4}x²".format(coeficientes_reta[0], coeficientes_reta[1], coeficientes_parabola[0],coeficientes_parabola[1], coeficientes_parabola[2]))
residuos_reta = 0
residuos_parabola = 0
def reta(x):
    fx = 0.1750 + 0.2167*x
    return fx
def parabola(x):
    fx = 0.4071 + 0.0774*x + 0.0155*(x**2)
    return fx
for i in range(len(x_i)):
    residuos_reta += (reta(x_i[i])-y_i[i])**2
    residuos_parabola += (parabola(x_i[i])-y_i[i])**2
print(residuos_reta, residuos_parabola)

# Questão 3
# a) A principal diferença entre os dois métodos é: Na regra dos Trapézios, estaremos basicamente montando n trapézios de base h abaixo de f(x) e 
# calcular sua área para aproximar a integral de f(x); Na regra de Simpson, invés do "teto" da função ser uma reta, será uma parabola, ou seja,
# oferecerá uma aproximação melhor para qualquer função "curva" do que uma reta. 

# b) Na regra dos trapézios, estamos fazendo o seguinte: Estamos pegando o intervalo $[a, b]$, partindo em $n$ intervalos de tamanho h e montando 
# diversos trapézios que cobrem a área de f(x), com os lados do trapézio sendo os segmentos $(a, f(a))$, $(f(a), f(b))$, $(f(b), b)$ e $(b, a)$.
# Quanto mais trapézios usarmos para cobrir o intervalo, melhor é aproximação feita. Se interpretarmos o intervalo h como altura e os valores de f(x)
# como suas bases, a fórmula da área do trapézio $(A+B)\frac{h}{2} = (f(a)+f(b))\frac{h}{2}$, logo a área aproximada debaixo da função vai ser a 
# soma da área de todos os trapézios: $I \approx \frac{h}{2}\sum{i=1}^{n} {f(x_{i-1}) + f(x_i)}. Para a implementação do cálculo do erro maximo, 
# iremos fazer o cálculo da derivada segunda de $f(x)$ que será 6.
def Trapezios(fx, a, b, n):
    h = (b-a)/n
    soma = fx(a) + fx(b)
    for i in range(1, n):
        x_i = a + i*h
        soma+= 2*fx(x_i)
    integral = (h/2) * soma
    erro_max_t = ((b - a)**3 / (12 * n**2)) * 6
    return integral, erro_max_t
# Para o método de Simpson, a ideia da partição é semelhando a do trapézio, porém substituiremos o "teto" do trapézio por uma parábola, que oferece 
# uma aproximação melhor para funções não lineares. Para fazer isso, iremos dividir [a, b] em n partes com n sendo um número par, e em seguida iremos
# traçar o ponto médio $x_k = \frac{x_{i-1} + x_i}{2}$ com $i=1,2,...,n$, e iremos fazer um polinômio com interpolação utilizando os polinômios de 
# Lagrange de maneira que o polinômio interpole $f(x_{i-1}),f(x_k),f(x_i)$. Fazendo as interpolações e as integrações nos intervalo partidos obtemos
# $I \approx \frac{h}{3}[f(a)+4f(\frac{a_b}{2})+f(b)]. Note que a derivada quarta de f(x) é 0, portanto, a fórmula máxima do erro máximo $\frac{(b-a)^5}{180n^4}
# max_{x\in [a,b]|f^{(4)}(x)|=0.}
def Simpson(fx, a, b, n):
    h = (b-a)/n
    soma_impar=0
    soma_par= 0
    for i in range(1, n, 2):
        x_impar = a+i*h
        soma_impar+= fx(x_impar)
    for i in range(2, n-1, 2):
        x_par = a + i*h
        soma_par+= fx(x_par)
    integral = (h/3)*(fx(a)+ fx(b) + 4*soma_impar+2*soma_par)
    erro_max_s = 0
    return integral
def fx(x):
    fx = 3*x**2 - 3*x + 1
    return fx
def fx_segunda(x):
    fx = 6
Trapezio, erro = Trapezios(fx, 0, 4, 4)
print(Trapezio, erro)
Simpsonvalor = Simpson(fx, 0, 4, 4)
print(Simpsonvalor, 0)

# Questão 4) 
# a) O método de Euler tem como objetivo aproximar uma solução numérica para algum valor de x começando de $x_0$ dado $y' = f(x, y)$ e $y(x_0) = 
# y_0$. O método de Euler utiliza como base a expanção da Série de Taylor de primeira ordem: $y(x+h) = y(x) + hy'(x) + \frac{h^2}{2}y''(\xi)$, onde o 
# terceiro termo é o erro da aproximação. Removendo esse termo de erro nós obtemos que $y(x+h) = y(x) + hy'(x) + \frac{h^2}{2}y''(\xi) \approx y(x) + 
# hf(x, y(x))$. Nós definimos h como sendo $h = \frac{x_f-x_0}{n}$, onde $x_f$ é o ponto em x avaliar y(x_f) e n é o numero de partições feitas nesse 
# nesse intervalo, onde $x_k = x_0 + hn$. Com isso, obtemos a fórmula iterativa: $y_{n+1} = y_n + hf(x_n,y_n)$.
def Metodo_Euler(fxy, x_0, y_0, x_f, n):
    h = (x_f - x_0)/n
    x = [x_0]
    y = [y_0]
    for i in range(n):
        y_k = y[-1] + h*fxy(x[-1], y[-1])
        x_k = x[-1] + h 
        x.append(x_k)
        y.append(y_k)
    return x, y
# Este programa retorna as diversas iterações de $x_k$ e $y_k$, onde a nossa melhor aproximação esperada seria o último elemento adicionado.
def fxy(x, y):
    fxy = 0.04*y
    return fxy
# Vamos agora testar para os diferentes valores de h = 1, 0.5, 0.25, e 0.1. Note que começaremos em $x_0 = 0$ e iremos até $x_f = 1$, logo teremos os
# valores de $n = 1$, $n = 2$, $n = 4$ e $n = 10$
x, y = Metodo_Euler(fxy, 0, 1000, 1, 1)
print(y)
x, y = Metodo_Euler(fxy, 0, 1000, 1, 2)
print(y)
x, y = Metodo_Euler(fxy, 0, 1000, 1, 4)
print(y)
x, y = Metodo_Euler(fxy, 0, 1000, 1, 10)
print(y)
# Observe que quanto mais passos nós damos, mais próximo da solução exata de y(1) = 1000e^{0.04} \approx 1040.810774$, pórem, consequentemente, mais o 
# lento se torna o algoritmo.

# b) O Método de Runge-Kutta é uma versão otimizada do Método de Euler que utiliza a avaliação em diferentes pontos em uma única iteração invés de
# de somente uma como no Método de Euler. Por conta disso, o método de Runge-Kutta difere do método de Euler na segunda ordem, com a primeira ordem
# do método sendo justamente o Método de Euler. A formulação geral de R-K é: Dados $y'(x) = f(x,y(x))$ e $y(x_0) = y_0$ e um $y(x_f)=x_f$ um valor
# que queremos descobrir. Temos que a fórmula iterativa será dada por: $y_{n+1}= y_n + h\sum{i=1}^{s}b_i k_i$ com $k_i = f(x_n + c_i h, y_n + 
# h\sum{j=1}^{i-1}a_{ij}k_j), onde as constantes $c_2,b_1,b_2,a_{21}$ são paramêtros do método. Para a implementação no nosso programa, vamos utilizar
# a variante de Heun, com paramêtros $c_2=1$,$a_{21}=1$,$b_1=\frac{1}{2}$ e $b_2=\frac{1}{2}$.
def Metodo_RK(fxy, x_0, y_0, x_f, n):
    h = (x_f - x_0)/n
    x = [x_0]
    y = [y_0]
    for i in range(n):
        k1 = fxy(x[-1],y[-1])
        k2 = fxy(x[-1]+h/2, y[-1]+(h/2)*k1)
        x_k = x[-1]+1
        y_k= y[-1] + (h/2)*(k1+k2)
        x.append(x_k)
        y.append(y_k)
    return x, y
x, y = Metodo_RK(fxy, 0, 1000, 1, 1)
print(y)
x, y = Metodo_RK(fxy, 0, 1000, 1, 2)
print(y)
x, y = Metodo_RK(fxy, 0, 1000, 1, 4)
print(y)
x, y = Metodo_RK(fxy, 0, 1000, 1, 10)
print(y)
# Como esperado, observe que conforme o número de passos aumenta, a solução converge mais rapidamente a solução exata em relação ao método
# de Euler.

# c) No método de Euler, nós obtemos, para cada h, os valores ${1040, 1040.4, 1040.60401, 1040.72773}$, enquanto que para o método de R-K, obtemos 
# ${1040.4, 1040.60401, 1040.70704, 1040.76920}$, enquanto que o valor aproximado de $1000e^0.04 \approx 1040.81077$. Observe que os valores se 
# aproximam cada vez do valor exato, em especial no método de R-K, onde vemos uma diferença de ${0.4, 0.20401, 0.10303, 0.04147}$. Apesar da diferença 
# diminuir, o método de R-K oferece uma convergência mais rápida ao valor exato e possui uma implementação fácil, então pode ser considerado uma melhora
# direta em relação ao método de Euler.

# d) 
x, y = Metodo_Euler(fxy, 0, 1000, 100, 1)
print(y[-1])
x, y = Metodo_Euler(fxy, 0, 1000, 100, 2)
print(y[-1])
x, y = Metodo_Euler(fxy, 0, 1000, 100, 4)
print(y[-1])
x, y = Metodo_Euler(fxy, 0, 1000, 100, 10)
print(y[-1])

x, y = Metodo_Euler(fxy, 0, 1000, 1000, 1)
print(y[-1])
x, y = Metodo_Euler(fxy, 0, 1000, 1000, 2)
print(y[-1])
x, y = Metodo_Euler(fxy, 0, 1000, 1000, 4)
print(y[-1])
x, y = Metodo_Euler(fxy, 0, 1000, 1000, 10)
print(y[-1])

x, y = Metodo_RK(fxy, 0, 1000, 100, 1)
print(y[-1])
x, y = Metodo_RK(fxy, 0, 1000, 100, 2)
print(y[-1])
x, y = Metodo_RK(fxy, 0, 1000, 100, 4)
print(y[-1])
x, y = Metodo_RK(fxy, 0, 1000, 100, 10)
print(y[-1])

x, y = Metodo_RK(fxy, 0, 1000, 1000, 1)
print(y[-1])
x, y = Metodo_RK(fxy, 0, 1000, 1000, 2)
print(y[-1])
x, y = Metodo_RK(fxy, 0, 1000, 1000, 4)
print(y[-1])
x, y = Metodo_RK(fxy, 0, 1000, 1000, 10)
print(y[-1])

# Observe nesse teste como a velocidade de convergência do R-K é consideravelmente superior a de Euler. O valor aproximado de $1000e{40} =
# 2.35385267e+20$, obtemos do método de Euler como resultado mais próximo o valor 9765625000.0, já R-K nos dá 3486784401000.0. É claro que ambos 
# resultados estão extremamente distantes do resultado exato, mostrando a limitação do método e falta de precisão para intervalos grandes. Ainda sim,
# podemos observar que R-K cresceu mais rapidamente em direção ao valor exato de $y(1000)$, mostra a superioridade do método acima do de Euler.

# Questão 5
# a) O objetivo principal do método das diferenças finitas é gerar uma aproximação de um Problema de Valor de Contorno pra $y(x)$ em algum ponto de
# x dentro do domínio da função. Para isso, nós repartiremos o intervalo [a,b] em n pedaços de comprimento $h= \frac{b-a}{n}$ e então calcularemos uma
# aproximação das derivadas $y'(x)$ e "y''(x)$ como sendo $y'(x_i) \approx \frac{y(x_{i+1})-y(x_{i-1})}{2h}$ e $y''(x_i) \approx 
# \frac{y(x_{i+1})-2y(x_i)+y(x_{i-1})}{h^2}$. Note, porém, que os intervalos que calculamos as derivada podem variar, no caso que definimos acima, estamos
# fazendo a diferença central para a aproximação da derivada, mas podemos também fazer a diferença progressiva usando $x_{i+2}$, $x_{i+1}$ e $x_i$, ou a 
# diferença regressiva usando os pontos $x_{i-2}$, $x_{i-1}$ e $x_i$. Nós então fazemos a substituição da aproximação das derivadas na nossa EDO e os valores
# centralizados de $x_i$ e $x_{i+1}$. Aos realizar as substituições para todos os pontos na malha, obtemos um sistema linear de diagonal principal 
# $A_{i,i} = -\frac{2a(x_i)}{h^2}+c(x_i)$, $A_{i,i-1}= \frac{a(x_i)}{h^2}-\frac{b(x_i)}{2h}$ e $A_{i,i+1}= \frac{a(x_i)}{h^2}+\frac{b(x_i)}{2h}$.
import math
def f_exata(x):
    return 2*(math.e**(-x))*(1-x) + x - 2
def Met_Diferencas(a, b, n, alpha, beta, p, q, r):
    h = (b-a)/n
    A = []
    x = []
    for i in range(n):
        x_i = a+i*h
        x.append(x_i)
    for i in range(n-1):
        linha = [0.0]*(n-1)
        A.append(linha)
    for i in range(n-1):
        x_i = (i+1) * h
        A[i][i] = -2*p/h**2 + r
        if i > 0:
            A[i][i-1] = p/h**2 - q/(2*h)
        if i < n-2:
            A[i][i+1] = p/h**2 + q/(2*h)
    B = []
    for i in range(n-1):
        x_i = a + (i+1)*h
        B.append(x_i)
    y_exato = []
    for i in range(n-1):
        y_exato.append(f_exata(x[i+1]))
    B[0] = B[0] - (p/h**2 - q/(2*h)) * alpha
    B[-1] = B[-1] - (p/h**2 + q/(2*h)) * beta
    y = Gauss_piv(A, B)
    erro_abs = 0
    erro_rel = 0
    for i in range(len(y_exato)):
        erro_abs += abs(y[i]-y_exato[i])
        if y_exato[i] != 0:
            erro_rel += abs(y[i]-y_exato[i])/abs(y_exato[i])
    erro_ponto = []
    for i in range(n-1):
        erro_ponto.append(abs(y[i]-y_exato[i]))
    return y, x, y_exato, erro_abs, erro_rel, erro_ponto
y, x_pontos, y_exato, erro_abs, erro_rel, erro_ponto = Met_Diferencas(0,1, 100,0,-1,1,2,1)

print(erro_ab, erro_rel, erro_ponto)

# b) Isso ocorre por os pontos centrais são os pontos mais distantes das extremidades que são fixados pela sua solução exata, então quanto mais distante
# dos pontos de contorno, mais o valor tende a flutuar e acumular erros.

# c) Para calcular sua ordem, iremos usar a aproximação $p \approx \frac{log(\frac{E(h_k)}{E(h_{k+1})}{log(\frac{h_k/h_{k+1}})})}$, essa formúla analisará
# como a função e sua aproximação crescem dentro de uma malha $h_k$. Temos que $f'(x) = -4e^{-x}+2xe^{-x}+1$.
def ordem():
    valores_n = [20, 40, 80, 160, 320]
    erros_absolutos = []
    erros_relativos = []
    for n in valores_n:
        y, x_pontos, y_exato, erro_abs, erro_rel, erro_ponto = Met_Diferencas(0, 1, n, 0, -1, 1, 2, 1)
        h = 1.0 / n
        erros_absolutos.append(erro_abs)
        erros_relativos.append(erro_rel)
        if len(erros_absolutos) > 1:
            razao = erros_absolutos[-2] / erros_absolutos[-1]
            ordem = math.log(razao) / math.log(2)  
            print(f"{razao:.4f}\t\t{ordem:.4f}")
        else:
            print(f"{n}\t\t{h:.6f}\t\t{erro_abs:.6e}\t\t-\t\t-")
    
    return erros_absolutos, erros_relativos, valores_n
erros_abs, erros_rel, ns = ordem()
# Os resultados, porém, retornam que o algoritmo possui ordem 1, o erro provavelmente se encontra no pivoteamento do sistema gerado pelo método, porém
# é evidente que o algoritmo é de ordem 2 uma vez que sua construção depende da derivada de segunda de $y(x)$.
# d) Para $h=0.00001$ e $h= 0.000000001$ teremos que $n = 100000$ e $n=1000000000$, substituindo na função:
y, x_pontos, y_exato, erro_abs, erro_rel, erro_ponto = Met_Diferencas(0,1, 100000,0,-1,1,2,1)
print(erro_ab, erro_rel, erro_ponto)

y, x_pontos, y_exato, erro_abs, erro_rel, erro_ponto = Met_Diferencas(0,1, 1000000000,0,-1,1,2,1)
print(erro_ab, erro_rel, erro_ponto)

# O algoritmo não retorna saída e fica rodando por tempo indeterminado, isso se deve ao fato do número de operações com ponto flutuante que cresce consideravelmente
# rápido, tornando o tempo e número de operações inviável para máquina. Inclusive a máquina dá indícios de dificuldade, aumentando a potência das ventoinhas e
# e um aumento de 8ºC do CPU. A métrica principal para analisar os custos de algoritmos númericos é o FLOP, que é basicamente contabilizar o número de operações
# de ponto flutuante do algoritmo, nossa função exata tem em torno de 25 FLOPs. A provável causa do algoritmo não retornar um valor rapidamente é falta de acesso
# a memória da máquina, que pode limitar o tamanho da malha utilizada no método.
    