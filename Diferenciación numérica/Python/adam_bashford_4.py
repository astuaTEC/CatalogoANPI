import numpy as np
import matplotlib.pyplot as plt
from sympy import *



# f: Funcion del tipo f(x,y)
#     #Entradas:
#	 interv: Lista con los valores frontera en x con la forma: [xi, xf]
#    h: tamaño de paso 
#    yo: Valor inicial de la solucion para y(x=0) 

#Salidas:
# 	 xk: coordenadas x evaluadas
# 	 yk: coodernadas y obtenidas
# 	 poly: polinomio de interpolacion en formato string
def adam_bashford_4(f, interv, pasoh, yo):

    
    a = interv[0]
    b = interv[1]

    n = pasoh 

    h = (b-a)/(n-1)

    xk = np.linspace(a,b,n)
    yk = np.zeros(n)

    x = symbols('x')
    y = symbols('y')
    f_s = sympify(funcion) #simbolica
    f = lambdify([x, y], f_s)

    yk[0] = yo
    yk[1] = yk[0] + h*f(xk[0], yk[0])
    yk[2] = yk[1] + h/2*(3*f(xk[1], yk[1]) - f(xk[0], yk[0]))
    yk[3] = yk[2] + h/12*(23*f(xk[2], yk[2]) - 16*f(xk[1], yk[1]) + 5*f(xk[0], yk[0]))
    for i in range(3, n-1):
        yk[i+1] = yk[i] + h/24*(55*f(xk[i], yk[i]) - 59*f(xk[i-1], yk[i-1]) + 37*f(xk[i-2], yk[i-2]) - 9*f(xk[i-3], yk[i-3]))
    
    p = lagrange(xk, yk) #Crea polinomio de interpolación
    p_n = lambdify(x, p)


    #Formato de la línea y los puntos
    plt.scatter(xk, p_n(xk), color='r',zorder=1) 
    plt.plot(xk, p_n(xk), color='b',zorder=2)

    #Especifiaciones de gráfica
    plt.title("Gráfica polinomio de interpolación")
    plt.xlabel("x")
    plt.ylabel("p(x)")

    #Muestra gráfica
    plt.show()

    return (xk, yk), p
    


"""
Función para interpolar una función mediante
el método de la Interpolación de Lagrange
Entradas: 
    xv: Puntos en el eje x usados para interpolar
    yv: Puntos en el eje y usados para interpolar
Salidas:
    p: Polinomio en formato simbólico
        que interpola a los puntos ingresados
"""
def lagrange(xv, yv):
    n = len(xv) #Largo del vector

    if (len(xv) != len(yv)): #Comprueba que los vectores sean iguales
        return ('Los venctores no cumplen con la condicion de tamano')

    p = 0
    for k in range(n):
        p = p + yv[k]*fun_Lk(xv, k)

    p = expand(p)
    return p

"""
Función utilizada para calcular el valor
de Lk según se indica en el método de Lagrange
Entradas:
    xv: vector de puntos en x
    k: iteración correspondiente
Salidas:
    Lk: polinomio Lk en formato simbólico
"""
def fun_Lk(xv, k):
    x = Symbol('x')
    n = len(xv)
    Lk = 1
    for j in range(n):
        if j != k:
            Lk = Lk*(x - xv[j])/(xv[k] - xv[j])

    Lk = expand(Lk)
    return Lk



if __name__ == '__main__':
    intervalo = [2, 4]
    num_pt = 11
    funcion = "1+(x-y)**2"
    yinicial = 1

    res = adam_bashford_4(funcion, intervalo, num_pt, yinicial)

    print("Pares ordenados: ", res[0])
    print("Polinomio de interpolacion: ", res[1])