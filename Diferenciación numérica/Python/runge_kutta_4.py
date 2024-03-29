import numpy as np
from sympy import *
import matplotlib.pyplot as plt


"""
Función del metodo de runge_kutta para diferenciacion numerica
Entradas:
    funcion: función a trabajar
    intervalo: intervalo a trabajar
    pasoh: número de puntos a utilizar
    yinicial: imagen inicial para ejecutar el algoritmo
Salidas:
    (xv, yv): Pares ordenados
    polinomio: Respectivo polinomio de interpolacion
    
"""
def runge_kutta_4(funcion, intervalo, pasoh, yinicial):
    x = Symbol('x')
    y  = Symbol('y')
    
    f_s = sympify(funcion) #simbolica

    n = pasoh

    f = lambdify([x, y], f_s)

    a = intervalo[0]
    b = intervalo[1]

    h = (b-a)/(n-1)
    xv = np.linspace(a, b, n)
    
    yv = np.empty(n)
    yv[0] = yinicial

    for i in range(n-1): 
        k1 = f(xv[i],yv[i])
        k2 = f(xv[i] + h/2, yv[i] + h*k1/2)
        k3 = f(xv[i] + h/2, yv[i] + h*k2/2)
        k4 = f(xv[i]+ h, yv[i] + h*k3)
        yv[i+1] = yv[i] + h/6 * (k1 + 2*k2 + 2*k3 + k4)

    p = lagrange(xv, yv) #Crea polinomio de interpolación
    p_n = lambdify(x, p)



    #Formato de la línea y los puntos
    plt.scatter(xv, p_n(xv), color='r',zorder=1) 
    plt.plot(xv, p_n(xv), color='b',zorder=2)

    #Especifiaciones de gráfica
    plt.title("Gráfica polinomio de interpolación")
    plt.xlabel("x")
    plt.ylabel("p(x)")

    #Muestra gráfica
    plt.show()

    return (xv, yv), p
    

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
    intervalo = [0, 1]
    num_pt = 11
    funcion = "-x*y+4*x/y"
    yinicial = 1

    res = runge_kutta_4(funcion, intervalo, num_pt, yinicial)

    print("Pares ordenados: ", res[0])
    print("Polinomio de interpolacion: ", res[1])