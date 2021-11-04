import numpy as np
from sympy import *

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
    #se definen los puntos en x y ya utilizar
    xv = [-2, 0, 1]
    yv = [0, 1, -1]

    print(lagrange(xv, yv))