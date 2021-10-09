import numpy as np
from sympy import *

def lagrange(xv, yv):
    x = Symbol('x')

    n = len(xv) #Largo del vector

    if (len(xv) != len(yv)): #Comprueba que los vectores sean iguales
        return ('Los venctores no cumplen con la condicion de tamano')

    p = 0
    for k in range(n):
        p = p + yv[k]*fun_Lk(xv, k)

    p = expand(p)
    return p

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
    xv = [-2, 0, 1]
    yv = [0, 1, -1]

    print(lagrange(xv, yv))