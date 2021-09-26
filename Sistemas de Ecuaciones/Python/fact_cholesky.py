from sustituciones import *
import numpy as np
import math

def verificar_positiva_definida(A, n):
    tol = np.finfo(float).eps
    for i in range(n):
        Ak = A[0:i, 0:i]
        deter = np.linalg.det(Ak)
        if deter < tol:
            return False
    return True

def verificar_simetria(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)


def fact_cholesky(A, b):
    n = len(A)
    if (not(verificar_simetria(A)) or not(verificar_positiva_definida(A, n))):
        return "Debe ser una matriz positiva definida y simetrica"

    L = np.zeros((n, n))
    j = 0
    L[j, j] = math.sqrt(A[j, j])

    for i in range(1, n):
        L[i, j] = A[i, j]/L[j, j]

    for j in range(1, n):
        aux1 = 0
        for k in range(0, j):
            aux1 += (L[j, k])**2
        
        L[j, j] = math.sqrt(A[j, j] - aux1)

        for i in range(j+1, n):
            aux2 = 0
            for k in range(0, j):
                aux2 += L[i, k] * L[j, k]
            
            L[i, j] = (A[i, j] - aux2)/L[j, j]


    y = sust_adelante(L, b, n) #se resuelve Ly = b
    x = sust_atras(L.T, y, n)
    return x

if __name__ == '__main__':
    A = np.eye(250,250,k=-1) + np.eye(250,250)*4 + np.eye(250,250,k=1)
    b =  np.full((250, 1), 3, dtype=float)
    b[0] = 2.5
    b[-1] = 2.5

    #A = np.matrix('25 15 -5 -10; 15 10 1 -7; -5 1 21 4; -10 -7 4 18')
    #b = np.matrix('-25; -19; -21; -5')

    print(fact_cholesky(A, b))