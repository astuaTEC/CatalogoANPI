import numpy as np
import math

def matriz_d(A, n):
    d = np.zeros_like(A)
    for i in range(0, n):
        d[i, i] = A[i, i]
    return d

def matriz_l(A, n):
    d = np.zeros_like(A)
    for i in range(0, n):
        for j in range(0, i):
            d[i, j] = A[i, j]
    return d

def matriz_u(A, n):
    u = np.zeros_like(A)
    for i in range(0, n):
        for j in range(i + 1, n):
            u[i, j] = A[i, j]
    return u

def sust_adelante(A, b, n):
    soluciones = np.zeros((n, 1))
    soluciones[0] = b[0] / A[0, 0]#Primera solucion

    for i in range(n):
        sumatoria = 0
        for j in range(i):#Mueve columnas
            sumatoria = sumatoria + (A[i, j] * soluciones[j])

        x = (b[i] - sumatoria) / A[i, i]
        soluciones[i] = x

    return soluciones


def diagonalmenteDominate(A):
    for f in range(0, len(A)):
        d = A[f,f] #se obtiene el valor de dicha fila

        f_ = A[f,:] #resto de valores de la fila
        f_[0, f] = 0

        if abs(d) <= np.sum(f_):
            return False
    
    return True
    

def gauss_seidel(A, b, x0, tol, iterMax):
    n = len(A)
    D = matriz_d(A, n)
    L = matriz_l(A, n)
    U = matriz_u(A, n)

    print(L)
    print(D)
    print(U)

    y = sust_adelante(L+D, b, n)
    xk = x0

    for i in range(iterMax):
        zk = sust_adelante(L+D, U*xk, n)
        xk = -zk + y

        error = np.linalg.norm(A*xk - b)

        if error < tol:
            break
    
    return xk, error, 

if __name__ == '__main__':
    tol = 1e-15
    iterMax = 2000
    x0 = np.matrix('1; 1; 1; 1')
    A = np.matrix('25 15 -5 -10; 15 10 1 -7; -5 1 21 4; -10 -7 4 18')
    b = np.matrix('-25; -19; -21; -5')

    print(gauss_seidel(A, b, x0, tol, iterMax))

#A=np.matrix('-11 2 3; 4 12 6; 7 8 19')
#print(diagonalmenteDominate(A))