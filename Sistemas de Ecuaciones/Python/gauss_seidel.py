import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

"""
Método para calcular la matriz D en el método de Gauss-seidel
Entradas: A es una matriz
          n es un número (el tamanio de la matriz)
Salidas: la matriz D calculada
"""
def matriz_d(A, n):
    d = np.zeros_like(A)
    for i in range(0, n):
        d[i, i] = A[i, i]
    return d

"""
Método para calcular la matriz L en el método de Gauss-seidel
Entradas: A es una matriz
          n es un número (el tamanio de la matriz)
Salidas: la matriz L calculada
"""
def matriz_l(A, n):
    d = np.zeros_like(A)
    for i in range(0, n):
        for j in range(0, i):
            d[i, j] = A[i, j]
    return d

"""
Método para calcular la matriz U en el método de Gauss-seidel
Entradas: A es una matriz
          n es un número (el tamanio de la matriz)
Salidas: la matriz U calculada
"""
def matriz_u(A, n):
    u = np.zeros_like(A)
    for i in range(0, n):
        for j in range(i + 1, n):
            u[i, j] = A[i, j]
    return u

"""
Método que realiza la sustitución hacia adelante
en un sistema de ecuaciones
Entradas: A es una matriz triangular inferior
          b es el vector de terminos idependientes
Salida: soluciones un vector con las soluciones del sistema
"""
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

"""
Método para verificar si una matriz es
diagonalmente dominante o no
Entradas: A es una matriz cuadrada
Salidas: true si es diagonalmente dominante
         false si no es diagonalmente dominante
"""
def diagonalmenteDominate(A):
    for f in range(0, len(A)):
        d = A[f,f] #se obtiene el valor de dicha fila

        f_ = A[f,:] #resto de valores de la fila
        f_[0, f] = 0

        if abs(d) <= np.sum(f_):
            return False
    
    return True
    
"""
Función para resolver sistemas de ecuaciones
mediante el método iterativo de Gauss-seidel
Entradas: A es una matriz cuadrada
          b es el vector de terminos idependientes
          x0 es el vector inicial
          tol es la tolerancia
          iterMax son las iteraciones máximas
Salidas: xk es el valor de la aproximación
         error es el error asociado
"""
def gauss_seidel(A, b, x0, tol, iterMax):
    #se revisa el teorema de convergencia
    if (not(diagonalmenteDominate(A.copy()))):
        return "La matriz tiene que ser diagonalmente dominante"

    n = len(A)
    D = matriz_d(A, n)
    L = matriz_l(A, n)
    U = matriz_u(A, n)

    y = sust_adelante(L+D, b, n)
    xk = x0

    er = []
    k = 0

    for i in range(iterMax):
        k += 1
        zk = sust_adelante(L+D, U*xk, n)
        xk = -zk + y

        error = np.linalg.norm(A*xk - b)
        er.append(sp.N(error))

        if error < tol:
            break
    
    fig, graf = plt.subplots()
    ejeX = np.arange(1, k+1, 1)
    graf.plot(ejeX, er)
    plt.show()

    return xk, error


if __name__ == '__main__':
    tol = 1e-5
    iterMax = 2000
    #x0 = np.matrix('0; 0; 0; 0')
    #A = np.matrix('25 15 -5 -10; 15 10 1 -7; -5 1 21 4; -10 -7 4 18')
    #b = np.matrix('-25; -19; -21; -5')

    A2 = np.matrix('5 1 1; 1 5 1; 1 1 5')
    b2 = np.matrix('7; 7; 7')
    x0 = np.matrix('0; 0; 0')

    print(gauss_seidel(A2, b2, x0, tol, iterMax))