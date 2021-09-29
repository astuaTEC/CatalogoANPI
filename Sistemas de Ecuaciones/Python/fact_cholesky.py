import numpy as np
import math

"""
Método que realiza la sustitución hacia adelante
en un sistema de ecuaciones
Entradas: A es una matriz triangular inferior
          b es el vector de terminos idependientes
Salida: soluciones es un vector con las soluciones del sistema
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
Método que realiza la sustitución hacia atrás
en un sistema de ecuaciones
Entradas: A es una matriz triangular superior
          b es el vector de terminos idependientes
Salida: soluciones un vector con las soluciones del sistema
"""
def sust_atras(A, b, n):
    soluciones = np.zeros((n, 1))
    for i in range(n-1, -1, -1):
        suma = 0
        for j in range(i+1, n):
            suma = suma + (A[i, j] * soluciones[j])
        
        x = (1 / A[i, i]) * (b[i] - suma)
        soluciones[i] = x
    
    return soluciones

"""
Método para verificar si una matriz
es positiva definida
Entradas: A es una matriz cuadrada
          n es el tamanio de la matriz
Salidas: true si cumple con la condición
         false si NO cumple con la condición
"""
def verificar_positiva_definida(A, n):
    tol = np.finfo(float).eps
    for i in range(n):
        Ak = A[0:i, 0:i]
        deter = np.linalg.det(Ak)
        if deter < tol:
            return False
    return True

"""
Método para verificar la simetría de 
una matriz
Entradas: A es una matriz cuadrada
          rtol es la tolerancia de filas
          atol es la tolerancia en la matriz
Salidas: true si cumple con la condición
         false si NO cumple con la condición
"""
def verificar_simetria(A, rtol=1e-05, atol=1e-08):
    return np.allclose(A, A.T, rtol=rtol, atol=atol)


"""
Función para resolver sistemas de ecuaciones
mediante el método directo de factorización de Cholesky
Entradas: A es una matriz cuadrada
          b es el vector de terminos idependientes
Salidas: x es el valor de la solución
"""
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
    A = np.matrix('25 15 -5 -10; 15 10 1 -7; -5 1 21 4; -10 -7 4 18')
    b = np.matrix('-25; -19; -21; -5')

    print(fact_cholesky(A, b))