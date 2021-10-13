import numpy as np
from sympy import *

"""
Función para calcular el trazador cubico que pasa por
los pares ordenados (xv, yv)
Entradas:
    xv: ector de puntos de tamano n
    yv: vector de imagenes de tamano n
Salida:
    S: vector de polinomios del trazador cubico
"""
def traz_cubico(xv, yv):
    n = len(xv) #Largo del vector

    if (len(xv) != len(yv)): #Comprueba que los vectores sean iguales
        return ('Los venctores no cumplen con la condicion de tamano')
    h = np.zeros(n-1)

    for i in range(0, n-1): #Calculo del vector h (distancia entre cada punto
                            # del trazador)
        h[i] = xv[i + 1] - xv[i]

    A = np.zeros((n-2,n-2)) #Matriz tridiagonal

    A[0][0] = 2*(h[0]+h[1]) #Calculo de la primer fila de la matriz n-1xn-1
    A[0][1] = h[1]

    for i in range(1, n-3): #Calculo del la de 1 a n-2 de la matriz
        A[i][i] = h[i]
        A[i][i+1] = 2*(h[i] + h[i+1])
        A[i][i+2] = h[i+1]
    
    A[n-3][n-4] = h[n-3]  #Calculo de la posicion n-1 de la matriz
    A[n-3][n-3] = 2*(h[n-3] + h[n-2])

    u = np.zeros(n-2) # Vector de variables u

    for i in range(0, n-2): #Calculo del vector u
        u[i] = 6*(((yv[i+2]-yv[i+1])/h[i+1])-((yv[i+1]-yv[i])/h[i]))

    M_temp = thomas(A,u) #Implementacion del metodo de Thomas para resolver
                         #El sistema Ax=U
    M = np.zeros(n)
    M[0] = 0 #Matriz M con M0 = 0 y Mn-1 = 0
    M[n-1] = 0

    for i in range(1, n-1): #Construccion de la matriz M apartir de la matriz 
                            #M_temp
        M[i] = M_temp[i-1]
    
    S = [] #Vector de polinimios de solucion S
    
    for i in range(0, n-1): #Calculo de las variables a,b,c,d
        a = (M[i+1]-M[i])/(6*h[i])
        b = M[i]/2
        c = (yv[i+1]-yv[i])/h[i] - (h[i]/6)*(M[i+1]+2*M[i])
        d = yv[i]
        S.append(crear_polinomio(a,b,c,d,xv[i])) #Creacion de polinomio Si 

    return S


"""
Método para crear un polinomio simbólico
Entradas:
    a, b, c, d: valores del polinomio
    xv: valor numerico de un punto del vector xv
Salida:
    p: polinomio simbólico
"""
def crear_polinomio(a, b ,c, d, xv):
    x = Symbol('x')
    polinomio = 0
    polinomio = sympify(polinomio)
    polinomio = a*(x - xv)**3 + b*(x - xv)**2 + c*(x - xv) + d 
    return polinomio

"""
Método para verificar si una matriz
es tridiagonal
Entradas:
    A: Matriz cuadrada mxm
    n: Tamanio de la matriz
Salidas:
    True: si cumple con la condición
    False: si no cumnple la condición
"""
def verificarTridiagonal(A, n):
    for i in range(n):
        for j in range(n):
            elem = A[i, j]
            if (i == j) or (i-1 == j) or (i+1 == j):
                if elem == 0:
                   return False
            else:
                if elem != 0:
                    return False
    return True

"""
Método para obtener los vectores
a, b y c de la matriz tridiagonal
Entradas:
    A: la matriz tridiagonal
    n: el tamanio de la matriz
Salidas:
    vector con los vectores a, b, c
"""
def obtener_abc(A, n):
    a = np.zeros((n, 1))#Matriz lleno de zeros nx1
    b = np.zeros((n, 1))
    c = np.zeros((n, 1))

    for i in range(n):#Mueve filas
        for j in range(n):#Mueve columnas
            if j == i: #Diagonal
                b[i] = A[i, j]#Guarda valor en la matriz
            elif i+1 == j:
                c[i] = A[i, j]
            elif i-1 == j:
                a[i] = A[i, j]
            else:
                if A[i, j] != 0:
                    raise ValueError("Esta matriz no es tridiagonal")
    return a, b, c

"""
Función para resolver sistemas de ecuaciones
mediante el método directo de Thomas
Entradas:
    A: una matriz cuadrada mxm
    b: vector de términos independientes
Salida:
    Solución del sistema
"""
def thomas(A, b):
    n = len(A)
    if not verificarTridiagonal(A, n):
        return "La matriz debe ser tridiagonal"
    
    ABC = obtener_abc(A, n)

    return thomas_aux(ABC[0], ABC[1], ABC[2], b, n)

"""
Método auxiliar de thomas para trabajar
con los vectores a, b, c y d
Entradas:
    a, b, c, d: vectores a trabajar
    n: tamanio de los vectores
Salida:
    xv: la solución del sistema
"""
def thomas_aux(a, b, c, d, n):
    q  =np.zeros((n, 1))#Matriz lleno de zeros nx1
    p = np.zeros((n, 1))
    xv = np.zeros((n, 1))

    p[0] = c[0]/b[0] #Primer coeficiente

    for i in range(1, n-1):
        if(b[i]- p[i-1]*a[i]) == 0:#Comprueba que el divisor no sea 0
            raise ValueError("No se puede dividir entre 0")
        else:
            p[i]= c[i]/(b[i]- p[i-1]*a[i])#Calcula los nuevos coeficientes
    
    q[0] = d[0]/b[0]
    for i in range(1, n):
        if(b[i] - p[i-1]*a[i]) == 0:#Comprueba que el divisor no sea 0
            raise ValueError("No se puede dividir entre 0")
        else:
            q[i] = (d[i] - q[i-1]*a[i])/(b[i] - p[i-1]*a[i])

    xv[n-1] = q[n-1]

    for i in range(n-2, -1, -1):
        xv[i] = q[i] - p[i]*xv[i+1]

    return xv

if __name__ == "__main__":
    xv = [1, 1.05, 1.07, 1.1]
    yv = [2.718282, 3.286299, 3.527609, 3.905416]

    print(traz_cubico(xv, yv))