import numpy as np
from sympy import *

#Parametros de entrada
#XK vector de puntos de tamano n, #YK vector de imagenes de tamano n
#Salida
#S vector de polinomios del trazador cubico
#Trazador cubico
def traz_cubico(Xk, Yk):
    n = len(Xk) #Largo del vector

    if (len(Xk) != len(Yk)): #Comprueba que los vectores sean iguales
        return ('Los venctores no cumplen con la condicion de tamano')
    h = np.zeros(n-1)

    for i in range(0, n-1): #Calculo del vector h (distancia entre cada punto
                            # del trazador)
        h[i] = Xk[i + 1] - Xk[i]

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
        u[i] = 6*(((Yk[i+2]-Yk[i+1])/h[i+1])-((Yk[i+1]-Yk[i])/h[i]))

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
        c = (Yk[i+1]-Yk[i])/h[i] - (h[i]/6)*(M[i+1]+2*M[i])
        d = Yk[i]
        S.append(crear_funcion(a,b,c,d,Xk[i])) #Creacion de polinomio Si 

    print(S)
    return([S, h])

#Parametros de entrada
#a, b, c, d valores del polinomio
#Xk constante X0
#Salida
#polinomio Si
#Polinomio Si
def crear_funcion(a,b,c,d,xk):
    x = Symbol('x')
    polinomio = 0
    polinomio = sympify(polinomio)
    polinomio = a*(x - xk)**3 + b*(x - xk)**2 + c*(x - xk) + d 
    return (polinomio)

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

def thomas(A, b):
    n = len(A)
    if not verificarTridiagonal(A, n):
        return "La matriz debe ser tridiagonal"
    
    ABC = obtener_abc(A, n)

    return thomas_aux(ABC[0], ABC[1], ABC[2], b, n)

def thomas_aux(a, b, c, d, n):
    q  =np.zeros((n, 1))#Matriz lleno de zeros nx1
    p = np.zeros((n, 1))
    xk = np.zeros((n, 1))

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

    xk[n-1] = q[n-1]

    for i in range(n-2, -1, -1):
        xk[i] = q[i] - p[i]*xk[i+1]

    return xk