from gaussiana import *
import numpy as np

def verificar_fact_lu(A, n):
    tol = np.finfo(float).eps
    for i in range(n):
        Ak = A[0:i, 0:i]
        deter = np.linalg.det(Ak)
        if abs(deter) < tol:
            return False
    return True


def obtener_lu(A, n):
    U = A
    L = np.identity(n)

    for k in range(n-1): #recorre las columnas
        for i in range(k+1, n): #recorre las filas
            m_ik = U[i, k] / U[k, k]
            L[i,k] = m_ik #se agrega el multiplicador
            
            for j in range(k, n): #actualiza los valores de U
                U[i, j] = U[i, j] - (m_ik * U[k, j])
    
    return L, U

def fact_lu(A, b):
    n = len(A)
    if(not(verificar_fact_lu(A, n))):
        return "NA"
    
    Matriz = obtener_lu(A, n)

    y = sust_adelante(Matriz[0], b, n) #resuelve Ly = b
    x = sust_atras(Matriz[1], y, n) #Resuelve Ux = y

    return x

        
if __name__ == '__main__':
    A = np.matrix('4 -2 1; 20 -7 12; -8 13 17')
    b = np.matrix('11; 70; 17')

    print(fact_lu(A, b))