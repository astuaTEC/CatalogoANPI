import numpy as np


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

def sust_atras(A, b, n):
    soluciones = np.zeros((n, 1))
    for i in range(n-1, -1, -1):
        suma = 0
        for j in range(i+1, n):
            suma = suma + (A[i, j] * soluciones[j])
        
        x = (1 / A[i, i]) * (b[i] - suma)
        soluciones[i] = x
    
    return soluciones

def gaussiana(A, b):
    n = len(A)
    res_triangular = obtener_triangular_superior(A, b, n)
    A_t = res_triangular[0]
    b_t = res_triangular[1]

    res = sust_atras(A_t, b_t, n)
    return res

def obtener_triangular_superior(A, b, n):
    A_aum = np.append(A, b, axis=1)
    for k in range(n): #recorre las columnas
        for i in range(k+1, n): #recorre las filas
            m_ik = A_aum[i, k] / A_aum[k, k]
            for j in range(k, n+1):
                A_aum[i, j] = A_aum[i, j] - (m_ik * A_aum[k, j])
    
    A_res = A_aum[:, 0:n]
    b_res = A_aum[:, n]
    return [A_res, b_res]


if __name__ == '__main__':
    A = np.matrix('1 1 -1 3;0 -1 -1 -5; 0 0 3 13; 0 0 0 -13', float)
    b = np.matrix('4; -7; 13; -13')
    n = 4

    A2 = np.matrix('5 1 1; 1 5 1; 1 1 5')
    b2 = np.matrix('7; 7; 7')

    #print(obtener_triangular_superior(A2, b2, 3))
    #print(sust_atras(A, b, n))
    print(gaussiana(A2, b2))

