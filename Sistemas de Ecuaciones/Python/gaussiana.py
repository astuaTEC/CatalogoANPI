import numpy as np

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
Método para obtener la matriz triangular superior
a partir de la matriz de coeficientes y el
vector de términos independientes
Entradas: A es la matriz de coeficientes
          b es el vector de términos independientes
          n es el tamanio de la matriz
Salidas: Una lista que contiene la nueva
         matriz de coeficientes y el nuevo vector
         de términos independientes
"""
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

"""
Función para resolver sistemas de ecuaciones
mediante el método directo de eliminación Gaussiana
Entradas: A es una matriz cuadrada
          b es el vector de terminos idependientes
Salidas: res es el valor de la solución
"""
def gaussiana(A, b):
    n = len(A)
    res_triangular = obtener_triangular_superior(A, b, n)
    A_t = res_triangular[0]
    b_t = res_triangular[1]

    res = sust_atras(A_t, b_t, n)
    return res

if __name__ == '__main__':
    A = np.matrix('25 15 -5 -10; 15 10 1 -7; -5 1 21 4; -10 -7 4 18')
    b = np.matrix('-25; -19; -21; -5')

    print(gaussiana(A, b))

