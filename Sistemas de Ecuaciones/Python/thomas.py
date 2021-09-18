import numpy as np


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
                a[i] = A[i, j]
            elif i-1 == j:
                c[i] = A[i, j]
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

if __name__ == '__main__':
    #Ejemplo de como se ingresan los datos
    A = np.matrix('5 1 0 0; 1 5 1 0; 0 1 5 1; 0 0 1 5')
    b = np.matrix('-12; -14; -14; -12')
    print("Solucion x de ejemplo:")
    print(thomas(A, b))

#np.eye(3,3,k=-1) + np.eye(3,3)*2 + np.eye(3,3,k=1)*3

