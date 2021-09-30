import numpy as np
import matplotlib.pyplot as plt


'''
Funcion para determinar si una matriz es diagonalmente dominante
Entradas:   matriz = matriz de coeficientes
Salida:     booleano (True=es diagonalmente dominante) (False=no es diagonalmente dominante)
'''
def verificarDiagonalmenteDominante(matriz):
    for fila in range(0, len(matriz)):
        valorEnLaDiagonal = matriz[fila,fila]
        filaActual = matriz[fila,:]
        filaActual[0, fila] = 0
        if abs(valorEnLaDiagonal) <= np.sum(filaActual):
            return False
    return True

'''
Funcion para calcular la matriz U a partir de una matriz de coeficientes
Entradas:   matriz = matriz de coeficientes
            tamanoMatriz = tamaño de la matriz de coeficientes
Salida:     matrizU = matriz U de la matriz de coeficientes
'''
def calcularMatrizU(matriz, tamanoMatriz):
    matrizU = np.zeros_like(matriz)
    for i in range(0, tamanoMatriz):
        for j in range(i + 1, tamanoMatriz):
            matrizU[i, j] = matriz[i, j]
    return matrizU

'''
Funcion para calcular la matriz D a partir de una matriz de coeficientes
Entradas:   matriz = matriz de coeficientes
            tamanoMatriz = tamano de la matriz de coeficientes
Salida:     matrizD = matriz D de la matriz de coeficientes
'''
def calcularMatrizD(matriz, tamanoMatriz):
    matrizD = np.zeros_like(matriz)
    for i in range(0, tamanoMatriz):
        matrizD[i, i] = matriz[i, i]
    return matrizD

'''
Funcion para calcular la matriz L a partir de una matriz de coeficientes
Entradas:   matriz = matriz de coeficientes
            tamanoMatriz = tamano de la matriz de coeficientes
Salida:     matrizL = matriz L de la matriz de coeficientes
'''
def calcularMatrizL(matriz, tamanoMatriz):
    matrizL = np.zeros_like(matriz)
    for i in range(0, tamanoMatriz):
        for j in range(0, i):
            matrizL[i, j] = matriz[i, j]
    return matrizL

'''
Funcion para graficar el metodo de Jacobi
Entradas:   errores = lista de errores
            iteracionActual = iteracion final del metodo
'''
def graficarMetodoJacobi(errores, iteracionActual):
    fig, grafica = plt.subplots()
    ejeX = np.arange(1, iteracionActual+1, 1)
    grafica.plot(ejeX, errores)
    plt.title('Metodo de Jabobi')
    grafica.set_xlabel('Cantidad de iteraciones')
    grafica.set_ylabel('Error')
    plt.show()

'''
Funcion para calcular el determinante de una matriz mediante el metodo de Jacobi
Entradas:   matrizCoeficientes = matriz de coeficientes
            vectorTerminosIndependientes = vector de terminos independientes
            valorInicial = vector con los valores iniciales
            tolerancia = tolerancia esperada por el metodo
            iteracionesMaximas = cantidad maxima de iteraciones
Salida:     matrizL = matriz L de la matriz de coeficientes
'''
def jacobi(matrizCoeficientes, vectorTerminosIndependientes, valorInicial, tolerancia, iteracionesMaximas):
    if verificarDiagonalmenteDominante(matrizCoeficientes.copy()):
        tamanoMatriz = len(matrizCoeficientes)
        matrizL = calcularMatrizL(matrizCoeficientes,tamanoMatriz)
        matrizU = calcularMatrizU(matrizCoeficientes,tamanoMatriz)
        matrizD = calcularMatrizD(matrizCoeficientes,tamanoMatriz)
        inversaMatrizD = matrizD.I
        xk = valorInicial
        errores = []
        iteracionActual = 0
        while iteracionActual <= iteracionesMaximas:
            iteracionActual = iteracionActual + 1
            xk = -inversaMatrizD * (matrizL + matrizU) * xk + inversaMatrizD*vectorTerminosIndependientes
            error = np.linalg.norm(matrizCoeficientes*xk-vectorTerminosIndependientes)
            errores.append(error)
            if error < tolerancia:
                break
        graficarMetodoJacobi(errores, iteracionActual)
        return "xk (aproximacion): \n" + str(xk), "Error aproximado: " + str(error)
    else:
        return "Error: la matriz \n" + str(matrizCoeficientes) + "\n no es diagonalmente dominante",""

if __name__ == '__main__':
    matrizCoeficientes = np.matrix('9 1 5; 2 5 2; 6 2 30')
    vectorTerminosIndependientes = np.matrix('1; 1; 1')
    valorInicial = np.matrix('0; 0; 0')
    tolerancia = 1e-5
    iteracionesMaximas = 1000
    resultado = jacobi(matrizCoeficientes,vectorTerminosIndependientes,valorInicial,tolerancia,iteracionesMaximas)
    print(resultado[0])
    print(resultado[1])