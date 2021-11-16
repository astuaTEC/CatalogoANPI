import numpy as np
import matplotlib.pyplot as plt


#Calculo de valores propios y vectores propios de una matriz
# Entradas:
# A = matriz tamano nxm, iterMax:Cantidad Maxima de iteraciones, tol; toleracion del error
# Salidas:
# ValoresPropios y Vactores Propios
def metodo_qr(A, iterMax, tol):

    Ak = A

    m = len(A)

    Uk = np.eye(m)

    er = []

    Vp = np.diag(Ak)

    for i in range(iterMax):
        #Factorización QR
        Q, R = np.linalg.qr(Ak)

        Ak = np.dot(R, Q)

        Uk = np.dot(Uk, Q)

        Vp_n = np.diag(Ak)

        error = np.linalg.norm(Vp_n - Vp)

        Vp = Vp_n

        er.append(error)

        if error < tol:
            break

    fig, graf = plt.subplots()  # se crea la gráfica
    ejeX = np.arange(0, i+1, 1)  # se crea el eje X (son las iteraciones)
    graf.plot(ejeX, er)  # se grafican los datos
    plt.title('Metodo QR')
    plt.show()  # se muestra la gráfica

    return returnValoresVectores(Ak, Uk)

#Funcion para acomodar los datos de forma estetica en la salida
#Entradas:
#     Una Matriz AK: Vector valor propios, Matriz Triangular
#     Una Matriz UK: Matriz de vectores propios
#Salida:
# Valores y vectores prpopios correspondientes
def returnValoresVectores(Ak, Uk):
    valores = []
    vectores = []
    for i in range(len(Ak)):
        valores.append(Ak[i][i])
        tmp = []
        for j in range(len(Uk)):
            tmp.append(Uk[j][i])
        
        vectores.append(tmp)

    return valores, vectores


if __name__ == "__main__":
    A = [[0, 11, -5], [-2, 17, -7], [-4, 26, -10]]
    
    iterMax = 100

    tol = 10e-6

    res = metodo_qr(A, iterMax, tol)

    for i in range(len(res[0])):
        print("Valor propio: ", res[0][i])
        print("Vector propio: ", res[1][i])
        print("###############################")

