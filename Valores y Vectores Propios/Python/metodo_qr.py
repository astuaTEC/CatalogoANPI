import numpy as np

def metodo_qr(A, iterMax):

    Ak = A

    m = len(A)

    Uk = np.eye(m);

    for i in range(iterMax):
        Q, R = np.linalg.qr(Ak)

        Ak = np.dot(R, Q)

        Uk = np.dot(Uk, Q)

    return returnValoresVectores(Ak, Uk)


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

    res = metodo_qr(A, iterMax)

    for i in range(len(res[0])):
        print("Valor propio: ", res[0][i])
        print("Vector propio: ", res[1][i])
        print("###############################")

