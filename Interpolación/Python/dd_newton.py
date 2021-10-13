import sympy as sp

'''
Funcion para calcular el polinomio de interpolacion, mediante el metodo
de diferencias divididas de Newton.

Entradas:       xv = vector con los puntos en x.
                yv = vector con los puntos en y.

Salidas:        polinomioInterpolacion = polinomio de interpolacion
                tras ser calculado.
'''
def dd_newton(xv,yv):
    x = sp.Symbol('x')
    n = len(xv)

    if len(xv) != len(yv):
        return 'Error, el tamano de los vectores xv y yv es diferente.'
    
    # Se calcula el primer termino del polinomio de interpolacion.
    polinomioInterpolacion = yv[0]

    variable = 1

    # Cantidad de operaciones a realizar por cada iteracion.
    contador = n - 1

    for i in range(1,n):

        # Se van almacenando las variables (x-x0)(x-x1)...(x-xn-1).
        variable = variable * (x-xv[i-1])

        # Lista para almacenar los multiplicadores.
        multiplicadores = []

        # Se calculan los multiplicadores.
        for j in range(0,contador):
            numerador = yv[j+1] - yv[j]
            denominador = xv[j+i] - xv[j]
            multiplicadores.append(numerador/denominador)
        
        # Se reduce en uno la cantidad de operaciones a realizar para 
        # la siguiente iteracion.
        contador = contador - 1

        # Se le agrega el nuevo termino al polinomio de interpolacion.
        polinomioInterpolacion = polinomioInterpolacion + multiplicadores[0]*variable

        # Se actualiza la lista en y.
        yv = multiplicadores
    
    return sp.expand(polinomioInterpolacion)


if __name__ == '__main__':
    xv = [-2, 0, 1]
    yv = [0, 1, -1]
    print('El polinomio de interpolacion es ', dd_newton(xv,yv))