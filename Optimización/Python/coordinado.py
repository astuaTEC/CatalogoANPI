import numpy as np
import sympy as sp
from scipy import optimize


def variables_simbolicas(variables):
    n = len(variables)
    tam = np.arange(0, n, 2)
    v_var = []
    for i in tam:
        v_var.append(sp.Symbol(variables[i]))
    return v_var


def gradiente(f, v_var):
    n = len(v_var)
    g = []
    for i in np.arange(0, n, 1):
        g.append(sp.diff(f, v_var[i]))
    return g


def coordinado(funcion, variables, puntoInicial, iterMax, tol):
    #funciones
    f_s = sp.sympify(funcion)
    f_n = sp.lambdify(variables, funcion)

    #gradirente
    g = gradiente(f_s, variables)
    g_n = sp.lambdify(variables, g)

    #punto inicial
    xk = puntoInicial

    numVariables = len(variables)

    k = 0
    n = 0
    while k < iterMax:
        k += 1
        while n < numVariables:
            simbolica = f_s.subs(variables[n], xk[n+1])

            funcion_num = sp.lambdify(variables[n+1], simbolica)

            punto = optimize.fmin(funcion_num, puntoInicial[n])

            print(punto)

            #xk = [x , y]
            n += 1
        
        error = np.linalg.norm(g_n(xk[0], xk[1]))
        if error < tol:
            break


if __name__ == "__main__":
    #se definen las variables
    variables = 'x y'
    #se crean las variables simbólicas
    v_var = variables_simbolicas(variables)
    #se define la función a trabajar
    f = '(x-2)^2+(y+3)^2 + x*y'
    #se define un punto inicial
    puntoInicial = [1, 1]
    #se definen las iteraciones máximas y la tolerancia
    iterMax = 9
    tol = 10e-6
    #se llama a la función
    coordinado(f, v_var, puntoInicial, iterMax, tol)