import numpy as np
import sympy as sp
from scipy import optimize
import json

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

    #gradirente
    g = gradiente(f_s, variables)
    g_n = sp.lambdify(variables, g)

    #punto inicial
    xk = puntoInicial

    numVariables = len(variables)

    k = 0
    
    while k < iterMax:
        k += 1
        n = 0
        indexSuperior = numVariables - 1
        while n < numVariables:
            
            simbolica = f_s.subs(variables[indexSuperior], xk[indexSuperior])

            funcion_num = sp.lambdify(variables[n], simbolica)

            punto = optimize.fmin(funcion_num, puntoInicial[n])

            xk[n] = punto[0]

            n += 1

            indexSuperior -= 1
        
        data ={}
        g_evaluado = []
        for i in range(0, numVariables):
            data[variables[i]] = xk[i]

        for i in range(0, numVariables):
            g_evaluado.append(float(g[i].subs(data).evalf()))

        error = np.linalg.norm(g_evaluado)
        if error < tol:
            break

    print(xk, error)
    return [xk, error]


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
