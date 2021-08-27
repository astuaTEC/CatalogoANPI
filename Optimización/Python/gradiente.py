import sympy as sp
import numpy as np


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


def gradiante_conjugado_no_lineal(funcion, variables, x0, y0, iterMax, tol):
    f_s = sp.sympify(funcion)
    f_n = sp.lambdify(variables, funcion)
    g = gradiente(f_s, variables)
    g_n = sp.lambdify(variables, g)

    k = 0
    gk = np.array(g_n(x0, y0)).T
    dk = -1*gk

    xk = np.array([x0, y0]).T

    iterMax2 = 1000 #iteraciones maximas para el valor de alpha_k

    while k < iterMax:
        delta = 0.5; #se define delta
        alpha_k = 1.0; #se define alpha_k
        k += 1
        n = 0
        while n < iterMax2:
            vectorEval = xk + alpha_k*dk

            izq = f_n(vectorEval[0], vectorEval[1]) - f_n(xk[0], xk[1])
            der = delta*alpha_k* np.sum((gk.T*dk))

            n += 1
            if izq <= der:
                break
            else:
                alpha_k = alpha_k/2

        xk = xk + (alpha_k*dk)
        error = np.linalg.norm(np.array(g_n(xk[0], xk[1])))
        
        if error < tol:
            break

        gk_anterior = gk
        gk = np.array(g_n(xk[0], xk[1])).T
        
        beta_k = (np.linalg.norm(gk))**2 / (np.linalg.norm(gk_anterior))**2
        dk = -gk + beta_k*dk

    print(xk, error)


if __name__ == "__main__":
    #se definen las variables
    variables = 'x y'
    #se crean las variables simbólicas
    v_var = variables_simbolicas(variables)
    #se define la función a trabajar
    f = '(x-2)^4+(x-2*y)^2'
    #se define un punto inicial
    x0 = 0
    y0 = 3
    #se definen las iteraciones máximas y la tolerancia
    iterMax = 13
    tol = 10e-6
    #se llama a la función
    gradiante_conjugado_no_lineal(f, v_var, x0, y0, iterMax, tol)


