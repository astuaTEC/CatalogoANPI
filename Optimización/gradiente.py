from numpy.core.fromnumeric import transpose
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
    gk = g_n(x0, y0)
    dk = -1*gk

    xk = []
    xk.append(x0)
    xk.append(y0)

    iterMax2 = 1000 #iteraciones maximas para el valor de alpha_k

    while k < iterMax:
        delta = 0.5; #se define delta
        alpha_k = 1; #se define alpha_k
        k += 1
        n = 0
        while n < iterMax2:
            vectorEval = xk + alpha_k*dk

            izq = f_n(vectorEval[0], vectorEval[1]) - f_n(xk[0], xk[1])
            der = delta*alpha_k*dk*gk

            n += 1
            if izq <= der:
                break
            else:
                alpha_k = alpha_k/2

        xk = xk + (alpha_k*dk)
        error = np.linalg.norm(g_n(xk[0], xk[1]))
        
        if error < tol:
            break

        gk_anterior = gk
        gk = g_n(xk[0], xk[1])
        
        beta_k = (np.linalg.norm(gk)^2)/(np.linalg.norm(gk_anterior)^2)
        dk = -gk + beta_k*dk

    print(xk)




"""Ejemplo de Variables"""
variables = 'x y'
v_var = variables_simbolicas(variables)
print(v_var)

"""Ejemplo del Gradiente Simbolico"""
f = '(x-2)^4+(x-2*y)^2'
f_s = sp.sympify(f)
g = gradiente(f_s, v_var)
print(g)

gradiante_conjugado_no_lineal(f, v_var, 0, 3, 13, 10e-3)

"""Ejemplo de Sustitución Simbolica y Numérica"""
"""
x0 = 0
y0 = 3
f_n = sp.lambdify(v_var, f)
print(f_n(x0, y0))
print(f_s.subs([(v_var[0], x0), (v_var[1], y0)]))

g_n = sp.lambdify(v_var, g)
g_trans = np.transpose(g)
g_n_trans = sp.lambdify(v_var, g_trans)
"""

