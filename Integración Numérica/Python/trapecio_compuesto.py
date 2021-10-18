import sympy as sp
import numpy as np

def trapecio_compuesto(func, n, intervalo):
    x = sp.Symbol('x')
    f_simbolica = sp.sympify(func) # Se define de la funcion funcion
    f_eval = sp.lambdify(x , f_simbolica) # La funcion inicial se pasa a a una funcion evaluable

    h = (intervalo[1] - intervalo[0])/(n-1)

    sumatoria = 0
    xi = []
    xi.append(intervalo[0])

    #calcular xi
    for i in range(1, n):
        xi += [xi[i-1] + 1]

    print(xi)
    for i in range(n-1):
        sumatoria += f_eval(xi[i]) + f_eval(xi[i + 1])

    aproximacion = (h/2)*sumatoria

    print(aproximacion)

trapecio_compuesto('ln(x)', 4, [2, 5])