import sympy as sp
import numpy as np
from ceros_pesos_cuad_gauss import *
from scipy.optimize import fminbound

def cuad_gaussiana(func, n, intervalo):
    x = sp.Symbol('x')
    f_simbolica = sp.sympify(func) # Se define de la funcion funcion

    a = intervalo[0]
    b = intervalo[1]

    y = ((b-a)*x + (b+a))/2

    gs = (b-a)/2*f_simbolica.subs(x, y)

    gn = sp.lambdify(x, gs) # La funcion inicial se pasa a a una funcion evaluable

    xw = ceros_pesos_cuad_gauss(n)

    xi = xw[0]
    wi = xw[1]

    aproximacion = 0

    for i in range(n):
        aproximacion += wi[i]*gn(xi[i])

    if n == 2:
        # Para sacar un máximo de la función en
        # un intervalo
        # es equivalente a sacarle el mínimo
        # a la función negada en ese intervalo
        # max(f) = min(-1*f)
        #se deriva 4 veces
        f4 = sp.Abs(gs.diff(x, 4))
        faux = -1*f4
        faux_4 = sp.lambdify(x , faux)
        
        #se obtiene el máximo en ese intervalo
        x_max = fminbound(faux_4, -1, 1)

        fn4 = sp.lambdify(x, f4)

        alpha_max = fn4(x_max)

        #se calcula la cota con la fórmula
        cota_error = alpha_max/135

        return aproximacion, cota_error

    else:
        cota_error = "El orden es distinto de 2"

        return aproximacion, cota_error

if __name__ == "__main__":
    f = 'log(x)'
    intervalo = [2, 5]
    orden = 2

    resultado = cuad_gaussiana(f, orden, intervalo)

    print("Aproximacion: ", resultado[0])
    print("Cota de error: ", resultado[1])