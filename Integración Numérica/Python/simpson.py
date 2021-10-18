import sympy as sp
import numpy as np
from scipy.optimize import fminbound

def simpson(func, intervalo):
    x = sp.Symbol('x')
    f_simbolica = sp.sympify(func) # Se define de la funcion funcion
    f_eval = sp.lambdify(x , f_simbolica) # La funcion inicial se pasa a a una funcion evaluable

    #se obtienen los extremos
    a = intervalo[0]
    b = intervalo[-1]

    aproximacion = ((b-a)/6)*(f_eval(a) + 4*f_eval((a+b)/2) + f_eval(b))

    # Para sacar un máximo de la función en
    # un intervalo
    # es equivalente a sacarle el mínimo
    # a la función negada en ese intervalo
    # max(f) = min(-1*f)
    #se deriva 4 veces
    f4 = abs(f_simbolica.diff(x, 4))
    faux = -1*f4
    faux_n = sp.lambdify(x , faux)
    
    #se obtiene el máximo en ese intervalo
    x_max = fminbound(faux_n, a, b)

    fn4 = sp.lambdify(x, f4)

    alpha_max = fn4(x_max)

    #se calcula la cota con la fórmula
    cota_error = (((b-a)**5) * alpha_max)/2880

    return aproximacion, cota_error

if __name__ == "__main__":
    resultado = simpson('ln(x)', [2, 5])

    print("Aproximacion: ", resultado[0])
    print("Cota de error: ", resultado[1])