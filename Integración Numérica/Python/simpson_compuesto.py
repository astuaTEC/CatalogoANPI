import sympy as sp
import numpy as np
from scipy.optimize import fminbound

def simpson_compuesto(func, n, intervalo):
    x = sp.Symbol('x')
    f_simbolica = sp.sympify(func) # Se define de la funcion funcion
    f_eval = sp.lambdify(x , f_simbolica) # La funcion inicial se pasa a a una funcion evaluable

    a = intervalo[0]
    b = intervalo[1]
    h = (b - a)/(n-1)
    
    sumatoria_pares = 0
    sumatoria_impares = 0
    xi = np.linspace(a, b, n)

    for i in range(2, n-1, 2):
        sumatoria_pares += f_eval(xi[i])
    
    for i in range(1, n-1, 2):
        sumatoria_impares += f_eval(xi[i])

    aproximacion = (h/3)*(f_eval(xi[0]) + 2*sumatoria_pares + 4*sumatoria_impares + f_eval(xi[n-1]))

    # Para sacar un máximo de la función en
    # un intervalo
    # es equivalente a sacarle el mínimo
    # a la función negada en ese intervalo
    # max(f) = min(-1*f)
    #se deriva 4 veces
    f4 = sp.Abs(f_simbolica.diff(x, 4))
    faux = -1*f4
    faux_4 = sp.lambdify(x , faux)
    
    #se obtiene el máximo en ese intervalo
    x_max = fminbound(faux_4, a, b)

    fn4 = sp.lambdify(x, f4)

    alpha_max = fn4(x_max)

    #se calcula la cota con la fórmula
    cota_error = (((b-a)*h**4) * alpha_max)/180

    return aproximacion, cota_error


if __name__ == "__main__":
    resultado = simpson_compuesto('ln(x)', 7, [2, 5])

    print("Aproximacion: ", resultado[0])
    print("Cota de error: ", resultado[1])