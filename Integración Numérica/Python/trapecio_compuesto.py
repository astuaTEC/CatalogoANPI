import sympy as sp
import numpy as np
from scipy.optimize import fminbound

def trapecio_compuesto(func, n, intervalo):
    x = sp.Symbol('x')
    f_simbolica = sp.sympify(func) # Se define de la funcion funcion
    f_eval = sp.lambdify(x , f_simbolica) # La funcion inicial se pasa a a una funcion evaluable

    h = (intervalo[1] - intervalo[0])/(n-1)
    a = intervalo[0]
    b = intervalo[1]
    
    sumatoria = 0
    xi = np.linspace(a, b, n)

    print(xi)
    for i in range(n-1):
        sumatoria += f_eval(xi[i]) + f_eval(xi[i + 1])

    aproximacion = (h/2)*sumatoria

    # Para sacar un máximo de la función en
    # un intervalo
    # es equivalente a sacarle el mínimo
    # a la función negada en ese intervalo
    # max(f) = min(-1*f)
    #se deriva 4 veces
    f2 = sp.Abs(f_simbolica.diff(x, 2))
    faux = -1*f2
    faux_2 = sp.lambdify(x , faux)
    
    #se obtiene el máximo en ese intervalo
    x_max = fminbound(faux_2, a, b)

    fn2 = sp.lambdify(x, f2)

    alpha_max = fn2(x_max)

    #se calcula la cota con la fórmula
    cota_error = (((b-a)*h**2) * alpha_max)/12

    return aproximacion, cota_error


if __name__ == "__main__":
    resultado = trapecio_compuesto('ln(x)', 4, [2, 5])

    print("Aproximacion: ", resultado[0])
    print("Cota de error: ", resultado[1])