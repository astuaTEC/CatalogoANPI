f = 'log(x) - exp(-x) - cos(x)'  # se define la función a utilizar.
a = 1  # se define el límite inferior del rango [a, b]
b = 2 # se define el límite superior del rang0 [a, b]
K = 0 # Valor inicial de la iteraciones
tolerancia = 10 ** -5  # se define la tolerancia máxima
iterMax = 500  # se definen las iteracioens máximas antes de detenerse


def falsaPosicion(f, a , b,  tolerancia , iterMax):
    x = sp.Symbol('x') # Asi se define la varible simbolica
    f_simbolica = sp.sympify(f) # Se define de la funcion funcion
    f_eval = sp.lambdify(x , f_simbolica) # La funcion inicial se pasa a a una funcion evaluable
    er = []  # se define el vector de los errores para graficar posteriormente
    xk = a  # Se define el primer valor o el valor inicial de
    xk_1 = b
    if( f_eval(a) * f_eval(b) < 0): # Se verficia el teorema de bolzano
        n = f(xk_1) * (xk_1 - xk)
         #numerador formula de la secante
        d = f(xk_1) - f(xk) #denominador de la formula de la secante
        xk = xk_1
        xk_1 = xk_1 - n/d
        while(k < iterMax):
            k = k + 1



    xk = a  # Se define el primer valor o el valor inicial de x0
    xk1 = b
    error = tolerancia + 1  # se define un valor inicial para el error
    k = 0  # se inicializa las iteraciones en cero
