

def ceros_pesos_cuad_gauss(n):
    if 2 <= n <= 10:
        x = []
        w = []
        if n == 2:
            x += [-0.5773502692, 0.5773502692] #x0 y x1
            w += [1, 1] #w0 y w1
            return x, w
        elif n == 3:
            x += [-0.7745966692, 0, 0.7745966692] #x0, x1 y x2
            w += [0.5555555556, 0.888888889, 0.5555555556] #w0. w1 y w2
            return x, w
        elif n == 4:
            x += [-0.86113631165, -0.3399810436, 0.3399810436, 0.8611363116] #x0, x1, x2 y x3
            w += [0.3478548451, 0.6521451549, 0.3478548451, 0.6521451549] #w0, w1, w2 y w3
            return x, w
        elif n == 5:
            x += [-0.9061798459, -0.5384693101, 0, 0.9061798459, 0.5384693101] #x0, x1, x2, x3 y x4
            w += [0.2369268851, 0.4786286705, 0.5688888889, 0.2369268851, 0.4786286705] #w0, w1, w2, w3 y x4
            return x, w
        elif n == 6:
            x += [-0.9324695142, -0.6612093865, -0.2386191861] #x0, x1 y x2
            x += [0.9324695142, 0.6612093865, 0.2386191861] #x3, x4 y x5
            w += [0.1713244924, 0.3607615730, 0.4679133946] #w0, w1 y w2
            w += [0.1713244924, 0.3607615730, 0.4679133946] #w3, w4 y w5
            return x, w
        elif n == 7:
            x += [-0.9491079123, -0.7415311856, -0.4058451514, 0] #x0, x1, x2 y x3
            x += [0.9491079123, 0.7415311856, 0.4058451514] #x4, x5 y x6
            w += [0.1294849662, 0.2797053915, 0.3818300505, 0.4179591837] #w0, w1, w2 y w3
            w += [0.1294849662, 0.2797053915, 0.3818300505] #w4, w5 y w6
            return x, w
        elif n == 8:
            x += [-0.9602898565, -0.7966664774, -0.5255324099, -0.1834346425] #x0, x1, x2 y x3
            x += [0.9602898565, 0.7966664774, 0.5255324099, 0.1834346425] #x4, x5, x6 y x7
            w += [0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834] #w0, w1, w2 y w3
            w += [0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834] #w4, w5, w6 y w7
            return x, w
        elif n == 9:
            x += [-0.9681602395, -0.8360311073, -0.6133714327, -0.3242534234, 0] #x0, x1, x2, x3 y x4
            x += [0.9681602395, 0.8360311073, 0.6133714327, 0.3242534234] #x5, x6, x7 y x8
            w += [0.0812743883, 0.1806481607, 0.2606106964, 0.3123470770, 0.3302393550] #w0, w1, w2, w3 y w4
            w += [0.0812743883, 0.1806481607, 0.2606106964, 0.3123470770] #w5, w6, w7 y w8
            return x, w
        elif n == 10:
            x += [-0.9739065285, -0.8650633667, -0.6794095683, -0.4333953941, -0.1488743390] #x0, x1, x2, x3 y x4
            x += [0.9739065285, 0.8650633667, 0.6794095683, 0.4333953941, 0.1488743390] #x5, x6, x7, x8 y x9
            w += [0.0666713443, 0.1494513492, 0.2190863625, 0.2692667193, 0.2955242247] #w0, w1, w2, w3 y w4
            w += [0.0666713443, 0.1494513492, 0.2190863625, 0.2692667193, 0.2955242247] #w5, w6, w7, w8 y w9
            return x, w
        else:
            return x, w
    else:
        print("Debe ingresar un n entre 2 y 10")
        return Exception