
function cuad_gausiana_aux
  pkg load symbolic;
  warning off;
  clc; clear;
  format long;
  
  f = 'log(x)';
  intervalo = [2 5];
  orden = 2;
  
  [aproximacion, cota_error] = cuad_gausiana(f, orden, intervalo)
  
end

function [aproximacion cota_error] = cuad_gausiana(func, n, intervalo)
    syms x;
    f_simbolica = sym(func);
    
    a = intervalo(1);
    b = intervalo(2);
    
    y = ((b-a)*x + (b+a))/2;
    
    gs = (b-a)/2*subs(f_simbolica, x, y);
    
    gn = matlabFunction(gs);
    
    [xi, wi] = ceros_pesos_cuad_gauss(n);
    
    aproximacion = 0;
    
    for i = 1:n
      aproximacion += wi(i)*gn(xi(i));
    endfor

    if n == 2
      f4 = abs(diff(gs, 4));
      faux = -1*f4;
      faux_4 = matlabFunction(faux);
        
      #se obtiene el máximo en ese intervalo
      x_max = fminbnd(faux_4, -1, 1);

      fn4 = matlabFunction(f4);

      alpha_max = fn4(x_max);

      #se calcula la cota con la fórmula
      cota_error = alpha_max/135;
      
    else
      cota_error = "El orden es distinto de 2";
    endif
end