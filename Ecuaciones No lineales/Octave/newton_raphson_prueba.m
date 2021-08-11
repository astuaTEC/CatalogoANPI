#################################
### METODO DE NEWTON RAPHSON ####
###                          ####
#################################


function newton_raphson_prueba
  clc; clear;
  format long;
  
  f = 'exp(x) + x-2';
  x0 = 5;
  tol=10^-9;s
  iterMax=100;
  [xk error] = newton_raphson(f, x0, tol, iterMax)
  
end


function [xk error] = newton_raphson(f, x0, tol, iterMax)
  %cargar el paquete symbolic
  pkg load symbolic
  
  f1 = matlabFunction(sym(f));
  df1 = matlabFunction(diff(sym(f)));
  er = [];
  iter = [];
  error = tol + 1;
  k = 0;
  xk = x0;
  while (error > tol && k < iterMax)
    iter = [iter; k];
    k = k + 1;
    n = f1(xk);
    d = df1(xk);
    
    if d > tol
      xk = xk - n/d;
      error = abs(f1(xk));
      er = [er; error];
    else
      break
    endif
  endwhile
  
  plot(iter, er)
  grid on
  
end
