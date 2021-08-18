############################
### Metodo de punto fijo####
############################

function punto_fijo
  clc; clear;
  format long;
  f = 'log(2*x + 1)';
  x0 = 1.5;
  tol=10^-9;
  iterMax=150;
  [xk error] = punto_fijo_aux(f, x0, tol, iterMax)
end

function [xk error] = punto_fijo_aux(f, x0, tol, iterMax)
  er = [];
  error = tol;
  iter = [];
  xk = x0;
  f1 = matlabFunction(sym(f));
  for(k=0:iterMax)
    iter = [iter; k];
    k +=1;
    xk = f1(xk);
    error = abs(f1(xk));
    er = [er; error];
    if(xk < tol)
      break;
    endif
  endfor
  plot(iter, er)
  grid on
end


