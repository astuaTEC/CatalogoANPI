

function potencia_aux
  A = [3, -1 ,0; -1, 2, -1; 0, -1, 3];
  x0 = [1 ;1 ; 1];
  
  iterMax = 12;
  tol = 1*10e-10;
  
  [ck xk] = potencia(A, x0, iterMax, tol)
  
endfunction

% Funcion para aproximar proxima el modulo del valor propio de 
% menor magnitud deuna matriz A y el vector propio normalizadoasociado 
% a dicho valor propio.
% Entradas:
% A: Matriz tamano nxm x0: Valor inicial
% iterMax: Iteraciones Maximas tol:Toleracia maximas del eror
% Salidas:
% xk:Vector Propio  ck:Valor Propio
function [ck xk] = potencia(A, x0, iterMax, tol)
  
  xk = x0;
  er = [];
  iter = [];
  
  for i=1:iterMax
    yk = A*xk;
    ck = norm(yk, "inf");
    xk_n = (1/ck)*yk;
    
    error = norm(xk_n - xk);
    xk = xk_n;
    
    er = [er error];
    iter = [iter i];
    
    if error < tol;
      break
    endif
  endfor
  
  plot(iter, er);
  
endfunction