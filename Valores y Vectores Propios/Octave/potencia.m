
function potencia_aux
  A = [3, -1 ,0; -1, 2, -1; 0, -1, 3];
  x0 = [1 ;1 ; 1];
  
  iterMax = 12;
  tol = 1*10e-10;
  
  [ck xk] = potencia(A, x0, iterMax, tol)
  
endfunction


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