
function trapecio_aux
  clc clear
  format long;
  warning off;
  pkg load symbolic 
  f = 'ln(x)';
  intervalo = [2, 5]; 
  [aprox cota_error] = trapecio(f ,intervalo)
end

function [aprox cota_error] = trapecio(f ,intervalo)
  syms x;
  f_simbolica = sym(f);
  f_evaluable = matlabFunction(f_simbolica);
  
  a = intervalo(1);
  b = intervalo(2); 
  
  aprox = ((b-a)*(f_evaluable(a)+f_evaluable(b)))/2;
  
  f_derivada_2 = diff(f_simbolica,2);
  f_derivada_2_evaluable = matlabFunction(f_derivada_2);
  f_aux = -1 * abs(f_derivada_2); 
  f_aux_num = matlabFunction(f_aux);
  x_max = fminbnd(f_aux_num,a,b);
  alpha_max = abs(f_derivada_2_evaluable(x_max));
  
  cota_error = (((b-a)^3)*alpha_max)/12;
end