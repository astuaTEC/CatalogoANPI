
function trapecio_compuesto_aux
  clc clear
  format long;
  warning off;
  pkg load symbolic
  
  f = 'ln(x)';
  n = 4;
  intervalo = [2, 5];
  
  [aprox cota_error] = trapecio_compuesto(f ,n ,intervalo)

end
%Se calcula la cota de error y el aproximado segun la formula del trapecio 
%compuesto
%Entradas: F => Funcion numerica , intervalo => intervalo a integras limites de 
%la integral n => numero de puntos.
%Salida: aprox => Valor Aproximado de la intergral
%cota_error => Error asociado al metodo

function [aprox cota_error] = trapecio_compuesto(f, n, intervalo)
  syms x;
  f_simbolica = sym(f);
  f_evaluable = matlabFunction(f_simbolica);
  
  h = (intervalo(2) - intervalo(1))/(n-1);
  a = intervalo(1);
  b = intervalo(2); 
  sumatoria = 0;
  xi = [];
  xi = [xi,intervalo(1)];
  
  %Calculo de la sumatoria de la integral
  for i=2:n
    xi = [xi, xi(i-1) + 1];
  endfor
  for i=1:(n-1);
    sumatoria += f_evaluable(xi(i))+f_evaluable(xi(i+1));
  endfor
  
  aprox = (h/2)*sumatoria;
  
  %Calculo dela cota de error
  f_derivada_2 = diff(f_simbolica,2);
  f_derivada_2_evaluable = matlabFunction(f_derivada_2);
  f_aux = -1 * abs(f_derivada_2); 
  f_aux_num = matlabFunction(f_aux);
  x_max = fminbnd(f_aux_num,a,b);
  alpha_max = abs(f_derivada_2_evaluable(x_max));
  
  cota_error = ((b-a)*(h^2)*alpha_max)/(12);
end
