%Predictor corrector

clc; clear; close all
a=0; b=5;
num_pt=11;
f=@(x,y) y-x.^2+1;
h=(b-a)/(num_pt-1);
xv=a:h:b;
yv=[0.5];
for n=1:num_pt-1 
  zv = yv(n)+h*f(xv(n),yv(n));
  yv(n+1)= yv(n) + h/2*(f(xv(n),yv(n)) + f(xv(n+1),zv));
end


%Grè´°fica de la aproximaciè´—n
hold on
plot(xv,yv,'r')

%Soluciè´—n analè´øtica
y_s=@(x) (x+1).^2-0.5*exp(x);
%Graficar la solucion
x_g=a:0.0001:b;
y_g=y_s(x_g);
plot(x_g,y_g,'b')