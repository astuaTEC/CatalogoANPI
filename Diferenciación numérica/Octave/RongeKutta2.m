

%Ronge-Kutta de orden 2

clc; clear; close all
a=0; b=1;
num_pt = 11;

f=@(x,y) -x*y + 4*x/y
h=(b-a)/(num_pt-1);
xv = a:h:b;
yv = [1];
for i=1:num_pt-1 
  k1 = f(xv(i),yv(i));
  k2 = f(xv(i) + h/2, yv(i) + h*k1/2);
  
  yv(i+1) = yv(i) + h*k2;
end

yv 

%Grè´°fica de la aproximaciè´—n
hold on
plot(xv,yv,'r')
