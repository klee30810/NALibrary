function [dXdt] = mckFunc(t, x1,x2)   % x1: displacement, x2: velocity
dXdt = zeros(2,1);

%example 2
 m=10;, k=4000;,c=20;
 FinDC=100; 
 %example 2
 Fin=FinDC*cos(2*pi*10*t);
 % y = x1, z = x2=ydot
 dXdt(1) = x2;
 dXdt(2) = 1/m * (Fin-c*x2-k*x1);
 