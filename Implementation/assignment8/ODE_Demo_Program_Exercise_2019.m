%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    ODE-IVP MATLAB Demo
%%%%
%%%%
%%%%    Numerical Methods 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Demo Problem 
% Example of solving a RC circuit with a sinusoidal input
% Defined in myRC.m
% tau=1; T=1/tau; f=10; Vm=1; w=2*pi*f;
% F=dvdt = @(x,y) -T*y + T*Vm*cos(w*x);
%
% function dydx = myRC(x,y)
% tau=1; T=1/tau; f=10; Vm=1;
% dydx =-T*y + 1*T*Vm*cos(2*pi*f*x);

%% Parameter Definition
clear all
close all
a=0; b=0.1; 
h=0.0001; 
yINI = 0;
x(1) = a;  y(1) = yINI;
N = (b-a)/h;

%% True Solution of RC circuit: Sinusoidal Input  
x=a:h:b;
tau=0.01; T=1/tau; f=100; Vm=1; w=2*pi*f;
A=Vm/(sqrt(1+(w*tau)^2));
alpha=-atan(w*tau);
y=-A*cos(alpha+pi/4)*exp(-x/tau)+A*cos(w*x+alpha+pi/4);


%% Option 1: Euler's explicit method
yE(1) = yINI;
for i = 1:N
    x(i+1) = x(i) + h;
    yE(i+1) = yE(i) + myRC(x(i),yE(i))*h;   % f(xi,yi) => myRC.m
end
figure(1)
plot(x,yE,'--b',x,y,'k')
xlabel('x'); ylabel('y')



%% Option 2: Euler's Modified method
yEu = 0;
Slope2=0; 
yEM(1) = yINI;
for i = 1:N
    x(i+1) = x(i) + h;
    Slope1 = myRC(x(i),yEM(i));
    yEu = yEM(i) + Slope1 * h ;
    Slope2 = myRC(x(i+1), yEu);
    yEM(i+1) = yEM(i) + (Slope1+Slope2)*h/2;
end
figure(1)
plot(x,yEM,'--b',x,y,'k')
xlabel('x'); ylabel('y')



%% Option 3: Mid Point method    
SlopeMid=0;
yMid(1) = yINI;
for i = 1:N
    Xm(i) = x(i) + h/2;
    Ym(i) = yMid(i) + myRC(x(i), yMid(i)) * h / 2;
    SlopeMid = myRC(Xm(i), Ym(i));
    yMid(i+1) = yMid(i) + SlopeMid*h;
    x(i+1) = x(i) + h;    
end

figure(1)
plot(x,yMid,'--b',x,y,'k')
xlabel('x'); ylabel('y')


%% Parameter Definition
clear all
close all
a=0; b=0.1; 
h=0.0001; 
yINI = 0;
x(1) = a;  y(1) = yINI;
N = (b-a)/h;

%% True Solution of RC circuit: Sinusoidal Input  
x=a:h:b;
tau=0.01; T=1/tau; f=100; Vm=1; w=2*pi*f;
A=Vm/(sqrt(1+(w*tau)^2));
alpha=-atan(w*tau);
y=-A*cos(alpha+pi/4)*exp(-x/tau)+A*cos(w*x+alpha+pi/4);



%% Option 4: RK Example

t0=a; tf=b; 
h=0.01; 
[xrk2, yrk2] = odeRK2(@myRC,t0,tf,h,yINI); % 2nd RK : 

figure(1)
plot(xrk2,yrk2,'--b',x,y,'k')
xlabel('x'); ylabel('y')


%% Option 5: MATLAB's function ODE45
[xmat,ymat] = ode45(@myRC, [a b], yINI);  % Fourth/Fifth RK
figure(1)
plot(xmat,ymat,'--b',x,y,'k')
xlabel('x'); ylabel('y')




