%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ODE-IVP MATLAB Demo
%%%%
%%%%
%%%% Numerical Methods 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demo Problem
% Example of solving an m-c-k system with a sinusoidal input
% Defined in mckFun.m
% m=10, k=4000,c=20;
% FinDC=100; Fin=FinDC*cos(2*pi*10*t);
%% Parameter Definition
clc
clear all
close all

xINI = 0;
vINI = 0.2;
a = 0;
b = 5;
h = 0.001;
N = (b-a)/h;

%% Actual Solution- MATLAB
tspan = [a:h:b];
xini = [xINI vINI];     % vector�� �־�� ��
[Time X] = ode45(@mckFunc1,tspan,xini);     % �����ؾ� ��

%% Option 1: RK 2nd order
[Time, Xrk, Vrk] = sys2RK2(@mckFunc,a,b,h,xINI,vINI);

%% Option 2: Euler Method
t(1) = a; xE(1) = xINI; vE(1) = vINI;
N = (b - a)/h;
for i = 1:N
    t(i+1) = t(i) + h;
    dXdt=mckFunc(t(i),xE(i),vE(i));     % x ���� & v ���� �ΰ��� ����
    xE(i+1) = xE(i) + dXdt(1)*h;
    vE(i+1) = vE(i) + dXdt(2)*h;
end

%% Plot
figure
subplot(2,1,1)
 plot(Time,X(:,1),'--b', Time,Xrk,'k')
 xlabel('Time (s)')
 ylabel('Position (m)')
 title('ode45 vs RK2')
subplot(2,1,2)
 plot(Time,X(:,2),'--b', Time,Vrk,'k')
 xlabel('Time (s)')
 ylabel('Velocity (m/s)')
figure
subplot(2,1,1)
 plot(Time,X(:,1),'--b', Time,xE,'k')
 xlabel('Time (s)')
 ylabel('Position (m)')
 title('ode45 vs Euler')
subplot(2,1,2)
 plot(Time,X(:,2),'--b', Time,vE,'k')
 xlabel('Time (s)')
 ylabel('Velocity (m/s)')