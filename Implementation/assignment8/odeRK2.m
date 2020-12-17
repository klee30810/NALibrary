%% 2nd ODE solver
function [xrk2, yrk2] = odeRK2(func,t0,tf,h,yINI) % 2nd RK

N = (tf-t0)/h;
alpha = 1/2;
beta = alpha;
c1 = 1-1/(2*alpha);    c2 = 1 / (2*alpha);
xrk2(1) = t0;
yrk2(1) = yINI;
for i = 1:N
    k1 = func(xrk2(i), yrk2(i));
    k2 = func(xrk2(i)+alpha*h, yrk2(i) + beta * k1 * h);
    xrk2(i+1) = xrk2(i) + alpha *h;
    yrk2(i+1) = yrk2(i) + c1*h*k1 + c2*h*k2;
end
