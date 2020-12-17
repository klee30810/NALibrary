function dydx = myRC(x,y)       % 1st order ODE, RC circuit
tau=0.01; T=1/tau; f=100; Vm=1;
dydx =-T*y + 1*T*Vm*cos(2*pi*f*x+pi/4);