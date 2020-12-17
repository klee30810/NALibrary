%% Initialization
clear;  clc;   %close all;

%% Function and Constants

% Constants
m = 2500;       %[]
k = 300000;     %[]
c = 36000;

syms w;

% Function and Derivation
F = @(w) w^-1 - 2;
%4*(m^2)*(w^4) - (8*k*m + 21*(c^2))*(w^2) - 21*(k^2);
dF = @(w) -1*w^-2;
%16 * (m^2) * (w^3) + ((-16 * k * m) - (42*(c^2))) * w ;


% plotting the function
 w = 0:0.01:1;
 %y = sqrt( ((w.^2*c^2) + k^2) / ( (k-m*(w.^2)).^2 + (w.*c).^2 )) - 0.4;
 y = w.^-1 - 2;
 dy = -w.^-2 ;
 z = 0*w;
 hold on
 plot(w, y)
 %plot(w, dy)
 plot(w,z)
 hold off

% Found that the local optimum is in (-40,-30), (30, 40)
X = fzero(F, 1.4);       % Finding the roots with fzero function near x0 = 37 => 37.1926
%% Bisection Method

% initial interval & Maxit
a = 0;
b = 10;
MaxIt = 100;
epsilon = 0.00001;

% Function results for Bisection Methods
Fa = F(a); Fb = F(b);
dFa = dF(a); dFb = dF(b);

% Algorithms for Bisection Method
if Fa * Fb > 0
    disp("a & b have the same sign!");
else
    disp("iteration     a       b         (xNS)        Tolerance");
    for i = 1:MaxIt
        x_ns = (a+b)/2;
        Fx_ns = F(x_ns);
        tol = (b-a)/2;
        
        fprintf("%3i  %11.6f %11.6f %11.6f %11.6f \n", i, a, b, x_ns, tol);
        if Fx_ns == 0
            fprintf("An exact solution x = %11.6f was not found.", x_ns);
            break;
        end
        
        if tol < epsilon
            disp("Tolerance condition satisfied!");
            break;
        end
        
        if F(a) * Fx_ns < 0
            b = x_ns;
        else
            a = x_ns;
        end
    end
end


%% Newton Rapson Method

x_c = 1.4;        % Initial condition
MaxIt = 100;
epsilon = 0.00001;
i=1;


while i <= MaxIt
    if (dF(x_c) == 0)
        fprintf("dF(%f) is zero!!",x_n);
        break
    end
    
    x_n = x_c - (F(x_c)/dF(x_c));                                 % Newton-Raphson method 
    tol = abs((x_n - x_c)/(x_c));
    fprintf("%3i x(n):%11.6f tolerance:%11.6f \n", i, x_n, tol);
    
    if tol < epsilon                     % stopping criterion when difference between iterations is below tolerance
        
        fprintf('Solution is %f \n', (x_n));
        break
    end
 
    i = i + 1;
    x_c = x_n;             %update 
end

if i == (MaxIt+1)
    fprintf('Solution did not coverge within %d iterations at a required precision of %d \n', MaxIt, epsilon)     %error for non-convergence within N iterations
end


%% Advanced Problem

a = 0;
b = 0.5;
x_c = 1.4;
MaxIt = 30;
epsilon = 0.00001;
i=0;


% Algorithms for Bisection Method

if F(a) * F(b) > 0
    disp("a & b have the same sign!");
end

while i <= MaxIt
    if (dF(x_c) == 0)
        fprintf("dF(%f) is zero!!",x_n);
        break
    end
    x_n = x_c - (F(x_c)/dF(x_c));                                 % Newton-Raphson method 
    tol = abs((x_n - x_c)/(x_c));
    fprintf("%3i x(n):%11.6f tolerance:%11.6f \n", i+1, x_n, tol);
    
    if (x_n>b)                                      % applying bisection method When NR gives the solution out of bounds
       x_n = (x_n - (b-a)) / 2;
    elseif (x_n<a)
       x_n = (x_n + (b-a)) / 2; 
    end
    
    if tol < epsilon                                % stopping criterion when difference between iterations is below tolerance
        
        fprintf('Solution is %f \n', (x_n));
        break
    end
 
    i = i + 1;
    x_c = x_n;             %update 
end

if i == (MaxIt+1)
    fprintf('Solution did not coverge within %d iterations at a required precision of %d \n', MaxIt, epsilon)     %error for non-convergence within N iterations
end

