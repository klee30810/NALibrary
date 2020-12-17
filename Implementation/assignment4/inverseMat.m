%% Assignment 3
% Initialization

clear; close all; clc;


%% Gauss eilimination when A is square matrix


A = [45 30 -25;
    30 -24 68;
    -25 68 80];
 

%% Gauss elimination inverse matrix

row_size = size(A,1); 
col_size = size(A,2);
I = eye(row_size, col_size);
A_temp = [A I];

for k = 1:row_size-1
   % upper elimination
   for i = k+1 : row_size
       m(i-1) = (A_temp(i,k)/A_temp(k,k));      % m => m-1 로 바꿈
     %  a(k,1) = 0;
       for j = k : col_size + col_size           % gausselim 에도 업데이트 시켜줄 것
           A_temp(i,j) = A_temp(i,j) - m(i-1) * A_temp(k, j); % m => m-1 로 바꿈
       end
   end
end

m = zeros(3,1);

  % lower elimination
for k = row_size:-1:1  
   for i = (k-1):-1:1
       m(i) = (A_temp(i,k)/A_temp(k,k));                    % m(i-1) => m(i)
     %  a(k,1) = 0;
       for j = k:col_size + col_size
           A_temp(i,j) = A_temp(i,j) - m(i) * A_temp(k, j);
       end
   end 
end

    % making jordan form
for i = 1:row_size
    m(i) = A_temp(i,i);
    for j = 1: col_size*2               % j = i:col_size*2   ???
        A_temp(i,j) = A_temp(i,j)/ m(i);  
    end
end

for i = 1:row_size
    for j = col_size+1 : col_size * 2
        I(i,j-col_size) = A_temp(i,j);
    end
end




%% LU decompostion

A_temp = A;
row_size = size(A,1); 
col_size = size(A,2);
L = eye(row_size, col_size);
P = eye(row_size, col_size);
P_temp = P;
U = zeros(row_size, col_size);

for k = 1:row_size-1
   % elimination
   
   for i = k+1 : row_size
       m(i) = (A_temp(i,k)/A_temp(k,k));
       L(i,k) = m(i);
     %  a(k,1) = 0;
       for j = k : col_size                 % +1 빠짐, gausselim에도 고쳐놓을 것
           A_temp(i,j) = A_temp(i,j) - m(i) * A_temp(k, j);
           P_temp(i,j) = P_temp(i,j) - m(i) * P_temp(k,j);
       end
   end
end


for i = 1:row_size
    for j = 1:col_size
        U(i,j) = A_temp(i,j);
    end
end


%% Forward substitution Ly = b

Y = zeros(row_size,1);
m = zeros(row_size,1);
b = [1;0;0];

for i = 1:row_size
   m(i) = L(i,i:row_size) * Y(i:row_size);
   Y(i) = (b(i) - m(i) ) / L(i,i);          
end


%% Back substitution for Ux = Y

X = zeros(row_size,1);
m = zeros(row_size,1);

for i = row_size:-1:1
   m(i) =  U(i,i:row_size) * X(i:row_size);
   X(i) = (Y(i) - m(i) )/U(i,i);
end




