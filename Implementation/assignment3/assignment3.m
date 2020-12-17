%% Assignment 3
% Initialization

clear; close all; clc;


%% Gauss eilimination when A is square matrix

addpath('assignment3');

A = [1	3	-2	4;
    2	-3	3	-1;
    -1	7	-4	2;
    3	-1	6	2];
 
b = [-11;6;-9;15];

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
       
       for j = k : col_size                
           A_temp(i,j) = A_temp(i,j) - m(i) * A_temp(k, j);
           P_temp(i,j) = P_temp(i,j) - m(i) * P_temp(k,j);
       end
       %bn(i) = bn(i) - m(i) * bn(k);
   end
end


for i = 1:row_size
    for j = 1:col_size
        U(i,j) = A_temp(i,j);
        P(i,j) = P_temp(i,j);
    end
end

%% Forward substitution Ly = b

Y = zeros(row_size,1);
m = zeros(row_size,1);

%bn = P * b;

Y(1,1) = b(1) / A_temp(1,1);
for i = 2:row_size
   m(i) = L(i,1:i-1) * Y(1:i-1);
   Y(i) = (b(i) - m(i) ) / L(i,i);          
end


%% Back substitution for Ux = Y

X = zeros(row_size,1);
m = zeros(row_size,1);

for i = row_size:-1:1
   m(i) =  U(i,i+1:row_size) * X(i+1:row_size);
   X(i) = (Y(i) - m(i) )/U(i,i);
end




%% Gauss elimination inverse matrix

row_size = size(A,1); 
col_size = size(A,2);
I = eye(row_size, col_size);
A_temp = [A I];

for k = 1:row_size-1
   % upper elimination
   for i = k+1 : row_size
       m(i-1) = (A_temp(i,k)/A_temp(k,k));
     %  a(k,1) = 0;
       for j = k : col_size + col_size           % gausselim 에도 업데이트 시켜줄 것
           A_temp(i,j) = A_temp(i,j) - m(i-1) * A_temp(k, j);
       end
   end
end



  % lower elimination
for k = row_size:-1:1  
   for i = (k-1):-1:1
       m(i-1) = (A_temp(i,k)/A_temp(k,k));
     %  a(k,1) = 0;
       for j = k:col_size + col_size
           A_temp(i,j) = A_temp(i,j) - m(i-1) * A_temp(k, j);
       end
   end 
end

    % making jordan form
for i = 1:row_size
    for j = col_size+1:col_size +col_size
        A_temp(i,j) = A_temp(i,j)/ A_temp(i,i);
        A_temp(i,i) = 1;
        I(i,j-col_size) = A_temp(i,j);
    end
end

 

%% Gauss elimination
A_temp = [A b];
row_size = size(A,1); 
col_size = size(A,2);
L = eye(row_size, col_size);
P = eye(row_size, col_size);
U = zeros(row_size, col_size);

for k = 1:row_size-1

   % pivoting
   for j = k:row_size
       max = A_temp(j,k);
       max_idx = j;
       for i = j+1 : row_size
           if (A_temp(i,k)) > max
               max = A_temp(i,k);
               max_idx = i;
           end
           temp = P(j,:);
           P(j,:) = P(max_idx,:);
           P(max_idx,:) = temp;
       end
       % switching
       temp = A_temp(j,:);
       A_temp(j,:) = A_temp(max_idx,:);
       A_temp(max_idx,:) = temp;       
   end
   % elimination
   
   for i = k+1 : row_size
       m(i) = (A_temp(i,k)/A_temp(k,k));
       L(i,k) = m(i);
     %  a(k,1) = 0;
       for j = k : col_size + 1
           A_temp(i,j) = A_temp(i,j) - m(i) * A_temp(k, j);
       end
       %bn(i) = bn(i) - m(i) * bn(k);
   end
end

for i=1:row_size
    bn(i,1) = A_temp(i,5);
end

for i = 1:row_size
    for j = 1:col_size
        U(i,j) = A_temp(i,j);
    end
end    
 

%% Back substitution for gauss_elimination

X = zeros(row_size,1);
m = zeros(row_size,1);

for i = row_size:-1:1
   m(i) =  U(i,i:row_size) * X(i:row_size);
   X(i) = (bn(i) - m(i) )/U(i,i);
end


