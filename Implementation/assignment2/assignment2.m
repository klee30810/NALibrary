%% Assignment 2
% Initialization

clear; close all; clc;

%% Gauss eilimination when A is square matrix

A = [1,3,-2,4;
     2,-3,3,-1;
     -1,7,-4,2;
     3,-1,6,2;];
 
b = [-11;6;-9;15];
bn = b;
%% LU decomposition
row_size = size(A,1); 
col_size = size(A,2);
L = eye(row_size, col_size);

for k = 1:row_size-1
   for i = k+1 : row_size
       m(i) = (A(i,k)/A(k,k));
       L(i,k) = m(i);
     %  a(k,1) = 0;
       for j = k : col_size
           A(i,j) = A(i,j) - m(i) * A(k, j);
       end
       bn(i) = bn(i) - m(i) * bn(k);
   end
end

for i = 1:row_size
    for j = 1:col_size
        U(i,j) = A(i,j);
    end
end

%% Back substitution

X = zeros(row_size,1);
%X(row_size) = A(row_size, col_size)/A(row_size,row_size);
m = zeros(row_size,1);

for i = row_size:-1:1
   m(i) =  U(i,i+1:row_size) * X(i+1:row_size);
   X(i) = (bn(i) - m(i) )/U(i,i);
end
