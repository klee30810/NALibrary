%% Assignment 4
% Initialization

clear; close all; clc;


%% Given Data

A = [45 30 -25;
    30 -24 68;
    -25 68 80];
 


%% Execution
U=A; % initialize A
for i = 66
 % A=QR using HouseHold Matrix
 fprintf("%d'th trial!!\n",i);
 [Q R] = QRFactorization(U);
 U = R*Q; % A_new=RQ
end
lamda = diag(U)

%% QR decomposition function
function [Q R] = QRFactorization(R)

% The function factors a matrix [A] into an orthogonal matrix [Q]
% and an upper-triangular matrix [R].
% Input variables:
% A The (square) matrix to be factored.
% Output variables:
% Q Orthogonal matrix.
% R Upper-triangular matrix.
n = size(R,1);
I = eye(n);
Q = I;
    for j = 1:n-1
         c = R(:,j);
         c(1:j-1) = 0;
         e(1:n,1)=0;
         if c(j) > 0
            e(j) = 1;
         else
             e(j) = -1;
         end
         v = c + norm(c)*e;
         H = I - 2/(v'*v)*(v*v');

         Q = Q*H
         R = H*R
    end
end

