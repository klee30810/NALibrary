/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 01-04-2019
Modified        : 01-04-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 4
/------------------------------------------------------------------------------------------*/


#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H 

#include "myMatrix.h"

/*--------------------------------------------------------------------------*/
/*							Non-linear Solver 								*/
/*--------------------------------------------------------------------------*/

/* Function
	_x		: variable */
double func(double _x);

/* Differential Function
	_x		: variable */
double dfunc(double _x);

/* Bisection Method
	_a		: initial value #1
	_b		: initial value #2
	_tol	: tolerance	*/
double bisectionNL(float _a, float _b, float _tol);

/* Newton Raphson Method
	_x0		: Initial value
	_tol	: tolerance	*/
double newtonRaphson(float _x0, float _tol);


/*--------------------------------------------------------------------------*/
/*							Linear Equation Solver 							*/
/*--------------------------------------------------------------------------*/


/* Gauss Elimination with Scaled Partial Pivoting
	_A		: input matrix _A
	_b		: input vector _b
	_U		: output matrix _U (upper triangular matrix)
	_bn		: output vector _bn
	_P		: output matrix _P (permutation matrix for both A & b))*/
void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _bn, Matrix _P);

/* LU Decomposition with Scaled Partial Pivoting
	_A		: input matrix _A
	_L		: output matrix _L (lower triangular matrix)
	_U		: output matrix _U (upper triangular matrix)
	_P		: output matrix _P (permutation matrix for both A & b) */
void LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P);


/* Back-substitution
	_U		: input matrix _U(upper triangular matrix)
	_bn		: input vector _bn
	_x		: output vector _x*/
void backsub(Matrix _U, Matrix _bn, Matrix _x);

/* Forward substitution
	_L		: input matrix _L(lower triangular matrix)
	_bn		: input vector _bn
	_y		: output vector _y*/
void fwdsub(Matrix _L, Matrix _bn, Matrix _y);

/* Solving linear system problem
	_A			: input matrix _A
	_b			: input vector _b
	_method		: input character _mehtod("g" or "l")
	_x			: output vector _x	*/
void solveLinear(Matrix _A, Matrix _b, const char* _method, Matrix _x);


/*--------------------------------------------------------------------------*/
/*							Eigenvalues & Eigenvectors						*/
/*--------------------------------------------------------------------------*/


/* 2nd Norm
	_A		: input matrix _A */
double  norm2(Matrix _A);

/* QR Decomposition
	_A		: input matrix _A
	_Q		: output matrix _Q (orthogonal matrix)
	_R		: output matrix _R (upper triangular matrix) */
void QRdecomp(Matrix _A, Matrix _Q, Matrix _R);

/* Eigen values
	_A		: input matrix _A	
	_Output : 3x1 output matrix(3 values) */
Matrix eig(Matrix _A);

/* Eigen vecotors
	_A		: input matrix _A
	_Output : 3x3 output matrix(3 vectors) */
Matrix eigvec(Matrix _A);

#endif