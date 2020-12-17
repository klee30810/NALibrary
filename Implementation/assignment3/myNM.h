/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Kyung-Min Lee
Created         : 12-04-2019
Modified        : 17-05-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 7 - myNM.h
------------------------------------------------------------------------------------------*/


#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H 

#include "myMatrix.h"
/*--------------------------------------------------------------------------*/
/*							Basic Operation 								*/
/*--------------------------------------------------------------------------*/
/* Matrix multiplication
	_A		: input matrix _A
	_B		: input matrix _B */
Matrix	multiMat(Matrix _A, Matrix _B);

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

/* Partial Pivoting
	_U		: input matrix _U (upper triangular matrix)
	_b		: input vector _b
	_P		: output matrix _P
	_k		: index integer _k */
void partialPivoting(Matrix _U, Matrix _b, Matrix _P, int _k);

/* Gauss Elimination with Partial Pivoting
	_A		: input matrix _A
	_b		: input vector _b
	_U		: output matrix _U (upper triangular matrix)
	_bn		: output vector _bn
	_P		: partial pivoted Matrix _P */
void gaussElim(Matrix _A, Matrix _b, Matrix &_U, Matrix _bn, Matrix & _P);

/* LU Decomposition with Scaled Partial Pivoting
	_A		: input matrix _A
	_L		: output matrix _L (lower triangular matrix)
	_U		: output matrix _U (upper triangular matrix)
	_P		: output matrix _P (permutation matrix) */
void LUdecomp(Matrix _A, Matrix &_L, Matrix &_U, Matrix _b, Matrix _P);


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
void	solveLinear(Matrix _A, Matrix _b, const char* _method, Matrix _x);


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

/*--------------------------------------------------------------------------*/
/*							Curve Fitting					*/
/*--------------------------------------------------------------------------*/

/* Curve fitting
	_x		: input vector _x
	_y		: input vector _y
	_n		: input value n (order of polynomial) 		
	_a		: output vector _a (coefficients)
	_e		: output value _e (sum of square error) */
void polyfit(Matrix _x, Matrix _y, int _n, Matrix &_a, double &_e);


/*--------------------------------------------------------------------------*/
/*							Differentiation					*/
/*--------------------------------------------------------------------------*/

/* Convert radian from degree
	_deg : input angle in degree _deg
	_rad : output angle in radian _rad */
double getRadian(double _deg);

/* 1st order derivation
	_v		: output vector _v
	_x		: input data vector _x
	_y		: input data vector _y
	_h		: the interval size
	_method : the method of derivation(2 forward, 2 bacward, 2 central, 3 backward, 3 forward*/
void deriv1(Matrix &_v, Matrix _x, Matrix _y, Matrix &_error, double _h, const char* _method);

/* 2nd order derivation
	_a		: output vector _a
	_x		: input data vector _x
	_y		: input data vector _y
	_h		: the interval size
	_method : the method of derivation(3 central, 3 backward, 3 forward*/
void deriv2(Matrix &_a, Matrix _x, Matrix _y, Matrix &_error, double _h, const char* _method);


/*--------------------------------------------------------------------------*/
/*							integration					*/
/*--------------------------------------------------------------------------*/

/* Function for the Problem 1
	_x		: variable */
double func1(double _x);

/* Function for the Problem 2
	_x		: variable */
double func2(double _x);

/* Trapezoidal Method
	_func	: input function _func
	_a		: input left end _a
	_b		: input right end _b
	_N		: the number of intervals */
double trapezoidal(double _func(double _x), double a, double b, int _N);

/* Simpson 1/3 Method
	_func	: input function _func (with input arguments _x)
	_a		: input left end _a
	_b		: input right end _b
	_N		: the number of intervals */
double simpson13(double _func(double _x), double a, double b, int _N);

/* Adaptive Simpson Method
	_func		: input function _func (with input arguments _x)
	_a			: input left end _a
	_b			: input right end _b
	_epsilon	: tolerance of the error _epsilon
	_N			: the number o intervals
	_itr		: the current number of iterations
	_itrMax		: the max number of iterations */
double adaptSimpson(double _func(double _x), double a, double b, double _epsilon, int _N, int _itr, int _itrMax);

#endif