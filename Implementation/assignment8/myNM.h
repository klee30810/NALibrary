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
void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _bn, Matrix _P);

/* LU Decomposition with Scaled Partial Pivoting
	_A		: input matrix _A
	_L		: output matrix _L (lower triangular matrix)
	_U		: output matrix _U (upper triangular matrix)
	_P		: output matrix _P (permutation matrix) */
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

/* Function for the Question 1, problem 1-c
	_x		: variable */
double _1cfunc(double _x);

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


/*--------------------------------------------------------------------------*/
/*							Test2_Preparation					*/
/*--------------------------------------------------------------------------*/


/* Interpolation using a Lagrange polynomial
	_x: x values of input data
	_y: y values of input data
	_value: output interpolating value*/
double LagrangeTNT(Matrix _x, Matrix _y, double _value);

/* Interpolation using a Lagrange polynomial
	_x: x values of input data
	_y: y values of input data
	_xint: The x coordinate of the at which y is to be interpolated */
double NewtonsTNT(Matrix _x, Matrix _y, double _xint);

/* Example function of 1st derivative ODE
	_x: input x data
	_y: input y data*/
double examplefunc(double _x, double _y);

/* Solving 1st order ODE with 2nd order Runge-Kutta
	_func: input 1st derivative function
	_a: left limit of the range
	_b: right limit of the range
	_h: size of interval
	_yinit: inital value of y	*/
void odeRK2(double func(double _x, double _y), double _a, double _b, double _h, double _yinit, Matrix _x, Matrix _y);


/* function for test2 Question 2 - b
_x: input variable */
//double _2bfunc(double _t, double _x);

/* Solving 1st order ODE with 3rd order Runge-Kutta (classical method)
	_func: input 1st derivative function
	_a: left limit of the range
	_b: right limit of the range
	_h: size of interval
	_yinit: inital value of y
	_X: output matrix x
	_Y output matrix y:*/
void odeRK3(double _func(double _x, double _y), double _a, double _b, double _h, double _yinit, Matrix _X, Matrix _Y);

/* mck function of dydt
	_x: input x data
	_y: input y data
	_v: input v data */
double mckdydt(double _x, double _y, double _v);

/* mck function of dvdt
	_x: input x data
	_y: input y data
	_v: input v data */
double mckdvdt(double _x, double _y, double _v);


/* Solving two 1st order ODEs with 3rd order Runge-Kutta (classical method)
	_mckdydt: input 1st derivative function
	_mckdvdt: input 2nd derivative function
	_a: left limit of the range
	_b: right limit of the range
	_h: size of interval
	_yinit: inital value of y
	_X: output matrix x
	_Y: output matrix y
	_V: output matrix v */
void ode2RK3(double _mckdydt(double _x, double _y, double _v), double _mckdvdt(double _x, double _y, double _v), double _a, double _b, double _h, double _yinit, double _vinit, Matrix &_X, Matrix &_Y, Matrix &_V);

/* velocity from two-point central method
	_y: input matrix which will be differentiated
	_h: the input length of interval	*/
Matrix derivative(Matrix _y, double _h);
#endif