/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Kyung-Min Lee
Created         : 10-03-2019
Modified        : 26-03-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 2
/------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myNM.h"
#include "myMatrix.h"


/*--------------------------------------------------------------------------*/
/*							Vector Arithmetic								*/
/*--------------------------------------------------------------------------*/

/* Add Vector (a + b)
	_a		: input vector #1
	_b		: input vector #2
	_out	: output vector (_a +_b)
	_vecLen	: length of input vector	*/
void addVec(double _a[], double _b[], double _out[], int _vecLen) {
	for (int i = 0; i < _vecLen; i++) {
		/*	 FILL IN THE BLANKS HERE 	*/
		_out[i] = _a[i] + _b[i];
	}
}

/* Subtract Vector (a - b)
	_a		: input vector #1
	_b		: input vector #2
	_out	: output vector (_a -_b)
	_vecLen	: length of input vector	*/
void subVec(double _a[], double _b[], double _out[], int _vecLen) {
	for (int i = 0; i < _vecLen; i++) {
		/*	 FILL IN THE BLANKS HERE (use pointer '*') 	*/
		*(_out + i) = *(_a + i) - *(_b + i);
	}
}

/* Dot product (a & b)
	_a		: input vector #1
	_b		: input vector #2
	_out	: output number (_a *_b) (pass by reference)
	_vecLen	: length of input vector	*/
void dotProduct(double _a[], double _b[], double* _out, int _vecLen) {
	*_out = 0;
	for (int i = 0; i < _vecLen; i++) {
		/*	 ADD YOUR LOGIC HERE (1 line)	*/
		*(_out) += *(_a + i) * *(_b + i);
	}
}

/* Print Vector Elements
	_a		: input vector
	_vecLen	: length of input vector	*/
void printVec(double _a[], int _vecLen) {
	for (int i = 0; i < _vecLen; i++)
		/*	 ADD YOUR LOGIC HERE (1 line)	*/
		printf("%lf\n", _a[i]);
	printf("\n");
}

/*--------------------------------------------------------------------------*/
/*							Matrix Arithmetic								*/
/*--------------------------------------------------------------------------*/

/* Add Matrix (A + B)
	_A		: input matrix #1
	_B		: input matrix #2
	_Out	: output matrix (_A +_B)
	_A_rows	: rows of input matrix _A
	_A_cols	: cols of input matrix _A	*/
void addMat(double _A[][N], double _B[][N], double _Out[][N], int _A_rows, int _A_cols) {
	for (int i = 0; i < _A_rows; i++)
		for (int j = 0; j < _A_cols; j++) {
			/*	 FILL IN THE BLANKS HERE 	*/
			_Out[i][j] = _A[i][j] + _B[i][j];
		}
}

/* Subtract Matrix (A - B)
	_A		: input matrix #1
	_B		: input matrix #2
	_Out	: output matrix (_A -_B)
	_A_rows	: rows of input matrix _A
	_A_cols	: cols of input matrix _A	*/
void subMat(double _A[][N], double _B[][N], double _Out[][N], int _A_rows, int _A_cols) {
	for (int i = 0; i < _A_rows; i++)
		for (int j = 0; j < _A_cols; j++) {
			/*	 FILL IN THE BLANKS HERE 	*/
			_Out[i][j] = _A[i][j] - _B[i][j];
		}
}

/* Multiply Matrix (A & B)
	_A		: input matrix #1
	_B		: input matrix #2
	_Out	: output matrix (_A *_B)
	_A_rows : rows of input matrix _A
	_A_cols : cols of input matrix _A
	_B_rows : rows of input matrix _B
	_B_cols : cols of input matrix _B	*/
void multMat(double _A[][N], double _B[][N], double _Out[][N], int _A_rows, int _A_cols, int _B_rows, int _B_cols) {
	for (int i = 0; i < _A_rows; i++) {
		for (int j = 0; j < _B_cols; j++) {
			_Out[i][j] = 0;
			if (_A_cols == _B_rows)
				/*	 ADD YOUR LOGIC HERE (2 lines)	*/
				for (int k = 0; k < _A_cols; k++)
					_Out[i][j] += _A[i][k] * _B[k][j];
			else printf("Matrices' sizes are different!");
		}
	}
}

/* Print Matrix Elements
	_A		: input matrix
	_rows	: rows of input matrix
	_cols	: cols of input matrix	*/
void printMat(double _A[][N], int _rows, int _cols) {
	for (int i = 0; i < _rows; i++) {
		for (int j = 0; j < _cols; j++) {
			/*	 ADD YOUR LOGIC HERE (1 lines)	*/
			printf("%lf\t", _A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

/*--------------------------------------------------------------------------*/
/*							Non-linear Solver 								*/
/*--------------------------------------------------------------------------*/

/* Function
	_x		: variable */
double func(double _x)
{
	double m = 2500;											// Setting the given constants
	double k = 300000;
	double c = 36000;

	// Function from the given formula
	return pow(_x, -1) - 2;													// 2nd function
	//4 * pow(m,2) *pow(_x, 4) - (8 * k * m + 21 * pow(c, 2)) * pow(_x, 2) - 21 * pow(k, 2);	=> 1st function
}

/* Differential Function
	_x		: variable */
double dfunc(double _x)
{
	double m = 2500;											// Setting the given constants
	double k = 300000;
	double c = 36000;

	// First derivation of the func.
	return -1 * pow(_x, -2);									// 2nd function
	//16 * pow(m,2) * pow(_x, 3) - (16 * k * m + 42 * pow(c, 2)) * _x;							=> 1st function
}

/* Bisection Method
	_a		: initial value #1
	_b		: initial value #2
	_tol	: tolerance	
	
	Cannot find true value when...
	1) f(x_true) == 0
	2) f(a) * f(b) > 0
	*/
double bisectionNL(float _a, float _b, float _tol)
{
	int i, MaxIt = 30;								// Set the maximum iteration number
	double x_ns;

	for (i = 0; i < MaxIt; i++)						// Iterate until MaxIteration.
	{
		
		if (func(_a)*func(_b) < 0)				// Check whether the interval is reasonable.
		{
			x_ns = (_a + _b) / 2;				// Find the midpoint of a and b
			double tol = fabs((_b - _a)/2);
			
			printf("Iteration:%3d	Left:%11.8lf   Right:%11.8lf   x_NS:%11.8lf   Tolerance: %11.8lf\n", i, _a, _b, x_ns, tol);
			

			if (tol < _tol)						// Check whether x_ns satisfies tolerance. 
			{
				break;
			}
			else
			{
				if (func(_a) * func(x_ns) < 0)
				{
					_b = x_ns;									// case 1. x_true is in between a & x_ns
				}
				else
				{
					_a = x_ns;									// case 2. x_true is in between x_ns & b
				}
			}
		}
		else
		{
			printf("True value(n_T) does not exist between 'a' and 'b' !!\n");
			break;
		}
		
	}
	
	if (i == MaxIt)												// When the solution is not found within the given number of iterations
	{
		printf("The solution is not found within the given number of iterations!! \n");  // Show the n_true is found within the max iteration.
		return x_ns;
	}
	else 
	{
		printf("The solution is found within the given number of iteration!! \n");
		return x_ns;
	}
	
}

/* Newton Raphson Method
	_x0		: Initial value
	_tol	: tolerance	*/
double newtonRaphson(float _x0, float _tol)
{
	double tol, x_NS, x_c = _x0;
	int i = 0, MaxIt = 30;								// Set the maximum iteration number

	while (i <= MaxIt)
	{
		x_NS = x_c - (func(x_c) / dfunc(x_c));			// Calculating derivative and applying 
		tol = (fabs(x_NS - x_c) / fabs(x_c));			// tolerance for NR Method
		printf("Iteration:%3d	x_NS: %11.8lf	tolerance:%11.8lf \n", i, x_NS, tol);

		if (dfunc(x_c) == 0)							// When dfunc = 0 so that x_NS diverses
		{
			printf("dFunc(%11.8lf) is zero!!", x_NS);
			break;
		}
		if (tol < _tol)									// When x_NS satisfies tolerance
		{
			printf("Solution is %11.8lf. \n", x_NS);
			return x_NS;
			break;
		}

		i += 1;
		x_c = x_NS;										// updating the value
	}
	if (i == (MaxIt + 1))								// When the solution is not found within the given number of iterationss
	{
		printf("Solution did not coverge within %d iterations at a required precision of %lf \n", MaxIt, _tol);
		return x_NS;
	}
}


/* Hybrid Method
	_a		: initial value #1 for bisection method
	_b		: initial value #2 for bisection method
	_x0		: Initial value for newton raphson method
	_tol	: tolerance	*/
double NR_hybrid(float _a, float _b, float _x0, float _tol)
{
	int MaxIt = 30, i = 0;					// Set the maximum iteration number
	float x_c = _x0, x_NS, tol;

	if (func(_a)*func(_b) > 0)				// Check whether the root is inside the interval
	{
		printf("True value(n_T) does not exist between 'a' and 'b' !!\n");
		printf("a= %f\t%f\t\tb = %f\t%f\n", _a, func(_a), _b, func(_b));
		return 0;
	}

	while (i <= MaxIt)						// Iteration 
	{
		if (dfunc(x_c) == 0)				// When dfunc = 0 so that x_NS diverses
		{
			printf("dFunc(%11.8lf) is zero!!\n", x_NS);
			break;
		}

		x_NS = x_c - (func(x_c) / dfunc(x_c));		// Calculating derivative and applying
		tol = (fabs(x_NS - x_c) / fabs(x_c));		// tolerance for hybrid Method
		printf("Iteration:%3d	x_NS: %11.8lf	tolerance:%11.8lf \n", i, x_NS, tol);

		if (x_NS > _b)								// When x_NS is in outside right of the interval
			x_NS = (x_NS - (_b - _a)) / 2;
		else if (x_NS < _a)							// When x_NS is in outside left of the interval
			x_NS = (x_NS + (_b - _a)) / 2;


		if (tol < _tol)								// When x_NS satisfies tolerance
		{
			printf("Solution is %11.8lf. \n", x_NS);
			return x_NS;
			break;
		}

		i += 1;
		x_c = x_NS;									// updating the value
	}
	if (i == (MaxIt + 1))							// When the solution is not found within the given number of iterationss
	{
		printf("Solution did not coverge within %d iterations at a required precision of %lf \n", MaxIt, _tol);
		return x_NS;
	}

}

/* Gauss Elimination
	_A		: input matrix _A
	_b		: input vector _b
	_U		: output matrix _U (upper triangular matrix)
	_bn		: output vector _bn */
void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _bn)
{
	double temp = 0;
	Matrix L = eye(_A.rows, _A.cols);						// Initialization of variables
	Matrix P = eye(_A.rows, _A.cols);						// P,L should be identity matrix for the calculation
	Matrix P_temp = copyMat(P);								// make temperary matrices for calculation and preserving original values
	Matrix m = createMat(_b.rows, 1);
	Matrix A_temp = copyMat(_A);
	Matrix _bn_temp = copyMat(_b);

	for (int k = 0; k < (A_temp.rows - 1); k++)
	{
		for (int i = k + 1; i < A_temp.rows; i++)
		{
			m.at[i][0] = (A_temp.at[i][k] / A_temp.at[k][k]);		// make m matrix for the efficient calculation
			if (A_temp.at[k][k] == 0)								// error case: divide by zero
			{
				printf("A[%d][%d] is zero that divides by zero!!\n", k, k);
				break;
			}
			else if (A_temp.rows != A_temp.cols)					// error case: not square matrix
			{
				printf("The matrix A is not square matrix!!\n");
				break;
			}
			else if (A_temp.rows != _b.rows)						// error case: A, b size disagreement
			{
				printf("The sizes of matrix A and vector b are not appropriate!!");
				break;
			}
			L.at[i][k] = m.at[i][0];
			for (int j = k; j < A_temp.cols; j++)					// gauss-elimination
			{
				A_temp.at[i][j] = A_temp.at[i][j] - m.at[i][0] * A_temp.at[k][j];
			}
			_bn_temp.at[i][0] = _bn_temp.at[i][0] - m.at[i][0] * _bn_temp.at[k][0];		// changing b is saved in bn
		}
	}
	
	for (int i = 0; i < _bn.rows; i++)
	{
		_bn.at[i][0] = _bn_temp.at[i][0];							// putting results into real matrix
	}

	for (int i = 0; i < _A.rows; i++)
	{
		
		for (int j = 0; j < _A.cols; j++)
		{
			_U.at[i][j] = A_temp.at[i][j];							// putting results into real matrix
		}
	}

	for (int i = 0; i < _U.rows; i++)								// Building permutation matrix
	{
		if (_U.at[i][i] < _U.at[i + 1][i + 1])						// compare the pivot values
		{
			for (int j = 0; j < _U.cols; j++)						// row exchange according to the pivot values
			{
				temp = P.at[i + 1][j];
				P.at[i + 1][j] = P.at[i][j];
				P.at[i][j] = temp;
			}
		}
		else break;
		
	}

	printf("\n[matrix P]\n\n");										// printing permutation matrix
	printMat(P);
}

/* Back-substitution
	_U		: input matrix _U(upper triangular matrix)
	_bn		: input vector _bn
	_x		: output vector _x*/
void backsub(Matrix _U, Matrix _bn, Matrix _x)
{
	Matrix m = zeros(_x.rows, 1);
	Matrix _x_temp = copyMat(_x);									// make temperary matrices for calculation and preserving original values
	

	for (int i = _U.rows-1; i >= 0; i--)						
	{
		for (int j = i+1; j < _U.cols; j++)
		{
			 m.at[i][0] += _U.at[i][j] * _x_temp.at[j][0];			// pre-calculation of multiplication for efficiency
		}
		_x_temp.at[i][0] = (_bn.at[i][0] - m.at[i][0]) / _U.at[i][i];	// back-substitution
	}

	for (int i = 0; i < _x.rows; i++)
	{
		_x.at[i][0] = _x_temp.at[i][0];								// putting results into real matrix
	}

}
