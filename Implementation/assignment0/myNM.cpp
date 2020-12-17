/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Kyung-Min Lee
Created         : 26-02-2019
Modified        : 26-02-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 0
/------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "myNM.h"

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
		_out[i]= _a[i] + _b[i];
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
			if(_A_cols == _B_rows)
			/*	 ADD YOUR LOGIC HERE (2 lines)	*/
				for(int k = 0; k < _A_cols; k++) 
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