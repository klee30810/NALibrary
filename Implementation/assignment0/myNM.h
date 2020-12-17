/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Kyung-Min Lee
Created         : 26-02-2019
Modified        : 26-02-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 0
/------------------------------------------------------------------------------------------*/

#ifndef		_myNM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_myNM_H

#define M 3			// # of rows  ¡é
#define N 3			// # of cols  ¡æ


/*--------------------------------------------------------------------------*/
/*							Vector Arithmetic								*/
/*--------------------------------------------------------------------------*/
/* Add Vector (a + b)
	_a		: input vector #1
	_b		: input vector #2
	_out	: output vector (_a +_b)
	_vecLen	: length of input vector	*/
void addVec(double _a[], double _b[], double _out[], int _vecLen);

/* Subtract Vector (a - b)
	_a		: input vector #1
	_b		: input vector #2
	_out	: output vector (_a -_b)
	_vecLen	: length of input vector	*/
void subVec(double _a[], double _b[], double _out[], int _vecLen);

/* Dot product (a & b)
	_a		: input vector #1
	_b		: input vector #2
	_out	: output number (_a *_b) (pass by reference)
	_vecLen	: length of input vector	*/
void dotProduct(double _a[], double _b[], double* _out, int _vecLen);


/* Print Vector Elements
	_a		: input vector
	_vecLen	: length of input vector	*/
void printVec(double _a[], int _vecLen);


/*--------------------------------------------------------------------------*/
/*							Matrix Arithmetic								*/
/*--------------------------------------------------------------------------*/
/* Add Matrix (A + B)
	_A		: input matrix #1
	_B		: input matrix #2
	_Out	: output matrix (_A +_B)
	_A_rows	: rows of input matrix _A
	_A_cols	: cols of input matrix _A	*/
void addMat(double _A[][N], double _B[][N], double _Out[][N], int _A_rows, int _A_cols);

/* Subtract Matrix (A - B)
	_A		: input matrix #1
	_B		: input matrix #2
	_Out	: output matrix (_A -_B)
	_A_rows	: rows of input matrix _A
	_A_cols	: cols of input matrix _A	*/
void subMat(double _A[][N], double _B[][N], double _Out[][N], int _A_rows, int _A_cols);

/* Multiply Matrix (A & B)
	_A		: input matrix #1
	_B		: input matrix #2
	_Out	: output matrix (_A *_B)
	_A_rows : rows of input matrix _A
	_A_cols : cols of input matrix _A
	_B_rows : rows of input matrix _B
	_B_cols : cols of input matrix _B	*/
void multMat(double _A[][N], double _B[][N], double _Out[][N], int _A_rows, int _A_cols, int _B_rows, int _B_cols);

/* Print Matrix Elements
	_A		: input matrix
	_rows	: rows of input matrix
	_cols	: cols of input matrix	*/
void printMat(double _A[][N], int _rows, int _cols);



#endif

