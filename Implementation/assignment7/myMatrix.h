/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Kyung-Min Lee
Created         : 01-04-2019
Modified        : 17-05-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 7 - myMatrix.h
/------------------------------------------------------------------------------------------*/

#ifndef		_MY_MATRIX_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_MATRIX_H

#include <iostream>
#include <string>
#include <fstream>
#define _USE_MATH_DEFINES
#define PI 3.1415926535897
#include "math.h"

typedef struct {
	double** at;
	int rows, cols;
}Matrix;

using namespace std;

// Create Matrix with specified size
extern	Matrix	createMat(int _rows, int _cols);

// Free a memory allocated matrix
extern	void	freeMat(Matrix _A);

// Create a matrix from a text file
extern	Matrix	txt2Mat(string _filePath, string _fileName);


/// It is recommended to create the following functions.

// initialization of Matrix elements
extern	void	initMat(Matrix _A, double _val);

// Create matrix of all zeros
extern	Matrix	zeros(int _rows, int _cols);

// Create matrix of all ones
extern	Matrix	ones(int _rows, int _cols);

// Create identity 
extern	Matrix eye(int _rows, int _cols);;

// Create Transpose matrix
extern	Matrix	transpose(Matrix _A);

// Copy matrix
extern	Matrix	copyMat(Matrix _A);

// Copy matrix Elements from A to B
extern	void	copyVal(Matrix _A, Matrix _B);

//// Print matrix
extern	void	printMat(Matrix _A);

//// Add column matrix
extern	Matrix	addColumn(Matrix _A, Matrix _B);

//// Add row matrix
extern Matrix	addRow(Matrix _A, Matrix _B);


/*--------------------------------------------------------------------------*/
/*							Vector Arithmetic								*/
/*--------------------------------------------------------------------------*/

/* Add Vector (a + b)
	_a		: input vector #1
	_b		: input vector #2
	_out	: output vector (_a +_b)
	_vecLen	: length of input vector	*/
void addVec(Matrix _A, Matrix _B, Matrix _out);

/* Subtract Vector (a - b)
	_a		: input vector #1
	_b		: input vector #2
	_out	: output vector (_a -_b)
	_vecLen	: length of input vector	*/
void subVec(Matrix _A, Matrix _B, Matrix _out);

/* Dot product (a & b)
	_a		: input vector #1
	_b		: input vector #2
	_out	: output number (_a *_b) (pass by reference)
	_vecLen	: length of input vector	*/
void dotProduct(Matrix _A, Matrix _B, double _out);

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
	_Out	: output matrix (_A +_B)	*/
void addMat(Matrix _A, Matrix _B, Matrix _Out);

/* Subtract Matrix (A - B)
	_A		: input matrix #1
	_B		: input matrix #2
	_Out	: output matrix (_A -_B)
	_A_rows	: rows of input matrix _A
	_A_cols	: cols of input matrix _A	*/
void subMat(Matrix _A, Matrix _B, Matrix _Out);

/* Multiply Matrix (A & B)
	_A		: input matrix #1
	_B		: input matrix #2
	_Out	: output matrix (_A *_B)	*/
void multiMat(Matrix _A, Matrix _B, Matrix _Out);

/* Multiply Matrix with Return (A & B)
	_A		: input matrix #1
	_B		: input matrix #2
	_Out	: output matrix (_A *_B)	*/
extern Matrix multiMat(Matrix _A, Matrix _B);

/* Deriving Inverse Matrix
	_A		: input matrix _A
	_Ainv	: output matrix => inverse matrix of _A*/
	//Matrix invMat(Matrix _A);
void invMat(Matrix _A, Matrix _Out);


#endif