/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Kyung-Min Lee
Created         : 13-03-2019
Modified        : 26-03-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 2
/------------------------------------------------------------------------------------------*/
#ifndef		_MY_MATRIX_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_MATRIX_H

#include <iostream>
#include <string>
#include <fstream>
#define _USE_MATH_DEFINES
#include "math.h"

typedef struct {
	double** at;
	int rows, cols;
}Matrix;

using namespace std;

// Create Matrix with specified size
extern	Matrix createMat(int _rows, int _cols);

// Free a memory allocated matrix
extern	void freeMat(Matrix _A);

// Create a matrix from a text file
extern	Matrix txt2Mat(string _filePath, string _fileName);


/// It is recommended to create the following functions.

// initialization of Matrix elements
extern void initMat(Matrix _A, double _val);

// Create matrix of all zeros
extern Matrix zeros(int _rows, int _cols);

// Create matrix of all ones
extern Matrix ones(int _rows, int _cols);

// Create identity 
extern	Matrix eye(int _rows, int _cols);;

// Create Transpose matrix
extern	Matrix transpose(Matrix _A);

// Copy matrix
extern	Matrix copyMat(Matrix _A);

// Copy matrix Elements from A to B
extern	void copyVal(Matrix _A, Matrix _B);

//// Print matrix
extern	void printMat(Matrix _A);

//// Add column matrix
extern Matrix addColumn(Matrix _A, Matrix _B);

//// Add row matrix
extern Matrix addRow(Matrix _A, Matrix _B);


#endif