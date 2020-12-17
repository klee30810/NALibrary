/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Kyung-Min Lee
Created          : 26-03-2018
Modified         : 26-03-2019
Language/ver     : C++ in MSVS2015

Description      : myMatrix.cpp
/----------------------------------------------------------------*/

#include "myMatrix.h"

// Create Matrix with specified size
Matrix	createMat(int _rows, int _cols)
{
	// check matrix dimension
	if (_rows < 0 || _cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'createMat' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}		

	Matrix Output;
	// 1. Allocate row array first
	Output.at = (double**)malloc(sizeof(double*) * _rows);
	// 2. Then, allocate column 
	for (int i = 0; i < _rows; i++)
		Output.at[i] = (double*)malloc(sizeof(double) * _cols);
	// 3. Initialize row & column values of a matrix
	Output.rows = _rows;
	Output.cols = _cols;

	return Output;
}

// Free a memory allocated matrix
void freeMat(Matrix _A)
{
	// 1. Free allocated column memory
	for (int i = 0; i < _A.rows; i++)
		free(_A.at[i]);
	// 2. Free allocated row memory
	free(_A.at);
}

// Create a matrix from a text file
Matrix	txt2Mat(string _filePath, string _fileName)
{
	ifstream file;
	string temp_string, objFile = _filePath + _fileName + ".txt";
	int temp_int = 0, nRows = 0;

	file.open(objFile);
	if (!file.is_open()) {
		printf("\n*********************************************");
		printf("\n  Could not access file: 'txt2Mat' function");
		printf("\n*********************************************\n");
		return createMat(0, 0);
	}
	while (getline(file, temp_string, '\t'))
		temp_int++;
	file.close();

	file.open(objFile);
	while (getline(file, temp_string, '\n'))
		nRows++;
	file.close();

	int nCols = (temp_int - 1) / nRows + 1;
	Matrix Output = createMat(nRows, nCols);

	file.open(objFile);
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++) {
			file >> temp_string;
			Output.at[i][j] = stof(temp_string);
		}
	file.close();

	return Output;
}

void initMat(Matrix _A, double _val)
{
	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols; j++)
		{
			_A.at[i][j] = _val;
		}
	}
}

void printMat(Matrix _A)
{
	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols; j++)
		{
			printf("%f\t", _A.at[i][j]);
		}
		printf("\n");
	}
}

Matrix zeros(int _rows, int _cols)
{
	Matrix _Mat = createMat(_rows, _cols);

	for (int i = 0; i < _Mat.rows; i++)
	{
		for (int j = 0; j < _Mat.cols; j++)
		{
			_Mat.at[i][j] = 0;
		}
	}

	return _Mat;
}

Matrix ones(int _rows, int _cols)
{
	Matrix _Mat = createMat(_rows, _cols);

	for (int i = 0; i < _Mat.rows; i++)
	{
		for (int j = 0; j < _Mat.cols; j++)
		{
			if 
				(i == j) _Mat.at[i][j] = 1;
			else 
				_Mat.at[i][j] = 0;
		}
	}

	return _Mat;
}

Matrix	eye(int _rows, int _cols)
{
	Matrix _Mat = createMat(_rows, _cols);

	for (int i = 0; i < _Mat.rows; i++)
	{
		for (int j = 0; j < _Mat.cols; j++)
		{
			if (i == j)
				_Mat.at[i][j] = 1;
			else
				_Mat.at[i][j] = 0;
		}
	}

	return _Mat;
}

Matrix	transpose(Matrix _A)
{
	Matrix _Mat = createMat(_A.rows, _A.cols);
	
	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols; j++)
		{
			_Mat.at[i][j] = _A.at[j][i];
		}
	}
	return _Mat;
}

Matrix	copyMat(Matrix _A)
{
	Matrix _Mat = createMat(_A.rows, _A.cols);

	for (int i = 0; i < _Mat.rows; i++)
	{
		for (int j = 0; j < _Mat.cols; j++)
		{
			_Mat.at[i][j] = _A.at[i][j];
		}
	}

	return _Mat;
}

void copyVal(Matrix _A, Matrix _B) 
{
	if ((_A.rows != _B.rows) || (_A.cols != _B.cols))
		printf("The size of Mat A and Mat B are different!\n");
	else
	{
		for (int i = 0; i < _A.rows; i++)
		{
			for (int j = 0; j < _A.cols; j++)
			{
				_B.at[i][j] = _A.at[i][j];
			}
		}
	}
}