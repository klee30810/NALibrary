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

// print Matrix
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
	Matrix temp_mat = zeros(_A.cols, _A.rows);
	
	for (int i = 0; i < _A.cols; i++) {
		for (int j = i; j < _A.rows; j++) {
			temp_mat.at[i][j] = _A.at[j][i];
		}
	}

	return temp_mat;
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

Matrix addColumn(Matrix _A, Matrix _B)
{
	Matrix matC = createMat(_A.rows, (_A.cols + _B.cols));

	if (_A.rows != _B.rows)
	{
		printf("The rows sizes of two matrice are not equal!!\n");
		return matC;
	}

	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols ; j++)
		{
			matC.at[i][j] = _A.at[i][j];
			matC.at[i][_A.cols + j] = _B.at[i][j];
		}
	}

	return matC;
}


Matrix addRow(Matrix _A, Matrix _B)
{
	Matrix matC = createMat(_A.rows + _B.rows, _A.cols);

	if (_A.cols != _B.cols)
	{
		printf("The column sizes of two matrice are not equal!!\n");
		return matC;
	}

	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols; j++)
		{
			matC.at[i][j] = _A.at[i][j];
			matC.at[_A.rows + i][j] = _B.at[i][j];
		}
	}

	return matC;
}

/*--------------------------------------------------------------------------*/
/*							Vector Arithmetic								*/
/*--------------------------------------------------------------------------*/

/* Add Vector (a + b)
	_a		: input vector #1
	_b		: input vector #2
	_out	: output vector (_a +_b)
	_vecLen	: length of input vector	*/
void addVec(Matrix _A, Matrix _B, Matrix _out) 
{
	for (int i = 0; i < _B.rows; i++) {
		if (_A.cols == 1 && _B.cols == 1 && _A.rows == _B.rows)
		_out.at[i][0] = _A.at[i][0] + _B.at[i][0];
		else
		{
			printf("_A & _B are not vectors!!\n");
			break;
		}
	}
}

/* Subtract Vector (a - b)
	_a		: input vector #1
	_b		: input vector #2
	_out	: output vector (_a -_b)
	_vecLen	: length of input vector	*/
void subVec(Matrix _A, Matrix _B, Matrix _out) 
{
	for (int i = 0; i < _B.rows; i++) {
		if (_A.cols == 1 && _B.cols == 1 && _A.rows == _B.rows)
			_out.at[i][0] = _A.at[i][0] - _B.at[i][0];
		else
		{
			printf("_A & _B are not vectors!!\n");
			break;
		}
	}
}

/* Dot product (a & b): element-wise multiplication of 2 vectors
	_a		: input vector #1
	_b		: input vector #2
	_out	: output number (_a *_b) 
	_vecLen	: length of input vector	*/
void dotProduct(Matrix _A, Matrix _B, double _out) 
{
	for (int i = 0; i < _B.rows; i++) {
		if (_A.cols == _B.cols && _A.rows == _B.rows)
			_out += _A.at[i][0] * _B.at[i][0];
		else
		{
			printf("_A & _B are not vectors!!\n");
			break;
		}
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
void addMat(Matrix _A, Matrix _B, Matrix _Out) {

	for (int i = 0; i < _A.rows; i++)
	{
		if ((_A.rows != _B.rows) || (_A.cols != _B.cols))
		{
			printf("The size of _A & _B are different!!\n");
			break;
		}

		for (int j = 0; j < _A.cols; j++) {
			_Out.at[i][j] = _A.at[i][j] + _B.at[i][j];
		}
	}
		
}

/* Subtract Matrix (A - B)
	_A		: input matrix #1
	_B		: input matrix #2
	_Out	: output matrix (_A -_B)
	_A_rows	: rows of input matrix _A
	_A_cols	: cols of input matrix _A	*/
void subMat(Matrix _A, Matrix _B, Matrix _Out) {
	for (int i = 0; i < _A.rows; i++)
	{
		if ((_A.rows != _B.rows) || (_A.cols != _B.cols))
		{
			printf("The size of _A & _B are different!!\n");
			break;
		}

		for (int j = 0; j < _A.cols; j++) {
			_Out.at[i][j] = _A.at[i][j] - _B.at[i][j];
		}
	}
}
/* Multiply Matrix (A & B) without return
	_A		: input matrix #1
	_B		: input matrix #2
	_Out	: output matrix (_A *_B)	*/
extern void multiMat(Matrix _A, Matrix _B, Matrix _Out)
{
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _B.cols; j++) {
			if (_A.cols == _B.rows)
				for (int k = 0; k < _A.cols; k++)
					_Out.at[i][j] += _A.at[i][k] * _B.at[k][j];
			else printf("Matrices' sizes are different!");
		}
	}
}

/* Multiply Matrix (A & B) with return
	_A		: input matrix #1
	_B		: input matrix #2	
	_Out	: output matrix (_A *_B)  */
extern Matrix multiMat(Matrix _A, Matrix _B)
{
	Matrix _Out = zeros(_A.rows, _B.cols);

	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _B.cols; j++) {
			if (_A.cols == _B.rows)
				for (int k = 0; k < _A.cols; k++)
					_Out.at[i][j] += _A.at[i][k] * _B.at[k][j];
			else printf("Matrices' sizes are different!\n");
		}
	}

	return _Out;
}

/* Print Matrix Elements
	_A		: input matrix
	_rows	: rows of input matrix
	_cols	: cols of input matrix	*/
//void printMat(double _A[][N], int _rows, int _cols) {
//	for (int i = 0; i < _rows; i++) {
//		for (int j = 0; j < _cols; j++) {
//			/*	 ADD YOUR LOGIC HERE (1 lines)	*/
//			printf("%lf\t", _A[i][j]);
//		}
//		printf("\n");
//	}
//	printf("\n");
//}



/* Deriving Inverse Matrix with allocation error
	_A		: input matrix _A
	_Out	: output matrix => inverse matrix of _A	
	*** Using [A | I] => [I | A^-1] Method ***	*/
void invMat(Matrix _A, Matrix _Out)
{
	int max_idx;
	double temp = 0, max = 0;

	Matrix row_temp = zeros(1, _A.cols);					// make temperary matrices for calculation and preserving original values
	Matrix m = createMat(_A.rows, 1);
	Matrix I = eye(_A.rows, _A.cols);
	Matrix A_aug = addColumn(_A, I);
	

	for (int k = 0; k < (A_aug.rows - 1); k++)
	{
		
		// upper elimination
		for (int i = k + 1; i < A_aug.rows; i++)
		{
			m.at[i-1][0] = (A_aug.at[i][k] / A_aug.at[k][k]);		// make m matrix for the efficient calculation

			if (_A.at[k][k] == 0)								// error case: divide by zero
			{
				printf("A[%d][%d] is zero that divides by zero!!\n", k, k);
				break;
			}
			else if (_A.rows != _A.cols)					// error case: not square matrix
			{
				printf("The matrix A is not square matrix!!\n");
				break;
			}
			for (int j = k; j < A_aug.cols; j++)					// gauss-elimination
			{
				A_aug.at[i][j] = A_aug.at[i][j] - m.at[i-1][0] * A_aug.at[k][j];
			}
		}
	}

	// lower elimination
	for (int k = A_aug.rows - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			m.at[i][0] = (A_aug.at[i][k] / A_aug.at[k][k]);		// make m matrix for the efficient calculation

			if (_A.at[k][k] == 0)								// error case: divide by zero
			{
				printf("A[%d][%d] is zero that divides by zero!!\n", k, k);
				break;
			}
			else if (_A.rows != _A.cols)					// error case: not square matrix
			{
				printf("The matrix A is not square matrix!!\n");
				break;
			}

			for (int j = k; j < A_aug.cols; j++)					// gauss-elimination
			{
				A_aug.at[i][j] = A_aug.at[i][j] - m.at[i][0] * A_aug.at[k][j];
			}
		}

	}

	for (int i = 0; i < _A.rows; i++)			// extracting inverse matrix from augmented matrices
	{
		m.at[i][0] = A_aug.at[i][i];
		for (int j = _A.cols; j < _A.cols + I.cols; j++)
		{
			A_aug.at[i][j] = A_aug.at[i][j] / m.at[i][0];
			_Out.at[i][j - _A.cols] = A_aug.at[i][j];
		}
	}

}
