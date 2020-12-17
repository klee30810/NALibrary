/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Kyung-Min Lee
Created         : 26-02-2019
Modified        : 26-02-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assigment0
/------------------------------------------------------------------------------------------*/

#include "myNM.h"
#include <iostream>

void main() {

	/*==========================================================================*/
	/*					Assignment 0 - Part 1									*/
	/*==========================================================================*/

	/************		Variables declaration & initialization		************/
	int vecLen = 5;
	double a[5] = { 10, 20, 30, 40, 50 };
	double b[5] = { 5, 10, 15, 20, 25 };
	double out[5] = { 0 };
	double dotVal = 0;

	/************		Test NM Functions & Show Output				************/
	printf("¦£----------------------------------¦¤\n");
	printf("¦¢           Part1 Results          ¦¢\n");
	printf("¦¦----------------------------------¦¥\n");
	printf("input vector a =\n");
	printVec(a, vecLen);

	printf("input vector b =\n");
	printVec(b, vecLen);

	addVec(a, b, out, vecLen);
	printf("\naddVec result\t (a + b) =\n");
	printVec(out, vecLen);

	subVec(a, b, out, vecLen);
	printf("subVec result\t (a - b) =\n");
	printVec(out, vecLen);

	dotProduct(a, b, &dotVal, vecLen);
	printf("dotProduct result(a * b) =\n%lf\n\n", dotVal);


	/*==========================================================================*/
	/*					Assignment 0 - Part 2									*/
	/*==========================================================================*/


	/************		Variables declaration & initialization		************/
	int A_rows = 4, A_cols = 3, B_rows = 4, B_cols = 3, C_rows = 3, C_cols = 3;
	double A[4][3] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120 };
	double B[4][3] = { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
	double C[3][3] = { 1, 0, 1, 0, 1, 0, 1, 0, 1 };
	double Out[4][3] = { 0 };


	/************		Test NM Functions & Show Output				************/
	printf("¦£----------------------------------¦¤\n");
	printf("¦¢           Part2 Results          ¦¢\n");
	printf("¦¦----------------------------------¦¥\n");

	printf("Input Matrix A =\n");
	printMat(A, A_rows, A_cols);
	printf("Input Matrix B =\n");
	printMat(B, B_rows, B_cols);
	printf("Input Matrix C =\n");
	printMat(C, C_rows, C_cols);

	addMat(A, B, Out, A_rows, A_cols);
	printf("\naddMat result\t(A + B)\t=\n");
	printMat(Out, A_rows, A_cols);

	subMat(A, B, Out, A_rows, A_cols);
	printf("subMat result\t(A - B)\t=\n");
	printMat(Out, A_rows, A_cols);

	multMat(A, C, Out, A_rows, A_cols, C_rows, C_cols);
	printf("multMat result\t(A * C)\t=\n");
	printMat(Out, A_rows, C_cols);

	system("pause");
}
