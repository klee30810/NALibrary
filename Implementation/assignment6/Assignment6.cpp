#pragma once
/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 24-04-2019
Modified        : 24-04-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 6
------------------------------------------------------------------------------------------*/

#define Assignment	6		// enter your assignment number
#define eval		0		// set 0

#include "myNM.h"
#include "myMatrix.h"

int main(int argc, char *argv[])
{
	/*	 [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
	string path = "C:/NM_data_2019/Assignment" + to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif

	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified file name						*/
	/*	   : For each assignment, the file name will be notified on HISNET		*/
	/*==========================================================================*/
	int n = 4;  

	// Problem 1
	const char* method1 = "3f";
	Matrix vecx1 = txt2Mat(path, "prob1_vecx");
	Matrix vecy1 = txt2Mat(path, "prob1_vecy");
	Matrix error1_1 = zeros(vecx1.rows - 1, 1);
	Matrix error2_1 = zeros(vecx1.rows - 2, 1);
	Matrix diff_11 = zeros(vecy1.rows, 1);
	Matrix diff_21 = zeros(vecy1.rows, 1);

	// Problem 2
	const char* method2 = "2f";
	Matrix vecx2 = zeros(2, 1);
	Matrix vecy2 = zeros(2, 1);
	Matrix error1_2 = zeros(1, 1);
	Matrix error2_2 = zeros(2, 1);
	Matrix diff_12 = zeros(1, 1);
	Matrix diff_22 = zeros(1, 1);
	


	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	
	// Problem 1
	deriv1(diff_11, vecx1, vecy1, error1_1, (vecx1.at[1][0] - vecx1.at[0][0]) , method1);
	deriv2(diff_21, vecx1, vecy1, error2_1, (vecx1.at[1][0] - vecx1.at[0][0]), method1);

	// Problem 2
	double h = 0.1;
	double error = 0;
	vecx2.at[0][0] = PI / 2;					// input in radian
	vecx2.at[1][0] = vecx2.at[0][0] + h;
	vecy2.at[0][0] = cos(vecx2.at[0][0]);		// cosine value in vecx
	vecy2.at[1][0] = cos(vecx2.at[1][0]);		// cosine value in vecx
	deriv1(diff_12, vecx2, vecy2, error1_2, h, method2);


	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	printf("----------------------------------------------------------------------------------------------\n");
	printf("			            Differentiation Results					\n");
	printf("----------------------------------------------------------------------------------------------\n");
	cout << "The problem 1" << endl;
	printf("\n[vector x]\n\n");
	printMat(vecx1); printf("\n");
	printf("\n[vector y]\n\n");
	printMat(vecy1); printf("\n");
	printf("\n[1st orderd differentiation between the adjacents points with three points forward method]\n\n");
	printMat(diff_11); printf("\n");
	printf("\n[2nd orderd differentiation between the adjacents points with three points forward method]\n\n");
	printMat(diff_21); printf("\n");
	printf("\n\n");


	cout << "The problem 2" << endl;
	cout << "\n[cos(PI/2)]: " << cos(PI/2)<< endl;
	cout << "\n[1st derivation of cos(PI/2) when h = 0.1]\n" << endl;
	printMat(diff_12); printf("\n");
	cout << "Error when h = 0.1: " << fabs((diff_12.at[0][0] + sin(PI / 2)) / sin(PI / 2)) << endl;	cout << "\n[1st derivation of cos(PI/2) when h = 0.01]\n" << endl;
	h = 0.01;
	vecx2.at[0][0] = PI / 2;					// input in radian
	vecx2.at[1][0] = vecx2.at[0][0] + h;
	vecy2.at[0][0] = cos(vecx2.at[0][0]);		// cosine value in vecx
	vecy2.at[1][0] = cos(vecx2.at[1][0]);
	deriv1(diff_12, vecx2, vecy2, error1_2, h, method2);
	printMat(diff_12); printf("\n");
	cout << "Error when h = 0.01: " << fabs((diff_12.at[0][0] + sin(PI/2))/sin(PI/2)) << endl;
	cout << "\n[1st derivation of cos(PI/2) when h = 0.001]\n" << endl;
	h = 0.001;
	vecx2.at[0][0] = PI / 2;					// input in radian
	vecx2.at[1][0] = vecx2.at[0][0] + h;
	vecy2.at[0][0] = cos(vecx2.at[0][0]);		// cosine value in vecx
	vecy2.at[1][0] = cos(vecx2.at[1][0]);
	deriv1(diff_12, vecx2, vecy2, error1_2, h, method2);
	printMat(diff_12); printf("\n");
	cout << "Error when h = 0.001: " << fabs((diff_12.at[0][0] + sin(PI / 2)) / sin(PI / 2)) << endl;
	
	printf("\n\n");
	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(vecx1); freeMat(vecy1);	
	freeMat(diff_11); freeMat(diff_21);


	system("pause");
	return 0;
}
