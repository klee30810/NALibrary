#pragma once
/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 24-04-2019
Modified        : 24-04-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 5
------------------------------------------------------------------------------------------*/

#define Assignment	5		// enter your assignment number
#define eval		0		// set 0

#include "myNM.h"

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
	double valSSE1 = 0, valSSE2 = 0;
	Matrix vecx1 = txt2Mat(path, "prob1_vecx");
	Matrix vecy1 = txt2Mat(path, "prob1_vecy");
	Matrix vecx2 = txt2Mat(path, "prob2_vecx");
	Matrix vecy2 = txt2Mat(path, "prob2_vecy");
	Matrix vecz1 = createMat(n + 1, 1);
	Matrix vecz2 = createMat(n + 1, 1);


	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	polyfit(vecx1, vecy1, n, vecz1, valSSE1);
	polyfit(vecx2, vecy2, n, vecz2, valSSE2);

	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	printf("----------------------------------------------------------------------------------------------\n");
	printf("			            Curve Fitting Results					\n");
	printf("----------------------------------------------------------------------------------------------\n");
	cout << "The problem 1" << endl;
	printf("\n[vector x]\n\n");
	printMat(vecx1); printf("\n");
	printf("\n[vector y]\n\n");
	printMat(vecy1); printf("\n");
	printf("\n[vector z]\n\n");
	printMat(vecz1); printf("\n");
	

	for (int i = n; i >= 0; i--)
	{
		if (i == n)
			printf("\n[f(x)]\n\n%f*x^%d", vecz1.at[n][0], n);
		else if (i > 0)
		{
			if (vecz1.at[i][0] >= 0)
				printf("+");
			printf("%f*x^%d", vecz1.at[i][0], i);
		}
		else
		{
			if (vecz1.at[i][0] >= 0)
				printf("+");
			printf("%f", vecz1.at[i][0]);
		}
	}

	printf("\n\n");
	printf("\n[Sum of Square Error]\n\n"); printf("%f", valSSE1); printf("\n\n");
	printf("\n");

	cout << "The problem 2" << endl;
	printf("\n[vector x]\n\n");
	printMat(vecx2); printf("\n");
	printf("\n[vector y]\n\n");
	printMat(vecy2); printf("\n");
	printf("\n[vector z]\n\n");
	printMat(vecz2); printf("\n");


	for (int i = n; i >= 0; i--)
	{
		if (i == n)
			printf("\n[f(x)]\n\n%f*x^%d", vecz2.at[n][0], n);
		else if (i > 0)
		{
			if (vecz2.at[i][0] >= 0)
				printf("+");
			printf("%f*x^%d", vecz2.at[i][0], i);
		}
		else
		{
			if (vecz2.at[i][0] >= 0)
				printf("+");
			printf("%f", vecz2.at[i][0]);
		}
	}

	printf("\n\n");
	printf("\n[Sum of Square Error]\n\n"); printf("%f", valSSE2); printf("\n\n");
	printf("\n");


	printf("\n\n");
	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(vecx1);	freeMat(vecx2);
	freeMat(vecy1);	freeMat(vecy2);
	freeMat(vecz1); freeMat(vecz2);


	system("pause");
	return 0;
}
