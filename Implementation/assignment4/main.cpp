/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Kyung-Min Lee
Created         : 25-03-2019
Modified        : 31-03-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 3
/------------------------------------------------------------------------------------------*/

#define Assignment	3		// enter your assignment number
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
	Matrix matA = txt2Mat(path, "prob3_matA");
	Matrix vecb = txt2Mat(path, "prob3_vecb");
	//Matrix matAinv = zeros(matA.rows, matA.cols);
	//invMat(matA, matAinv);

	Matrix matU_G = createMat(matA.rows, matA.cols);
	Matrix vecbn = createMat(vecb.rows, 1);
	Matrix matP_G = createMat(matA.rows, matA.cols);
	Matrix vecx_G = createMat(vecb.rows, 1);

	Matrix matL = createMat(matA.rows, matA.cols);
	Matrix matU_LU = createMat(matA.rows, matA.cols);
	Matrix matP_LU = createMat(matA.rows, matA.cols);
	Matrix vecx_LU = createMat(vecb.rows, 1);

	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	gaussElim(matA, vecb, matU_G, vecbn, matP_G);
	LUdecomp(matA, matL, matU_LU, matP_LU);
	solveLinear(matA, vecb, "g", vecx_G);
	solveLinear(matA, vecb, "l", vecx_LU);
	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	printf("----------------------------------------------------------------------------------------------\n");
	printf("			Gauss Elimination with scaled partial pivoting Results					\n");
	printf("----------------------------------------------------------------------------------------------\n");
	printf("\n[matrix A]\n\n");
	printMat(matA);
	printf("\n[vector b]\n\n");
	printMat(vecb);
	printf("\n[matrix U]\n\n");
	printMat(matU_G);
	printf("\n[matrix P]\n\n");
	printMat(matP_G);

	printf("----------------------------------------------------------------------------------------------\n");
	printf("			LU Decomposition with scaled partail pivoting Results					\n");
	printf("----------------------------------------------------------------------------------------------\n");
	printf("\n[matrix A]\n\n");
	printMat(matA);
	printf("\n[matrix L]\n\n");
	printMat(matL);
	printf("\n[matrix U]\n\n");
	printMat(matU_LU);
	printf("\n[matrix P]\n\n");
	printMat(matP_LU);

	printf("------------------------------------------------------------------------------------\n");
	printf("			Solving Linear Equation					\n");
	printf("------------------------------------------------------------------------------------\n");

	printf("{Using Gauss Elimination}:\n");
	printf("\n[matrix A]\n\n");
	printMat(matA);
	printf("\n[vector b]\n\n");
	printMat(vecb);
	printf("\n[vector x]\n\n");
	printMat(vecx_G);
	printf("\n");
	printf("{Using LU Decomposition}:\n");
	printf("\n[matrix A]\n\n");
	printMat(matA);
	printf("\n[vector b]\n\n");
	printMat(vecb);
	printf("\n[vector x]\n\n");
	printMat(vecx_LU);

	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(matA);
	freeMat(vecb);

	freeMat(matU_G);
	freeMat(vecbn);
	freeMat(matP_G);
	freeMat(vecx_G);

	freeMat(matL);
	freeMat(matU_LU);
	freeMat(matP_LU);
	freeMat(vecx_LU);

	system("pause");
	return 0;
}
