/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Kyung-Min Lee
Created         : 13-03-2019
Modified        : 26-03-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 2
/------------------------------------------------------------------------------------------*/

#define Assignment	2		// enter your assignment number
#define eval		0		// set 0

#include "myMatrix.h"
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
	Matrix matA = txt2Mat(path, "prob1_matA");
	Matrix vecb = txt2Mat(path, "prob1_vecb");
	Matrix matU = createMat(matA.rows, matA.cols);
	Matrix vecbn = createMat(vecb.rows, 1);
	Matrix vecx = createMat(vecb.rows, 1);
	

	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	gaussElim(matA, vecb, matU, vecbn);
	backsub(matU, vecbn, vecx);

	
	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	printf("------------------------------------------------------------------------------------\n");
	printf("			Gauss Elimination Results            \n");
	printf("------------------------------------------------------------------------------------\n");
	printf("\n[matrix A]\n\n");
	printMat(matA);
	printf("\n[vector b]\n\n");
	printMat(vecb);
	printf("\n[matrix U]\n\n");
	printMat(matU);
	printf("\n[vector bn]\n\n");
	printMat(vecbn);
	printf("\n[vector x]\n\n");
	printMat(vecx);
	printf("\n");
	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(matA);	freeMat(vecb);	freeMat(matU);	freeMat(vecbn);	freeMat(vecx);


	system("pause");
	return 0;
}
