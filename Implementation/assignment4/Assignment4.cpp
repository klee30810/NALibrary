/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Kyung-Min Lee
Created         : 01-04-2019
Modified        : 05-04-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 4
------------------------------------------------------------------------------------------*/

#define Assignment	4		// enter your assignment number
#define eval		0		// set 0

#include "myNM.h"

int main(int argc, char *argv[])
{
	/*	 [※ DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
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
	Matrix matQ = createMat(matA.rows, matA.cols);
	Matrix matB = createMat(3, 3);
	Matrix matR = createMat(matA.rows, matA.cols);
	
	Matrix valEig;
	Matrix vecEig;
	
	
	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	QRdecomp(matA, matQ, matR);
	valEig = eig(matA);
	vecEig = eigvec(matA); //4x4 행렬 이상에서도 가능

	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	printf("----------------------------------------------------------------------------------------------\n");
	printf("			            QR Decomposition Results					\n");
	printf("----------------------------------------------------------------------------------------------\n");
	printf("\n[matrix A]\n\n");
	printMat(matA);
	printf("\n[matrix Q]\n\n");
	printMat(matQ);
	printf("\n[matrix R]\n\n");
	printMat(matR);

	printf("\n\n");

	printf("----------------------------------------------------------------------------------------------\n");
	printf("			       Eigenvalue & Eigenvector Results					\n");
	printf("----------------------------------------------------------------------------------------------\n");
	printf("\n[Eigen value]\n\n");
	printMat(valEig);
	printf("\n[Eigen vector]\n\n");
	printMat(vecEig);
	
	printf("\n\n");
	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(matA);	freeMat(matQ);	freeMat(matR);
	freeMat(valEig);	freeMat(vecEig);

	system("pause");
	return 0;
}

