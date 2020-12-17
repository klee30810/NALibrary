/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Kyung-Min Lee
Created         : 08-03-2019
Modified        : 15-03-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 1
/------------------------------------------------------------------------------------------*/

#include "myNM.h"					// Including header file

void main() {

	/*==========================================================================*/
	/*					Assignment 1 - Bisection Method							*/
	/*==========================================================================*/

	/************		Variables declaration & initialization		************/
	float tol = 0.00001;			// the tolerance 
	//float a = 0.3;				// initial left value for bisection method, checking 1st function
	//float a = 0.6;				// initial right value for bisection method, checking 1st function
	float a = 0;					// initial left value for bisection method, checking 2nd function
	float b = 0.5;					// initial right value for bisection method, checking 2nd function
	double BM_result;

	/************		Test NM Functions & Show Output				************/
	printf("------------------------------------------------------------------------------------\n");
	printf("			Bisection Method Results             \n");
	printf("------------------------------------------------------------------------------------\n");

	printf("Bisection Method:\n");
	BM_result = bisectionNL(a, b, tol);					// Bisection Method
	
	printf("Final Solution for Bisection Method: %f \t", BM_result);
	printf("\n");

	/*==========================================================================*/
	/*					Assignment 1 - Newton-Raphson Method					*/
	/*==========================================================================*/

	/************		Variables declaration & initialization		************/
	
	float x0 = 1.4;					// Initial value for Newton-Raphson Method
	double NR_result, HB_result;

	/************		Test NM Functions & Show Output				************/
	printf("------------------------------------------------------------------------------------\n");
	printf("			Newton-Raphson Method Results             \n");
	printf("------------------------------------------------------------------------------------\n");

	printf("Newton-Raphson Method Result:\n");
	NR_result = newtonRaphson(x0, tol);					// Newton-Raphson Method

	printf("Final Solution for Newton-Raphson Method: %f \t", NR_result);
	printf("\n");

	/************		Test NM Functions & Show Output				************/
	printf("------------------------------------------------------------------------------------\n");
	printf("			Hybrid Method Results             \n");
	printf("------------------------------------------------------------------------------------\n");
	
	printf("Hybrid Method Result:\n");
	HB_result = NR_hybrid(a, b, x0, tol);				// Hybrid Method

	printf("Final Solution for NR_hybrid Method: %lf \t", HB_result);
	printf("\n");

	system("pause");




}
