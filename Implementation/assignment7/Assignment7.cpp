#pragma once
/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Kyung-Min Lee
Created         : 24-04-2019
Modified        : 17-05-2019
Language/ver	: C in MSVS2017
Course			: Numerical method 2019-Spring

Description     : Assignment 7
------------------------------------------------------------------------------------------*/

#define Assignment	7		// enter your assignment number
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
	
	// Problem 1
	double a1 = 0, b1 = PI;					// Integration limits
	int N1 = 10, N2 = 4;						// the number of intervals
	double epsilon1 = pow(10, -5);		// error tolerance
	int itr1 = 0, itrMax1 = 8;
	double lr_1 = 0;				// Analytic solution
	double lt_1 = 0, err1 = 0;		// trapezoidal method
	double ls_1 = 0, err2 = 0;		// simpson 1/3 method
	double la_1 = 0, err3 = 0;		// adaptive simpson method

	// Problem 2
	double a2 = 0, b2 = 2*PI;					// Integration limits
	int N3 = 10, N4 = 4;						// the number of intervals
	double epsilon2 = pow(10, -5);		// error tolerance
	int itr2 = 0, itrMax2 = 8;
	double lr_2 = 0;				// Analytic solution
	double lt_2 = 0, err4 = 0;		// trapezoidal method
	double ls_2 = 0, err5 = 0;		// simpson 1/3 method
	double la_2 = 0, err6 = 0;		// adaptive simpson method



	/*==========================================================================*/
	/*				Apply your numerical method algorithm of N = 10				*/
	/*==========================================================================*/
	
	// Problem 1
	lr_1 = PI / 2;
	lt_1 = trapezoidal(func1, a1, b1, N1);
	ls_1 = simpson13(func1, a1, b1, N1);
	la_1 = adaptSimpson(func1, a1, b1, epsilon1, N1, itr1, itrMax1);
	err1 = fabs(lt_1 - lr_1);
	err2 = fabs(ls_1 - lr_1);
	err3 = fabs(la_1 - lr_1);

	// Problem 2
	lr_2 = (double)1/5 - exp(-2*PI)/5;
	lt_2 = trapezoidal(func2, a2, b2, N3);
	ls_2 = simpson13(func2, a2, b2, N3);
	la_2 = adaptSimpson(func2, a2, b2, epsilon2, N3, itr2, itrMax2);
	err4 = fabs(lt_2 - lr_2);
	err5 = fabs(ls_2 - lr_2);
	err6 = fabs(la_2 - lr_2);


	/*==========================================================================*/
	/*					  Print your results of N = 10							*/
	/*==========================================================================*/
	printf("----------------------------------------------------------------------------------------------\n");
	printf("			           Integration Results	N = 10				\n");
	printf("----------------------------------------------------------------------------------------------\n");
	cout << "The problem 1 when N = 10" << endl;
	
	cout << "fx=sin(x)^2, [0,pi]\n" << endl;			
	cout << "Exact result of integration of fx : " << lr_1 << endl; // Analytic Solution
	cout << "Result from Trapezoidal Method: " << lt_1 << endl;
	cout << "Error from Trapezoidal Method: " << err1 << endl;
	cout << "Result from Simpson 1/3 Method: " << ls_1 << endl;
	cout << "Error from Simpson 1/3 Method: " << err2 << endl;
	cout << "Result from Adaptive Simpson Method: " << la_1 << endl;
	cout << "Error from Adaptive Simpson Method: " << err3 << endl;
	printf("\n\n");
	
	cout << "The problem 2 when N = 10" << endl;

	cout << "fx= cos(2x) * exp(-x), [0,2pi]\n" << endl;
	cout << "Exact result of integration of fx : " << lr_2 << endl; // Analytic Solution
	cout << "Result from Trapezoidal Method: " << lt_2 << endl;
	cout << "Error from Trapezoidal Method: " << err4 << endl;
	cout << "Result from Simpson 1/3 Method: " << ls_2 << endl;
	cout << "Error from Simpson 1/3 Method: " << err5 << endl;
	cout << "Result from Adaptive Simpson Method: " << la_2 << endl;
	cout << "Error from Adaptive Simpson Method: " << err6 << endl;
	printf("\n\n");
	

	/*==========================================================================*/
	/*				Apply your numerical method algorithm of N = 4				*/
	/*==========================================================================*/

	// Problem 1
	lr_1 = PI / 2;
	lt_1 = trapezoidal(func1, a1, b1, N2);
	ls_1 = simpson13(func1, a1, b1, N2);
	la_1 = adaptSimpson(func1, a1, b1, epsilon1, N2, itr1, itrMax1);
	err1 = fabs(lt_1 - lr_1);
	err2 = fabs(ls_1 - lr_1);
	err3 = fabs(la_1 - lr_1);

	// Problem 2
	lr_2 = (double)1 / 5 - exp(-2 * PI) / 5;
	lt_2 = trapezoidal(func2, a2, b2, N4);
	ls_2 = simpson13(func2, a2, b2, N4);
	la_2 = adaptSimpson(func2, a2, b2, epsilon2, N4, itr2, itrMax2);
	err4 = fabs(lt_2 - lr_2);
	err5 = fabs(ls_2 - lr_2);
	err6 = fabs(la_2 - lr_2);


	/*==========================================================================*/
	/*						  Print your results of N = 4						*/
	/*==========================================================================*/
	printf("----------------------------------------------------------------------------------------------\n");
	printf("			           Integration Results	N = 4				\n");
	printf("----------------------------------------------------------------------------------------------\n");
	cout << "The problem 1 when N = 4" << endl;

	cout << "fx=sin(x)^2, [0,pi]\n" << endl;
	cout << "Exact result of integration of fx : " << lr_1 << endl; // Analytic Solution
	cout << "Result from Trapezoidal Method: " << lt_1 << endl;
	cout << "Error from Trapezoidal Method: " << err1 << endl;
	cout << "Result from Simpson 1/3 Method: " << ls_1 << endl;
	cout << "Error from Simpson 1/3 Method: " << err2 << endl;
	cout << "Result from Adaptive Simpson Method: " << la_1 << endl;
	cout << "Error from Adaptive Simpson Method: " << err3 << endl;
	printf("\n\n");

	cout << "The problem 2 when N = 4" << endl;

	cout << "fx= cos(2x) * exp(-x), [0,2pi]\n" << endl;
	cout << "Exact result of integration of fx : " << lr_2 << endl; // Analytic Solution
	cout << "Result from Trapezoidal Method: " << lt_2 << endl;
	cout << "Error from Trapezoidal Method: " << err4 << endl;
	cout << "Result from Simpson 1/3 Method: " << ls_2 << endl;
	cout << "Error from Simpson 1/3 Method: " << err5 << endl;
	cout << "Result from Adaptive Simpson Method: " << la_2 << endl;
	cout << "Error from Adaptive Simpson Method: " << err6 << endl;
	printf("\n\n");


	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/


	system("pause");
	return 0;
}
