#include "myNM.h"
#include "myMatrix.h"


int main(int argc, char *argv[])
{
	
	Matrix x = zeros(4, 1);
	Matrix y = zeros(4, 1);

	// Input data
	x.at[0][0] = 1; x.at[1][0] = 1.5; x.at[2][0] = 2; x.at[3][0] = 2.5;
	y.at[0][0] = 3.8; y.at[1][0] = 4.9; y.at[2][0] = 5.2; y.at[3][0] = 5.4;


	// Question 1
	printf("----------------------------------------------------------------------------------------------\n");
	printf("			           Question 1										\n");
	printf("----------------------------------------------------------------------------------------------\n");
	cout << "The problem 1-a" << endl;
	Matrix Pcoef = zeros(3, 1);
	double RMSE = 0;
	polyfit(x, y, 2, Pcoef, RMSE);
	
	cout << "	The values of a2, a1, a0: " << Pcoef.at[0][0] << ',' << Pcoef.at[1][0] << ',' << Pcoef.at[2][0] << endl;
	printf("\n");
	cout << "	The fitted curve: "<< Pcoef.at[0][0] << "x^2 + " << Pcoef.at[1][0] << "x + " << Pcoef.at[2][0] << endl;
	printf("\n");

	cout << "The problem 1-b" << endl;
	cout << "	RMSE value: " << RMSE << endl;
	printf("\n");

	cout << "The problem 1-c" << endl;
	double integration = 0;
	int itr = 0;
	integration = adaptSimpson(_1cfunc, 1.25, 25.25,0.00001, 10, itr, 8);
	printf("	The area underneath the curve: %.6lf", integration);
	printf("\n"); printf("\n");

	// Question 2
	printf("----------------------------------------------------------------------------------------------\n");
	printf("			           Question 2										\n");
	printf("----------------------------------------------------------------------------------------------\n");
	cout << "The problem 2-a" << endl;
	Matrix matA = zeros(2, 2);
	matA.at[0][0] = 0; matA.at[0][1] = 1; matA.at[1][0] = -3; matA.at[1][1] = -20;
	Matrix eigval = zeros(2, 1);
	eigval = eig(matA);
	cout << "	The 1st eigenvalue of the given system: " << eigval.at[0][0] << endl;
	cout << "	The 2nd eigenvalue of the given system: " << eigval.at[1][0] << endl;
	printf("\n"); printf("\n");
	
	cout << "The problem 2-b" << endl;
	printf("\n");
	double a1 = 0, b1 = 0.2, a2 = 1.8, b2 = 2, h = 0.01, xinit = 0.05, vinit = 0.2; // all given values
	Matrix time1 = zeros((int)((b1 - a1)/h), 1);
	Matrix time2 = zeros((int)((b2 - a2) / h), 1);
	Matrix displacement1 = zeros((int)((b1 - a1) / h), 1);
	Matrix displacement2 = zeros((int)((b2 - a2) / h), 1);
	Matrix velocity1 = zeros((int)((b1 - a1) / h), 1);
	Matrix velocity2 = zeros((int)((b2 - a2) / h), 1);

	ode2RK3(mckdydt, mckdvdt, a1, b1, h, xinit, vinit, time1, displacement1, velocity1);
	ode2RK3(mckdydt, mckdvdt, a2, b2, h, xinit, vinit, time2, displacement2, velocity2);
	
	cout << "The x(t) from t=0 to t=0.2" << endl;
	printMat(displacement1); printf("\n");
	cout << "The x(t) from t=1.8 to t=2" << endl;
	printMat(displacement2); printf("\n");
	printf("\n"); printf("\n");
	
	cout << "The problem 2-c" << endl;
	printf("\n");
	cout << "The v(t) from t=0 to t=0.2 from ODE" << endl;
	printMat(velocity1); printf("\n");
	cout << "The v(t) from t=1.8 to t=2 from ODE" << endl;
	printMat(velocity2); printf("\n");
	
	Matrix velocity1_2 = derivative(displacement1, h);
	Matrix velocity2_2 = derivative(displacement2, h);

	cout << "The v(t) from t=0 to t=0.2 from derivative()" << endl;
	printMat(velocity1_2); printf("\n");
	cout << "The v(t) from t=1.8 to t=2 from derivative()" << endl;
	printMat(velocity2_2); printf("\n");
	
	
	system("pause");
	return 0;
}