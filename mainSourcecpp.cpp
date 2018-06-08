#include <iostream>

using namespace std;

void differentialQuotient(double**R, double *tabx, double *taby, int k)
{
	R[0] = taby;
	
	for (int j = 1; j < k; j++)
	{
		R[j] = new double[k - j];

		for (int i = 0; i < (k - j); i++)
		{
			double nominator = (R[j-1][i+1]-R[j-1][i]);
			double denominator=(tabx[i+j]-tabx[i]);
			R[j][i] = nominator / denominator;
		}
	}
}

void showQuotientList(double **R, int k)
{
	for (int j = 1; j < k; j++)
	{
		cout << j << endl;
		for (int i = 0; i < (k - j); i++)
			cout << R[j][i] << "\t";
		cout << endl;
	}
}

double newtonPolynominal(double *tabx, int k, double x)
{
	double p = 1;
	if (k == 0)
		return p;
	for (int i = 0; i < k; i++)
		p *=(x - tabx[i]);
	return p;
}

double newtonInterpolation(double *tabx, double *taby, double**R, int k, double x)
{
	double output = taby[0];
	if (k == 1)
		return output;
	double p;
	for (int j = 1; j < k; j++)
	{
		p = newtonPolynominal(tabx, j, x);
		output += R[j][0] * p;
	}
	return output;
}

double lagrangeInterpolation(double *tabx, double *taby, int k, double x)
{
	double output = 0;
	double l;
	for (int i = 0; i < k; i++)
	{
		l = 1;
		for (int j = 0; j < k; j++)
		{
			if (j != i)
				l *= (x - tabx[j]) / (tabx[i] - tabx[j]);
		}
		output += l*taby[i];
	}
	return output;
}

int main()
{
	int k1 = 5;
	double tabx1[] = { 0,2,3,5,6 };
	double taby1[] = { 0,8,27,125,216 };
	double **R1 = new double*[k1];

	int k2 = 3;
	double tabx2[] = { 0,2,4 };
	double taby2[] = { 1,5,17 };
	double **R2 = new double*[k2];

	
	//===============differential quotient function test:===================================
	differentialQuotient(R1, tabx1, taby1, k1);
	differentialQuotient(R2, tabx2, taby2, k2);

	cout << "quotient function test:" << endl;
	cout << "tab1 - quotient list: " << endl;
	showQuotientList(R1, k1);
	cout << "tab2 - quotient list: " << endl;
	showQuotientList(R2, k2);

	cout << "finished" << endl;

	//================Newton Interpolation test:=============================================
	double iN27 = newtonInterpolation(tabx1, taby1, R1, k1, 3);
	cout << "Newton Interpolation - the resut should equals 27: " << iN27 << endl;

	double iN216 = newtonInterpolation(tabx1, taby1, R1, k1, 6);
	cout << "Newton Interpolation - result should equals 216: " << iN216 << endl;

	double iN17 = newtonInterpolation(tabx2, taby2, R2, k2, 4);
	cout << "Newton Interpolation - result should equals 17: " << iN17 << endl;

	//==============Lagrange Interpolation test:==============================================
	double iL27 = lagrangeInterpolation(tabx1, taby1,k1, 3);
	cout << "Lagrange Interpolation - the resut should equals 27: " << iL27 << endl;
	cout << "error estimation: " << endl;

	double iL216 = lagrangeInterpolation(tabx1, taby1, k1, 6);
	cout << "Lagrange Interpolation - result should equals 216: " << iN216 << endl;
	cout << "error estimation: " << endl;

	double iL17 = lagrangeInterpolation(tabx2, taby2, k2, 4);
	cout << "Lagrange Interpolation - result should equals 17: " << iN17 << endl;

	system("PAUSE");
	return 0;
}