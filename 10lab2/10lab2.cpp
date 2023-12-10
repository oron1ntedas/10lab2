#include <iostream>
#include <iomanip>

using namespace std;

double E1 = 1e-9;
double E2 = 1e-9;

double f1(double x1, double x2)
{
	return sin( x1 +1) - x2 - 1;
}
double f2(double x1, double x2)
{
	return 2*x1+cos(x2)-2;
}
double df1dx1(double x1, double x2)
{
	return cos(x1 +1);
}
double df1dx2(double x1, double x2)
{
	return (-1);
}
double df2dx1(double x1, double x2)
{
	return 2;
}
double df2dx2(double x1, double x2)
{
	return (-sin(x2));
}

double* Gauss(double** Array_A, double* Array_B, int n);
void Jcobian(double x1, double x2, double** J)
{
	J[0][0] = df1dx1(x1, x2);
	J[0][1] = df1dx2(x1, x2);
	J[1][0] = df2dx1(x1, x2);
	J[1][1] = df2dx2(x1, x2);
}
void Jcobian_M(double x1, double x2, double** J, double M)
{
	J[0][0] = (f1(x1 + M * x1, x2) - f1(x1, x2)) / (M * x1);
	J[0][1] = (f1(x1, x2 + M * x2) - f1(x1, x2)) / (M * x2);
	J[1][0] = (f2(x1 + M * x1, x2) - f2(x1, x2)) / (M * x1);
	J[1][1] = (f2(x1, x2 + M * x2) - f2(x1, x2)) / (M * x2);
}
double* Newton(double x1, double x2, int NIT, double E1, double E2);
double* Newton_M(double x1, double x2, int NIT, double E1, double E2, double M);

int main()
{
	setlocale(LC_ALL, "rus");
	int NIT = 20;
	cout << "Аналитический метод для точки (1,1): \n";
	double* x = new double[2];

	x = Newton(1, 1, NIT, E1, E2);
	cout << "\nЧерез матрицу Якоби для точки (1,1) и М=0.01: \n";
	double* x3 = new double[2], * x4 = new double[2], * x5 = new double[2];
	x3 = Newton_M(1, 1, NIT, E1, E2, 0.01);
	cout << "\nЧерез матрицу Якоби для точки (1,1) и М=0.05: \n";
	x4 = Newton_M(1, 1, NIT, E1, E2, 0.05);
	cout << "\nЧерез матрицу Якоби для точки (1,1) и М=0.1: \n";
	x5 = Newton_M(1, 1, NIT, E1, E2, 0.1);
	system("pause");
	return 0;
}

double* Gauss(double** Array_A, double* Array_B, int n)
{
	double* X = new double[n];
	for (int k = 0; k < n; k++)
	{
		for (int i = k + 1; i < n; i++)
		{
			if (abs(Array_A[i][k]) > abs(Array_A[k][k]))
			{
				swap(Array_A[i], Array_A[k]); 
				swap(Array_B[i], Array_B[k]);
			}
		}
		double A_Main = Array_A[k][k];
		if (A_Main == 0)
		{
			cout << "error\n";
			system("pause");
			exit(0);
		}
		for (int i = k; i < n; i++)
		{
			Array_A[k][i] /= A_Main;
		}
		Array_B[k] /= A_Main;
		for (int i = k + 1; i < n; i++)
		{
			double s = Array_A[i][k];
			for (int j = k; j < n; j++)
			{
				Array_A[i][j] -= s * Array_A[k][j];
			}
			Array_B[i] -= s * Array_B[k];
		}
	}
	for (int k = n - 1; k >= 0; k--)
	{
		X[k] = Array_B[k];
		for (int i = n - 1; i > k; i--)
		{
			X[k] -= Array_A[k][i] * X[i];
		}
	}
	return X;
}

double* Newton(double x1, double x2, int NIT, double E1, double E2)
{
	int k = 1;
	cout << left << setw(4) << "k" << setw(20) << "d1" << setw(20) << "d2" << setw(20) << "x1" << setw(20) << "x2" << endl;
	double* F = new double[2];
	double** J = new double* [2];
	for (int i = 0; i < 2; i++)
	{
		J[i] = new double[2];
	}
	double* dX = new double[2];
	double x1k, x2k;
	double d1, d2;
	double tmp;
	do {
		F[0] = -f1(x1, x2);
		F[1] = -f2(x1, x2);
		Jcobian(x1, x2, J);

		dX = Gauss(J, F, 2);
		x1k = x1 + dX[0];
		x2k = x2 + dX[1];
		d1 = abs(f1(x1, x2));
		tmp = abs(f2(x1, x2));
		if (tmp > d1)
		{
			d1 = tmp;
		}
		d2 = abs(x1k - x1) / (x1k >= 1 ? x1k : 1);
		tmp = abs(x2k - x2) / (x2k >= 1 ? x2k : 1);
		if (tmp > d2)
		{
			d2 = tmp;
		}
		x1 = x1k;
		x2 = x2k;
		cout << left << setw(4) << k << setw(20) << d1 << setw(20) << d2 << setw(20) << x1 << setw(20) << x2 << endl;
		if (k >= NIT)
		{
			cout << "IER=2\n";
			system("pause");
			exit(2);
		}
		k++;
	} while (d1 > E1 && d2 > E2);
	dX[0] = x1;
	dX[1] = x2;
	return dX;
}

double* Newton_M(double x1, double x2, int NIT, double E1, double E2, double M)
{
	int k = 1;
	cout << left << setw(4) << "k" << setw(20) << "d1" << setw(20) << "d2" << setw(20) << "x1" << setw(20) << "x2" << endl;
	double* F = new double[2];
	double** J = new double* [2];
	for (int i = 0; i < 2; i++)
	{
		J[i] = new double[2];
	}
	double* dX = new double[2];
	double x1k, x2k;
	double d1, d2;
	double tmp;
	do {
		F[0] = -f1(x1, x2);
		F[1] = -f2(x1, x2);
		Jcobian_M(x1, x2, J, M);

		dX = Gauss(J, F, 2);
		x1k = x1 + dX[0];
		x2k = x2 + dX[1];
		d1 = abs(f1(x1, x2));
		tmp = abs(f2(x1, x2));
		if (tmp > d1)
		{
			d1 = tmp;
		}
		d2 = abs(x1k - x1) / (x1k >= 1 ? x1k : 1);
		tmp = abs(x2k - x2) / (x2k >= 1 ? x2k : 1);
		if (tmp > d2)
		{
			d2 = tmp;
		}
		x1 = x1k;
		x2 = x2k;
		cout << left << setw(4) << k << setw(20) << d1 << setw(20) << d2 << setw(20) << x1 << setw(20) << x2 << endl;
		if (k >= NIT)
		{
			cout << "IER=2\n";
			exit(2);
		}
		k++;
	} while (d1 > E1 && d2 > E2);
	dX[0] = x1;
	dX[1] = x2;
	return dX;
}
