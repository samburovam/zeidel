#include "numerical.h"
#include <math.h>
#include <iomanip>
#include <stdio.h> 
#include <iostream>
using namespace std;

double Uxy(double x, double y)
{

	return x * x * x + y * y + 3;
}

double M1(double y)
{
	return (y * y + 3);
}

double M2(double y)
{
	return (y * y + 4);
}

double M3(double x)
{
	return (x * x * x + 3);
}

double M4(double x)
{
	return (x * x * x + 4);
}

double f(double x, double y)
{
	return (6 * x + 2);
}



double** CreateMatrix(int n, int m)
{
	double** Matrix = NULL;

	Matrix = new double* [n];
	for (int i = 0; i < n; i++)
		Matrix[i] = new double[m];

	return Matrix;
}



void ShowSolution(double** V, int n, int m, int s)
{
	cout << "iteration number  " << s << endl;;
	for (int j = m; j >= 0; j--)
	{
		for (int i = 0; i <= m; i++) 
		{
			printf("%5.3lf   ", V[i][j]);
		}
		printf("\n");
	}
}

double** FillRightSide(double** F, int n, int m, double a, double c, double h, double k)
{
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			double xi, yj, sum = 0;
			xi = a + i * h;
			yj = c + j * k;

			if (j == 1)
				sum += (1 / (k * k)) * M3(xi);
			else
				if (j == m - 1)
					sum += (1 / (k * k)) * M4(xi);
			if (i == 1)
				sum += (1 / (h * h)) * M1(yj);
			else
				if (i == n - 1)
					sum += (1 / (h * h)) * M2(yj);

			F[i][j] = f(xi, yj) - sum;
		}
	return F;
}

void FillStartSolution(double** V, int n, int m, double a, double b, double c, double d)
{
	double h, k;
	h = (b - a) / n;
	k = (d - c) / m;

	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
		{
			if (i == 0 || j == 0 || i == n || j == m)
			{
				double xi, yj, sum = 0;
				xi = a + i * h;
				yj = c + j * k;
				if (j == 0)
					V[i][j] = M3(xi);
				else
					V[i][j] = M4(xi);
				if (i == 0)
					V[i][j] = M1(yj);
				else
					if (i == n)
						V[i][j] = M2(yj);
			}
			else
				V[i][j] = 0;
		}
}

double ZeidelsMethod(double** V, int n, int m, double a, double b, double c, double d, double eps, int Nmax, double& epsMax, int& S)
{
	double	epsCur = 0;
	double	a2, k2, h2;
	double	v_old;
	double	v_new;
	double h;
	double k;

	h = (b - a) / n;
	k = (c - d) / m;

	h2 = 1 / (h * h);
	k2 = 1 / (k * k);
	a2 = -2 * (h2 + k2);

	while (true)
	{
		epsMax = 0;
		for (int j = 1; j < m; j++)
			for (int i = 1; i < n; i++)
			{
				double xi, yj;
				xi = a + i * h;
				yj = c + j * k;

				v_old = V[i][j];
				v_new = -(h2 * (V[i + 1][j] + V[i - 1][j]) + k2 * (V[i][j + 1] + V[i][j - 1]));
				v_new = v_new + f(xi, yj);
				v_new = v_new / a2;

				epsCur = fabs(v_old - v_new);
				if (epsCur > epsMax)
					epsMax = epsCur;

				V[i][j] = v_new;
			}
		if (S == 1 || S == 2)
		{
			ShowSolution(V, n, m, S);
		}
		++S;

		if ((epsMax < eps) || (S >= Nmax))
			break;
	}
	ShowSolution(V, n, m, S);
	return epsMax;
}

double DiscrepancyOfSolution(double** V, int n, int m, double a, double b, double c, double d)
{
	double	a2, k2, h2;
	double  h, k;
	double** F;
	double rs = 0;

	h = (b - a) / n;
	k = (d - c) / m;

	h2 = 1 / (h * h);
	k2 = 1 / (k * k);
	a2 = -2 * (h2 + k2);

	F = CreateMatrix(n + 1, m + 1);
	double** D = FillRightSide(F, n, m, a, c, h, k);

	for (int j = 1; j < m; j++)
	{
		for (int i = 1; i < n; i++)
		{
			double r = 0;
			double value = 0;
			r = a2 * V[i][j] + h2 * V[i - 1][j] + h2 * V[i + 1][j] + k2 * V[i][j - 1] + k2 * V[i][j + 1] - f(a + h * i, c + k * j);
			rs += r * r;
		}
	}
	return sqrt(rs);
}


double CheckError(double** V, int n, int m, double a, double b, double c, double d)
{
	double** U = CreateMatrix(n + 1, m + 1);
	double h, k;
	double zs = 0;

	h = (b - a) / n;
	k = (d - c) / m;

	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
		{
			double xi, yj;
			xi = a + i * h;
			yj = c + j * k;

			U[i][j] = Uxy(xi, yj);
		}

	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			double z = abs(U[i][j] - V[i][j]);
			if (z > zs)
			{
				zs = z;
			}
		}

	return zs;
}
