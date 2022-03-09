#include "numerical.h"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <stdio.h> 

using namespace std;


void clrscr()
{
	system("@cls||clear");
}


int main(void)
{
	setlocale(LC_ALL, "rus");

	int		Nmax = 10000;
	int		S = 1;
	double	eps = 0.0000001;
	double	epsMax = 0;
	int		n = 3, m = 3;
	double** V = NULL;
	double** F = NULL;
	double	a, b, c, d;
	int Exit = 1, Show = 0;

	a = 0;
	b = 1;
	c = 0;
	d = 1;

	V = CreateMatrix(n + 1, m + 1);
	FillStartSolution(V, n, m, a, b, c, d);
	ZeidelsMethod(V, n, m, a, b, c, d, eps, Nmax, epsMax, S);

	cout << endl;
	cout << "Справка: " << endl;
	cout << endl;
	cout << "Выполнено шагов: " << S << endl;
	cout << "Точность: " << eps << endl;
	cout << "Точность на выходе: " << epsMax << endl;
	cout << "Общая погрешность: " << CheckError(V, n, m, a, b, c, d) << endl;
	cout << "Невязка: " << DiscrepancyOfSolution(V, n, m, a, b, c, d) << endl;

	cout << endl;

	return 0;
}