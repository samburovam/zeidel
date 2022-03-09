#include "Numerical.h"
#include "Numerical.cpp"
#include "Numerical.h"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <stdio.h> 
#include <time.h>
#include <chrono>
#include <thread>
using namespace std;


void clrscr()
{
	system("@cls||clear");
}


// U"xx + U"yy = -f(x, y)
// a <= x <= b, c <= y <= d
// U(a, y) = M1(y), U(b, y) = M2(y)
// U(x, c) = M3(x), U(x, d) = M3(x)

int main(void)
{
	setlocale(LC_ALL, "rus");

	std::cout << std::endl << std::endl;

	int		Nmax = 10000;	
	int		S = 0;				
	double	eps = 0.0000001;	
	double	epsMax = 0;			
	int		n = 3, m = 3;		
	double** V = NULL;			
	double** F = NULL;			
	double	a, b, c, d;			
	int Exit = 1, Show = 0;

	a = 0;
	b = 2;
	c = 0;
	d = 1;

	V = MemoryAllocator(n + 1, m + 1);
	FillStartSolution(V, n, m, a, b, c, d);
	ZeidelsMethod(V, n, m, a, b, c, d, eps, Nmax, epsMax, S);
	ShowSolution(V, n, m);

	cout << endl << endl;
	cout << n << "x" << m << endl;
	cout.precision(15);
	cout << "S is" << S << endl;
	cout << "eps is " << eps << endl;
	cout << "epsMax is " << epsMax << endl;
	cout << "Error is " << CheckError(V, n, m, a, b, c, d) << endl;
	cout << "Disc is  " << DiscrepancyOfSolution(V, n, m, a, b, c, d) << endl;

	MemoryCleaner(V, n);

	cout << endl;
	cin >> Nmax;

	return 0;
}