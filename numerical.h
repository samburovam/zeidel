#pragma once


double Uxy(double x, double y);
double f(double x, double y);
double M1(double y);
double M2(double y);
double M3(double x);
double M4(double x);
double** CreateMatrix(int n, int m);
void ShowSolution(double** V, int n, int m, int s);
double** FillRightSide(double** F, int n, int m, double a, double c, double h, double k);
void FillStartSolution(double** V, int n, int m, double a, double b, double c, double d);
double ZeidelsMethod(double** V, int n, int m, double a, double b, double c, double d, double eps, int Nmax, double& epsMax, int& S);
double DiscrepancyOfSolution(double** V, int n, int m, double a, double b, double c, double d);
double CheckError(double** V, int n, int m, double a, double b, double c, double d);
