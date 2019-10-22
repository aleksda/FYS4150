#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <string>

using namespace std;

void legendre(double a, double b, int steps, string file);
double legendre_func_to_integrate(double x_1, double y_1, double z_1, double x_2, double y_2, double z_2);
void gauleg(double x1, double x2, double x[], double w[], int n);

#endif // LEGENDRE_H
