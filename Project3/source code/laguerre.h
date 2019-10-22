#ifndef LAGUERRE_H
#define LAGUERRE_H

#include <string>

double function_integrand(double r_1, double r_2, double theta_1, double theta_2, double phi_1, double phi_2);
void laguerre(double a, double b, int steps, std::string file);
void gauleg(double x1, double x2, double x[], double w[], int n);
void gauss_laguerre(double *x, double *w, int n, double alf);
double gammln( double xx);

#endif // LAGUERRE_H
