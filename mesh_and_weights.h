#ifndef MESH_AND_WEIGHTS_H
#define MESH_AND_WEIGHTS_H


double legendre(int order, double x);
double dLegendre(int order, double x);
void r_legendre(int order, double *roots);
void weights_legendre(int order, double *weights, double *roots, double a, double b);
double gammln( double xx);
void gauss_laguerre(double *x, double *w, int n, double alf);

#endif // MESH_AND_WEIGHTS_H
