#ifndef JACOBI_H
#define JACOBI_H

#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include<time.h>
#include<algorithm>
#include<vector>
#include<armadillo>

using namespace arma;

mat makematrix(int dim);

/* Algoritm for finding the maximum matrix element
 * From lecture-notes Ch.7 */
void maxoffdiag(mat &A, int& k, int& l, double& max, int dim);

/* Rotation algoritm
 * For finding sine and cosine
 * From lecture-notes Ch. 7 */
void rotation(mat &A, int &k, int &l, int dim);

mat jacobimethod(mat &A, mat &R, int dim);

void writefile(mat &R, std::string fn);

vec analytic_eigen(int dim);
vec eigensolv(mat &F, int dim);
void eigenvector(mat &A);

#endif // JACOBI_H
