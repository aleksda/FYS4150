#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

#include "tridiagtest.h"

using namespace std; using namespace arma;

vec testForEigens(int dim) {

	mat A = zeros<mat>(dim, dim);
	//mat A = mat(dim, dim, fill::eye);

    vec eigen(dim);

    double step, diag, nondiag, r_max;
    r_max = 1.0;

	// Step-lengths
	step 	= r_max / dim;
	diag 	= 2.0 / (step * step);
	nondiag = -1.0 / (step * step);

	// Setting up the matrix
	A(0, 0) = diag;
	A(0, 0) = nondiag;

	for(int i = 1; i < dim - 1; i++) {
		A(i, i - 1)	= nondiag;
		A(i, i)		= diag;
		A(i, i + 1)	= nondiag;
	}

	/* Try writing it here
	A(0, 0) = diag;
	A(0, 0) = nondiag;
	*/

    A(dim - 1, dim - 1) = diag;
    A(dim - 1, dim - 2) = nondiag;

	// Getting the values, setting them in a variable
	//sol = eig_sym(eigen, A);
    eig_sym(eigen, A);

	return eigen;
}
