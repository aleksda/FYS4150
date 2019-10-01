#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

#include "catch.hpp"
#include "tridiagtest.h"

//using namespace std; using namespace arma;


// Defining pi ($cos_{-1}(-1.0)$)
#define PI acos(-1.0)
//#define PI 3.14159265

TEST_CASE("Testing for eigenvalues for a tridiagonal matrix") {

	// Getting exact values
	int dim = 10;

    vec eigen(dim), exac(dim);

    double step, diag, nondiag;
    double r_max;
    r_max = 1.0;

	// Step-lengths
	step 	= r_max / dim;
	diag 	= 2.0 / (step * step);
	nondiag = -1.0 / (step * step);

	for(int i = 0; i < dim; i++) {

        exac(i) = diag + 2 * nondiag * cos ((i + 1) * PI / (dim + 1));
	}

	// Nummerical values
    eigen = testForEigens(dim);

	// Defining tolerance
    double tol = 1e14;	// or 12 ?

	REQUIRE(eigen(0) == Approx(exac(0)).epsilon(tol));
	REQUIRE(eigen[1] == Approx(exac(1)).epsilon(tol));
	REQUIRE(eigen[2] == Approx(exac(2)).epsilon(tol));

}

#undef PI

