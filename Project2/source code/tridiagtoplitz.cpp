#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <armadillo>

using namespace  std;
using namespace  arma;

#define PI acos(-1.0)
#define POTENTIAL(x) ( (x) * (x) )

/* This is the first part of 2.b)
 * Where we are going to set up
 * the matrix and diagonalizing
 * using Armadillo */
mat tridiagtoplitz(string method, int dim) {


    /* Backround
    * Defining variables */
    //int dim; // N, mesh points

    double rhomin, rhomax, step, d, a;
    rhomin = 0.0; rhomax = 1.0; //dim = 4;

    /* Defining the step-length
    * Defining the main diagonal for the matrix
    * Defining the upper/lower diagonal */
    step    = (rhomax - rhomin) / dim;
    d       = 2.0 / (step * step);	// Diagonal Constant
    a       = -1.0 / (step * step);	// Upper/Lower Diagonal Constant

    // Empty dim x dim matrix
    mat A = zeros<mat>(dim,dim);
    // Eigenvalue vecotr
    vec eigenval(dim);
    // Step-Size vector
    vec stepvec = zeros<vec>(5);
    // Potential vector
    vec potentialvec = zeros<vec>(dim);
    //
    mat eigenmatrix = zeros(5, 20);

    for(int i = 0; i < dim; i++) {
        double rho = rhomin + (i + 1) * step;
        potentialvec[i] = POTENTIAL(rho); // POTENTIAL?
    }

    /* This is the first part of 2.b)
     * Where we are going to set up
     * the matrix and diagonalizing
     * using Armadillo */
    if (method == "buckling" || method == "0") {
        /* Setting up the tridiagonal matrix
           Diagonalizing using Armadillo */

        cout << "Running part 2b) " << endl;

        A(0, 0) = d;
        A(0, 1) = a;

        for(int i = 1; i < dim - 1; i++) {

            A(i, i - 1)    = a;
            A(i, i)        = d;
            A(i, i + 1)    = a;
        }

        A(dim - 1, dim - 2) = a;
        A(dim - 1, dim - 1) = d;

        /* Diagonalization
           Obtaining eigenvalues */
        vec eigenval(dim);
        eig_sym(eigenval, A);


        cout << "RESULTS: " << endl;
        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout <<"Number of Eigenvalues = " << setw(15) << dim << endl;
        cout << "Exact versus numerical eigenvalues:" << endl;
        for(int i = 0; i < dim; i++) {
        double Exact = d + 2 * a * cos((i + 1) * PI / (dim + 1));	//pi
        cout << setw(15) << setprecision(8) << fabs(eigenval[i] - Exact) << endl;

        // /*
        cout << setw(15) << setprecision(8) << "Exact: " << Exact << endl;
        cout << setw(15) << setprecision(8) << "Numerical: " << eigenval[i] << endl;
        cout << "------------------------" << endl;
        // */

        }

    /* This is the solution to 2.d)
    * Where we are going to set up
    * the TDM for the two dimensional
    * quantum harmonic oscillator TODO */
    } else if (method == "quantumone" || method == "1") {

        cout << "Running part 2d) " << endl;

        A(0, 0) = d;
        A(0, 1) = a;

        for(int i = 1; i < dim - 1; i++) {

            A(i, i - 1)    = a;
            A(i, i)        = d + POTENTIAL(i); // def this!!
            A(i, i + 1)    = a;
        }

        A(dim - 1, dim - 2) = a;
        A(dim - 1, dim - 1) = d + POTENTIAL(dim - 1);

        /* Diagonalization
           Obtaining eigenvalues */
        vec eigenval(dim);
        eig_sym(eigenval, A);

    /* TODO */
    } else if (method == "quantumtwo" || method == "2") {

        cout << "Running part 2e) " << endl;

        double omega = 1.0; // No less than 1.0
        omega *= omega;

        A(0, 0) = d + omega + 1.0 / potentialvec[1];
        A(0, 1) = a;

        for(int i = 1; i < dim - 1; i++) {

            A(i, i - 1)    = a;
            A(i, i)        = d + omega * POTENTIAL(i) + 1.0 / potentialvec[1];
            A(i, i + 1)    = a;
        }

        A(dim - 1, dim - 2) = a;
        A(dim - 1, dim - 1) = d + omega * POTENTIAL(dim - 1) + 1.0 / potentialvec[dim - 1];

        /* Diagonalization
           Obtaining eigenvalues */
        vec eigenval(dim);
        eig_sym(eigenval, A);

    } else {

        cout << "Call for tridiagtoplitz() must have one argument " << endl;
        cout << "Either 0: buckling, 1: quantumone, 2: quantumtwo " << endl;
        exit(1);
    }

    //cout << A << endl;
    //cout <<  "VALUES " << eig_sym(A) << endl;

    return A;
} //  end of function

#undef PI
