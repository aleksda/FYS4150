
#include "jacobi.h"

#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<iomanip>
#include<time.h>
#include<algorithm>
#include<vector>
#include<armadillo>

//using namespace std;
using namespace arma;

std::ofstream ofile;

mat makematrix(int dim) {

    //int dim = 10; // N, mesh points
    double rhomin, rhomax, step, d, a;
    rhomin = 0.0; rhomax = 1.0;

    /* Defining the step-length
    * Defining the main diagonal for the matrix
    * Defining the upper/lower diagonal */
    step  = (rhomax - rhomin) / dim;
    d     = 2.0 / (step * step);   // Diagonal Constant
    a     = -1.0 / (step * step); // Upper/Lower Diagonal Constant

    // Empty dim x dim matrix
    mat A = zeros<mat>(dim, dim);

    /* Setting up the tridiagonal matrix
     * Diagonalizing using Armadillo */
    A(0, 0) = d;
    A(0, 1) = a;

    for(int i = 1; i < dim-1; i++) {

      A(i, i - 1)    = a;
      A(i, i)        = d;
      A(i, i + 1)    = a;

    }

    A(dim - 1, dim - 2) = a;
    A(dim - 1, dim - 1) = d;

    // Returns the matrix
    return A;
}

/* Algoritm for finding the maximum matrix element
 * From lecture-notes Ch.7 */
void maxoffdiag(mat &A, int& k, int& l, double& max, int dim){	// &A
    //int dim = 10;

    //double max = 0.0;
    //std::cout<<"maxoffdiag"<<std::endl;
    //std::cout << A << std::endl;
    for(int i = 0; i < dim; i++){
        for(int j = i + 1; j < dim; j++){
            if (fabs(A(i, j)) > max){
                //cout << A(i, j) << " " << max << " " << (fabs(A(i,j)) > max) << endl;
                //cout << '------------------------------------' << '--' << endl;
                //cout <<  "j = " << j << endl;
                max = fabs(A(i, j));
                l = i; k = j;
            }
        }
    }
    cout << "l = " << l << "k = " << k << endl;

}

/* Rotation algoritm
 * For finding sine and cosine
 * From lecture-notes Ch. 7 */
void rotation(mat &A, int &k, int &l, int dim) {	// &A

    // A.shape, A.cols.size
    //int dim = 10;
    double c, s;
    mat R = eye(dim, dim);
    //cout << k << l <<endl << A;
    //double theta;

    //s=sin(theta), c=cos(theta);

    // Defining the variables
    if(A(k, l) != 0.0){
        // Tau
        double tau;
        tau = (A(l, l) - A(k, k)) / (2 * A(k, l));
        //tau =cot(2 * theta);

        double t;
        // t = tan(theta);
        // Value 't' can be either positive or negative
        if (tau >= 0) {

            t = 1.0 / (tau + sqrt(1 + pow(tau, 2)));

        } else {

            t = -1.0 / (-tau + sqrt(1 + pow(tau, 2)));
        }
        // Sine and cosine
        c = 1.0 / (sqrt(1 + pow(t,2)));
        s = c * t;

    // Division by zero
    } else {
        s = 1.0; c = 0.0;
    }


    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k, k);
    a_ll = A(l, l);

    A(k, k) = a_kk * c * c - 2 * A(k,l) * c * s + a_ll * s * s;
    A(l, l) = a_ll * c * c + 2 * A(k,l) * c * s + a_kk * s * s;
    A(l, k) = 0.0; // Non-diagonal elements
    A(k, l) = 0.0;

    for (int i = 0; i < dim; i++) {
        if (i != k && i != l) {

            a_ik = A(i, k);
            a_il = A(i, l);
            A(i, k) = a_ik * c - a_il * s;
            A(k, i) = A(i, k);
            A(i, l) = a_il * c + a_ik * s;
            A(l, i) = A(i, l);
            //cout <<"it54f = "<< i;
        }

	    // Theese are the new eigenvectors
        r_ik = R(i, k);
        r_il = R(i, l);
        R(i, k) = c * r_ik - s * r_il;
        R(i, l) = c * r_il + s * r_ik;

    }

    //return A, R;

}

mat jacobimethod(mat &A, mat &R, int dim) {

    //int dim = 10;
    int k, l;
    double eps = 1e-8;

    // Maximum number of iterations
    double maxiter = 4;//(double) dim * (double) dim * (double) dim;
    int iteration = 0;
    double max_offdiag;
    maxoffdiag(A, k, l, max_offdiag, dim);

    // The rotation algoritm
    while (max_offdiag > eps && (double) iteration <= maxiter) {
        max_offdiag = 0;
        rotation(A, k, l, dim);
        //max_offdiag = maxoffdiag(A, k, l, max_offdiag);
        maxoffdiag(A, k, l, max_offdiag, dim);
        iteration++;
        //cout << A <<endl;

    }

    std::cout << "Iterations: " << iteration << std::endl;
    return A;
}
/*
void writefile(mat &R, std::string fn) {

    //int dim = (int) A.size() / 2;
    int dim = 2;
    int index;

    std::ofstream ofile;
    ofile.open(fn.c_str());

    for (int i = 0; i < dim; i++) {
        ofile << R(i, index) << std::endl << std::endl;

     }

     ofile.close();
}
*/
/* Not a part of the jacobi algorithm */

vec analytic_eigen(int dim) {

    double step = 1.0 / dim;

    double diag = 2.0 / (step * step);
    double notdiag = -1.0 / (step * step);

    vec v(dim);
    const double pi = acos(-1.0);

    for (int i = 1; i <= dim; i++){
        v(i - 1) = diag + 2 * notdiag * cos(i * pi / (dim + 1));
    }

    return v;
}

vec eigensolv(mat &F, int dim) {

    int i;
    vec e(dim);

    for (i = 0; i < dim; i++) {
        e(i) = F(i, i);
    }

    e = sort(e);

    return e;   // e.sort() ??
}


void eigenvector(mat &A) {

    int dim = A.n_cols;
    mat B = A.t() * A;

    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, A);
    eigval.print("Eigenvalues for A: ");
    eigvec.print("Eigenvectors for A: ");

    //cout << eigvec << endl;

    //return eigvec;
}
