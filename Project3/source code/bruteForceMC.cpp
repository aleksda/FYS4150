
#include "bruteForceMC.h"
#include "lib.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

ofstream ofile; ifstream ifile;

double const PI = acos(-1); // Defining pi
double ANALYTICAL = 5.0 * PI * PI / (16.0*16.0); // Setting the analytical solution

// Integrand function in cartesian coordinates
double func_to_integrate(double *x) {

    double alpha = 2;

    double x1 = x[0]; double y1 = x[1]; double z1 = x[2];
    double x2 = x[3]; double y2 = x[4]; double z2 = x[5];

    double exponential1 = -2 * alpha * sqrt(x1*x1 + y1*y1 + z1*z1);
    double exponential2 = -2 * alpha * sqrt(x2*x2 + y2*y2 + z2*z2);

    double denominator = sqrt(pow((x1-x2), 2) + pow((y1-y2), 2) + pow((z1-z2), 2));

    if (denominator < pow(10., -6.)) return 0;
    else return exp(exponential1 + exponential2) / denominator;


} // End func_to_integrate()


double bruteForceMC(double a, double b, int n, string file) {

    double func;

    int dim      = 6;
    double cmc 	 = 0.0;
    double sigma = 0.0;

    double jdet     = pow((b-a), dim);	 // Jacobi-determinant
    double var      = 0.0;

    high_resolution_clock::time_point t0 = high_resolution_clock::now(); // Start clock

    // Evaluate integral with importance sampling
    double x[dim]; long idum = -1;
    for (int i = 0; i < n; i++) {
        // x[] contains the random numbers
        for (int j = 0; j < dim; j++) {
            x[j] = a + (b-a) * ran0(&idum);	// See 'lib.h' for the ran0() function
            //x[j] = double(rand()) * (1./RAND_MAX) * (b-a) + a;
        } // End for j

        func = func_to_integrate(x);

        cmc 	+= func;
        sigma   += func * func;

    } // End for i

    //cmc   /= n; cmc   *= jdet;
    //sigma /= n; sigma *= jdet;

    cmc = cmc / (double) n; sigma = cmc / (double) n;

    var = sigma - cmc * cmc; // Varance

    cmc *= jdet; sigma *= jdet;

    high_resolution_clock::time_point t1 = high_resolution_clock::now(); // Stop clock
    auto time = duration_cast<nanoseconds>(t1 - t0).count() / 1e9; // time difference in nano-sec

	double error = abs(cmc - ANALYTICAL);

    // The file to open should be a .csv file, that way, it will be easy to analyze it with Python
    ofile.open(file, ios::app); ifile.open(file); // File is read as input argument

    ofile << setiosflags(ios::showpoint | ios::uppercase);

    /* This if-statement checks if the file is empty
     * That way, it will not not have to print out the labels each time */
    if (ifile.peek() == EOF) {
        ofile << "steps, exact, numeric, error, variance, time" << endl;
    } // End if

    ofile << setw(1) << setprecision(3) <<  n     << ", ";
    ofile << setw(1) <<   ANALYTICAL 	<<           ", ";
    ofile << setw(1) << setprecision(3) <<  cmc   << ", ";
    ofile << setw(1) << setprecision(3) <<  error << ", ";
    ofile << setw(1) << setprecision(3) <<  var   << ", ";
    ofile << setw(1) << setprecision(3) <<  time  << endl;

    cout << time << endl;

    ofile.close();

    return 0.0;

} // End bruteForceMC()
