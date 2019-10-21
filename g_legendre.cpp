#include "g_legendre.h"
#include "g_laguerre.h"
#include "mesh_and_weights.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>


using namespace std;
using namespace std::chrono;

ofstream ofile; ifstream ifile;
double const pi = acos(-1);
double ANALYTICAL = 5.0 * pi * pi / (16.0*16.0); // Setting the analytical solution

double legendre_func_to_integrate(double x_1, double y_1, double z_1, double x_2, double y_2, double z_2) {
    double alpha = 2.0;

    double exponent1 = -2*alpha*sqrt(x_1*x_1 + y_1*y_1 + z_1*z_1);
    double exponent2 = -2*alpha*sqrt(x_2*x_2 + y_2*y_2 + z_2*z_2);
    double denominator = sqrt(pow(x_1 - x_2, 2) + pow(y_1 - y_2, 2) + pow(z_1 - z_2, 2));

    if (denominator < 1e-7) {
        denominator = 0;
    }

    double result = exp(exponent1 + exponent2)*denominator;

    return result;
}


double guass_legendre(double a, double b, int steps, string file) {

    double *weights = new double [steps];
    double *mesh_p = new double [steps];

    r_legendre(steps, mesh_p);
    weights_legendre(steps, weights, mesh_p, a, b);

    double integral_segments = 0;

    high_resolution_clock::time_point t0 = high_resolution_clock::now(); // Stop clock

    for (int i = 0; i < steps; i++) {
        for (int j = 0; j < steps; j++) {
            for (int k = 0; k < steps; k++) {
                for (int l = 0; l < steps; l++) {
                    for (int m = 0; m < steps; m++) {
                        for (int n = 0; n < steps; n++) {

                            integral_segments += weights[i]*weights[j]*weights[k]*weights[l]*weights[m]*weights[n]*
                                                 legendre_func_to_integrate(mesh_p[i], mesh_p[j], mesh_p[k], mesh_p[l], mesh_p[m], mesh_p[n]);
                        }
                    }
                }
            }
        }
    }


    delete[] weights;
    delete[] mesh_p;

    high_resolution_clock::time_point t1 = high_resolution_clock::now(); // Stop clock
    auto time = duration_cast<nanoseconds>(t1 - t0).count() / 1e9; // time difference in nano-sec

    double error = abs(integral_segments - ANALYTICAL);

    // The file to open should be a .csv file, that way, it will be easy to analyze it with Python
    ofile.open("results.csv", ios::app); ifile.open("results.csv"); // File is read as input argument

    ofile << setiosflags(ios::showpoint | ios::uppercase);

    /* This if-statement checks if the file is empty
     * That way, it will not not have to print out the labels each time */
    if (ifile.peek() == EOF) {
        ofile << "steps, exact, numeric, error, variance, time" << endl;
    } // End if

    ofile << setw(1) << setprecision(3) <<  steps     << ", ";
    ofile << setw(1) <<   ANALYTICAL 	<<           ", ";
    ofile << setw(1) << setprecision(3) <<  integral_segments   << ", ";
    ofile << setw(1) << setprecision(3) <<  error << ", ";
    ofile << setw(1) << setprecision(3) <<  time  << endl;


}


