#include "samplingMC.h"
#include "lib.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

ofstream outfile; ifstream infile;

double const PI = acos(-1); // Defining pi
double const ANALYTICAL = 5.0 * PI * PI / (16.0*16.0); // Setting the analytical solution

// Integrand function in polar coordinates
double func_to_integrate(double u, double v, double theta1, double theta2, double phi1, double phi2) {

    double alpha = 2.0;

    double cos_beta = cos(theta1) * cos(theta2) + sin(theta1) * sin(theta2) * cos(phi1-phi2);
    double r12 = sqrt(u*u + v*v - 2*u*v*cos_beta);

    if (r12 < 1e-6 || isnan(r12)) return 0;
    else return 1.0 / (pow(2*alpha, 5)) * sin(theta1) * sin(theta2) * u * u * v * v / r12;

} // End func_to_integrate()

double samplingMC(int a, int b, int n, string file) {

    double func;

    double cmc 	  = 0.0;
    double sigma  = 0.0;

    double jdet 	= 4 * pow(PI, 4);	//pow((b-a), dim);
    //double var 	= sigma - cmc * cmc;

    high_resolution_clock::time_point t0 = high_resolution_clock::now(); // Start clock

    long idum       = -1;
   	double u, v, theta1, theta2, phi1, phi2;

    for(int i = 0; i < n; i++) {

        u 		= -log(1-ran0(&idum));
        v 		= -log(1-ran0(&idum));		// v = u
        theta1 	= PI * ran0(&idum);
        theta2 	= PI * ran0(&idum);			// theta2 = theta1
        phi1 	= 2 * PI * ran0(&idum);
        phi2 	= 2 * PI * ran0(&idum);		// phi1 = phi2

        func = func_to_integrate(u, v, theta1, theta2, phi1, phi2);

        cmc   += func;
        sigma += func * func;

    } // End for

    //cmc   /= N; cmc   *= jdet;
    //sigma /= N; sigma *= jdet;

    cmc /= (double) n; sigma /= (double) n;

	double var 	= sigma - cmc * cmc;

    cmc *= jdet; //sigma *= jdet;

    high_resolution_clock::time_point t1 = high_resolution_clock::now(); // Stop clock
    auto time = duration_cast<nanoseconds>(t1 - t0).count() / 1e9; // time difference in nano-sec

    double error = abs(cmc - ANALYTICAL);

    // The file to open should be a .csv file, that way, it will be easy to analyze it with Python
    outfile.open(file, ios::app); infile.open(file); // File is read as input argument

    outfile << setiosflags(ios::showpoint | ios::uppercase);

    /* This if-statement checks if the file is empty
     * That way, it will not not have to print out the labels each time */
    if (infile.peek() == EOF) {
        outfile << "steps, exact, numeric, error, variance, time" << endl;
    } // End if

    outfile << setw(1) << setprecision(3) <<  n     << ", ";
    outfile << setw(1) <<   ANALYTICAL 	<<             ", ";
    outfile << setw(1) << setprecision(3) <<  cmc   << ", ";
    outfile << setw(1) << setprecision(3) <<  error << ", ";
    outfile << setw(1) << setprecision(3) <<  var   << ", ";
    outfile << setw(1) << setprecision(3) <<  time  << endl;

    cout << time << endl;

    outfile.close();

    return 0.0;

} // End samplingMC()
