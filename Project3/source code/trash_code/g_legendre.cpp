#include "g_legendre.h"
#include "g_laguerre.h"
#include "mesh_and_weights.h"


#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>



#define   ZERO       1.0E-10

using namespace std;
using namespace std::chrono;

ofstream ofile; ifstream ifile;
double const pi = acos(-1);
double ANALYTICAL = 5.0 * pi * pi / (16.0*16.0); // Setting the analytical solution

double legendre_func_to_integrate(double x_1, double y_1, double z_1, double x_2, double y_2, double z_2) {
    double alpha = 2.0;

    double exponent1 = exp(-2*alpha*sqrt(x_1*x_1 + y_1*y_1 + z_1*z_1));
    double exponent2 = exp(-2*alpha*sqrt(x_2*x_2 + y_2*y_2 + z_2*z_2));
    double denominator = sqrt(pow(x_1 - x_2, 2) + pow(y_1 - y_2, 2) + pow(z_1 - z_2, 2));

    if (denominator < 1e-7) {
        denominator = 0;
    }

    double result = (exponent1 + exponent2)*denominator;

    return result;
}


double guass_legendre(double a, double b, int steps, string file) {

    double *weights = new double [steps];
    double *mesh_p = new double [steps];

    r_legendre(steps, mesh_p);
    weights_legendre(steps, weights, mesh_p, a, b);

    //gauleg(a, b, weights, mesh_p, steps);


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
    ofile.open("file_legend.csv", ios::app); ifile.open("file_legend.csv"); // File is read as input argument

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


void gauleg(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
       ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.


      do {
         p1 =1.0;
     p2 =0.0;

       /*
       ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x


     for(j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
     }


       ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.


     pp = n * (z * p1 - p2)/(z * z - 1.0);
     z1 = z;
     z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
      ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart


      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()

