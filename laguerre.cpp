#include "laguerre.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <chrono>

#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10


using namespace std;
using namespace std::chrono;

ofstream ofile;
ifstream ifile;
double const pi = acos(-1);
double ANALYTICAL = 5.0 * pi * pi / (16.0*16.0); // Setting the analytical solution


using namespace std;



double function_integrand(double r_1, double r_2, double theta_1, double theta_2, double phi_1, double phi_2){
    /*
     * This function does...
     * ...
     * ...
     * @param r_1
     * this parameter is the ...
     * @param r_2
     * ...
     * @return 0.0 if ..
     * @return exp...
     *
    */
    double alpha = 2.0;

    double cos_beta = cos(theta_1)*cos(theta_2) + sin(theta_1)*sin(theta_2)*cos(phi_1 - phi_2);
    double r_12 = r_1*r_1 + r_2*r_2 - 2*r_1*r_2*cos_beta;

    if (r_12 < 1E-9){
        return 0.0;
    }
    else{
        return (exp(-3.0*(r_1 + r_2))*sin(theta_1)*sin(theta_2))/(sqrt(r_12));
    }
}

void laguerre(double a, double b, int steps, string file){
    double alfa = 2.0;

    double angle_start = 0.0;
    double angle_end = 3.141592;

    double *weight_r = new double[steps+1];
    double *weight_theta = new double[steps + 1];
    double *weight_phi = new double[steps + 1];

    double *x_r = new double[steps + 1];
    double *x_theta = new double[steps + 1];
    double *x_phi = new double[steps + 1];

    gauss_laguerre(x_r, weight_r, steps, alfa);
    gauleg(angle_start, angle_end, x_theta, weight_theta, steps);
    gauleg(angle_start, 2*angle_end, x_phi, weight_phi, steps);

    double sum = 0.0;

    high_resolution_clock::time_point t0 = high_resolution_clock::now(); // Stop clock

    for(int i = 0; i <= steps; i++){
        for(int j = 0; j <= steps; j++){
            for(int k = 0; k <= steps; k++){
                for(int l = 0; l <= steps; l++){
                    for(int m = 0; m <= steps; m++){
                        for(int n = 0; n <= steps; n++){

                           sum += weight_r[i]*weight_r[j]*weight_theta[k]*weight_theta[l]*weight_phi[m]*weight_phi[n]*(
                                        function_integrand(x_r[i],x_r[j], x_theta[k], x_theta[l], x_phi[m], x_phi[n]));

                        }
                    }
                }
            }
        }
    }

    high_resolution_clock::time_point t1 = high_resolution_clock::now(); // Stop clock
    auto time = duration_cast<nanoseconds>(t1 - t0).count() / 1e9; // time difference in nano-sec

    cout << time << ", ";

    double error = abs(sum - ANALYTICAL);

    // The file to open should be a .csv file, that way, it will be easy to analyze it with Python
    ofile.open(file, ios::app); ifile.open(file); // File is read as input argument

    ofile << setiosflags(ios::showpoint | ios::uppercase);

    /* This if-statement checks if the file is empty
     * That way, it will not not have to print out the labels each time */
    if (ifile.peek() == EOF) {
        ofile << "steps, exact, numeric, error, variance, time" << endl;
    } // End if

    ofile << setw(1) << setprecision(3) <<  steps << ", ";
    ofile << setw(1) <<   ANALYTICAL 	<<           ", ";
    ofile << setw(1) << setprecision(3) <<  sum   << ", ";
    ofile << setw(1) << setprecision(3) <<  error << ", ";
    ofile << setw(1) << setprecision(3) <<  time  << endl;


    cout << sum << endl;
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
           */

      do {
         p1 =1.0;
     p2 =0.0;

       /*
       ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

     for(j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
     }

       /*
       ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */

     pp = n * (z * p1 - p2)/(z * z - 1.0);
     z1 = z;
     z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
      ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()


void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}
// end function gaulag

double gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

