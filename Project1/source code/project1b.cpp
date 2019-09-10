//#include "project1b.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

ofstream static ofile;

// Numeric solution
inline double f(double x) {
    return 100.0 * exp(-10.0 * x);
}
// Analytic solution
inline double u(double x) {
    return 1.0 - (1 - exp(-10)) * x - exp(-10 * x);
}


int project1b(string fn, int ex) {
    /** Short comment: Part 1b
     *
     * Algorithm for solving for solving part b in [1]
     *  ...
     *  ...
     *  ...
     * 
     * @param fn
     * Takes filename as argument
     * 
     * @param ex
     * Exponent (i.e 10, 100, 1000)
     *
     * @return 0
     *
     * Source
     * ------
     * [1] https://github.com/CompPhysics/ComputationalPhysics/tree/master/doc/Projects/2019/Project1
     *
    */
    int exponent = ex;
    string filename = fn;

    // Loop over powers of 10
    for (int i = 1; i <= exponent; i++) {
        int  n = int(pow(10.0, i));

        string fileout = filename;

        // Giving the text file a number corresponding to the exponent so that it will be easier to distinguish from the others
        string fileNr = to_string(i);
        fileout.append(fileNr);

        // Step-size
        double h = 1.0 / (n + 1);
        double hh = h * h;

        double *diag = new double [n+1]; double *a = new double [n+1]; double *c = new double [n+1]; double *f_tilde = new double [n+1]; double *v = new double [n+1];
        double *x = new double[n+1]; double *b_tilde = new double [n+1]; double *d = new double [n+1]; double *u_exact = new double [n+1];

        // zeros at index 0 for every array

//        diag[0] = 0;
//        b_tilde[0] = 0;
//        v[0] = 0;
//        a[0] = 0;
//        c[0] = 0;
//        f_tilde[0] = 0;
//        u_exact[0] = 0; but i found out that i dont need to do that!

        //Declaring x and u (the exact function)
        for (int i = 0; i <= n; i++){
            x[i]= i * h;
            u_exact[i] = u(x[i]);
        }

        //Filling with corresponding values
        for (int i = 1; i <= n; i++){
            a[i] = -1;
            c[i] = -1;
            d[i] = 2;
            b_tilde[i] = hh * f(x[i]);

        }
        // Boundary conditions
        a[1] = 0;
        c[n] = 0;

        // Start timing
        clock_t start, finish;
        start = clock();

        // Forward substitution
        //Boundary conditions
        diag[1] = d[1];
        f_tilde[1] = b_tilde[1];
        for (int i = 2; i <= n; i++) {
            diag[i] = d[i] - ((a[i] * c[i - 1]) / diag[i - 1]);
            f_tilde[i] = b_tilde[i] - ((f_tilde[i - 1] * a[i]) / diag[i - 1]);
        }

        //Backward substitution
        //Boundary condition
        v[n] = f_tilde[n] / diag[n];
        for (int i = n - 1; i >= 1; i--){
            v[i] = (f_tilde[i] - c[i] * v[i + 1]) / diag[i];
        }
        //End timing
        finish = clock();

        double timeused = double(finish - start) / (double(CLOCKS_PER_SEC));
        //cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << setprecision(10) << setw(20) << "Time used  for  computation=" << timeused  << endl;

        ofile.open(fileout);
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << "     x:              approx:         exact:      " << endl;
        for (int i = 1; i <= n; i++) {
            ofile << setw(15) << setprecision(8) << x[i];
            ofile << setw(15) << setprecision(8) << v[i];
            ofile << setw(15) << setprecision(8) << u_exact[i] << endl;
        }
        ofile.close();
        //Memory Allocation: Deleting pointers
        delete [] x; delete [] diag; delete [] f_tilde; delete [] v; delete [] a; delete [] c; delete [] b_tilde; delete [] d; delete [] u_exact;
    }
    return 0;
}
