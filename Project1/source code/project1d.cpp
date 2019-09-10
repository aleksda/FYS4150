#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

ofstream static ofile;

inline double f(double x) {
    return 100.0*exp(-10.0*x);
}
inline double exact(double x) {
    return 1.0-(1-exp(-10))*x-exp(-10*x);
}


int project1d(string fn, int ex) {
    /** Short comment: This is for ...
     * 
     * Describing the code
     *  ...
     *  ...
     *  ...
     * 
     * @param fn
     * String that does something ...
     * 
     * @param ex
     * Integer that does something
     *
     * @return 0
    */
    int exponent = ex;
    string filename = fn;

    // Loop over powers of 10
    for (int i = 1; i <= exponent; i++) {
        int  n = int(pow(10.0, i));

        string fileout = filename;
        // Giving the text file a number corresponding to the exponent so that it will be easier to distinguish from the others
        string argument = to_string(i);
        fileout.append(argument);

        double h = 1.0 / (n + 1);
        double hh = h * h;

        double *d = new double [n+1]; double *b = new double [n+1]; double *v = new double [n+1]; double *diag = new double [n+1];
        double *x = new double[n+2]; double *u = new double[n+1]; double *f_tilde = new double[n+1];
        // boundary conditions
        v[1] = v[n] = 0.0;

        x[n + 1] = 1;
        for (int i = 1; i <= n; i++){
            d[i] = 2.0;

        }
        for (int i = 0; i <= n; i++){
            x[i] = i * h;
            u[i] = exact(x[i]);
            b[i] = hh * f(x[i]);
        }

      // Start timing
        clock_t start, finish;
        start = clock();

      // Forward substitution
        f_tilde[1] = b[1];
        diag[1] = d[1];
        for (int i = 2; i <= n; i++) {
            diag[i] = d[i] - (1.0 / diag[i - 1]);
            f_tilde[i] = b[i] + f_tilde[i - 1] / diag[i - 1];
        }
      // Backward substitution
        v[n] = f_tilde[n] / diag[n];
        for (int i = n - 1; i >= 1; i--){
            v[i] = (f_tilde[i] + v[i + 1]) / diag[i];
        }

        finish = clock();
        double timeused = double(finish - start) / (double(CLOCKS_PER_SEC) );
        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << setprecision(10) << setw(20) << "Time used  for  computation=" << timeused  << endl;

        ofile.open(fileout);
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << "     x:             approx:        exact:       relativeerror:" << endl;
        for (int i = 1; i <=n; i++) {
            double RelativeError = fabs((v[i]-u[i]) / u[i]);
            ofile << setw(15) << setprecision(8) << x[i];
            ofile << setw(15) << setprecision(8) << v[i];
            ofile << setw(15) << setprecision(8) << u[i];
            ofile << setw(15) << setprecision(8) << RelativeError << endl;
        }
        ofile.close();
        //Memory allocation: Deleting pointers
        delete [] x; delete [] d; delete [] b; delete [] v; delete [] u; delete [] f_tilde; delete [] diag;
    }
    return 0;
}
