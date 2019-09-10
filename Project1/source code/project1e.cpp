#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

ofstream static ofile;


double f(double x){
    return 100.0 * exp(-10.0 * x);
}
double exact(double x) {
    return 1.0 - (1 - exp(-10)) * x-exp(-10 * x);
}


int project1e(string fn, int ex ){
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

    // Loop over powers of 10
    for (int i = 1; i <= exponent; i++){
        unsigned int  n = int(pow(10.0, i));

        string fileout = fn;

        string argument = to_string(i);

        fileout.append(argument);
        double h = 1.0 / (n + 1);
        double hh = h * h;
        // Inorder for solve to work A and b has to have the same dimentions
        mat A = zeros(n, n);
        vec b(n);
        // x(0) = 0 and x(n + 1)=1 thats why x(n + 2)
        vec x(n + 2);

        A(0, 0) = 2.0;  A(0, 1) = -1; x(0) = 0;  x(1) = h;  b(0) =  hh * f(x(1));
        x(n) = (n) * h; b(n - 1) = hh * f(x(n));
        x(n + 1) = 1;

 	    // Start timing
        clock_t start, finish;
        start = clock();

        for (unsigned int i = 1; i < n - 1; i++){
            x(i + 1)     = (i + 1) * h;
            b(i)         = hh * f(x(i + 1));
            A(i, i - 1)  = -1.0;
            A(i, i)      = 2.0;
            A(i, i + 1)  = -1.0;
        }
        A(n-1, n-1) = 2;
        A(n-1, n-2) = -1;

        vec solution  = solve(A, b);
        finish = clock();
        double timeused = double(finish - start) / (double(CLOCKS_PER_SEC) );
        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << setprecision(10) << setw(20) << "Time used  for  computation= " << timeused  << endl;

        ofile.open(fileout);
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        for (unsigned int i = 0; i < n ; i++) {
            double RelativeError = fabs((exact(x(i + 1)) - solution(i)) / exact(x(i + 1)));
            ofile << setw(15) << setprecision(8) << x(i + 1);
            ofile << setw(15) << setprecision(8) << solution(i);
            ofile << setw(15) << setprecision(8) << exact(x(i + 1));
            ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
        }
        ofile.close();
    }
    return 0;
}
