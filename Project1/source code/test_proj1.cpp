#include "project1b.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;
/*
inline double f(double x) {
    return 100.0 * exp(-10.0 * x);
}
// Analytic solution
inline double u(double x) {
    return 1.0 - (1 - exp(-10)) * x - exp(-10 * x);
}*/

string test_proj1()
{/*
    unsigned int n = 3;

    double *x = new double[n];
    double *y = new double[n];
    double *z = new double[n-1];
    double *w = new double[n-1];

    vector<int> vec1 {1, 2, 3};
    vector<int> vec2 {10, 20, 30};
    vector<int> vec3 {4, 5};
    vector<int> vec4 {6, 7};

    for (unsigned int i = 0; i < n; i++) {
        x[i] = vec1;
        y[i] = vec2;
    }
    for (unsigned int i = 0; i < n-1; i++) {
        z[i] = vec3;
        w[i] = vec4;
    }

    string str = "";
    int e  = 1;
    project1b(str, e);

    for (unsigned int i = 0; i < n; i++) {
        cout << func[i] << endl;
    }

    //delete vec1, vec2, vec3, vec4, *x, *y, *z, *w*/

    int ex;
    int exponent = ex;

    // Loop over powers of 10
    for (int i = 1; i <= exponent; i++) {
        int  n = int(pow(10.0, i));

        // Step-size
        double h = 1.0 / (n + 1);
        double hh = h * h;

        double *diag = new double [n+1]; double *a = new double [n+1]; double *c = new double [n+1]; double *f_tilde = new double [n+1]; double *v = new double [n+1];
        double *x = new double[n+1]; double *b_tilde = new double [n+1]; double *d = new double [n+1]; double *u_exact = new double [n+1];

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

        delete [] x; delete [] diag; delete [] f_tilde; delete [] v; delete [] a; delete [] c; delete [] b_tilde; delete [] d; delete [] u_exact;
    }
    return "fwergw";
}
