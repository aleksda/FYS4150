#include "jacobi.h"
#include "tridiagtoplitz.h"
#include "quantumdot_two.h"
#include <ctime>
#include <string>
#include <armadillo>
#include <iomanip>

using namespace std;

int main() {

    ofstream ofile;

    int dim = 100;
    double omega = 5.0;
    double rhomax;

    string method = "quantumtwo";
    mat A = tridiagtoplitz(method, dim);
    //cout << A << endl;

    eigenvector(A);

    //mat A = makematrix(dim);
    //cout <<  A << endl;

    //quantumone();

    //std::cout << quantumone(4, 1.0, "qdot2", 1.0) << std::endl;

    int n = dim / 2;
    mat R = eye(n, n);

    int k, l;
    double max;
    maxoffdiag(A, k, l, max, dim);
    //std::cout << "kl= " << k << std::endl << l << std::endl;
    //mat B, R; = rotation(A, k, l);
    rotation(A, k, l, dim);
    mat F = jacobimethod(A, R, dim);
    //cout << F << endl;

    string file1 = method + "_"; string file2;
    if (method == "buckling" || method == "0") {
        file2 = to_string(dim) + ".txt";

    } else if (method == "quantumone" || method == "1") {
        file2 = to_string(dim) + ".txt";

    } else if (method == "quantumtwo" || method == "2") {
        file2 = to_string(dim) + "_rhomax_" + "0" + to_string(int(rhomax)) + "_w_" + to_string(omega) + ".txt";

    } else {cout << "Wrong input argument " << endl;}

    file1.append(file2);

/* ----------------------------------------------------------------------------------------------------------------------- */

    // Not used
    ofile.open(file1);

    vec eigen = eigensolv(F, dim);
    if (method == "buckling" || method == "0") {
        vec analytic = analytic_eigen(dim);

        //ofile << "Analytical" << "\t" << "Computed" << "\t" << "Relative Error" << endl;
        for (int i = 0; i < dim; i++) {
            double rError = (fabs((analytic(i) - eigen(i)) / analytic(i)));
            ofile << setw(20) << setprecision(8) << analytic(i);
            ofile << setw(20) << setprecision(8) << eigen(i);
            ofile << setw(20) << setprecision(4) << rError << endl;
        }


    } else if (method == "quantumone" || method == "1") {
        vec analytic4x4; analytic4x4 << 3 << 7 << 11 << 15 << endr; // for 4x4 matrix
        //ofile << "Analytical" << "\t" << "Computed" << "\t" << "Relative Error" << endl;
        for (int i = 0; i < dim; i++) {
            double rError = fabs(analytic4x4(i) - eigen(i) / analytic4x4(i));
            ofile << setw(20) << setprecision(8) << analytic4x4(i);
            ofile << setw(20) << setprecision(8) << eigen(i);
            ofile << setw(20) << setprecision(4) << rError << endl;
        }

    }











    clock_t begin = clock();

    //cout << eig_sym(A) << endl;
    //std::cout << A << std::endl;

    //std::cout << F << std::endl;

//    std::string file = "test.txt";
    //int dummy; writefile(R, file);

    //clock_t end = clock();
    //double time = (end - begin) / (double) CLOCKS_PER_SEC;

//    std::cout << time << std::endl;
//    std::cout << file.c_str() << endl;
    return 0;
}
