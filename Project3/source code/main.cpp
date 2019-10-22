#include "legendre.h"
#include "laguerre.h"
#include "bruteForceMC.h"
#include "samplingMC.h"
#include "samplingMC_omp.h"

#include <iostream>
#include <string>

using namespace std;

string const dir = "/home/aleksandar/Desktop/FYS4150/Project3/";

int main() {

    double a = -3.0; // lower integration limit
    double b = 3.0;  // upper integration limit
    int n    = 5;

    string file_legendre = dir + "data/legendre.csv";
    string file_laguerre = dir + "data/laguerre.csv";

    laguerre(a, b, n, file_laguerre);
    legendre(a, b, n, file_legendre);

    string file_bfmc    = dir + "data/bruteForceMC.csv";
    string file_mc      = dir + "data/samplingMC.csv";
    string file2_omp_mc = dir + "data/parallel_O2_samplingMC.csv";


    //bruteForceMC(a, b, n, file_bfmc);
    //samplingMC(a, b, n, file_mc);
    //samplingMC_omp(a, b, n, file2_omp_mc);

    return 0;
} // End main()
