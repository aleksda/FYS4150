/* This file contains the algorithm for Brute-Force MonteCarlo integration
 * All files in the project can be compiled normally
 * You can also compile the entire project directly using an IDE that supports .pro files (i.e QTCreator)
 */
#ifndef BRUTEFORCEMC_H
#define BRUTEFORCEMC_H

#include <string>

double func_to_integrate(double *x);
double bruteForceMC(double a, double b, int n, std::string file);

#endif // BRUTEFORCEMC_H
