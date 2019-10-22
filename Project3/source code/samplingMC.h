#ifndef SAMPLINGMC_H
#define SAMPLINGMC_H

#include <string>

double func_to_integrate(double u, double v, double theta1, double theta2, double phi1, double phi2);

double samplingMC(int a, int b, int n, std::string file);

#endif // SAMPLINGMC_H
