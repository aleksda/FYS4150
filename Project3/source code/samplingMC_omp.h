#ifndef SAMPLINGMC_OMP_H
#define SAMPLINGMC_OMP_H

#include <string>

double func_to_integrate(double u, double v, double theta1, double theta2, double phi1, double phi2);

inline double ran();

double samplingMC_omp(int a, int b, int n, std::string file);

#endif // SAMPLINGMC_OMP_H
