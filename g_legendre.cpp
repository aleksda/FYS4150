#include "g_legendre.h"

#include <cmath>


using namespace std;


double guass_legendre(double x_1, double y_1, double z_1, double x_2, double y_2, double z_2, int steps) {
    double integral_segments = 0;
    double alpha = 2.0;

    double *weights = new double [steps];
    double *mesh_p = new double [steps];

    double exponent1 = -2*alpha*sqrt(x_1*x_1 + y_1*y_1 + z_1*z_1);
    double exponent2 = -2*alpha*sqrt(x_2*x_2 + y_2*y_2 + z_2*z_2);
    double denominator = sqrt(pow(x_1 - x_2, 2) + pow(y_1 - y_2, 2) + pow(z_1 - z_2, 2));

    for (int i = 0; i < steps; i++) {
        for (int j = 0; j < steps; j++) {
            for (int k = 0; k < steps; k++) {
                for (int l = 0; l < steps; l++) {
                    for (int m = 0; m < steps; m++) {
                        for (int n = 0; n < steps; n++) {
                            if (denominator < 1e-7) {
                                denominator = 0;
                            }
                            integral_segments += weights[i]*weights[j]*weights[k]*weights[l]*weights[m]*weights[n]*
                                                 exp(exponent1 + exponent2)*denominator;
                        }
                    }
                }
            }
        }
    }
    delete[] weights;
    delete[] mesh_p;

    return integral_segments;
}
