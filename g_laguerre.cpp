#include <Integrands.h>
#include <cmath>


using namespace std;


double const pi = acos(-1);
double alpha = 2;
double a = -1.0;
double b = 1.0;






double guass_laguerre(double r_1, double r_2, double theta_1, double theta_2, double phi_1 double phi_2, int steps) {
    double integral_segments = 0;

    double *theta = new double[steps];
    double *phi = new double[steps];
    double *r = new double[steps + 1];

    double *w_theta = new double[steps];
    double *w_phi = new double[steps];
    double *w_r = new double[steps + 1];

    double theta_1, theta_2,    delta_theta_1, delta_theta_2;
    double phi_1, phi_2,    delta_phi_1, delta_phi_2;
    double r_1, r_2,    delta_r_1, delta_r_2;


    double cos_beta = cos(theta_1)*cos(theta_2) + sin(theta_1)*sin(theta_2)*cos(phi_1 - phi_2);
    double r_12 = 1/sqrt(r_1*r_1 + r_2*r_2 - 2*r_1*r_2*cos_beta);

    for (int i = 0; i < steps; i++) {
        delta_r_1 = w_r[i + 1];
        r_1 = r[i + 1];
        for (int j = 0; j < steps; j++) {
            delta_r_2 = w_r[j + 1];
            r_2 = r[j + 1];
            for (int k = 0; k < steps; k++) {
                delta_theta_1 = w_theta[k];
                theta_1 = theta[k];
                for (int l = 0; l < steps; l++) {
                    delta_theta_2 = w_theta[l];
                    theta_2 = theta[l];
                    for (int m = 0; m < steps; m++) {
                        delta_phi_1 = w_phi[m];
                        phi_1 = phi[m];
                        for (int n = 0; n < steps; n++) {
                            delta_phi_2 = w_phi[n];
                            phi_2 = phi[n];
                            if (r_12 < 1e-7) {
                                r_12 = 0;
                            }
                            integral_segments += delta_r_1*delta_r_2*delta_theta_1*delta_theta_2*delta_phi_1*delta_phi_2*
                                                 (1.0 / (32*pow(alpha, 5))*sin(theta_1)*sin(theta_2)*r_12);
                        }
                    }
                }
            }
        }
    }
    delete[] theta;
    delete[] phi;
    delete[] r;

    delete[] w_theta;
    delete[] w_phi;
    delete[] w_r;

    return integral_segments;
}


