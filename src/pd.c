#include "include/pd.h"

void PD_loop(double* r_e, double* r_omega, double* tau)
/* From the error (r_e) in axis-angle form, the angular velocity (r_omega) in axis-angle form, returns the torque wanted */
{

    double placeholder = 1.0;
    //TODO: find these constants
    double Kp[3] = {-placeholder, -placeholder, -placeholder};
    double Kd[3] = {-placeholder, -placeholder, -placeholder};

    for (int i = 0; i < 3; i++) {

        double Pi = Kp[i] * r_e[3] * r_e[i];
        double Di = Kd[i] * r_omega[3] * r_omega[i];

        tau[i] = Pi + Di;
    }
}
