#include "pd.h"

void PD_loop(float* r_e, float* r_omega, float* tau)
/* From the error (r_e) in axis-angle form, the angular velocity (r_omega) in axis-angle form, returns the torque wanted */
{

    float placeholder = 1.0;
    //TODO: find these constants
    float Kp[3] = {-placeholder, -placeholder, -placeholder};
    float Kd[3] = {-placeholder, -placeholder, -placeholder};

    for (int i = 0; i < 3; i++) {

        float Pi = Kp[i] * r_e[3] * r_e[i];
        float Di = Kd[i] * r_omega[3] * r_omega[i];

        tau[i] = Pi + Di;
    }
}