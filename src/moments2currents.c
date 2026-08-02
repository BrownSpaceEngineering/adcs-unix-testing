#include "include/moments2currents.h"
#include <math.h>

void moment2current3axis(double* m, double* Imax, double* I_out)
{
    double n[3]    = {1.0, 1.0, 1.0};
    double A[3]    = {1.0, 1.0, 1.0};
    double L[3]    = {1.0, 1.0, 1.0};
    double r[3]    = {1.0, 1.0, 1.0};
    double mu_r[3] = {1.0, 1.0, 1.0};

    double G[3];
    for (int i = 0; i < 3; i++) {
        if (mu_r[i] > 1.0) {
            double ratio = L[i] / r[i];
            double Nd    = 4.0 * (log(ratio) - 1.0) / (ratio * ratio - 4.0 * log(ratio));
            G[i] = 1.0 + (mu_r[i] - 1.0) / (1.0 + (mu_r[i] - 1.0) * Nd);
        } else {
            G[i] = 1.0;
        }
    }

    for (int i = 0; i < 3; i++) {
        I_out[i] = m[i] / (n[i] * A[i] * G[i]);
    }

    for (int i = 0; i < 3; i++) {
        if(Imax[i] == -1.0){
            continue;
        } else if (I_out[i] >  Imax[i]) {
            I_out[i] =  Imax[i];
        } else if (I_out[i] < -Imax[i]) {
            I_out[i] = -Imax[i];
        }
    }
}
