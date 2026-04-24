#include <math.h>

/**
 * moment2current3axis
 * Converts desired moment into currents for 3 orthogonal magnetorquers
 *
 * Inputs:
 *   m       : 3-element array, magnetic moment [mx, my, mz] [A·m2]
 *   Imax    : 3-element array of current limits [A] (use -1 at an axis if no current limit on that axis)
 *   I_out   : 3-element output array for coil currents [Ix, Iy, Iz] [A]
 */
void moment2current3axis(float* m, float* Imax, float* I_out)
{
    /* Hard-coded parameters (to be measured/calculated)
     *   n    : turns per coil
     *   A    : coil areas [m2]
     *   L, r : rod geometry (ignored if mu_r == 1)
     *   mu_r : relative permeability (1 for air-core)
     */
    float n[3]    = {1.0, 1.0, 1.0};
    float A[3]    = {1.0, 1.0, 1.0};
    float L[3]    = {1.0, 1.0, 1.0};
    float r[3]    = {1.0, 1.0, 1.0};
    float mu_r[3] = {1.0, 1.0, 1.0};

    /* Compute core amplification factors G */
    float G[3];
    for (int i = 0; i < 3; i++) {
        if (mu_r[i] > 1.0) {
            float ratio = L[i] / r[i];
            float Nd    = 4.0 * (log(ratio) - 1.0) / (ratio * ratio - 4.0 * logf(ratio));

            /* Core amplification */
            G[i] = 1.0 + (mu_r[i] - 1.0) / (1.0 + (mu_r[i] - 1.0) * Nd);
        } else {
            G[i] = 1.0;
        }
    }

    /* Convert moments to currents: I = m / (n * A * G) */
    for (int i = 0; i < 3; i++) {
        I_out[i] = m[i] / (n[i] * A[i] * G[i]);
    }

    /* Saturation clamp (only applied when Imax is provided) */
    for (int i = 0; i < 3; i++) {
        if(Imax[i] == -1){
            continue;
        }
        else if (I_out[i] >  Imax[i]) {
            I_out[i] =  Imax[i];
        }
        else if (I_out[i] < -Imax[i]) {
            I_out[i] = -Imax[i];
        }
    }
}