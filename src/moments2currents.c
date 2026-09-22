#include <math.h>
#include "arm_math.h"

/**
 * \fn moment2current3axis
 * \brief Converts desired moment into currents for 3 orthogonal magnetorquers
 *
 * Inputs:
 *  \param[in] m       : 3-element array, magnetic moment [mx, my, mz] [A·m2]
 *  \param[in] Imax    : 3-element array of current limits [A] (use -1 at an axis if no current limit on that axis)
 *  \param[out] I_out   : 3-element output array for coil currents [Ix, Iy, Iz] [A]
 */
void moment2current3axis(float32_t *m, float32_t *Imax, float32_t *I_out) {

    /* Hard-coded parameters (to be measured/calculated)
     *   n    : turns per coil
     *   A    : coil areas [m2]
     *   L, r : rod geometry (ignored if mu_r == 1)
     *   mu_r : relative permeability (1 for air-core)
     */
    float32_t n[3]    = {1.0, 1.0, 1.0};
    float32_t A[3]    = {1.0, 1.0, 1.0};
    float32_t L[3]    = {1.0, 1.0, 1.0};
    float32_t r[3]    = {1.0, 1.0, 1.0};
    float32_t mu_r[3] = {1.0, 1.0, 1.0};

    /* Compute core amplification factors G */
    float32_t G[3];
    float32_t ratio, Nd;

    for (int i = 0; i < 3; i++) {
        if (mu_r[i] > 1.0) {
            ratio = L[i] / r[i];
            Nd    = 4.0 * (logf(ratio) - 1.0) / (ratio * ratio - 4.0 * logf(ratio));

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