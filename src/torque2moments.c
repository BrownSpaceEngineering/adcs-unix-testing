#include "include/torque2moments.h"
#include "include/laextension.h"
#include "arm_math.h"

/**
 * \fn torque_2_moments
 *
 * \brief Magnetorquer dipole for a desired torque: m = (B x tau) / |B|^2.
 *
 * A magnetorquer produces tau = m x B, and (B x tau) x B / |B|^2 is the part of tau
 * perpendicular to B (the only part a magnetorquer can make).
 *
 * NOTE: torque2moment3axis.m computes pinv(skew(B)) * tau = -(B x tau) / |B|^2, i.e. it solves
 * B x m = tau, which gives the opposite torque. This C version keeps the physically correct
 * sign; the MATLAB should probably be fixed to match.
 */
void torque_2_moments(const float* B, const float* torques, float* moments) {
    float32_t squared_norm = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];
    if (!(squared_norm > 1e-30f)) {
        moments[0] = 0.0f;
        moments[1] = 0.0f;
        moments[2] = 0.0f;
        return;
    }
    float32_t Bxt[3];
    cross(B, torques, Bxt);
    for (int i = 0; i < 3; i++) {
        moments[i] = Bxt[i] / squared_norm;
    }
}
