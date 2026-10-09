#include "include/pd.h"
#include "include/quat.h"
#include "arm_math.h"

const float32_t PD_INERTIA[9] = {3.0054115e-02f,  -5.1674000e-05f, 2.5075000e-05f,
                                 -5.1674000e-05f, 1.1430611e-02f,  -3.6481690e-03f,
                                 2.5075000e-05f,  -3.6481690e-03f, 2.3052295e-02f};

/**
 * \fn pd_loop
 *
 * \brief PD attitude controller, port of PD_loop.m.
 *
 * The previous version read r_e[3] / r_omega[3] past the end of 3-element arrays and used
 * those garbage values as gains.
 *
 * \param[in] q_error Body-frame error quaternion (rotation from current to desired attitude)
 * \param[in] omega Body angular rate (rad/s)
 * \param[out] tau Commanded torque (N m, body frame)
 */
void pd_loop(const float32_t* q_error, const float32_t* omega, float32_t* tau) {
    // quat2rotationvec already takes the shortest rotation (flips q if w < 0), which gives
    // theta * axis = r_e(4) * r_e(1:3) from the MATLAB code
    float32_t r_e[3];
    quat2rotationvec(q_error, r_e);

    float32_t t[3];
    for (int i = 0; i < 3; i++) {
        // omega_mag * omega_axis == omega
        t[i] = PD_KP * r_e[i] - PD_KD * omega[i];
    }
    for (int i = 0; i < 3; i++) {
        float32_t v = PD_INERTIA[i * 3 + 0] * t[0] + PD_INERTIA[i * 3 + 1] * t[1]
                      + PD_INERTIA[i * 3 + 2] * t[2];
        if (v > PD_MAX_TAU) {
            v = PD_MAX_TAU;
        } else if (v < -PD_MAX_TAU) {
            v = -PD_MAX_TAU;
        }
        tau[i] = v;
    }
}
