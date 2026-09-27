#include "include/pd.h"
#include "arm_math.h"
#include <math.h>

#define PD_EPS 1e-12f

/**
 * \fn pd_loop
 *
 * \brief PD control loop for attitude control.
 *
 * \param[in]  omega    Angular velocity w_eci2b in body coordinates, deg/s (3 elements).
 * \param[in]  q_error  Error quaternion from pointing_error, scalar first, WXYZ (4 elements).
 *                       May be negated in place to pick the shortest-path rotation.
 * \param[out] tau      Output torque in N-m, body coordinates (3 elements).
 */
void pd_loop(float32_t *omega, float32_t *q_error, float32_t *tau) {
    // Pick the shortest possible rotation
    if (q_error[0] < 0.0f) {
        arm_negate_f32(q_error, q_error, 4);
    }

    float32_t angle = 2.0f * acosf(q_error[0]);

    // Convert error quaternion to axis-angle: r_e = [axis(3); angle]
    float32_t r_e[4];
    if (angle > PD_EPS) {
        float32_t half_sin = arm_sin_f32(angle / 2.0f);
        arm_scale_f32(&q_error[1], 1.0f / half_sin, r_e, 3);
    } else {
        r_e[0] = 0.0f;
        r_e[1] = 0.0f;
        r_e[2] = 0.0f;
    }
    r_e[3] = angle;

    // Convert omega (deg/s) to axis-angle form: r_omega = [axis(3); magnitude]
    float32_t omega_rad[3];
    arm_scale_f32(omega, PI / 180.0f, omega_rad, 3);

    float32_t omega_mag_sq;
    arm_power_f32(omega_rad, 3, &omega_mag_sq);
    float32_t omega_mag;
    arm_sqrt_f32(omega_mag_sq, &omega_mag);

    float32_t r_omega[4];
    if (omega_mag > PD_EPS) {
        arm_scale_f32(omega_rad, 1.0f / omega_mag, r_omega, 3);
    } else {
        r_omega[0] = 0.0f;
        r_omega[1] = 0.0f;
        r_omega[2] = 0.0f;
    }
    r_omega[3] = omega_mag;

    // PD controller gains
    const float32_t Kp0 = 0.034f;
    const float32_t Kd0 = 0.42f;
    const float32_t max_tau = 5.0f;
    const float32_t min_tau = -5.0f; //THIS NEEDS TUNING

    //THIS ALSO NEEDS TUNING
    float32_t I_body_data[9] = {
        3.0054115e-02f, -5.1674000e-05f,  2.5075000e-05f,
       -5.1674000e-05f,  1.1430611e-02f, -3.6481690e-03f,
        2.5075000e-05f, -3.6481690e-03f,  2.3052295e-02f
    };
    arm_matrix_instance_f32 I_body;
    arm_mat_init_f32(&I_body, 3, 3, I_body_data);

    float32_t tau_raw[3];
    for (int i = 0; i < 3; i++) {
        float32_t Pi = Kp0 * r_e[3] * r_e[i];
        float32_t Di = Kd0 * r_omega[3] * r_omega[i];
        tau_raw[i] = Pi - Di;
    }

    arm_mat_vec_mult_f32(&I_body, tau_raw, tau);

    for (int i = 0; i < 3; i++) {
        if (tau[i] > max_tau) tau[i] = max_tau;
        if (tau[i] < min_tau) tau[i] = min_tau;
    }
}
