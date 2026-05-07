#include "include/pd.h"
#include "Include/arm_math_types.h"
#include "Include/dsp/basic_math_functions.h"
#include "arm_math.h"
/**
 * \fn pd_loop
 * 
 * \brief A simple PD control loop for attitude control.
 * 
 * From the error (r_e) in axis-angle form, the angular velocity (r_omega) in axis-angle form, returns the torque wanted 
 * 
 * \param[in] r_e The attitude error in axis-angle form (r_e[0], r_e[1], r_e[2]) is the rotation vector, and r_e[3] is the angle of rotation in radians. 
 * \param[in] r_omega The angular velocity error in axis-angle form (r_omega[0], r_omega[1], r_omega[2]) is the angular velocity vector, and r_omega[3] is the magnitude of the angular velocity in radians per second.
 * \param[out] tau The output torque vector to apply to the satellite in order to correct the attitude error and angular velocity error.
 */
void pd_loop(float32_t *r_e, float32_t *r_omega, float32_t *tau) {

    float32_t placeholder = 1.0;
    float32_t placeholder_kp = placeholder * r_e[3]; 
    float32_t placeholder_kd = placeholder * r_omega[3]; 
    //TODO: find these constants
    float32_t Kp[3] = {-placeholder_kp, -placeholder_kp, -placeholder_kp};
    float32_t Kd[3] = {-placeholder_kd, -placeholder_kd, -placeholder_kd};

    float32_t p[3]; 
    float32_t d[3];

    arm_mult_f32(Kp, r_e, p, 3); 
    arm_mult_f32(Kd, r_omega, d, 3);

    arm_add_f32(p, d, tau, 3);

}