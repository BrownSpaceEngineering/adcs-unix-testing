#include "include/ecef2eci.h"
#include "Include/dsp/matrix_functions.h"
#include "declareFunctions.h"
#include "arm_math.h"
#include <math.h>

float32_t unix_2_jd(int unix){
    return ((float32_t)unix / (float32_t)86400.0) + (float32_t)2440587.5;
}

float32_t jd_2_gmst(float32_t jd){
    float32_t t = (jd - 2451545.0)/36525.0;
    float32_t theta_deg = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * t * t - (t*t*t)/ 38710000.0;
    return fmodf(theta_deg, 360.0);
}
void ecef_2_eci(arm_matrix_instance_f32* ecef, arm_matrix_instance_f32* eci, int unix){
    float32_t jd = unix_2_jd(unix);
    float32_t theta_deg = jd_2_gmst(jd);
    float32_t theta = theta_deg * M_PI / 180.0;
    float32_t R_data[] = {cosf(theta), -sinf(theta), 0, 
                  sinf(theta), cosf(theta), 0,
                   0, 0, 1};
    const arm_matrix_instance_f32 R = {3, 3, (float32_t *)R_data};


    arm_mat_mult_f32(&R, ecef, eci);
}