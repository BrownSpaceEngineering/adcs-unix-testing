#include "include/torque2moments.h"
#include "include/laextension.h"
#include "Include/dsp/basic_math_functions.h"
#include "arm_math.h"
#include <string.h>
void torque_2_moments(float* B, float* torques, float* moments){
    float32_t Bxt[3] = {0,0,0};
    cross(B, torques, Bxt);
    float32_t squared_norm = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
    arm_scale_f32(Bxt, 1.0f/squared_norm, Bxt, 3);
    memcpy(moments, Bxt, 3 * sizeof(float32_t));
}
