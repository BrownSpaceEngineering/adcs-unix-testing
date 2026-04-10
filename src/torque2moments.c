#include "include/torque2moments.h"
#include "declareFunctions.h"
#include "include/laextension.h"
void torque_2_moments(float* B, float* torques, float* moments){
    float Bxt[3] = {0,0,0};
    cross(B, torques, Bxt);
    float squared_norm = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
    scalar_multiply(Bxt, 1.0f/squared_norm, 3);
    memcpy(moments, Bxt, 3 * sizeof(float));
}
