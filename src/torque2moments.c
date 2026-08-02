#include "include/torque2moments.h"
#include "declareFunctions.h"
#include "include/laextension.h"
void torque_2_moments(double* B, double* torques, double* moments){
    double Bxt[3] = {0, 0, 0};
    cross(B, torques, Bxt);
    double squared_norm = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
    scale(Bxt, 1.0/squared_norm, 1, 3);
    memcpy(moments, Bxt, 3 * sizeof(double));
}

