#include "include/bdot.h"

void Bdot(float* M_T, float* M_TMINUS1, float dT, float* moments){
    float kx = 1;//Adjustable
    float ky = 1;
    float kz = 1;

    float bdot[3] = {-kx / dT * (M_T[0] - M_TMINUS1[0]),
        -ky / dT * (M_T[1] - M_TMINUS1[2]),
        -kz / dT * (M_T[1] - M_TMINUS1[2])};
    memcpy(moments, bdot, 3 * sizeof(float));
}