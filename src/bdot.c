#include "include/bdot.h"
#include "string.h"

void Bdot(double* M_T, double* M_TMINUS1, double dT, double* moments){
    double kx = 1;//Adjustable
    double ky = 1;
    double kz = 1;

    double bdot[3] = {-kx / dT * (M_T[0] - M_TMINUS1[0]),
        -ky / dT * (M_T[1] - M_TMINUS1[2]),
        -kz / dT * (M_T[1] - M_TMINUS1[2])};
    memcpy(moments, bdot, 3 * sizeof(double));
}
