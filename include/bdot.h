#ifndef BDOT_H
#define BDOT_H

#include "arm_math.h"

void Bdot(float32_t* M_T, float32_t* M_TMINUS1, float32_t dT, float32_t* moments); 

#endif // BDOT_H