#ifndef BDOT_H
#define BDOT_H

#include "arm_math.h"

// m = -k * (M_T - M_TMINUS1) / dT   (port of bDot.m). Writes zeros if dT <= 0.
void Bdot(const float32_t* M_T, const float32_t* M_TMINUS1, float32_t k, float32_t dT,
          float32_t* moments);

#endif // BDOT_H
