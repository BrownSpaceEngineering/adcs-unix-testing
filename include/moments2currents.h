#ifndef MOMENTS2CURRENTS
#define MOMENTS2CURRENTS
#include "arm_math.h"
// I = m / (n A G) per axis, clamped to +-Imax[i]. A negative Imax[i] means "no limit".
void moment2current3axis(const float32_t* m, const float32_t* Imax, float32_t* I_out);
#endif
