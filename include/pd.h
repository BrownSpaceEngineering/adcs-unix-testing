#ifndef PD_H
#define PD_H
#include "arm_math.h"
void pd_loop(float32_t *omega, float32_t *q_error, float32_t *tau);
#endif
