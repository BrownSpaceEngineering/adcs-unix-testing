#ifndef PD_H
#define PD_H
#include "arm_math.h"
void pd_loop(float32_t *r_e, float32_t *r_omega, float32_t *tau); 
#endif