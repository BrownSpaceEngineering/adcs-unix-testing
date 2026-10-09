#ifndef PD_H
#define PD_H
#include "arm_math.h"

// Gains/limits from PD_loop.m
#define PD_KP 0.034f
#define PD_KD 0.42f
#define PD_MAX_TAU 5.0f
extern const float32_t PD_INERTIA[9]; // body inertia (kg m^2), row-major

// Port of PD_loop.m: tau = I * (Kp * theta * axis - Kd * omega), clamped to +-PD_MAX_TAU.
// q_error: body-frame error quaternion (q_b2eci^-1 * q_want, scalar first).
// omega: body rates in RAD/S (the MATLAB version takes deg/s and converts).
void pd_loop(const float32_t* q_error, const float32_t* omega, float32_t* tau);
#endif
