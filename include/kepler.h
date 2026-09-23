#ifndef KEPLER
#define KEPLER
#include "arm_math.h"
// Two-body propagation (port of "Two Body Propogation.m").
// a in metres, angles in DEGREES; output = [a, e, i, RAAN, argp, nu_new], nu_new in (-180, 180].
void propogateOrbitalElements(float32_t semi_major_axis,
                              float32_t eccentricity,
                              float32_t inclination,
                              float32_t ascending_node,
                              float32_t periapsis,
                              float32_t true_anomaly,
                              float32_t time_delta,
                              float32_t *output);
#endif
