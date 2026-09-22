#ifndef KEPLER
#define KEPLER
#include "arm_math.h"
void propogateOrbitalElements(float32_t semi_major_axis,
                              float32_t eccentricity,
                              float32_t inclination,
                              float32_t ascending_node,
                              float32_t periapsis,
                              float32_t true_anomaly,
                              float32_t time_delta,
                              float32_t *output);
#endif