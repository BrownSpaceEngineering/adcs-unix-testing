#ifndef SGP4_PROPAGATE_H
#define SGP4_PROPAGATE_H
#include "include/sgp4_init.h"

int sgp4_propagate(const sgp4_sat_t *sat, float dt,
                   float *r_gcrf, float *v_gcrf, float *oe_osc, float *oe_mean,
                   float *r_teme, float *v_teme);
#endif
