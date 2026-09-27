#ifndef SGP4_PROPAGATE_H
#define SGP4_PROPAGATE_H
#include "include/sgp4_init.h"

int sgp4_propagate(const sgp4_sat_t *sat, double dt,
                   double *r_gcrf, double *v_gcrf, double *oe_osc, double *oe_mean,
                   double *r_teme, double *v_teme);
#endif
