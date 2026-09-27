#ifndef SGP4_STEP_H
#define SGP4_STEP_H

int sgp4_step(const double *oe_epoch, double bstar, double epoch_jd, double dt,
              double *r_gcrf, double *v_gcrf, double *oe_osc);
#endif
