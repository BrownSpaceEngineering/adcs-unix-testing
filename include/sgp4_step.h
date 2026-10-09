#ifndef SGP4_STEP_H
#define SGP4_STEP_H

int sgp4_step(const float *oe_epoch, float bstar, double epoch_jd, float dt,
              float *r_gcrf, float *v_gcrf, float *oe_osc);
#endif
