#ifndef PHOTODIODE_DETERMINATION_H
#define PHOTODIODE_DETERMINATION_H

#include <stdbool.h>

extern const float MAX_READING;
extern const float PHOTODIODES[18][3];
extern const float MIN_READING;
bool get_vec_from_photodiode_readings(float* photodiode_readings, float* estimated_sun_vector);//returns whether we're in the sun or not

#endif