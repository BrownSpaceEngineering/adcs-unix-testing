#ifndef PHOTODIODE_DETERMINATION_H
#define PHOTODIODE_DETERMINATION_H

#include <stdbool.h>

extern const double MAX_READING;
extern const double PHOTODIODES[18][3];
extern const double MIN_READING;
bool get_vec_from_photodiode_readings(double* photodiode_readings, double* estimated_sun_vector);//returns whether we're in the sun or not

#endif