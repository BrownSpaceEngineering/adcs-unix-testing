#ifndef PHOTODIODE_DETERMINATION_H
#define PHOTODIODE_DETERMINATION_H

#include <stdbool.h>

extern const float MAX_READING;
extern const float PHOTODIODES[18][3];
extern const float MIN_READING;
#define NUM_DIODES 18
// Least-squares sun vector (body frame, unit) from the 5 brightest pair-winners.
// Returns whether we're in the sun (>= 3 pairs above MIN_READING and a solvable system).
bool get_vec_from_photodiode_readings(const float* photodiode_readings, float* estimated_sun_vector);//returns whether we're in the sun or not

#endif
