#ifndef PHOTODIODE_DETERMINATION_H
#define PHOTODIODE_DETERMINATION_H

extern const float MAX_READING;
extern const float PHOTODIODES[18][3];

void get_vec_from_photodiode_readings(float photodiode_readings[],
                                      float estimated_sun_vector[]);

#endif