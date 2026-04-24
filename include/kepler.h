#ifndef KEPLER
#define KEPLER
void propogateOrbitalElements(float semi_major_axis, float eccentricity, float inclination, float ascending_node, float periapsis, float true_anomaly, float time_delta, float *array);
#endif