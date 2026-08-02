#ifndef KEPLER
#define KEPLER
void propogateOrbitalElements(double semi_major_axis, double eccentricity, double inclination, double ascending_node, double periapsis, double true_anomaly, double time_delta, double *array);
#endif