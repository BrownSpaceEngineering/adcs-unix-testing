#include "include/kepler.h"
#include "math.h"
void propogateOrbitalElements(double semi_major_axis, double eccentricity, double inclination, double ascending_node, double periapsis, double true_anomaly, double time_delta, double *output) {

double u = 3.986004418e14;
double a = semi_major_axis;
double e = eccentricity;
double f = true_anomaly * M_PI / 180.0;
double dt = time_delta;

double E = 2 * atan2(sqrt(1 - e) * sin(f/2), sqrt(1 + e) * cos(f/2));
double M = E - e * sin(E);
double n = sqrt(u/pow(a, 3.0));
double M_new = M + n * dt;

double E_current = M_new;
double tolerance = 1e-12;
int MAXIMUM_ITERATIONS = 10;

for (int i = 1; i < MAXIMUM_ITERATIONS; ++i) {
    double E_change = ((E_current - e * sin(E_current) - M_new)/(1 - e * cos(E_current)));
    E_current = E_current - E_change;
    if(fabs(E_change) < tolerance)
        break;
}
output[0] = a;
output[1] = e;
output[2] = inclination;
output[3] = ascending_node;
output[4] = periapsis;

double new_true_anomaly = 2 * atan2(sqrt(1 + e) * sin(E_current/2), sqrt(1 - e) * cos(E_current/2));
output[5] = new_true_anomaly * 180.0 / M_PI;
}
