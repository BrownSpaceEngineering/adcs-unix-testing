#include "include/ecef2eci.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Port of eceftoeci.m.
 *
 * The Julian date and GMST are computed in double: a JD (~2.46e6) held in a float32 only
 * resolves 0.25 days, and 360.98 * (JD - 2451545) (~3e6 degrees) only resolves ~0.25 deg, so
 * the float32 version was off by up to ~90 degrees of Earth rotation. Only a handful of
 * double operations are needed per call, so this is cheap even with soft-float doubles.
 */

double unix_2_jd(double unix_time) { return unix_time / 86400.0 + 2440587.5; }

double jd_2_gmst_deg(double jd) {
    double d = jd - 2451545.0;
    double t = d / 36525.0;
    double theta_deg = 280.46061837 + 360.98564736629 * d + 0.000387933 * t * t - (t * t * t) / 38710000.0;
    theta_deg = fmod(theta_deg, 360.0);
    if (theta_deg < 0.0) {
        theta_deg += 360.0;
    }
    return theta_deg;
}

void ecef_2_eci(const float* ecef, float* eci, int unix_time) {
    double theta = jd_2_gmst_deg(unix_2_jd((double)unix_time)) * (M_PI / 180.0);
    float c = (float)cos(theta);
    float s = (float)sin(theta);
    float x = c * ecef[0] - s * ecef[1];
    float y = s * ecef[0] + c * ecef[1];
    eci[0] = x;
    eci[1] = y;
    eci[2] = ecef[2];
}
