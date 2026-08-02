#include "include/ecef2eci.h"
#include "declareFunctions.h"
#include <math.h>

static double unix_2_jd(int unix){
    return unix / 86400.0 + 2440587.5;
}

static double jd_2_gmst(double jd){
    double t = (jd - 2451545.0) / 36525.0;
    double theta_deg = 280.46061837 + 360.98564736629*(jd - 2451545.0) + 0.000387933*t*t - (t*t*t)/38710000.0;
    return fmod(theta_deg, 360.0);
}
void ecef_2_eci(double* ecef, double* eci, int unix){
    double jd = unix_2_jd(unix);
    double theta_deg = jd_2_gmst(jd);
    double theta = theta_deg * M_PI / 180.0;
    double R[] = {cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1};
    mul(R, ecef, false, eci, 3, 3, 1);
}
