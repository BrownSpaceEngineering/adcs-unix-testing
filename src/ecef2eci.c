#include "include/ecef2eci.h"
#include "declareFunctions.h"
#include <math.h>

float unix_2_jd(int unix){
    return unix / 86400.0 + 2440587.5;
}

float jd_2_gmst(float jd){
    float t = (jd - 2451545.0)/36525.0;
    float theta_deg = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * t * t - (t*t*t)/ 38710000.0;
    return fmodf(theta_deg, 360.0);
}
void ecef_2_eci(float* ecef, float* eci, int unix){
    float jd = unix_2_jd(unix);
    float theta_deg = jd_2_gmst(jd);
    float theta = theta_deg * M_PI / 180.0;
    float R[] = {cosf(theta), -sinf(theta), 0, sinf(theta), cosf(theta), 0, 0, 0, 1};
    mulf(R, ecef, false, eci, 3, 3, 1);
}