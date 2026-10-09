#ifndef ECEF2ECI
#define ECEF2ECI
double unix_2_jd(double unix_time);
double jd_2_gmst_deg(double jd); // GMST in degrees, [0, 360)
void ecef_2_eci(const float* ecef, float* eci, int unix_time);
#endif
