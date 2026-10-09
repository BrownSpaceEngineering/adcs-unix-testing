#include "include/sunvec.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double sind(double in) { return sin(in * M_PI / 180.0); }
static double cosd(double in) { return cos(in * M_PI / 180.0); }

/*
 * Port of sunVectorECI.m. Time arithmetic is in double: a float32 Julian date only resolves
 * 0.25 days (~0.25 deg of solar longitude).
 */
void sun_vec(int unix_time, float* sun) {
    double julianOffset = (double)unix_time / 86400.0 + 2440587.5 - 2451545.0;
    double julianC = julianOffset / 36525.0;

    double meanAnomaly = fmod(357.529 + 35999.050 * julianC, 360.0);
    double meanLongitude = fmod(280.459 + 36000.770 * julianC, 360.0);

    double sunCenter = (1.914602 - 0.004817 * julianC - 0.000014 * julianC * julianC) * sind(meanAnomaly)
                       + (0.019993 - 0.000101 * julianC) * sind(2 * meanAnomaly)
                       + 0.000289 * sind(3 * meanAnomaly);

    double eclipticLongitude = fmod(meanLongitude + sunCenter, 360.0);

    double obliquityEcliptic = 23 + 26.0 / 60.0 + 21.448 / 3600
                               - (46.8150 * julianC + 0.00059 * julianC * julianC
                                  - 0.001813 * julianC * julianC * julianC)
                                     / 3600;

    sun[0] = (float)cosd(eclipticLongitude);
    sun[1] = (float)(cosd(obliquityEcliptic) * sind(eclipticLongitude));
    sun[2] = (float)(sind(obliquityEcliptic) * sind(eclipticLongitude));
}
