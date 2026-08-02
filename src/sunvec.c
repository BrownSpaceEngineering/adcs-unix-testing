#include "include/sunvec.h"
#include <math.h>
#include <string.h>
static double sind(double in){ return sin(in * M_PI / 180.0); }
static double cosd(double in){ return cos(in * M_PI / 180.0); }
void sun_vec(int unix, double* sun){
    double julianDate = (double)unix / 86400.0 + 2440587.5;
    double julianOffset = julianDate - 2451545.0;
    double julianC = julianOffset / 36525.0;

    double meanAnomaly  = fmod(357.529 + 35999.050*julianC, 360.0);
    double meanLongitude = fmod(280.459 + 36000.770*julianC, 360.0);

    double sunCenter = (1.914602 - 0.004817*julianC - 0.000014*julianC*julianC)*sind(meanAnomaly)
        + (0.019993 - 0.000101*julianC)*sind(2*meanAnomaly)
        + 0.000289*sind(3*meanAnomaly);

    double eclipticLongitude = fmod(meanLongitude + sunCenter, 360.0);

    double obliquityEcliptic = 23.0 + 26.0/60.0 + 21.448/3600.0
        - (46.8150*julianC + 0.00059*julianC*julianC - 0.001813*julianC*julianC*julianC)/3600.0;

    sun[0] = cosd(eclipticLongitude);
    sun[1] = cosd(obliquityEcliptic) * sind(eclipticLongitude);
    sun[2] = sind(obliquityEcliptic) * sind(eclipticLongitude);
}
