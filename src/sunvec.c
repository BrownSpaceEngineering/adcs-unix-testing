#include "include/sunvec.h"
#include <math.h>
#include <string.h>
float sind(float in){
    return sin(in * M_PI / 180.0);
}
float cosd(float in){
    return cos(in * M_PI/180.0);
}
void sun_vec(int unix, float* sun){
    float julianDate = (float)unix/86400 + 2440587.5;
    float julianOffset = julianDate - 2451545.0;
    float julianC = julianOffset / 36525;

    float meanAnomaly = 357.529 + 35999.050*julianC;
    float meanLongitude = 280.459 + 36000.770*julianC;

    meanAnomaly = fmod(meanAnomaly, 360);
    meanLongitude = fmod(meanLongitude, 360);

    float sunCenter = (1.914602 - 0.004817*julianC - 0.000014*julianC*julianC)*sind(meanAnomaly) + (0.019993 - 0.000101*julianC)*sind(2*meanAnomaly) + 0.000289*sind(3*meanAnomaly);

    float eclipticLongitude = meanLongitude + sunCenter;
    eclipticLongitude = fmod(eclipticLongitude, 360);

    float obliquityEcliptic = 23 + 26.0f/60.0f + 21.448f/3600 - (46.8150f*julianC + 0.00059*julianC*julianC - 0.001813*julianC*julianC*julianC)/3600;

    float ans[3] = {cosd(eclipticLongitude), cosd(obliquityEcliptic) * sind(eclipticLongitude), sind(obliquityEcliptic) * sind(eclipticLongitude)};
    memcpy(sun, ans, 3 * sizeof(float));
}