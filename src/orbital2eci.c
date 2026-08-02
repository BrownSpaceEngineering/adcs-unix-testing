#include "include/orbital2eci.h"
#include "include/quat.h"
#include "math.h"

void orbitalToECI(double* kepler6, double* eci){
    double smAxis       = kepler6[0];
    double eccentricity = kepler6[1];
    double inclination  = kepler6[2];
    double aNodeLongitude = kepler6[3];
    double periapsisArg = kepler6[4];
    double trueAnomaly  = kepler6[5];

    double r0 = smAxis*(1-eccentricity*eccentricity)/(1+(eccentricity*cos(trueAnomaly)));
    double r1[3] = {r0*cos(trueAnomaly), r0*sin(trueAnomaly), 0};

    double qPerapsisArg[4] = {cos(periapsisArg/2), 0, 0, sin(periapsisArg/2)};
    double qInclination[4] = {cos(inclination/2), sin(inclination)/2, 0, 0};
    double qLongitude[4]   = {cos(aNodeLongitude/2), 0, 0, sin(aNodeLongitude/2)};

    double qInclinationPeriapsis[4];
    quat_multiply(qInclination, qPerapsisArg, qInclinationPeriapsis);
    double qComb[4];
    quat_multiply(qLongitude, qInclinationPeriapsis, qComb);

    quat_apply(qComb, r1, eci);
}
