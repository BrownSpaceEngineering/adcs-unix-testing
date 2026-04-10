#include "include/orbital2eci.h"
#include "include/quat.h"

void orbitalToECI(float* kepler6, float* eci){
    // % Converts the 6 orbital elements into ECI position vector
    // % 
    // % All angles must be in radians
    // % Units are consistent between semi-major axis and the resultant vector
    // % Now position in km, velocity in km/s
    
    // % r0 is the current radius in the perifocal frame (plane of orbit)
    // % since the orbit is not purely circular, it gives the current distance
    
    // % r1 is the position in the perifocal frame
    // % Adjusted for the offset of the true anomaly
    // % Axes: periapsis, 90 deg, normal

    float smAxis = kepler6[0];
    float eccentricity = kepler6[1];
    float inclination = kepler6[2];
    float aNodeLongitude = kepler6[3];
    float periapsisArg=  kepler6[4];
    float trueAnomaly = kepler6[5];

    float r0 = smAxis*(1-eccentricity*eccentricity)/(1+(eccentricity*cos(trueAnomaly)));
    float r1[3] ={r0*cosf(trueAnomaly), r0*sinf(trueAnomaly), 0};

    float gravConst = 398600.4418; 

    float h = sqrtf(gravConst * smAxis * (1 - eccentricity * eccentricity));

    float qPerapsisArg[4] = {cosf(periapsisArg / 2), 0, 0, sinf(periapsisArg / 2)};
    float qInclination[4] = {cosf(inclination/2), sinf(inclination)/2, 0, 0};
    float qLongitude[4] = {cosf(aNodeLongitude / 2), 0, 0, sinf(aNodeLongitude / 2)};

    float qInclinationPeriapsis[4];
    quat_multiply(qInclination, qPerapsisArg, qInclinationPeriapsis);
    float qComb[4];
    quat_multiply(qLongitude, qInclinationPeriapsis, qComb);

    quat_apply(qComb, r1, eci);
}