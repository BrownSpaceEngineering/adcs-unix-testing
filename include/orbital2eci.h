#ifndef ORBITAL2ECI
#define ORBITAL2ECI
#include "arm_math.h"

#define MU_EARTH_M3_S2 3.986004418e14f

// kepler6 = [a, e, i, RAAN, argp, nu], angles in RADIANS, a in any length unit.
// eci gets the position in the same unit as a.
void orbital_to_eci(const float32_t* kepler6, float32_t* eci);
// Same, plus velocity. mu must match a's unit (MU_EARTH_M3_S2 for metres -> m/s; the MATLAB
// orbitalToECI.m uses km and 398600.4418 km^3/s^2).
void orbital_to_eci_posvel(const float32_t* kepler6, float32_t mu, float32_t* r_eci, float32_t* v_eci);
#endif
