#include "include/orbital2eci.h"
#include "Include/dsp/quaternion_math_functions.h"
#include "include/quat.h"
#include "arm_math.h"
#include "math.h"

/**
 * \fn orbitalToECI
 * 
 * \brief Converts the 6 orbital elements into ECI position vector
 * 
 * \warning All angles must be in radians
 * \note Units are consistent between semi-major axis and the resultant vector
 *      position in km, velocity in km/s
 * 
 * \note r0 is the current radius in the perifocal frame (plane of orbit)
 *      since the orbit is not purely circular, it gives the current distance
 * 
 * \note r1 is the position in the perifocal frame adjusted for the offset of the true anomaly
 * \note Axes: periapsis, 90 deg, normal
 * 
 * Parametres:
 *  \param[in] kepler6 : 6-element array of orbital elements [a, e, i, Ω, ω, ν] [km, dimensionless, rad, rad, rad, rad]
 *  \param[out] eci     : 3-element output array for ECI position vector [x, y, z] [km]
 */
void orbital_to_eci(float32_t *kepler6, float32_t *eci){
    // % Axes: periapsis, 90 deg, normal

    float32_t sm_axis = kepler6[0];
    float32_t eccentricity = kepler6[1];
    float32_t inclination = kepler6[2];
    float32_t a_node_longitude = kepler6[3];
    float32_t periapsis_arg=  kepler6[4];
    float32_t true_anomaly = kepler6[5];
    
    float32_t r0 = sm_axis*(1-eccentricity*eccentricity)/(1+(eccentricity*arm_cos_f32(true_anomaly)));
    float32_t r1[3] ={r0*arm_cos_f32(true_anomaly), r0*arm_sin_f32(true_anomaly), 0};

    float32_t grav_const = 398600.4418; 

    float32_t h = sqrtf(grav_const * sm_axis * (1 - eccentricity * eccentricity));

    float32_t q_perapsis_arg[4] = {arm_cos_f32(periapsis_arg / 2), 0, 0, arm_sin_f32(periapsis_arg / 2)};
    float32_t q_inclination[4] = {arm_cos_f32(inclination/2), 0.5f * arm_sin_f32(inclination), 0, 0};
    float32_t q_longitude[4] = {arm_cos_f32(a_node_longitude / 2), 0, 0, arm_sin_f32(a_node_longitude / 2)};

    float32_t q_inclination_periapsis[4];
    float32_t q_comb[4];

    arm_quaternion_product_f32(q_inclination, q_perapsis_arg, q_inclination_periapsis, 1);
    arm_quaternion_product_f32(q_longitude, q_inclination_periapsis, q_comb, 1);

    quat_apply(q_comb, r1, eci);
}