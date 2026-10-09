#include "include/orbital2eci.h"
#include "include/quat.h"
#include "arm_math.h"
#include "math.h"

/**
 * \fn orbital_to_eci_posvel
 *
 * \brief Port of orbitalToECI.m: converts the 6 orbital elements into ECI position and velocity.
 *
 * Rotation perifocal->ECI is qZ(RAAN) * qX(i) * qZ(argp) (Hamilton, active).
 *
 * The previous version built the inclination quaternion as [cos(i/2), 0.5 sin(i), 0, 0]
 * instead of [cos(i/2), sin(i/2), 0, 0], which is not a unit quaternion and gives the wrong
 * position for any inclined orbit.
 *
 * \param[in] kepler6 [a, e, i, RAAN, argp, nu], angles in radians
 * \param[in] mu gravitational parameter in a's units
 * \param[out] r_eci position (a's unit)
 * \param[out] v_eci velocity (a's unit per second)
 */
void orbital_to_eci_posvel(const float32_t* kepler6, float32_t mu, float32_t* r_eci, float32_t* v_eci) {
    float32_t sm_axis = kepler6[0];
    float32_t eccentricity = kepler6[1];
    float32_t inclination = kepler6[2];
    float32_t a_node_longitude = kepler6[3];
    float32_t periapsis_arg = kepler6[4];
    float32_t true_anomaly = kepler6[5];

    float32_t cos_nu = cosf(true_anomaly);
    float32_t sin_nu = sinf(true_anomaly);
    float32_t p = sm_axis * (1 - eccentricity * eccentricity);

    // Perifocal frame. Axes: periapsis, 90 deg, normal
    float32_t r0 = p / (1 + eccentricity * cos_nu);
    float32_t r1[3] = {r0 * cos_nu, r0 * sin_nu, 0};

    float32_t h = sqrtf(mu * p);
    float32_t v1[3] = {-mu / h * sin_nu, mu / h * (eccentricity + cos_nu), 0};

    float32_t q_periapsis_arg[4] = {cosf(periapsis_arg / 2), 0, 0, sinf(periapsis_arg / 2)};
    float32_t q_inclination[4] = {cosf(inclination / 2), sinf(inclination / 2), 0, 0};
    float32_t q_longitude[4] = {cosf(a_node_longitude / 2), 0, 0, sinf(a_node_longitude / 2)};

    float32_t q_inclination_periapsis[4];
    float32_t q_comb[4];
    quat_multiply(q_inclination, q_periapsis_arg, q_inclination_periapsis);
    quat_multiply(q_longitude, q_inclination_periapsis, q_comb);

    if (r_eci != NULL) {
        quat_apply(q_comb, r1, r_eci);
    }
    if (v_eci != NULL) {
        quat_apply(q_comb, v1, v_eci);
    }
}

/**
 * \fn orbital_to_eci
 *
 * \brief Position-only version of orbital_to_eci_posvel (angles in radians).
 */
void orbital_to_eci(const float32_t* kepler6, float32_t* eci) {
    orbital_to_eci_posvel(kepler6, MU_EARTH_M3_S2, eci, NULL);
}
