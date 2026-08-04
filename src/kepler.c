#include "include/kepler.h"
#include "math.h"
#include "arm_math.h"

/**
 * \fn propogateOrbitalElements
 *
 * \brief Propogate orbital elements using Kepler's equation
 * 
 * \param[in] semi_major_axis: Semi major axis of the orbit in meters
 * \param[in] eccentricity: Eccentricity of the orbit   
 * \param[in] inclination: Inclination of the orbit in degrees
 * \param[in] ascending_node: Longitude of the ascending node in degrees
 * \param[in] periapsis: Argument of periapsis in degrees
 * \param[in] true_anomaly: True anomaly at time of epoch in degrees
 * \param[in] time_delta: Time delta in seconds to propogate the orbit
 * \param[out] output: Array of size 6 to hold the updated orbital elements
 *
 * \return void
 */
void propogateOrbitalElements(float32_t semi_major_axis, 
                                float32_t eccentricity, 
                                float32_t inclination, 
                                float32_t ascending_node, 
                                float32_t periapsis, 
                                float32_t true_anomaly, 
                                float32_t time_delta, 
                                float32_t *output) {
    //Mass of earth is extremely larger than the satellite, so:

    //gravitational  
    float32_t u = 3.986004418e14; //[m^3/s^2]
    float32_t f = true_anomaly * M_PI / 180.0; //[radians]

    
    //convert true anomaly to eccentric anamoly
    float32_t sqrt_1_minus_e, sqrt_1_plus_e;
    arm_sqrt_f32(1 - eccentricity, &sqrt_1_minus_e);
    arm_sqrt_f32(1 + eccentricity, &sqrt_1_plus_e);
    float32_t E_atan2;
    arm_atan2_f32(sqrt_1_minus_e * arm_sin_f32(f/2), sqrt_1_plus_e * arm_cos_f32(f/2), &E_atan2);
    float32_t E = 2 * E_atan2;

    // atan2( sqrt(1 - e) * sin(f/2), sqrt(1 + e) * cos(f/2));

    //Get current true anomaly from Kepler's equation
    float32_t M = E - eccentricity * arm_sin_f32(E);

    //define mean motion
    float32_t n;
    arm_sqrt_f32(u/pow(semi_major_axis, 3.0), &n);

    //Find mean anomaly at T + dt
    float32_t M_new = M + n * time_delta;

    //Use newton's method to iterate until we find the new_true_anomly. 
    float32_t E_current = M_new;
    float32_t tolerance = 1e-12;
    int MAXIMUM_ITERATIONS = 10;

    for (int i = 1; i < MAXIMUM_ITERATIONS; ++i) {
        //Calculate the derivative of kepler's equation
        float32_t E_change = ((E_current - eccentricity * arm_sin_f32(E_current) - M_new)
                        /(1 - eccentricity * arm_cos_f32(E_current)));
    
        //calculate next iteration
        E_current = E_current - E_change;

        //if tolerance is small enough, break out of the loop
        if(fabs(E_change) < tolerance)
            break;
    }
    output[0] = semi_major_axis;
    output[1] = eccentricity;
    output[2] = inclination;
    output[3] = ascending_node;
    output[4] = periapsis;

    //updated true anomaly from our updated eccentric anomaly
    float32_t new_atan2;
    arm_atan2_f32(sqrt_1_plus_e * arm_sin_f32(E_current/2), sqrt_1_minus_e * arm_cos_f32(E_current/2), &new_atan2);
    float new_true_anomaly = 2 * new_atan2;
    output[5] = new_true_anomaly * 180.0 / M_PI; //convert back to degrees

}