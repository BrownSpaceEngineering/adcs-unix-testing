#ifndef ADCS_ALGORITHMS_H
#define ADCS_ALGORITHMS_H

#include <stdbool.h>

/**
 * @brief B-Dot detumbling algorithm.
 * 
 * @param mag_t0 Magnetometer reading at current time [3].
 * @param mag_t1 Magnetometer reading at previous time [3].
 * @param dt Time step between readings.
 * @param k Gain.
 * @param moments Output magnetic moments [3].
 */
void bdot_detumbling(const float* mag_t0, const float* mag_t1, float dt, float k, float* moments);

/**
 * @brief Calculate the quaternion to point the satellite's z-axis towards a target (Nadir).
 * 
 * @param sat_pos_eci Satellite position in ECI [3].
 * @param q_world Current satellite orientation quaternion (WXYZ).
 * @param target_pos_eci Target position in ECI (e.g., Providence) [3].
 * @param q_down Output target quaternion (WXYZ).
 */
void down_quaternion(const float* sat_pos_eci, const float* q_world, const float* target_pos_eci, float* q_down);

/**
 * @brief Calculate the Sun vector in ECI frame for a given Julian Date.
 * 
 * @param jd Julian Date.
 * @param sun_vec Output Sun vector [3].
 */
void sun_vector_eci(double jd, float* sun_vec);

/**
 * @brief Convert orbital elements to ECI position and velocity.
 * 
 * @param smAxis Semi-major axis (km).
 * @param eccentricity Eccentricity.
 * @param inclination Inclination (rad).
 * @param aNodeLongitude Right Ascension of Ascending Node (rad).
 * @param periapsisArg Argument of Periapsis (rad).
 * @param trueAnomaly True Anomaly (rad).
 * @param pos Output position [3].
 * @param vel Output velocity [3].
 */
void orbital_to_eci(float smAxis, float eccentricity, float inclination, 
                    float aNodeLongitude, float periapsisArg, float trueAnomaly,
                    float* pos, float* vel);

/**
 * @brief PD control loop to calculate required magnetic moment.
 * 
 * @param q_error_aa Error vector in axis-angle format (axis*angle) [3].
 * @param omega_aa Angular velocity in axis-angle format [3].
 * @param B Magnetic field vector in body frame [3].
 * @param Kp Proportional gain [3].
 * @param Kd Derivative gain [3].
 * @param m Output magnetic moment [3].
 */
void pd_control_loop(const float* q_error_aa, const float* omega_aa, const float* B, 
                     const float* Kp, const float* Kd, float* m);

/**
 * @brief World Magnetic Model (WMM) to get magnetic field in ECI.
 * 
 * @param r_eci Position in ECI [3].
 * @param jd Julian Date.
 * @param B_eci Output magnetic field in ECI [3].
 */
void wmm_eci(const float* r_eci, double jd, float* B_eci);

/**
 * @brief Multi-sensor Unscented Kalman Filter (MUKF) iteration.
 * 
 * @param state Current state [6] (3 for error rotation, 3 for gyro bias).
 * @param q_est Current orientation estimate (WXYZ).
 * @param P Current covariance matrix [6x6].
 * @param measurement Sensor measurement (e.g., normalized magnetometer) [3].
 * @param gyro_meas Gyroscope measurement [3].
 * @param B_ref Reference vector in ECI (e.g., WMM prediction) [3].
 * @param dt Time step.
 * @param next_state Output next state [6].
 * @param next_q_est Output next orientation estimate (WXYZ).
 * @param next_P Output next covariance matrix [6x6].
 */
void mukf_iterate(float* state, float* q_est, float* P, 
                  const float* measurement, const float* gyro_meas, const float* B_ref, 
                  float dt, float* next_state, float* next_q_est, float* next_P);

/**
 * @brief Convert magnetic moment to coil currents.
 * 
 * @param m Magnetic moment [3].
 * @param Imax Maximum current.
 * @param I Output currents [3].
 */
void moment_to_current(const float* m, float Imax, float* I);

/**
 * @brief Convert ECEF coordinates to ECI.
 * 
 * @param r_ecef Position in ECEF [3].
 * @param unix_seconds Unix timestamp.
 * @param r_eci Output position in ECI [3].
 */
void ecef_to_eci(const float* r_ecef, double unix_seconds, float* r_eci);

/**
 * @brief Propagate orbital elements using two-body Keplerian model.
 * 
 * @param semi_major_axis Semi-major axis (m).
 * @param eccentricity Eccentricity.
 * @param true_anomaly_deg True anomaly (degrees).
 * @param time_delta Time step (seconds).
 * @param new_true_anomaly_deg Output new true anomaly (degrees).
 */
void propagate_orbital_elements(float semi_major_axis, float eccentricity, 
                                float true_anomaly_deg, float time_delta, 
                                float* new_true_anomaly_deg);

#endif // ADCS_ALGORITHMS_H
