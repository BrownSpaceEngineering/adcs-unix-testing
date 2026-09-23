#ifndef BODY
#define BODY
#include <stdbool.h>

#ifndef NULLPTR
#define NULLPTR 0x0
#endif

extern float PVD_ECEF[3];
extern const float DETUMBLING_MAG_THRESHOLD;
extern const float DETUMBLING_GYRO_THRESHOLD;
extern const float BDOT_GAIN;
extern float Imax[3];
extern const float INIT_P[6 * 6];
extern const float Q_RATE[6 * 6];
extern const float R_IN_SUN[6 * 6];
extern const float R_IN_SHADOW[3 * 3];

/**
 * One ADCS step.
 *
 * Units / conventions:
 *   magnetometer: body frame, TESLA (same units as the WMM model and torque_2_moments)
 *   gyro:         body frame, rad/s
 *   photodiodes:  NUM_DIODES raw readings in PHOTODIODES order
 *   posn_update:  [a (m), e, i, RAAN, argp, nu (deg)] or NULLPTR to propagate the last one
 *   unix_time / jd_scalar + jd_frac: the same instant, UTC
 *   output_currents: magnetorquer currents (A); always written (zeros when idle)
 */
void body(const float* last_magnetometer_measurements, // 1x3
          const float* magnetometer_measurements,      // 1x3
          const float* gyro_measurements,              // 1x3, radians
          const float* photodiode_measurements,        // 1xNUM_DIODES, raw readings
          const float* posn_update,                    // 1x6, may be NULLPTR
          int unix_time,                               // unix time
          int jd_scalar,
          float jd_frac,
          float dt,                                    // time since last update
          float* output_currents);                     // 1x3

// Clears all filter / mode state (power-on state)
void body_reset(void);
bool body_is_pointing(void);
// Current body->ECI estimate; returns false if the attitude hasn't been initialized
bool body_get_attitude(float* q_body_to_eci);
void body_get_gyro_bias(float* bias);
// Number of times the filter has been declared failed and reset to detumbling
int body_get_filter_resets(void);
#endif
