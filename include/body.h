#include "stdbool.h"
#include "magnetosphere.h"
#ifndef BODY
#define BODY

#ifndef NULLPTR
#define NULLPTR 0x0
#endif
extern float PVD_ECEF[3];
extern const float DETUMBLING_MAG_THRESHOLD;
extern const float DETUMBLING_GYRO_THRESHOLD;
extern float Imax[3];
extern const float INIT_P[6 * 6];
extern const float Q[6 * 6];
extern const float R_IN_SUN[6 * 6];
extern const float R_IN_SHADOW[6 * 6];

static float kepler_posn[6];
static float estimated_quat[4];
static float estimated_gyro_bias[3];
static float error_quat_state[6];
static float error_quat_cov[6 * 6];

static float last_ref_magnetosphere[3];

static bool last_ref_magnetosphere_initialized = false;
static bool pointing = false;
static bool posn_initialized = false;
static bool attitude_initialized = false;

void body(float* last_magnetometer_measurements, //1x3
    float* magnetometer_measurements,  //1x3
    float* gyro_measurements, //1x3, radians
    float* photodiode_measurements, //1xNUM_DIODES, raw readings
    float* posn_update,//1x6, may be NULLPTR
    int unix,//unix time
    int jd_scalar,
    float jd_frac,
    float dt,//time since last update
    float* output_currents);//1x3
#endif