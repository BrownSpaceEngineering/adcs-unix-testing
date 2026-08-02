#include "stdbool.h"
#include "magnetosphere.h"
#ifndef BODY
#define BODY

#ifndef NULLPTR
#define NULLPTR 0x0
#endif
extern double PVD_ECEF[3];
extern const double DETUMBLING_MAG_THRESHOLD;
extern const double DETUMBLING_GYRO_THRESHOLD;
extern double Imax[3];
extern const double INIT_P[6 * 6];
extern const double Q[6 * 6];
extern const double R_IN_SUN[6 * 6];
extern const double R_IN_SHADOW[3 * 3];

static double kepler_posn[6];
static double estimated_quat[4];
static double error_quat_state[6];
static double error_quat_cov[6 * 6];

static bool pointing = false;
static bool posn_initialized = false;
static int  filter_step = 0;

void body(double* last_magnetometer_measurements, //1x3
    double* magnetometer_measurements,  //1x3
    double* gyro_measurements, //1x3, radians
    double* photodiode_measurements, //1xNUM_DIODES, raw readings
    double* posn_update,//1x6, may be NULLPTR
    int unix,//unix time
    int jd_scalar,
    double jd_frac,
    double dt,//time since last update
    double* output_currents);//1x3
#endif