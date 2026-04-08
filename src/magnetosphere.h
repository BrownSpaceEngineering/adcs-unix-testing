#define WMM_NMAX        12
#define WMM_DIM        (WMM_NMAX + 2)
#define WMM_EPOCH       2025.0f     /* WMM2025 reference epoch              */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stdint.h>

#ifdef MAG_SPHERE_ORIGINAL
/* Full double implementation */
double jd2year(double JD);
double gmst_from_jd(double JD);
void   ecef_to_geodetic(const double r_ecef[3],
                     double *lat_rad, double *lon_rad, double *alt_m);
void   synthesize_mag_field(double lat_rad, double lon_rad, double alt_m,
                          const double g[WMM_DIM][WMM_DIM],
                          const double h[WMM_DIM][WMM_DIM],
                          double *B_N, double *B_E, double *B_D);
void   load_wmm2025(double JD,
                   double g[WMM_DIM][WMM_DIM],
                   double h[WMM_DIM][WMM_DIM]);
void   wmm_eci_embedded(const double r_ECI[3], double JD, double B_ECI[3]);
#elif defined(MAG_SPHERE_UPDATE)
/* Targeted double promotion — float interface, double JD arithmetic */
float jd2year(double JD);
float gmst_from_jd(double JD);
void  ecef_to_geodetic(const float r_ecef[3],
                    float *lat_rad, float *lon_rad, float *alt_m);
void  synthesize_mag_field(float lat_rad, float lon_rad, float alt_m,
                         const float g[WMM_DIM][WMM_DIM],
                         const float h[WMM_DIM][WMM_DIM],
                         float *B_N, float *B_E, float *B_D);
void  load_wmm2025(double JD,
                  float g[WMM_DIM][WMM_DIM],
                  float h[WMM_DIM][WMM_DIM]);
void  wmm_eci_embedded_v2(const float r_ECI[3],
                         int32_t JD_int, float JD_frac,
                         float B_ECI[3]);
#endif

void test_magnetosphere(void);
