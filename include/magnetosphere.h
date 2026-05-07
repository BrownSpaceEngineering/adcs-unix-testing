#define WMM_NMAX        12
#define WMM_DIM        (WMM_NMAX + 2)
#define WMM_EPOCH       2025.0f     /* WMM2025 reference epoch              */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stdint.h>

/* Targeted double promotion — float interface, double JD arithmetic */
float jd2year(float64_t JD);
float gmst_from_jd(float64_t JD);
void  ecef_to_geodetic(const float32_t r_ecef[3],
                    float32_t *lat_rad, float32_t *lon_rad, float32_t *alt_m);
void  synthesize_mag_field(float32_t lat_rad, float32_t lon_rad, float32_t alt_m,
                         const float32_t g[WMM_DIM][WMM_DIM],
                         const float32_t h[WMM_DIM][WMM_DIM],
                         float32_t *B_N, float32_t *B_E, float32_t *B_D);
void  load_wmm2025(float64_t JD,
                  float32_t g[WMM_DIM][WMM_DIM],
                  float32_t h[WMM_DIM][WMM_DIM]);
void wmm_eci_embedded_v2(const float32_t r_ECI[3],
                        int32_t JD_int, float32_t JD_frac,
                        float32_t B_ECI[3]);