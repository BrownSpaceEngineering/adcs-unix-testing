#define WMM_NMAX        12
#define WMM_DIM        (WMM_NMAX + 2)
#define WMM_EPOCH       2025.0f     /* WMM2025 reference epoch              */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stdint.h>

/* Targeted double promotion — double interface, double JD arithmetic */
double jd2year(double JD);
double gmst_from_jd(double JD);
void  ecef_to_geodetic(const double r_ecef[3],
                    double*lat_rad, double*lon_rad, double*alt_m);
void  synthesize_mag_field(double lat_rad, double lon_rad, double alt_m,
                         const double g[WMM_DIM][WMM_DIM],
                         const double h[WMM_DIM][WMM_DIM],
                         double*B_N, double*B_E, double*B_D);
void  load_wmm2025(double JD,
                  double g[WMM_DIM][WMM_DIM],
                  double h[WMM_DIM][WMM_DIM]);
void  wmm_eci_embedded_v2(const double r_ECI[3],
                         int32_t JD_int, double JD_frac,
                         double B_ECI[3]);