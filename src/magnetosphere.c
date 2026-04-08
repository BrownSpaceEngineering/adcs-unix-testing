/**
 * wmm_eci.c
 *
 * C translation of wmm_eci_embedded.m (WMM2025, ECI output).
 * gcc -O2 -Wall -Wextra -o magnetosphere_test magnetosphere.c -lm
 * Sample Inputs are from the matlab function wrldmagm that we're trying to mimic
 *
 * All internal computation uses double precision for comparative testing
 * against magnetosphere_update.c (which uses targeted double promotion).
 */

#ifdef MAG_SPHERE_ORIGINAL

#include "magnetosphere.h"
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

void wmm_eci_embedded(const double r_ECI[3], double JD, double B_ECI[3]);

/* =========================================================================
 * Internal constants
 * ========================================================================= */

/* WGS-84 / geomagnetic reference constants */
#define WGS84_A         6378137.0
#define WGS84_F         (1.0 / 298.257223563)
#define WGS84_B         (WGS84_A * (1.0 - WGS84_F))
#define WGS84_E2        (WGS84_F * (2.0 - WGS84_F))
#define WGS84_EP2       (WGS84_E2 / (1.0 - WGS84_E2))
#define GEO_REF_RADIUS  6371200.0

/* =========================================================================
 * WMM2025 embedded coefficients
 * Columns: n  m  g(nT)  h(nT)  gdot(nT/yr)  hdot(nT/yr)
 * ========================================================================= */

typedef struct {
    int32_t n;
    int32_t m;
    double  g;
    double  h;
    double  gdot;
    double  hdot;
} WmmRecord;

static const WmmRecord WMM2025_DATA[] = {
    { 1,  0, -29351.8,     0.0,    12.0,    0.0 },
    { 1,  1,  -1410.8,  4545.4,     9.7,  -21.5 },
    { 2,  0,  -2556.6,     0.0,   -11.6,    0.0 },
    { 2,  1,   2951.1, -3133.6,    -5.2,  -27.7 },
    { 2,  2,   1649.3,  -815.1,    -8.0,  -12.1 },
    { 3,  0,   1361.0,     0.0,    -1.3,    0.0 },
    { 3,  1,  -2404.1,   -56.6,    -4.2,    4.0 },
    { 3,  2,   1243.8,   237.5,     0.4,   -0.3 },
    { 3,  3,    453.6,  -549.5,   -15.6,   -4.1 },
    { 4,  0,    895.0,     0.0,    -1.6,    0.0 },
    { 4,  1,    799.5,   278.6,    -2.4,   -1.1 },
    { 4,  2,     55.7,  -133.9,    -6.0,    4.1 },
    { 4,  3,   -281.1,   212.0,     5.6,    1.6 },
    { 4,  4,     12.1,  -375.6,    -7.0,   -4.4 },
    { 5,  0,   -233.2,     0.0,     0.6,    0.0 },
    { 5,  1,    368.9,    45.4,     1.4,   -0.5 },
    { 5,  2,    187.2,   220.2,     0.0,    2.2 },
    { 5,  3,   -138.7,  -122.9,     0.6,    0.4 },
    { 5,  4,   -142.0,    43.0,     2.2,    1.7 },
    { 5,  5,     20.9,   106.1,     0.9,    1.9 },
    { 6,  0,     64.4,     0.0,    -0.2,    0.0 },
    { 6,  1,     63.8,   -18.4,    -0.4,    0.3 },
    { 6,  2,     76.9,    16.8,     0.9,   -1.6 },
    { 6,  3,   -115.7,    48.8,     1.2,   -0.4 },
    { 6,  4,    -40.9,   -59.8,    -0.9,    0.9 },
    { 6,  5,     14.9,    10.9,     0.3,    0.7 },
    { 6,  6,    -60.7,    72.7,     0.9,    0.9 },
    { 7,  0,     79.5,     0.0,    -0.0,    0.0 },
    { 7,  1,    -77.0,   -48.9,    -0.1,    0.6 },
    { 7,  2,     -8.8,   -14.4,    -0.1,    0.5 },
    { 7,  3,     59.3,    -1.0,     0.5,   -0.8 },
    { 7,  4,     15.8,    23.4,    -0.1,    0.0 },
    { 7,  5,      2.5,    -7.4,    -0.8,   -1.0 },
    { 7,  6,    -11.1,   -25.1,    -0.8,    0.6 },
    { 7,  7,     14.2,    -2.3,     0.8,   -0.2 },
    { 8,  0,     23.2,     0.0,    -0.1,    0.0 },
    { 8,  1,     10.8,     7.1,     0.2,   -0.2 },
    { 8,  2,    -17.5,   -12.6,     0.0,    0.5 },
    { 8,  3,      2.0,    11.4,     0.5,   -0.4 },
    { 8,  4,    -21.7,    -9.7,    -0.1,    0.4 },
    { 8,  5,     16.9,    12.7,     0.3,   -0.5 },
    { 8,  6,     15.0,     0.7,     0.2,   -0.6 },
    { 8,  7,    -16.8,    -5.2,    -0.0,    0.3 },
    { 8,  8,      0.9,     3.9,     0.2,    0.2 },
    { 9,  0,      4.6,     0.0,    -0.0,    0.0 },
    { 9,  1,      7.8,   -24.8,    -0.1,   -0.3 },
    { 9,  2,      3.0,    12.2,     0.1,    0.3 },
    { 9,  3,     -0.2,     8.3,     0.3,   -0.3 },
    { 9,  4,     -2.5,    -3.3,    -0.3,    0.3 },
    { 9,  5,    -13.1,    -5.2,     0.0,    0.2 },
    { 9,  6,      2.4,     7.2,     0.3,   -0.1 },
    { 9,  7,      8.6,    -0.6,    -0.1,   -0.2 },
    { 9,  8,     -8.7,     0.8,     0.1,    0.4 },
    { 9,  9,    -12.9,    10.0,    -0.1,    0.1 },
    {10,  0,     -1.3,     0.0,     0.1,    0.0 },
    {10,  1,     -6.4,     3.3,     0.0,    0.0 },
    {10,  2,      0.2,     0.0,     0.1,   -0.0 },
    {10,  3,      2.0,     2.4,     0.1,   -0.2 },
    {10,  4,     -1.0,     5.3,    -0.0,    0.1 },
    {10,  5,     -0.6,    -9.1,    -0.3,   -0.1 },
    {10,  6,     -0.9,     0.4,     0.0,    0.1 },
    {10,  7,      1.5,    -4.2,    -0.1,    0.0 },
    {10,  8,      0.9,    -3.8,    -0.1,   -0.1 },
    {10,  9,     -2.7,     0.9,    -0.0,    0.2 },
    {10, 10,     -3.9,    -9.1,    -0.0,   -0.0 },
    {11,  0,      2.9,     0.0,     0.0,    0.0 },
    {11,  1,     -1.5,     0.0,    -0.0,   -0.0 },
    {11,  2,     -2.5,     2.9,     0.0,    0.1 },
    {11,  3,      2.4,    -0.6,     0.0,   -0.0 },
    {11,  4,     -0.6,     0.2,     0.0,    0.1 },
    {11,  5,     -0.1,     0.5,    -0.1,   -0.0 },
    {11,  6,     -0.6,    -0.3,     0.0,   -0.0 },
    {11,  7,     -0.1,    -1.2,    -0.0,    0.1 },
    {11,  8,      1.1,    -1.7,    -0.1,   -0.0 },
    {11,  9,     -1.0,    -2.9,    -0.1,    0.0 },
    {11, 10,     -0.2,    -1.8,    -0.1,    0.0 },
    {11, 11,      2.6,    -2.3,    -0.1,    0.0 },
    {12,  0,     -2.0,     0.0,     0.0,    0.0 },
    {12,  1,     -0.2,    -1.3,     0.0,   -0.0 },
    {12,  2,      0.3,     0.7,    -0.0,    0.0 },
    {12,  3,      1.2,     1.0,    -0.0,   -0.1 },
    {12,  4,     -1.3,    -1.4,    -0.0,    0.1 },
    {12,  5,      0.6,    -0.0,    -0.0,   -0.0 },
    {12,  6,      0.6,     0.6,     0.1,   -0.0 },
    {12,  7,      0.5,    -0.1,    -0.0,   -0.0 },
    {12,  8,     -0.1,     0.8,     0.0,    0.0 },
    {12,  9,     -0.4,     0.1,     0.0,   -0.0 },
    {12, 10,     -0.2,    -1.0,    -0.1,   -0.0 },
    {12, 11,     -1.3,     0.1,    -0.0,    0.0 },
    {12, 12,     -0.7,     0.2,    -0.1,   -0.1 },
};

#define WMM_NRECORDS  (int32_t)(sizeof(WMM2025_DATA) / sizeof(WMM2025_DATA[0]))

/**
 * \brief Compute the Earth's magnetic field vector at a satellite position in ECI.
 *
 * Evaluates the WMM2025 spherical harmonic model at the given position and
 * time, returning the field in the Earth-Centred Inertial (ECI) frame.
 * Pipeline: ECI -> ECEF (GMST rotation) -> geodetic (WGS-84) ->
 * NED (spherical harmonic synthesis) -> ECEF -> ECI.
 *
 * \param r_ECI  Satellite position in ECI frame (metres), length-3 array.
 * \param JD     Julian Date of the observation epoch.
 * \param B_ECI  Output magnetic field vector in ECI frame (Tesla), length-3 array.
 */
void wmm_eci_embedded(const double r_ECI[3], double JD, double B_ECI[3])
{
    /* --- Build time-extrapolated Gauss coefficients (stack allocated) --- */
    double g[WMM_DIM][WMM_DIM];
    double h[WMM_DIM][WMM_DIM];
    load_wmm2025(JD, g, h);

    /* --- ECI -> ECEF via GMST rotation (Rz) --- */
    double theta = gmst_from_jd(JD);
    double cos_t = cos(theta);
    double sin_t = sin(theta);

    double r_ecef[3];
    r_ecef[0] =  cos_t * r_ECI[0] + sin_t * r_ECI[1];
    r_ecef[1] = -sin_t * r_ECI[0] + cos_t * r_ECI[1];
    r_ecef[2] =  r_ECI[2];

    /* --- ECEF -> geodetic (WGS-84) --- */
    double lat_rad, lon_rad, alt_m;
    ecef_to_geodetic(r_ecef, &lat_rad, &lon_rad, &alt_m);

    /* --- Synthesize field in NED frame (nT) --- */
    double B_N, B_E, B_D;
    synthesize_mag_field(lat_rad, lon_rad, alt_m, g, h, &B_N, &B_E, &B_D);

    /* --- NED -> ECEF --- */
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);
    double sin_lon = sin(lon_rad);
    double cos_lon = cos(lon_rad);

    double B_ecef[3];
    B_ecef[0] = (-sin_lat * cos_lon) * B_N + (-sin_lon) * B_E + (-cos_lat * cos_lon) * B_D;
    B_ecef[1] = (-sin_lat * sin_lon) * B_N + ( cos_lon) * B_E + (-cos_lat * sin_lon) * B_D;
    B_ecef[2] = ( cos_lat          ) * B_N +           0 * B_E + (-sin_lat           ) * B_D;

    /* --- ECEF -> ECI  (Rz^T = Rz(-theta)) --- */
    double B_eci_nT[3];
    B_eci_nT[0] =  cos_t * B_ecef[0] - sin_t * B_ecef[1];
    B_eci_nT[1] =  sin_t * B_ecef[0] + cos_t * B_ecef[1];
    B_eci_nT[2] =  B_ecef[2];

    /* --- nT -> T --- */
    B_ECI[0] = B_eci_nT[0] * 1.0e-9;
    B_ECI[1] = B_eci_nT[1] * 1.0e-9;
    B_ECI[2] = B_eci_nT[2] * 1.0e-9;
}

/**
 * \brief Load WMM2025 Gauss coefficients linearly extrapolated to a given epoch.
 *
 * Reads the embedded WMM2025_DATA table and applies secular-variation rates
 * to produce time-adjusted g[n][m] and h[n][m] arrays.
 *
 * \param JD  Julian Date of the target epoch.
 * \param g   Output cosine Gauss coefficients g[n][m] in nT (WMM_DIM x WMM_DIM).
 * \param h   Output sine Gauss coefficients h[n][m] in nT (WMM_DIM x WMM_DIM).
 */
void load_wmm2025(double JD,
                 double g[WMM_DIM][WMM_DIM],
                 double h[WMM_DIM][WMM_DIM])
{
    double decimalYear = jd2year(JD);
    double dt          = decimalYear - WMM_EPOCH;

    memset(g, 0, sizeof(double) * WMM_DIM * WMM_DIM);
    memset(h, 0, sizeof(double) * WMM_DIM * WMM_DIM);

    for (int32_t i = 0; i < WMM_NRECORDS; ++i) {
        int32_t n = WMM2025_DATA[i].n;
        int32_t m = WMM2025_DATA[i].m;

        g[n][m] = WMM2025_DATA[i].g + dt * WMM2025_DATA[i].gdot;
        h[n][m] = WMM2025_DATA[i].h + dt * WMM2025_DATA[i].hdot;
    }
}

/**
 * \brief Convert a Julian Date to a decimal year.
 *
 * Uses the J2000.0 epoch (JD 2451545.0 = 2000 Jan 1.5) as the reference.
 *
 * \param JD  Julian Date.
 * \return    Corresponding decimal year (e.g. 2025.5).
 */
double jd2year(double JD)
{
    return 2000.0 + (JD - 2451545.0) / 365.25;
}

/**
 * \brief Compute Greenwich Mean Sidereal Time (GMST) from a Julian Date.
 *
 * Evaluates the IAU polynomial in Julian centuries T referenced to J2000.0.
 * The result is the rotation angle between the ECI and ECEF frames about
 * the Z-axis.
 *
 * \param JD  Julian Date.
 * \return    GMST angle in radians, in [0, 2π).
 */
double gmst_from_jd(double JD)
{
    double T        = (JD - 2451545.0) / 36525.0;
    double gmst_sec = 67310.54841
                    + (876600.0 * 3600.0 + 8640184.812866) * T
                    + 0.093104 * T * T
                    - 6.2e-6   * T * T * T;
    double gmst_deg = fmod(gmst_sec / 240.0, 360.0);
    if (gmst_deg < 0.0) gmst_deg += 360.0;
    return gmst_deg * M_PI / 180.0;
}

/**
 * \brief Convert ECEF Cartesian coordinates to WGS-84 geodetic coordinates.
 *
 * Uses the Bowring iterative method (up to 10 iterations).  Handles the
 * polar singularity (p ≈ 0) as a special case.
 *
 * \param r_ecef   Input position in ECEF frame (metres), length-3 array [x, y, z].
 * \param lat_rad  Output geodetic latitude (radians), in [-π/2, π/2].
 * \param lon_rad  Output geodetic longitude (radians), in (-π, π].
 * \param alt_m    Output altitude above the WGS-84 ellipsoid (metres).
 */
void ecef_to_geodetic(const double r_ecef[3],
                   double *lat_rad, double *lon_rad, double *alt_m)
{
    const double a   = WGS84_A;
    const double b   = WGS84_B;
    const double e2  = WGS84_E2;
    const double ep2 = WGS84_EP2;

    double x = r_ecef[0];
    double y = r_ecef[1];
    double z = r_ecef[2];

    *lon_rad = atan2(y, x);

    double p = sqrt(x * x + y * y);

    /* Pole singularity */
    if (p < 1.0e-12) {
        *lon_rad = 0.0;
        if (z >= 0.0) {
            *lat_rad = M_PI / 2.0;
            *alt_m   = z - b;
        } else {
            *lat_rad = -M_PI / 2.0;
            *alt_m   = -z - b;
        }
        return;
    }

    /* Bowring initial estimate */
    double theta = atan2(z * a, p * b);

    double lat = 0.0, alt = 0.0;
    for (int32_t iter = 0; iter < 10; ++iter) {
        double sin_t = sin(theta);
        double cos_t = cos(theta);

        double num = z + ep2 * b * sin_t * sin_t * sin_t;
        double den = p - e2  * a * cos_t * cos_t * cos_t;
        lat = atan2(num, den);

        double sin_lat = sin(lat);
        double cos_lat = cos(lat);
        double N = a / sqrt(1.0 - e2 * sin_lat * sin_lat);

        alt = p * cos_lat + z * sin_lat - a * sqrt(1.0 - e2 * sin_lat * sin_lat);

        double theta_new = atan2(z * (1.0 - e2 * N / (N + alt)), p);
        if (fabs(theta_new - theta) < 1.0e-10) {
            theta = theta_new;
            break;
        }
        theta = theta_new;
    }

    *lat_rad = lat;
    *alt_m   = alt;
}

/**
 * \brief Evaluate the WMM spherical harmonic expansion to produce the magnetic field in NED.
 *
 * Converts geodetic coordinates to geocentric, builds Schmidt semi-normalised
 * associated Legendre polynomials P[n][m] and their colatitude derivatives
 * dP[n][m] up to degree/order WMM_NMAX, sums the series in geocentric
 * spherical components (Br, Bt, Bphi), then rotates to the geodetic NED frame.
 *
 * \param lat_rad  Geodetic latitude (radians).
 * \param lon_rad  Geodetic longitude (radians).
 * \param alt_m    Altitude above the WGS-84 ellipsoid (metres).
 * \param g        Cosine Gauss coefficients g[n][m] in nT (from load_wmm2025).
 * \param h        Sine Gauss coefficients h[n][m] in nT (from load_wmm2025).
 * \param B_N      Output North component of the magnetic field (nT).
 * \param B_E      Output East component of the magnetic field (nT).
 * \param B_D      Output Down component of the magnetic field (nT).
 */
void synthesize_mag_field(double lat_rad, double lon_rad, double alt_m,
                        const double g[WMM_DIM][WMM_DIM],
                        const double h[WMM_DIM][WMM_DIM],
                        double *B_N, double *B_E, double *B_D)
{
    /* ---- Geodetic -> geocentric ---- */
    const double Re  = WGS84_A;
    const double f   = WGS84_F;
    double sin_lat   = sin(lat_rad);
    double cos_lat   = cos(lat_rad);
    double rc        = Re / sqrt(1.0 - (2.0 * f - f * f) * sin_lat * sin_lat);

    double p_gc = (rc + alt_m) * cos_lat;
    double z_gc = (rc * (1.0 - f) * (1.0 - f) + alt_m) * sin_lat;
    double r    = sqrt(p_gc * p_gc + z_gc * z_gc);

    double phi_prime = asin(z_gc / r);          /* geocentric latitude    */
    double psi       = lat_rad - phi_prime;     /* geodetic - geocentric  */

    /* ---- Spherical harmonic synthesis ---- */
    double theta     = M_PI / 2.0 - phi_prime;  /* geocentric colatitude  */
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);

    /* Stack-allocated Legendre arrays — WMM_DIM x WMM_DIM */
    double P [WMM_DIM][WMM_DIM];
    double dP[WMM_DIM][WMM_DIM];
    memset(P,  0, sizeof(P));
    memset(dP, 0, sizeof(dP));

    /* Seed values (n=0,m=0) and (n=1) — stored 0-indexed in P/dP */
    P [0][0] = 1.0;
    P [1][0] = cos_theta;           /* P(1,0) */
    P [1][1] = sin_theta;           /* P(1,1) */
    dP[0][0] = 0.0;
    dP[1][0] = -sin_theta;          /* dP(1,0)/dtheta */
    dP[1][1] =  cos_theta;          /* dP(1,1)/dtheta */

    for (int32_t n = 2; n <= WMM_NMAX; ++n) {
        for (int32_t m = 0; m <= n; ++m) {

            if (m == n) {
                /* Diagonal recurrence */
                double factor = sqrt(1.0 - 1.0 / (2.0 * (double)n));
                P [n][m] = sin_theta * P [n-1][m-1] * factor;
                dP[n][m] = (sin_theta * dP[n-1][m-1] + cos_theta * P[n-1][m-1]) * factor;
            } else {
                /* Off-diagonal recurrence */
                double n2    = (double)(n * n);
                double m2    = (double)(m * m);
                double nm1_2 = (double)((n-1)*(n-1));
                double denom = sqrt(n2 - m2);
                double K     = (2.0 * (double)n - 1.0) / denom;
                double M_val = sqrt(nm1_2 - m2)        / denom;

                P [n][m] = K * cos_theta * P [n-1][m] - M_val * P [n-2][m];
                dP[n][m] = K * (cos_theta * dP[n-1][m] - sin_theta * P[n-1][m])
                         - M_val * dP[n-2][m];
            }
        }
    }

    /* ---- Sum up the series ---- */
    double Br   = 0.0;   /* radial    (outward +)   */
    double Bt   = 0.0;   /* theta     (southward +) */
    double Bphi = 0.0;   /* phi       (eastward +)  */

    const double a_ref = GEO_REF_RADIUS;

    for (int32_t n = 1; n <= WMM_NMAX; ++n) {
        double ar = a_ref / r;
        double ratio = ar;
        for (int32_t k = 0; k < n + 1; ++k) ratio *= ar;   /* ar^(n+2) */

        for (int32_t m = 0; m <= n; ++m) {
            double gm = g[n][m];
            double hm = h[n][m];

            double sin_mlon = sin((double)m * lon_rad);
            double cos_mlon = cos((double)m * lon_rad);

            double xy_comp = gm * cos_mlon + hm * sin_mlon;
            double z_comp  = gm * sin_mlon - hm * cos_mlon;

            Br += (double)(n + 1) * ratio * xy_comp * P[n][m];
            Bt += ratio * xy_comp * dP[n][m];

            if (m > 0 && fabs(sin_theta) > 1.0e-10) {
                Bphi += ratio * (double)m * z_comp * P[n][m] / sin_theta;
            }
        }
    }

    /* ---- Geocentric (r, theta, phi) -> geodetic NED ---- */
    double B_X_geo = -Bt;    /* North, geocentric */
    double B_Z_geo = -Br;    /* Down,  geocentric */

    double cos_psi = cos(psi);
    double sin_psi = sin(psi);

    *B_N = -(B_X_geo * cos_psi + B_Z_geo * sin_psi);
    *B_E =  Bphi;
    *B_D =  (B_X_geo * sin_psi + B_Z_geo * cos_psi);
}

#endif /* MAG_SPHERE_ORIGINAL */
