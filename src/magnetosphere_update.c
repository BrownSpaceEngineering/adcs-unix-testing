/**
 * magnetosphere_update.c
 *
 * Precision-improved version of magnetosphere.c.
 *
 * Changes vs magnetosphere.c:
 *   Fix 1 – gmst_from_jd() now computes entirely in double to avoid
 *            float precision loss in the large intermediate GMST value
 *            and the subsequent fmod(~3.5e6, 360) call.
 *   Fix 2 – jd2year() now computes in double to avoid catastrophic
 *            cancellation when subtracting 2451545.0 from a large JD.
 *   Fix 3 – wmm_eci_embedded() interface split into (JD_int, JD_frac)
 *            so the caller's Julian Date is never truncated to float
 *            precision (~7 sig figs) before arithmetic.  Both helpers
 *            now accept a double JD reconstructed from those two parts.
 *
 * gcc -O2 -Wall -Wextra -o magnetosphere_update_test magnetosphere_update.c -lm
 */

#ifdef MAG_SPHERE_UPDATE

#include "magnetosphere.h"
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Updated public signature — JD split into integer day + fractional day */
void wmm_eci_embedded_v2(const float r_ECI[3],
                        int32_t JD_int, float JD_frac,
                        float B_ECI[3]);

/* =========================================================================
 * Internal constants
 * ========================================================================= */

#define WGS84_A         6378137.0f
#define WGS84_F         (1.0f / 298.257223563f)
#define WGS84_B         (WGS84_A * (1.0f - WGS84_F))
#define WGS84_E2        (WGS84_F * (2.0f - WGS84_F))
#define WGS84_EP2       (WGS84_E2 / (1.0f - WGS84_E2))
#define GEO_REF_RADIUS  6371200.0f

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =========================================================================
 * WMM2025 embedded coefficients
 * Columns: n  m  g(nT)  h(nT)  gdot(nT/yr)  hdot(nT/yr)
 * ========================================================================= */

typedef struct {
    int32_t n;
    int32_t m;
    float   g;
    float   h;
    float   gdot;
    float   hdot;
} WmmRecord;

static const WmmRecord WMM2025_DATA[] = {
    { 1,  0, -29351.8f,     0.0f,    12.0f,    0.0f },
    { 1,  1,  -1410.8f,  4545.4f,     9.7f,  -21.5f },
    { 2,  0,  -2556.6f,     0.0f,   -11.6f,    0.0f },
    { 2,  1,   2951.1f, -3133.6f,    -5.2f,  -27.7f },
    { 2,  2,   1649.3f,  -815.1f,    -8.0f,  -12.1f },
    { 3,  0,   1361.0f,     0.0f,    -1.3f,    0.0f },
    { 3,  1,  -2404.1f,   -56.6f,    -4.2f,    4.0f },
    { 3,  2,   1243.8f,   237.5f,     0.4f,   -0.3f },
    { 3,  3,    453.6f,  -549.5f,   -15.6f,   -4.1f },
    { 4,  0,    895.0f,     0.0f,    -1.6f,    0.0f },
    { 4,  1,    799.5f,   278.6f,    -2.4f,   -1.1f },
    { 4,  2,     55.7f,  -133.9f,    -6.0f,    4.1f },
    { 4,  3,   -281.1f,   212.0f,     5.6f,    1.6f },
    { 4,  4,     12.1f,  -375.6f,    -7.0f,   -4.4f },
    { 5,  0,   -233.2f,     0.0f,     0.6f,    0.0f },
    { 5,  1,    368.9f,    45.4f,     1.4f,   -0.5f },
    { 5,  2,    187.2f,   220.2f,     0.0f,    2.2f },
    { 5,  3,   -138.7f,  -122.9f,     0.6f,    0.4f },
    { 5,  4,   -142.0f,    43.0f,     2.2f,    1.7f },
    { 5,  5,     20.9f,   106.1f,     0.9f,    1.9f },
    { 6,  0,     64.4f,     0.0f,    -0.2f,    0.0f },
    { 6,  1,     63.8f,   -18.4f,    -0.4f,    0.3f },
    { 6,  2,     76.9f,    16.8f,     0.9f,   -1.6f },
    { 6,  3,   -115.7f,    48.8f,     1.2f,   -0.4f },
    { 6,  4,    -40.9f,   -59.8f,    -0.9f,    0.9f },
    { 6,  5,     14.9f,    10.9f,     0.3f,    0.7f },
    { 6,  6,    -60.7f,    72.7f,     0.9f,    0.9f },
    { 7,  0,     79.5f,     0.0f,    -0.0f,    0.0f },
    { 7,  1,    -77.0f,   -48.9f,    -0.1f,    0.6f },
    { 7,  2,     -8.8f,   -14.4f,    -0.1f,    0.5f },
    { 7,  3,     59.3f,    -1.0f,     0.5f,   -0.8f },
    { 7,  4,     15.8f,    23.4f,    -0.1f,    0.0f },
    { 7,  5,      2.5f,    -7.4f,    -0.8f,   -1.0f },
    { 7,  6,    -11.1f,   -25.1f,    -0.8f,    0.6f },
    { 7,  7,     14.2f,    -2.3f,     0.8f,   -0.2f },
    { 8,  0,     23.2f,     0.0f,    -0.1f,    0.0f },
    { 8,  1,     10.8f,     7.1f,     0.2f,   -0.2f },
    { 8,  2,    -17.5f,   -12.6f,     0.0f,    0.5f },
    { 8,  3,      2.0f,    11.4f,     0.5f,   -0.4f },
    { 8,  4,    -21.7f,    -9.7f,    -0.1f,    0.4f },
    { 8,  5,     16.9f,    12.7f,     0.3f,   -0.5f },
    { 8,  6,     15.0f,     0.7f,     0.2f,   -0.6f },
    { 8,  7,    -16.8f,    -5.2f,    -0.0f,    0.3f },
    { 8,  8,      0.9f,     3.9f,     0.2f,    0.2f },
    { 9,  0,      4.6f,     0.0f,    -0.0f,    0.0f },
    { 9,  1,      7.8f,   -24.8f,    -0.1f,   -0.3f },
    { 9,  2,      3.0f,    12.2f,     0.1f,    0.3f },
    { 9,  3,     -0.2f,     8.3f,     0.3f,   -0.3f },
    { 9,  4,     -2.5f,    -3.3f,    -0.3f,    0.3f },
    { 9,  5,    -13.1f,    -5.2f,     0.0f,    0.2f },
    { 9,  6,      2.4f,     7.2f,     0.3f,   -0.1f },
    { 9,  7,      8.6f,    -0.6f,    -0.1f,   -0.2f },
    { 9,  8,     -8.7f,     0.8f,     0.1f,    0.4f },
    { 9,  9,    -12.9f,    10.0f,    -0.1f,    0.1f },
    {10,  0,     -1.3f,     0.0f,     0.1f,    0.0f },
    {10,  1,     -6.4f,     3.3f,     0.0f,    0.0f },
    {10,  2,      0.2f,     0.0f,     0.1f,   -0.0f },
    {10,  3,      2.0f,     2.4f,     0.1f,   -0.2f },
    {10,  4,     -1.0f,     5.3f,    -0.0f,    0.1f },
    {10,  5,     -0.6f,    -9.1f,    -0.3f,   -0.1f },
    {10,  6,     -0.9f,     0.4f,     0.0f,    0.1f },
    {10,  7,      1.5f,    -4.2f,    -0.1f,    0.0f },
    {10,  8,      0.9f,    -3.8f,    -0.1f,   -0.1f },
    {10,  9,     -2.7f,     0.9f,    -0.0f,    0.2f },
    {10, 10,     -3.9f,    -9.1f,    -0.0f,   -0.0f },
    {11,  0,      2.9f,     0.0f,     0.0f,    0.0f },
    {11,  1,     -1.5f,     0.0f,    -0.0f,   -0.0f },
    {11,  2,     -2.5f,     2.9f,     0.0f,    0.1f },
    {11,  3,      2.4f,    -0.6f,     0.0f,   -0.0f },
    {11,  4,     -0.6f,     0.2f,     0.0f,    0.1f },
    {11,  5,     -0.1f,     0.5f,    -0.1f,   -0.0f },
    {11,  6,     -0.6f,    -0.3f,     0.0f,   -0.0f },
    {11,  7,     -0.1f,    -1.2f,    -0.0f,    0.1f },
    {11,  8,      1.1f,    -1.7f,    -0.1f,   -0.0f },
    {11,  9,     -1.0f,    -2.9f,    -0.1f,    0.0f },
    {11, 10,     -0.2f,    -1.8f,    -0.1f,    0.0f },
    {11, 11,      2.6f,    -2.3f,    -0.1f,    0.0f },
    {12,  0,     -2.0f,     0.0f,     0.0f,    0.0f },
    {12,  1,     -0.2f,    -1.3f,     0.0f,   -0.0f },
    {12,  2,      0.3f,     0.7f,    -0.0f,    0.0f },
    {12,  3,      1.2f,     1.0f,    -0.0f,   -0.1f },
    {12,  4,     -1.3f,    -1.4f,    -0.0f,    0.1f },
    {12,  5,      0.6f,    -0.0f,    -0.0f,   -0.0f },
    {12,  6,      0.6f,     0.6f,     0.1f,   -0.0f },
    {12,  7,      0.5f,    -0.1f,    -0.0f,   -0.0f },
    {12,  8,     -0.1f,     0.8f,     0.0f,    0.0f },
    {12,  9,     -0.4f,     0.1f,     0.0f,   -0.0f },
    {12, 10,     -0.2f,    -1.0f,    -0.1f,   -0.0f },
    {12, 11,     -1.3f,     0.1f,    -0.0f,    0.0f },
    {12, 12,     -0.7f,     0.2f,    -0.1f,   -0.1f },
};

#define WMM_NRECORDS  (int32_t)(sizeof(WMM2025_DATA) / sizeof(WMM2025_DATA[0]))

/* =========================================================================
 * Internal helper prototypes
 * All helpers that need the Julian Date now accept double JD to avoid
 * re-introducing truncation after the split at the public interface.
 * ========================================================================= */
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

/**
 * \brief Compute the Earth's magnetic field vector at a satellite position in ECI (precision-improved).
 *
 * Precision-improved variant of wmm_eci_embedded().  The Julian Date is
 * accepted as a split (integer day + fractional day) so it is never truncated
 * to float precision before arithmetic.  GMST and jd2year computations are
 * performed in double internally; the remainder of the pipeline uses float.
 *
 * \param r_ECI   Satellite position in ECI frame (metres), length-3 array.
 * \param JD_int  Integer part of the Julian Date (e.g. 2461116 for JD 2461116.223974).
 * \param JD_frac Fractional part of the Julian Date in [0, 1) (e.g. 0.223974f).
 * \param B_ECI   Output magnetic field vector in ECI frame (Tesla), length-3 array.
 */
void wmm_eci_embedded_v2(const float r_ECI[3],
                        int32_t JD_int, float JD_frac,
                        float B_ECI[3])
{
    /* Reconstruct full JD in double — no cancellation, no truncation */
    double JD = (double)JD_int + (double)JD_frac;

    /* --- Build time-extrapolated Gauss coefficients (stack allocated) --- */
    float g[WMM_DIM][WMM_DIM];
    float h[WMM_DIM][WMM_DIM];
    load_wmm2025(JD, g, h);

    /* --- ECI -> ECEF via GMST rotation (Rz) --- */
    float theta = gmst_from_jd(JD);   /* Fix 1: computed in double internally */
    float cos_t = cosf(theta);
    float sin_t = sinf(theta);

    float r_ecef[3];
    r_ecef[0] =  cos_t * r_ECI[0] + sin_t * r_ECI[1];
    r_ecef[1] = -sin_t * r_ECI[0] + cos_t * r_ECI[1];
    r_ecef[2] =  r_ECI[2];

    /* --- ECEF -> geodetic (WGS-84) --- */
    float lat_rad, lon_rad, alt_m;
    ecef_to_geodetic(r_ecef, &lat_rad, &lon_rad, &alt_m);

    /* --- Synthesize field in NED frame (nT) --- */
    float B_N, B_E, B_D;
    synthesize_mag_field(lat_rad, lon_rad, alt_m, g, h, &B_N, &B_E, &B_D);

    /* --- NED -> ECEF --- */
    float sin_lat = sinf(lat_rad);
    float cos_lat = cosf(lat_rad);
    float sin_lon = sinf(lon_rad);
    float cos_lon = cosf(lon_rad);

    float B_ecef[3];
    B_ecef[0] = (-sin_lat * cos_lon) * B_N + (-sin_lon) * B_E + (-cos_lat * cos_lon) * B_D;
    B_ecef[1] = (-sin_lat * sin_lon) * B_N + ( cos_lon) * B_E + (-cos_lat * sin_lon) * B_D;
    B_ecef[2] = ( cos_lat          ) * B_N +           0 * B_E + (-sin_lat           ) * B_D;

    /* --- ECEF -> ECI  (Rz^T = Rz(-theta)) --- */
    float B_eci_nT[3];
    B_eci_nT[0] =  cos_t * B_ecef[0] - sin_t * B_ecef[1];
    B_eci_nT[1] =  sin_t * B_ecef[0] + cos_t * B_ecef[1];
    B_eci_nT[2] =  B_ecef[2];

    /* --- nT -> T --- */
    B_ECI[0] = B_eci_nT[0] * 1.0e-9f;
    B_ECI[1] = B_eci_nT[1] * 1.0e-9f;
    B_ECI[2] = B_eci_nT[2] * 1.0e-9f;
}

/**
 * \brief Load WMM2025 Gauss coefficients linearly extrapolated to a given epoch.
 *
 * Reads the embedded WMM2025_DATA table and applies secular-variation rates
 * to produce time-adjusted g[n][m] and h[n][m] arrays.  JD is accepted as
 * double to avoid precision loss from the split-JD interface.
 *
 * \param JD  Julian Date of the target epoch (double precision).
 * \param g   Output cosine Gauss coefficients g[n][m] in nT (WMM_DIM x WMM_DIM).
 * \param h   Output sine Gauss coefficients h[n][m] in nT (WMM_DIM x WMM_DIM).
 */
void load_wmm2025(double JD,
                        float g[WMM_DIM][WMM_DIM],
                        float h[WMM_DIM][WMM_DIM])
{
    float decimalYear = jd2year(JD);   /* Fix 2: double arithmetic inside */
    float dt          = decimalYear - WMM_EPOCH;

    memset(g, 0, sizeof(float) * WMM_DIM * WMM_DIM);
    memset(h, 0, sizeof(float) * WMM_DIM * WMM_DIM);

    for (int32_t i = 0; i < WMM_NRECORDS; ++i) {
        int32_t n = WMM2025_DATA[i].n;
        int32_t m = WMM2025_DATA[i].m;

        g[n][m] = WMM2025_DATA[i].g + dt * WMM2025_DATA[i].gdot;
        h[n][m] = WMM2025_DATA[i].h + dt * WMM2025_DATA[i].hdot;
    }
}

/**
 * \brief Convert a Julian Date to a decimal year (double precision internally).
 *
 * Performs the subtraction from J2000.0 in double to avoid the 3-4 digit
 * precision loss that occurs when a large JD is held in a float before
 * the cancellation step.
 *
 * \param JD  Julian Date (double precision).
 * \return    Corresponding decimal year as float (e.g. 2025.5f).
 */
float jd2year(double JD)
{
    return (float)(2000.0 + (JD - 2451545.0) / 365.25);
}

/**
 * \brief Compute Greenwich Mean Sidereal Time (GMST) from a Julian Date (double precision internally).
 *
 * The IAU polynomial produces an intermediate value of ~3.5e6 seconds.
 * In float this value has a precision of ~0.25, causing up to ~0.25° error
 * after fmod(gmst_sec, 360) — directly rotating the ECI<->ECEF frame.
 * Computing in double reduces this error to below 1e-9 degrees.
 *
 * \param JD  Julian Date (double precision).
 * \return    GMST angle in radians as float, in [0, 2π).
 */
float gmst_from_jd(double JD)
{
    double T        = (JD - 2451545.0) / 36525.0;
    double gmst_sec = 67310.54841
                    + (876600.0 * 3600.0 + 8640184.812866) * T
                    + 0.093104 * T * T
                    - 6.2e-6   * T * T * T;
    double gmst_deg = fmod(gmst_sec / 240.0, 360.0);
    if (gmst_deg < 0.0) gmst_deg += 360.0;
    return (float)(gmst_deg * M_PI / 180.0);
}

/**
 * \brief Convert ECEF Cartesian coordinates to WGS-84 geodetic coordinates.
 *
 * Uses the Bowring iterative method (up to 10 iterations).  Handles the
 * polar singularity (p ≈ 0) as a special case.  Float precision is
 * sufficient for this stage of the pipeline.
 *
 * \param r_ecef   Input position in ECEF frame (metres), length-3 array [x, y, z].
 * \param lat_rad  Output geodetic latitude (radians), in [-π/2, π/2].
 * \param lon_rad  Output geodetic longitude (radians), in (-π, π].
 * \param alt_m    Output altitude above the WGS-84 ellipsoid (metres).
 */
void ecef_to_geodetic(const float r_ecef[3],
                          float *lat_rad, float *lon_rad, float *alt_m)
{
    const float a   = WGS84_A;
    const float b   = WGS84_B;
    const float e2  = WGS84_E2;
    const float ep2 = WGS84_EP2;

    float x = r_ecef[0];
    float y = r_ecef[1];
    float z = r_ecef[2];

    *lon_rad = atan2f(y, x);

    float p = sqrtf(x * x + y * y);

    /* Pole singularity */
    if (p < 1.0e-12f) {
        *lon_rad = 0.0f;
        if (z >= 0.0f) {
            *lat_rad = (float)M_PI / 2.0f;
            *alt_m   = z - b;
        } else {
            *lat_rad = (float)(-M_PI / 2.0f);
            *alt_m   = -z - b;
        }
        return;
    }

    /* Bowring initial estimate */
    float theta = atan2f(z * a, p * b);

    float lat = 0.0f, alt = 0.0f;
    for (int32_t iter = 0; iter < 10; ++iter) {
        float sin_t = sinf(theta);
        float cos_t = cosf(theta);

        float num = z + ep2 * b * sin_t * sin_t * sin_t;
        float den = p - e2  * a * cos_t * cos_t * cos_t;
        lat = atan2f(num, den);

        float sin_lat = sinf(lat);
        float cos_lat = cosf(lat);
        float N = a / sqrtf(1.0f - e2 * sin_lat * sin_lat);

        alt = p * cos_lat + z * sin_lat - a * sqrtf(1.0f - e2 * sin_lat * sin_lat);

        float theta_new = atan2f(z * (1.0f - e2 * N / (N + alt)), p);
        if (fabsf(theta_new - theta) < 1.0e-10f) {
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
 * Float precision is sufficient for this stage of the pipeline.
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
void synthesize_mag_field(float lat_rad, float lon_rad, float alt_m,
                               const float g[WMM_DIM][WMM_DIM],
                               const float h[WMM_DIM][WMM_DIM],
                               float *B_N, float *B_E, float *B_D)
{
    /* ---- Geodetic -> geocentric ---- */
    const float Re  = WGS84_A;
    const float f   = WGS84_F;
    float sin_lat   = sinf(lat_rad);
    float cos_lat   = cosf(lat_rad);
    float rc        = Re / sqrtf(1.0f - (2.0f * f - f * f) * sin_lat * sin_lat);

    float p_gc = (rc + alt_m) * cos_lat;
    float z_gc = (rc * (1.0f - f) * (1.0f - f) + alt_m) * sin_lat;
    float r    = sqrtf(p_gc * p_gc + z_gc * z_gc);

    float phi_prime = asinf(z_gc / r);
    float psi       = lat_rad - phi_prime;

    /* ---- Spherical harmonic synthesis ---- */
    float theta     = (float)(M_PI / 2.0) - phi_prime;
    float cos_theta = cosf(theta);
    float sin_theta = sinf(theta);

    float P [WMM_DIM][WMM_DIM];
    float dP[WMM_DIM][WMM_DIM];
    memset(P,  0, sizeof(P));
    memset(dP, 0, sizeof(dP));

    P [0][0] = 1.0f;
    P [1][0] = cos_theta;
    P [1][1] = sin_theta;
    dP[0][0] = 0.0f;
    dP[1][0] = -sin_theta;
    dP[1][1] =  cos_theta;

    for (int32_t n = 2; n <= WMM_NMAX; ++n) {
        for (int32_t m = 0; m <= n; ++m) {
            if (m == n) {
                float factor = sqrtf(1.0f - 1.0f / (2.0f * (float)n));
                P [n][m] = sin_theta * P [n-1][m-1] * factor;
                dP[n][m] = (sin_theta * dP[n-1][m-1] + cos_theta * P[n-1][m-1]) * factor;
            } else {
                float n2    = (float)(n * n);
                float m2    = (float)(m * m);
                float nm1_2 = (float)((n-1)*(n-1));
                float denom = sqrtf(n2 - m2);
                float K     = (2.0f * (float)n - 1.0f) / denom;
                float M_val = sqrtf(nm1_2 - m2)        / denom;

                P [n][m] = K * cos_theta * P [n-1][m] - M_val * P [n-2][m];
                dP[n][m] = K * (cos_theta * dP[n-1][m] - sin_theta * P[n-1][m])
                         - M_val * dP[n-2][m];
            }
        }
    }

    /* ---- Sum up the series ---- */
    float Br   = 0.0f;
    float Bt   = 0.0f;
    float Bphi = 0.0f;

    const float a_ref = GEO_REF_RADIUS;

    for (int32_t n = 1; n <= WMM_NMAX; ++n) {
        float ar = a_ref / r;
        float ratio = ar;
        for (int32_t k = 0; k < n + 1; ++k) ratio *= ar;

        for (int32_t m = 0; m <= n; ++m) {
            float gm = g[n][m];
            float hm = h[n][m];

            float sin_mlon = sinf((float)m * lon_rad);
            float cos_mlon = cosf((float)m * lon_rad);

            float xy_comp = gm * cos_mlon + hm * sin_mlon;
            float z_comp  = gm * sin_mlon - hm * cos_mlon;

            Br += (float)(n + 1) * ratio * xy_comp * P[n][m];
            Bt += ratio * xy_comp * dP[n][m];

            if (m > 0 && fabsf(sin_theta) > 1.0e-10f) {
                Bphi += ratio * (float)m * z_comp * P[n][m] / sin_theta;
            }
        }
    }

    /* ---- Geocentric (r, theta, phi) -> geodetic NED ---- */
    float B_X_geo = -Bt;
    float B_Z_geo = -Br;

    float cos_psi = cosf(psi);
    float sin_psi = sinf(psi);

    *B_N = -(B_X_geo * cos_psi + B_Z_geo * sin_psi);
    *B_E =  Bphi;
    *B_D =  (B_X_geo * sin_psi + B_Z_geo * cos_psi);
}

#endif /* MAG_SPHERE_UPDATE */