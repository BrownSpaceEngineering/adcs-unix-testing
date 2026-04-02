/**
 * wmm_eci.c
 *
 * C translation of wmmECI_embedded.m (WMM2025, ECI output).
 * gcc -O2 -Wall -Wextra -o magnetosphere_test magnetosphere.c -lm
 * Sample Inputs are from the matlab function wrldmagm that we're trying to mimic
 */

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

void wmmECI_embedded(const float r_ECI[3], float JD, float B_ECI[3]);

/* =========================================================================
 * Internal constants
 * ========================================================================= */

#define WMM_NMAX        12          
#define WMM_DIM        (WMM_NMAX + 2)   
#define WMM_EPOCH       2025.0f     /* WMM2025 reference epoch              */

/* WGS-84 / geomagnetic reference constants */
#define WGS84_A         6378137.0f
#define WGS84_F         (1.0f / 298.257223563f)
#define WGS84_B         (WGS84_A * (1.0f - WGS84_F))
#define WGS84_E2        (WGS84_F * (2.0f - WGS84_F))
#define WGS84_EP2       (WGS84_E2 / (1.0f - WGS84_E2))
#define GEO_REF_RADIUS  6371200.0f

#ifndef M_PI
#define M_PI 3.14159265358979323846f
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
 * ========================================================================= */
static float jd2year(float JD);
static float gmst_from_jd(float JD);
static void  ecef2geodetic(const float r_ecef[3],
                           float *lat_rad, float *lon_rad, float *alt_m);
static void  synthesizeMagField(float lat_rad, float lon_rad, float alt_m,
                                const float g[WMM_DIM][WMM_DIM],
                                const float h[WMM_DIM][WMM_DIM],
                                float *B_N, float *B_E, float *B_D);
static void  loadWMM2025(float JD,
                         float g[WMM_DIM][WMM_DIM],
                         float h[WMM_DIM][WMM_DIM]);

/* =========================================================================
 * Public function
 * ========================================================================= */
void wmmECI_embedded(const float r_ECI[3], float JD, float B_ECI[3])
{
    /* --- Build time-extrapolated Gauss coefficients (stack allocated) --- */
    float g[WMM_DIM][WMM_DIM];
    float h[WMM_DIM][WMM_DIM];
    loadWMM2025(JD, g, h);

    /* --- ECI -> ECEF via GMST rotation (Rz) --- */
    float theta = gmst_from_jd(JD);
    float cos_t = cosf(theta);
    float sin_t = sinf(theta);

    float r_ecef[3];
    r_ecef[0] =  cos_t * r_ECI[0] + sin_t * r_ECI[1];
    r_ecef[1] = -sin_t * r_ECI[0] + cos_t * r_ECI[1];
    r_ecef[2] =  r_ECI[2];

    /* --- ECEF -> geodetic (WGS-84) --- */
    float lat_rad, lon_rad, alt_m;
    ecef2geodetic(r_ecef, &lat_rad, &lon_rad, &alt_m);

    /* --- Synthesize field in NED frame (nT) --- */
    float B_N, B_E, B_D;
    synthesizeMagField(lat_rad, lon_rad, alt_m, g, h, &B_N, &B_E, &B_D);

    /* --- NED -> ECEF --- */
    float sin_lat = sinf(lat_rad);
    float cos_lat = cosf(lat_rad);
    float sin_lon = sinf(lon_rad);
    float cos_lon = cosf(lon_rad);

    float B_ecef[3];
    B_ecef[0] = (-sin_lat * cos_lon) * B_N + (-sin_lon) * B_E + (-cos_lat * cos_lon) * B_D;
    B_ecef[1] = (-sin_lat * sin_lon) * B_N + ( cos_lon) * B_E + (-cos_lat * sin_lon) * B_D;
    B_ecef[2] = ( cos_lat          ) * B_N +          0 * B_E + (-sin_lat           ) * B_D;

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

/* =========================================================================
 * Load WMM2025 coefficients extrapolated to JD
 * ========================================================================= */
static void loadWMM2025(float JD,
                        float g[WMM_DIM][WMM_DIM],
                        float h[WMM_DIM][WMM_DIM])
{
    float decimalYear = jd2year(JD);
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

/* =========================================================================
 * Julian Date -> decimal year
 * ========================================================================= */
static float jd2year(float JD)
{
    return 2000.0f + (JD - 2451545.0f) / 365.25f;
}

/* =========================================================================
 * Greenwich Mean Sidereal Time (radians)
 * ========================================================================= */
static float gmst_from_jd(float JD)
{
    float T        = (JD - 2451545.0f) / 36525.0f;
    float gmst_sec = 67310.54841f
                   + (876600.0f * 3600.0f + 8640184.812866f) * T
                   + 0.093104f * T * T
                   - 6.2e-6f   * T * T * T;
    float gmst_deg = fmodf(gmst_sec / 240.0f, 360.0f);
    if (gmst_deg < 0.0f) gmst_deg += 360.0f;
    return gmst_deg * (float)M_PI / 180.0f;
}

/* =========================================================================
 * ECEF -> WGS-84 geodetic  (Bowring iterative)
 * ========================================================================= */
static void ecef2geodetic(const float r_ecef[3],
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

/* =========================================================================
 * Synthesize magnetic field in NED frame (nT)
 * using Schmidt semi-normalised associated Legendre polynomials
 * ========================================================================= */
static void synthesizeMagField(float lat_rad, float lon_rad, float alt_m,
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

    float phi_prime = asinf(z_gc / r);          /* geocentric latitude       */
    float psi       = lat_rad - phi_prime;      /* geodetic - geocentric    */

    /* ---- Spherical harmonic synthesis ---- */
    float theta     = (float)(M_PI / 2.0f) - phi_prime;  /* geocentric colatitude */
    float cos_theta = cosf(theta);
    float sin_theta = sinf(theta);

    /* Stack-allocated Legendre arrays — WMM_DIM x WMM_DIM */
    float P [WMM_DIM][WMM_DIM];
    float dP[WMM_DIM][WMM_DIM];
    memset(P,  0, sizeof(P));
    memset(dP, 0, sizeof(dP));

    /* Seed values (n=0,m=0) and (n=1) — stored 0-indexed in P/dP */
    P [0][0] = 1.0f;
    P [1][0] = cos_theta;           /* P(1,0) */
    P [1][1] = sin_theta;           /* P(1,1) */
    dP[0][0] = 0.0f;
    dP[1][0] = -sin_theta;          /* dP(1,0)/dtheta */
    dP[1][1] =  cos_theta;          /* dP(1,1)/dtheta */

    for (int32_t n = 2; n <= WMM_NMAX; ++n) {
        for (int32_t m = 0; m <= n; ++m) {

            if (m == n) {
                /* Diagonal recurrence */
                float factor = sqrtf(1.0f - 1.0f / (2.0f * (float)n));
                P [n][m] = sin_theta * P [n-1][m-1] * factor;
                dP[n][m] = (sin_theta * dP[n-1][m-1] + cos_theta * P[n-1][m-1]) * factor;
            } else {
                /* Off-diagonal recurrence */
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
    float Br   = 0.0f;   /* radial    (outward +)   */
    float Bt   = 0.0f;   /* theta     (southward +) */
    float Bphi = 0.0f;   /* phi       (eastward +)  */

    const float a_ref = GEO_REF_RADIUS;

    for (int32_t n = 1; n <= WMM_NMAX; ++n) {
        float ar = a_ref / r;
        float ratio = ar;
        for (int32_t k = 0; k < n + 1; ++k) ratio *= ar;   /* ar^(n+2) */

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
    float B_X_geo = -Bt;    /* North, geocentric */
    float B_Z_geo = -Br;    /* Down,  geocentric */

    float cos_psi = cosf(psi);
    float sin_psi = sinf(psi);

    *B_N = -(B_X_geo * cos_psi + B_Z_geo * sin_psi);
    *B_E =  Bphi;
    *B_D =  (B_X_geo * sin_psi + B_Z_geo * cos_psi);
}

/* =========================================================================
 * TEST CASES (Main Function)
 * ========================================================================= */
typedef struct {
    const char *name;
    float r[3];
    float expected_ned[3];
} TestCase;

int main(void)
{
    /* Passed back to standard 32-bit float */
    float jd = 2461116.223974f;
    
    /* Array of the 10 orbital test cases with expected Matlab NED vectors */
    TestCase tests[] = {
        {"SAA Center",  { 4.5e6f, -3.7e6f, -3.4e6f}, { 2.0598e+04f,  6.0145e+03f, -1.8521e+04f}},
        {"ISS Max Lat", { 4.2e6f,  0.0f,    5.3e6f}, { 1.2509e+04f, -2.9673e+03f,  4.3991e+04f}},
        {"Mag N Pole",  { 0.5e6f,  0.0f,    6.75e6f}, { 1.2450e+03f, -1.3110e+03f,  4.7655e+04f}},
        {"Mag S Pole",  {-2.0e6f,  2.0e6f, -6.1e6f}, { 6.2109e+03f, -1.3507e+04f, -3.8299e+04f}},
        {"Eq 90 deg E", { 0.0f,    6.77e6f, 0.0f},   { 2.3820e+04f, -2.0337e+02f, -1.1299e+04f}},
        {"Bermuda",     { 2.5e6f, -5.3e6f,  3.4e6f}, { 2.0728e+04f,  4.2250e+03f,  2.6296e+04f}},
        {"Pacific 180", {-6.77e6f, 0.0f,    0.0f},   { 3.3273e+04f, -9.4567e+01f, -9.5107e+03f}},
        {"45N 45E",     { 3.4e6f,  3.4e6f,  4.8e6f}, { 1.8653e+04f, -3.0303e+03f,  3.3657e+04f}},
        {"SSO Polar",   { 1.1e6f,  0.0f,   -6.6e6f}, { 1.3022e+04f,  8.2517e+03f, -3.9577e+04f}},
        {"Vandenberg",  {-2.8e6f, -4.8e6f,  3.8e6f}, { 2.3925e+04f,  2.3729e+02f,  2.5251e+04f}}
    };

    int num_tests = sizeof(tests) / sizeof(tests[0]);

    float g[WMM_DIM][WMM_DIM];
    float h[WMM_DIM][WMM_DIM];
    loadWMM2025(jd, g, h);

    printf("Test Name       | Input Mag    | Output B_ECI [Bx, By, Bz] (nT)\n");
    printf("----------------------------------------------------------------\n");

    for (int i = 0; i < num_tests; ++i) {
        float r_check[3];
        r_check[0] = tests[i].r[0];
        r_check[1] = tests[i].r[1];
        r_check[2] = tests[i].r[2];

        float B_vec[3];
        
        wmmECI_embedded(r_check, jd, B_vec);

        float theta = gmst_from_jd(jd);
        float cos_t = cosf(theta);
        float sin_t = sinf(theta);
        
        float r_ecef[3];
        r_ecef[0] =  cos_t * r_check[0] + sin_t * r_check[1];
        r_ecef[1] = -sin_t * r_check[0] + cos_t * r_check[1];
        r_ecef[2] =  r_check[2];
        
        float lat_rad, lon_rad, alt_m;
        ecef2geodetic(r_ecef, &lat_rad, &lon_rad, &alt_m);
        
        float B_N, B_E, B_D;
        synthesizeMagField(lat_rad, lon_rad, alt_m, g, h, &B_N, &B_E, &B_D);

        float mat_N = tests[i].expected_ned[0];
        float mat_E = tests[i].expected_ned[1];
        float mat_D = tests[i].expected_ned[2];

        float dot_prod = (mat_N * B_N) + (mat_E * B_E) + (mat_D * B_D);
        float norm_mat = sqrtf(mat_N*mat_N + mat_E*mat_E + mat_D*mat_D);
        float norm_our = sqrtf(B_N*B_N + B_E*B_E + B_D*B_D);
        
        float accuracy = dot_prod / (norm_mat * norm_our);
        if (accuracy > 1.0f) accuracy = 1.0f;
        if (accuracy < -1.0f) accuracy = -1.0f;
        float angle_error_deg = acosf(accuracy) * 180.0f / (float)M_PI;

        float r_mag = sqrtf(r_check[0]*r_check[0] + r_check[1]*r_check[1] + r_check[2]*r_check[2]) / 1000.0f;
        float bx_nt = B_vec[0] * 1.0e9f;
        float by_nt = B_vec[1] * 1.0e9f;
        float bz_nt = B_vec[2] * 1.0e9f;
        float b_mag = sqrtf(bx_nt*bx_nt + by_nt*by_nt + bz_nt*bz_nt);

        printf("Matlab NED Vector (Correct magnetic field): [%.4e, %.4e, %.4e]\n", mat_N, mat_E, mat_D);
        printf("Our NED Vector: [%.4e, %.4e, %.4e]\n", B_N, B_E, B_D);
        printf("Our accuracy: %f\n", accuracy);
        printf("Our angular error (degrees): %f\n", angle_error_deg);
        printf("Our Starting ECI position: [%.4e, %.4e, %.4e]\n", r_check[0], r_check[1], r_check[2]);
        printf("%-15s | %6.0f km   | [%8.1f, %8.1f, %8.1f] (Total: %.1f nT)\n",
               tests[i].name, r_mag, bx_nt, by_nt, bz_nt, b_mag);
        printf("----------------------------------------------------------------\n");
    }

    return 0;
}