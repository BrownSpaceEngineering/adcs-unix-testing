#include "include/sgp4_init.h"
#include "math.h"

// ================= TUNABLE CONSTANTS =================================
// Gravity model: 72 = WGS-72 (use this for NORAD/Space-Track TLEs --
// they are generated with WGS-72), 84 = WGS-84.
#define GRAV_MODEL 72
// TAI-UTC leap seconds (37 s since 2017-01-01; update if IERS adds one)
#define DELTA_AT 37.0
// =====================================================================

// C = A * B, all 3x3 row-major
static void matmul3(const double *A, const double *B, double *C) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            C[3 * i + j] = A[3 * i] * B[j] + A[3 * i + 1] * B[3 + j] + A[3 * i + 2] * B[6 + j];
        }
    }
}

// Rk(a) rotates the coordinate frame by +a about axis k (row-major)
static void rot1(double a, double *M) {
    double c = cos(a), s = sin(a);
    M[0] = 1.0; M[1] = 0.0; M[2] = 0.0;
    M[3] = 0.0; M[4] = c;   M[5] = s;
    M[6] = 0.0; M[7] = -s;  M[8] = c;
}

static void rot2(double a, double *M) {
    double c = cos(a), s = sin(a);
    M[0] = c;   M[1] = 0.0; M[2] = -s;
    M[3] = 0.0; M[4] = 1.0; M[5] = 0.0;
    M[6] = s;   M[7] = 0.0; M[8] = c;
}

static void rot3(double a, double *M) {
    double c = cos(a), s = sin(a);
    M[0] = c;   M[1] = s;   M[2] = 0.0;
    M[3] = -s;  M[4] = c;   M[5] = 0.0;
    M[6] = 0.0; M[7] = 0.0; M[8] = 1.0;
}

/**
 * \fn teme2gcrf_matrix
 *
 * \brief Rotation TEME -> GCRF using IAU-76 precession and IAU-80 nutation
 *        (Vallado teme2eci: r_gcrf = P * N * R3(-eqe) * r_teme).
 *
 * Nutation truncated to the 10 largest IAU-80 terms (omitted terms each
 * < 0.02 arcsec, < ~1 m at LEO); the ~23 mas FK5->GCRF frame bias is ignored.
 *
 * \param[in]  T  TT Julian centuries since J2000
 * \param[out] R  3x3 rotation matrix, row-major
 */
static void teme2gcrf_matrix(double T, double *R) {
    const double as2r = SGP4_PI / (180.0 * 3600.0);
    double T2 = T * T;
    double T3 = T2 * T;

    // IAU-76 precession angles
    double zeta  = (2306.2181 * T + 0.30188 * T2 + 0.017998 * T3) * as2r;
    double theta = (2004.3109 * T - 0.42665 * T2 - 0.041833 * T3) * as2r;
    double z     = (2306.2181 * T + 1.09468 * T2 + 0.018203 * T3) * as2r;

    // mean obliquity of the ecliptic
    double epsb = (84381.448 - 46.8150 * T - 0.00059 * T2 + 0.001813 * T3) * as2r;

    // IAU-80 fundamental arguments
    const double r = 1296000.0;
    double l  = fmod(485866.733  + (1325.0 * r +  715922.633) * T + 31.310 * T2 + 0.064 * T3, r) * as2r;
    double lp = fmod(1287099.804 + (  99.0 * r + 1292581.224) * T -  0.577 * T2 - 0.012 * T3, r) * as2r;
    double F  = fmod(335778.877  + (1342.0 * r +  295263.137) * T - 13.257 * T2 + 0.011 * T3, r) * as2r;
    double D  = fmod(1072261.307 + (1236.0 * r + 1105601.328) * T -  6.891 * T2 + 0.019 * T3, r) * as2r;
    double Om = fmod(450160.280  - (   5.0 * r +  482890.539) * T +  7.455 * T2 + 0.008 * T3, r) * as2r;

    // largest IAU-80 nutation terms, units 0.0001 arcsec
    //  [l lp F D Om], dpsi = (A + B T) sin(arg), deps = (C + D T) cos(arg)
    static const double nut[10][9] = {
        { 0,  0,  0,  0,  1,  -171996.0, -174.2,  92025.0,   8.9},
        { 0,  0,  2, -2,  2,   -13187.0,   -1.6,   5736.0,  -3.1},
        { 0,  0,  2,  0,  2,    -2274.0,   -0.2,    977.0,  -0.5},
        { 0,  0,  0,  0,  2,     2062.0,    0.2,   -895.0,   0.5},
        { 0,  1,  0,  0,  0,     1426.0,   -3.4,     54.0,  -0.1},
        { 1,  0,  0,  0,  0,      712.0,    0.1,     -7.0,   0.0},
        { 0,  1,  2, -2,  2,     -517.0,    1.2,    224.0,  -0.6},
        { 0,  0,  2,  0,  1,     -386.0,   -0.4,    200.0,   0.0},
        { 1,  0,  2,  0,  2,     -301.0,    0.0,    129.0,  -0.1},
        { 0, -1,  2, -2,  2,      217.0,   -0.5,    -95.0,   0.3}
    };
    double dpsi = 0.0, deps = 0.0;
    for (int k = 0; k < 10; k++) {
        double arg = nut[k][0] * l + nut[k][1] * lp + nut[k][2] * F + nut[k][3] * D + nut[k][4] * Om;
        dpsi += (nut[k][5] + nut[k][6] * T) * sin(arg);
        deps += (nut[k][7] + nut[k][8] * T) * cos(arg);
    }
    dpsi = dpsi * 1.0e-4 * as2r;
    deps = deps * 1.0e-4 * as2r;
    double epst = epsb + deps;

    // equation of the equinoxes (TEME -> TOD is a rotation about z)
    double eqe = dpsi * cos(epsb);

    //  TOD   = R3(-eqe) * TEME
    //  MOD   = N * TOD,  N = R1(-epsb) * R3(dpsi) * R1(epst)
    //  J2000 = P * MOD,  P = R3(zeta) * R2(-theta) * R3(z)
    double A[9], B[9], tmp[9], N[9], P[9], PN[9];

    rot1(-epsb, A);
    rot3(dpsi, B);
    matmul3(A, B, tmp);
    rot1(epst, A);
    matmul3(tmp, A, N);

    rot3(zeta, A);
    rot2(-theta, B);
    matmul3(A, B, tmp);
    rot3(z, A);
    matmul3(tmp, A, P);

    matmul3(P, N, PN);
    rot3(-eqe, A);
    matmul3(PN, A, R);
}

/**
 * \fn sgp4_init
 *
 * \brief One-time SGP4 initialization for a new set of uplinked elements.
 *
 * Call once each time new elements arrive (1/day), then call
 * sgp4_propagate(sat, dt, ...) at the control rate (10 Hz).
 * Near-Earth branch of Vallado, Crawford, Hujsak, Kelso, "Revisiting
 * Spacetrack Report #3", AIAA 2006-6753 (initl + sgp4init).
 *
 * \param[in]  oe_epoch  SGP4/TLE *mean* elements at epoch (TEME frame), 6 elements:
 *                       [n (rev/day) mean motion (Kozai, as in TLE line 2),
 *                        e (-) eccentricity,
 *                        i (rad) inclination,
 *                        RAAN (rad) right ascension of ascending node,
 *                        argp (rad) argument of perigee,
 *                        M (rad) mean anomaly]
 * \param[in]  bstar     B* drag term (1/earth radii), TLE line 1 cols 54-61
 * \param[in]  epoch_jd  Julian date (UTC) of the element epoch, double (a float JD only
 *                       resolves 0.25 day). Only used for the TEME->GCRF rotation, which is
 *                       evaluated once here at epoch; the propagation itself only sees dt.
 * \param[out] sat       Precomputed constants for sgp4_propagate.
 *                       sat->initErr = 0 ok, 7 = deep-space orbit (period >= 225 min),
 *                       not supported.
 */
void sgp4_init(const float *oe_epoch, float bstar_f, double epoch_jd, sgp4_sat_t *sat)
{
#if GRAV_MODEL == 84
    const double radiusearthkm = 6378.137;
    const double mu = 398600.5;
    const double j2 =  0.00108262998905;
    const double j3 = -0.00000253215306;
    const double j4 = -0.00000161098761;
#else // WGS-72
    const double radiusearthkm = 6378.135;
    const double mu = 398600.8;
    const double j2 =  0.001082616;
    const double j3 = -0.00000253881;
    const double j4 = -0.00000165597;
#endif
    const double bstar = (double)bstar_f;
    const double xke   = 60.0 / sqrt(radiusearthkm * radiusearthkm * radiusearthkm / mu); // sqrt(mu) in er^1.5/min
    const double j3oj2 = j3 / j2;
    const double twopi = 2.0 * SGP4_PI;
    const double x2o3  = 2.0 / 3.0;
    const double temp4 = 1.5e-12;

    double no_kozai = (double)oe_epoch[0] * twopi / 1440.0; // rev/day -> rad/min
    double ecco  = (double)oe_epoch[1];
    double inclo = (double)oe_epoch[2];
    double nodeo = (double)oe_epoch[3];
    double argpo = (double)oe_epoch[4];
    double mo    = (double)oe_epoch[5];

    // ------------------------- initl ------------------------------------
    double eccsq  = ecco * ecco;
    double omeosq = 1.0 - eccsq;
    double rteosq = sqrt(omeosq);
    double cosio  = cos(inclo);
    double cosio2 = cosio * cosio;

    // un-Kozai the mean motion (Brouwer mean motion)
    double ak   = pow(xke / no_kozai, x2o3);
    double d1   = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
    double del  = d1 / (ak * ak);
    double adel = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
    del = d1 / (adel * adel);
    double no_unkozai = no_kozai / (1.0 + del);

    double ao    = pow(xke / no_unkozai, x2o3);
    double sinio = sin(inclo);
    double po    = ao * omeosq;
    double con42 = 1.0 - 5.0 * cosio2;
    double con41 = -con42 - cosio2 - cosio2;
    double posq  = po * po;
    double rp    = ao * (1.0 - ecco);

    int initErr = 0;
    if (twopi / no_unkozai >= 225.0) {
        initErr = 7;
    }

    // ------------------------- sgp4init ---------------------------------
    double ss         = 78.0 / radiusearthkm + 1.0;
    double qzms2ttemp = (120.0 - 78.0) / radiusearthkm;
    double qzms2t     = pow(qzms2ttemp, 4);

    // isimp = 1 for perigee < 220 km: drop the higher-order drag terms
    int isimp = (rp < (220.0 / radiusearthkm + 1.0));

    double sfour  = ss;
    double qzms24 = qzms2t;
    double perige = (rp - 1.0) * radiusearthkm;
    if (perige < 156.0) {
        sfour = perige - 78.0;
        if (perige < 98.0) {
            sfour = 20.0;
        }
        qzms24 = pow((120.0 - sfour) / radiusearthkm, 4);
        sfour  = sfour / radiusearthkm + 1.0;
    }
    double pinvsq = 1.0 / posq;

    double tsi   = 1.0 / (ao - sfour);
    double eta   = ao * ecco * tsi;
    double etasq = eta * eta;
    double eeta  = ecco * eta;
    double psisq = fabs(1.0 - etasq);
    double coef  = qzms24 * pow(tsi, 4);
    double coef1 = coef / pow(psisq, 3.5);
    double cc2   = coef1 * no_unkozai * (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
                   0.375 * j2 * tsi / psisq * con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    double cc1   = bstar * cc2;
    double cc3   = 0.0;
    if (ecco > 1.0e-4) {
        cc3 = -2.0 * coef * tsi * j3oj2 * no_unkozai * sinio / ecco;
    }
    double x1mth2 = 1.0 - cosio2;
    double cc4    = 2.0 * no_unkozai * coef1 * ao * omeosq * (eta * (2.0 + 0.5 * etasq) +
                    ecco * (0.5 + 2.0 * etasq) - j2 * tsi / (ao * psisq) * (-3.0 * con41 *
                    (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) + 0.75 * x1mth2 *
                    (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * argpo)));
    double cc5    = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
    double cosio4 = cosio2 * cosio2;
    double temp1  = 1.5 * j2 * pinvsq * no_unkozai;
    double temp2  = 0.5 * temp1 * j2 * pinvsq;
    double temp3  = -0.46875 * j4 * pinvsq * pinvsq * no_unkozai;
    double mdot   = no_unkozai + 0.5 * temp1 * rteosq * con41 +
                    0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
    double argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                     temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
    double xhdot1  = -temp1 * cosio;
    double nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
    double omgcof  = bstar * cc3 * cos(argpo);
    double xmcof   = 0.0;
    if (ecco > 1.0e-4) {
        xmcof = -x2o3 * coef * bstar / eeta;
    }
    double nodecf = 3.5 * omeosq * xhdot1 * cc1;
    double t2cof  = 1.5 * cc1;
    double xlcof;
    if (fabs(cosio + 1.0) > 1.5e-12) {
        xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    } else {
        xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
    }
    double aycof  = -0.5 * j3oj2 * sinio;
    double delmo  = pow(1.0 + eta * cos(mo), 3);
    double sinmao = sin(mo);
    double x7thm1 = 7.0 * cosio2 - 1.0;

    double d2 = 0.0, d3 = 0.0, d4 = 0.0;
    double t3cof = 0.0, t4cof = 0.0, t5cof = 0.0;
    if (!isimp) {
        double cc1sq = cc1 * cc1;
        d2 = 4.0 * ao * tsi * cc1sq;
        double temp = d2 * tsi * cc1 / 3.0;
        d3 = (17.0 * ao + sfour) * temp;
        d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * cc1;
        t3cof = d2 + 2.0 * cc1sq;
        t4cof = 0.25 * (3.0 * d3 + cc1 * (12.0 * d2 + 10.0 * cc1sq));
        t5cof = 0.2 * (3.0 * d4 + 12.0 * cc1 * d3 + 6.0 * d2 * d2 + 15.0 * cc1sq * (2.0 * d2 + cc1sq));
    }

    // TEME -> GCRF at epoch (TT Julian centuries since J2000). Held fixed between uplinks:
    // it only drifts ~0.3 arcsec/day (< ~10 m at LEO over a day).
    double T0 = (epoch_jd + (DELTA_AT + 32.184) / 86400.0 - 2451545.0) / 36525.0;
    double R[9];
    teme2gcrf_matrix(T0, R);

    sat->initErr       = initErr;
    sat->isimp         = isimp;
    sat->radiusearthkm = (float)radiusearthkm;
    sat->xke           = (float)xke;
    sat->j2            = (float)j2;
    sat->vkmpersec     = (float)(radiusearthkm * xke / 60.0);
    sat->bstar         = (float)bstar;
    sat->ecco          = (float)ecco;
    sat->inclo         = (float)inclo;
    sat->nodeo         = (float)nodeo;
    sat->argpo         = (float)argpo;
    sat->mo            = (float)mo;
    sat->no_unkozai    = (float)no_unkozai;
    sat->mdot          = (float)mdot;
    sat->argpdot       = (float)argpdot;
    sat->nodedot       = (float)nodedot;
    sat->nodecf        = (float)nodecf;
    sat->cc1           = (float)cc1;
    sat->cc4           = (float)cc4;
    sat->cc5           = (float)cc5;
    sat->t2cof         = (float)t2cof;
    sat->omgcof        = (float)omgcof;
    sat->xmcof         = (float)xmcof;
    sat->eta           = (float)eta;
    sat->delmo         = (float)delmo;
    sat->sinmao        = (float)sinmao;
    sat->d2            = (float)d2;
    sat->d3            = (float)d3;
    sat->d4            = (float)d4;
    sat->t3cof         = (float)t3cof;
    sat->t4cof         = (float)t4cof;
    sat->t5cof         = (float)t5cof;
    sat->con41         = (float)con41;
    sat->x1mth2        = (float)x1mth2;
    sat->x7thm1        = (float)x7thm1;
    sat->xlcof         = (float)xlcof;
    sat->aycof         = (float)aycof;
    for (int k = 0; k < 9; k++) {
        sat->R_teme2gcrf[k] = (float)R[k];
    }
}
