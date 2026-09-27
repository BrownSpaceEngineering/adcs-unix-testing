#include "include/sgp4_propagate.h"
#include "math.h"
#include <string.h>

// MATLAB mod(x, y) for y > 0: result in [0, y)
static double mod_pos(double x, double y) {
    return x - floor(x / y) * y;
}

static double dot3d(const double *a, const double *b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void cross3d(const double *a, const double *b, double *c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

static double norm3d(const double *a) {
    return sqrt(dot3d(a, a));
}

// C = A * B, all 3x3 row-major
static void matmul3(const double *A, const double *B, double *C) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            C[3 * i + j] = A[3 * i] * B[j] + A[3 * i + 1] * B[3 + j] + A[3 * i + 2] * B[6 + j];
        }
    }
}

// y = A * x, A 3x3 row-major
static void matvec3(const double *A, const double *x, double *y) {
    for (int i = 0; i < 3; i++) {
        y[i] = A[3 * i] * x[0] + A[3 * i + 1] * x[1] + A[3 * i + 2] * x[2];
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
 * \fn rv2coe
 *
 * \brief Osculating classical elements from r (km), v (km/s).
 *
 * For near-circular or near-equatorial orbits the undefined angles are set to 0
 * and absorbed into nu (so orbital_to_eci(a,e,i,RAAN,argp,nu) still reproduces r, v).
 *
 * \param[in]  r   Position (km), 3 elements
 * \param[in]  v   Velocity (km/s), 3 elements
 * \param[in]  mu  Gravitational parameter (km^3/s^2)
 * \param[out] oe  [a; e; i; RAAN; argp; nu; M] (km, rad), 7 elements
 */
static void rv2coe(const double *r, const double *v, double mu, double *oe) {
    const double small = 1.0e-10;
    double rmag = norm3d(r);
    double vmag = norm3d(v);
    double h[3];
    cross3d(r, v, h);
    double hmag = norm3d(h);
    double nvec[3] = {-h[1], h[0], 0.0};
    double nmag = norm3d(nvec);
    double rdotv = dot3d(r, v);
    double evec[3];
    for (int k = 0; k < 3; k++) {
        evec[k] = ((vmag * vmag - mu / rmag) * r[k] - rdotv * v[k]) / mu;
    }
    double e  = norm3d(evec);
    double xi = 0.5 * vmag * vmag - mu / rmag;
    double a  = -mu / (2.0 * xi);
    double incl = acos(fmax(-1.0, fmin(1.0, h[2] / hmag)));

    int equatorial = nmag < small * hmag;
    int circular   = e < small;

    double RAAN, argp, nu;
    double c[3];

    if (equatorial) {
        RAAN = 0.0;
    } else {
        RAAN = atan2(nvec[1], nvec[0]);
    }

    if (circular) {
        argp = 0.0;
        if (equatorial) {
            nu = atan2(r[1], r[0]); // true longitude
            if (h[2] < 0) {
                nu = -nu;
            }
        } else {
            cross3d(nvec, r, c);
            nu = atan2(dot3d(c, h) / hmag, dot3d(nvec, r)); // arg of latitude
        }
    } else {
        if (equatorial) {
            argp = atan2(evec[1], evec[0]); // longitude of periapsis
            if (h[2] < 0) {
                argp = -argp;
            }
        } else {
            cross3d(nvec, evec, c);
            argp = atan2(dot3d(c, h) / hmag, dot3d(nvec, evec));
        }
        cross3d(evec, r, c);
        nu = atan2(dot3d(c, h) / hmag, dot3d(evec, r));
    }

    double E = 2.0 * atan2(sqrt(1.0 - e) * sin(nu / 2.0), sqrt(1.0 + e) * cos(nu / 2.0));
    double M = E - e * sin(E);

    const double twopi = 2.0 * SGP4_PI;
    oe[0] = a;
    oe[1] = e;
    oe[2] = incl;
    oe[3] = mod_pos(RAAN, twopi);
    oe[4] = mod_pos(argp, twopi);
    oe[5] = mod_pos(nu, twopi);
    oe[6] = mod_pos(M, twopi);
}

/**
 * \fn sgp4_propagate
 *
 * \brief SGP4 propagation step (near-Earth / LEO only).
 *
 * Near-Earth branch of Vallado, Crawford, Hujsak, Kelso, "Revisiting
 * Spacetrack Report #3", AIAA 2006-6753 (sgp4 routine).
 *
 * ALWAYS propagate from the uplinked epoch elements: dt is the total time
 * since the element epoch, NOT the time since the last call. Do not feed
 * the outputs back in as new elements.
 *
 * \param[in]  sat      Constants from sgp4_init (recompute only when new elements arrive)
 * \param[in]  dt       Time since element epoch (s), UTC-consistent clock
 * \param[out] r_gcrf   Position (km) in GCRF (~ICRF), 3 elements
 * \param[out] v_gcrf   Velocity (km/s) in GCRF, 3 elements
 * \param[out] oe_osc   Osculating Keplerian elements in GCRF, 7 elements
 *                      [a (km); e; i; RAAN; argp; nu (true anom); M] (rad)
 *                      (same convention and mu as orbital_to_eci)
 * \param[out] oe_mean  SGP4 mean elements at dt (TEME), 6 elements
 *                      [a (km); e; i; RAAN; argp; M] (rad) -- secular+drag only
 * \param[out] r_teme   Native SGP4 position output in TEME (km), 3 elements
 * \param[out] v_teme   Native SGP4 velocity output in TEME (km/s), 3 elements
 *
 * \return errCode: 0 ok
 *                  1 mean e out of range   2 mean motion < 0
 *                  4 semi-latus rectum < 0 6 satellite has decayed
 *                  7 deep-space orbit (from sgp4_init), not supported
 */
int sgp4_propagate(const sgp4_sat_t *sat, double dt,
                   double *r_gcrf, double *v_gcrf, double *oe_osc, double *oe_mean,
                   double *r_teme, double *v_teme)
{
    memset(r_gcrf, 0, sizeof(double) * 3);
    memset(v_gcrf, 0, sizeof(double) * 3);
    memset(r_teme, 0, sizeof(double) * 3);
    memset(v_teme, 0, sizeof(double) * 3);
    memset(oe_osc, 0, sizeof(double) * 7);
    memset(oe_mean, 0, sizeof(double) * 6);

    if (sat->initErr != 0) {
        return sat->initErr;
    }

    const double twopi = 2.0 * SGP4_PI;
    const double x2o3  = 2.0 / 3.0;
    const double xke   = sat->xke;
    double tsince = dt / 60.0; // minutes since epoch

    // secular gravity and atmospheric drag
    double xmdf   = sat->mo + sat->mdot * tsince;
    double argpdf = sat->argpo + sat->argpdot * tsince;
    double nodedf = sat->nodeo + sat->nodedot * tsince;
    double argpm  = argpdf;
    double mm     = xmdf;
    double t2     = tsince * tsince;
    double nodem  = nodedf + sat->nodecf * t2;
    double tempa  = 1.0 - sat->cc1 * tsince;
    double tempe  = sat->bstar * sat->cc4 * tsince;
    double templ  = sat->t2cof * t2;

    if (sat->isimp == 0) {
        double delomg   = sat->omgcof * tsince;
        double delmtemp = 1.0 + sat->eta * cos(xmdf);
        double delm     = sat->xmcof * (delmtemp * delmtemp * delmtemp - sat->delmo);
        double temp     = delomg + delm;
        mm    = xmdf + temp;
        argpm = argpdf - temp;
        double t3 = t2 * tsince;
        double t4 = t3 * tsince;
        tempa = tempa - sat->d2 * t2 - sat->d3 * t3 - sat->d4 * t4;
        tempe = tempe + sat->bstar * sat->cc5 * (sin(mm) - sat->sinmao);
        templ = templ + sat->t3cof * t3 + t4 * (sat->t4cof + tsince * sat->t5cof);
    }

    double nm    = sat->no_unkozai;
    double em    = sat->ecco;
    double inclm = sat->inclo;

    if (nm <= 0.0) {
        return 2;
    }
    double am = pow(xke / nm, x2o3) * tempa * tempa;
    nm = xke / pow(am, 1.5);
    em = em - tempe;

    if (em >= 1.0 || em < -0.001) {
        return 1;
    }
    if (em < 1.0e-6) {
        em = 1.0e-6;
    }
    mm = mm + sat->no_unkozai * templ;
    double xlm = mm + argpm + nodem;

    nodem = fmod(nodem, twopi);
    argpm = fmod(argpm, twopi);
    xlm   = fmod(xlm, twopi);
    mm    = fmod(xlm - argpm - nodem, twopi);

    oe_mean[0] = am * sat->radiusearthkm;
    oe_mean[1] = em;
    oe_mean[2] = inclm;
    oe_mean[3] = mod_pos(nodem, twopi);
    oe_mean[4] = mod_pos(argpm, twopi);
    oe_mean[5] = mod_pos(mm, twopi);

    double sinip = sin(inclm);
    double cosip = cos(inclm);

    // long period periodics (near-Earth: no lunar-solar terms)
    double axnl = em * cos(argpm);
    double temp = 1.0 / (am * (1.0 - em * em));
    double aynl = em * sin(argpm) + temp * sat->aycof;
    double xl   = mm + argpm + nodem + temp * sat->xlcof * axnl;

    // solve Kepler's equation
    double u    = fmod(xl - nodem, twopi);
    double eo1  = u;
    double tem5 = 9999.9;
    int ktr = 1;
    double sineo1 = 0.0, coseo1 = 1.0;
    while (fabs(tem5) >= 1.0e-12 && ktr <= 10) {
        sineo1 = sin(eo1);
        coseo1 = cos(eo1);
        tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
        tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
        if (fabs(tem5) >= 0.95) {
            if (tem5 > 0.0) {
                tem5 = 0.95;
            } else {
                tem5 = -0.95;
            }
        }
        eo1 = eo1 + tem5;
        ktr = ktr + 1;
    }

    // short period preliminary quantities
    double ecose = axnl * coseo1 + aynl * sineo1;
    double esine = axnl * sineo1 - aynl * coseo1;
    double el2   = axnl * axnl + aynl * aynl;
    double pl    = am * (1.0 - el2);
    if (pl < 0.0) {
        return 4;
    }
    double rl     = am * (1.0 - ecose);
    double rdotl  = sqrt(am) * esine / rl;
    double rvdotl = sqrt(pl) / rl;
    double betal  = sqrt(1.0 - el2);
    temp = esine / (1.0 + betal);
    double sinu  = am / rl * (sineo1 - aynl - axnl * temp);
    double cosu  = am / rl * (coseo1 - axnl + aynl * temp);
    double su    = atan2(sinu, cosu);
    double sin2u = (cosu + cosu) * sinu;
    double cos2u = 1.0 - 2.0 * sinu * sinu;
    temp = 1.0 / pl;
    double temp1 = 0.5 * sat->j2 * temp;
    double temp2 = temp1 * temp;

    // update for short period periodics
    double mrt   = rl * (1.0 - 1.5 * temp2 * betal * sat->con41) + 0.5 * temp1 * sat->x1mth2 * cos2u;
    su = su - 0.25 * temp2 * sat->x7thm1 * sin2u;
    double xnode = nodem + 1.5 * temp2 * cosip * sin2u;
    double xinc  = inclm + 1.5 * temp2 * cosip * sinip * cos2u;
    double mvt   = rdotl - nm * temp1 * sat->x1mth2 * sin2u / xke;
    double rvdot = rvdotl + nm * temp1 * (sat->x1mth2 * cos2u + 1.5 * sat->con41) / xke;

    // orientation vectors
    double sinsu = sin(su),    cossu = cos(su);
    double snod  = sin(xnode), cnod  = cos(xnode);
    double sini  = sin(xinc),  cosi  = cos(xinc);
    double xmx   = -snod * cosi;
    double xmy   =  cnod * cosi;
    double uvec[3] = {xmx * sinsu + cnod * cossu, xmy * sinsu + snod * cossu, sini * sinsu};
    double vvec[3] = {xmx * cossu - cnod * sinsu, xmy * cossu - snod * sinsu, sini * cossu};

    for (int k = 0; k < 3; k++) {
        r_teme[k] = mrt * uvec[k] * sat->radiusearthkm;
        v_teme[k] = (mvt * uvec[k] + rvdot * vvec[k]) * sat->vkmpersec;
    }

    if (mrt < 1.0) {
        return 6;
    }

    // TEME -> GCRF, evaluated at the current time (precession ~0.14"/day)
    double R[9];
    teme2gcrf_matrix(sat->T0 + dt / (86400.0 * 36525.0), R);
    matvec3(R, r_teme, r_gcrf);
    matvec3(R, v_teme, v_gcrf);

    rv2coe(r_gcrf, v_gcrf, 398600.4418, oe_osc);
    return 0;
}
