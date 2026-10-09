#include "include/sgp4_propagate.h"
#include "math.h"
#include <string.h>

// MATLAB mod(x, y) for y > 0: result in [0, y)
static float mod_pos(float x, float y) {
    return x - floorf(x / y) * y;
}

static float dot3d(const float *a, const float *b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void cross3d(const float *a, const float *b, float *c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

static float norm3d(const float *a) {
    return sqrtf(dot3d(a, a));
}

// y = A * x, A 3x3 row-major
static void matvec3(const float *A, const float *x, float *y) {
    for (int i = 0; i < 3; i++) {
        y[i] = A[3 * i] * x[0] + A[3 * i + 1] * x[1] + A[3 * i + 2] * x[2];
    }
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
static void rv2coe(const float *r, const float *v, float mu, float *oe) {
    const float small = 1.0e-10f;
    float rmag = norm3d(r);
    float vmag = norm3d(v);
    float h[3];
    cross3d(r, v, h);
    float hmag = norm3d(h);
    float nvec[3] = {-h[1], h[0], 0.0f};
    float nmag = norm3d(nvec);
    float rdotv = dot3d(r, v);
    float evec[3];
    for (int k = 0; k < 3; k++) {
        evec[k] = ((vmag * vmag - mu / rmag) * r[k] - rdotv * v[k]) / mu;
    }
    float e  = norm3d(evec);
    float xi = 0.5f * vmag * vmag - mu / rmag;
    float a  = -mu / (2.0f * xi);
    float incl = acosf(fmaxf(-1.0f, fminf(1.0f, h[2] / hmag)));

    int equatorial = nmag < small * hmag;
    int circular   = e < small;

    float RAAN, argp, nu;
    float c[3];

    if (equatorial) {
        RAAN = 0.0f;
    } else {
        RAAN = atan2f(nvec[1], nvec[0]);
    }

    if (circular) {
        argp = 0.0f;
        if (equatorial) {
            nu = atan2f(r[1], r[0]); // true longitude
            if (h[2] < 0) {
                nu = -nu;
            }
        } else {
            cross3d(nvec, r, c);
            nu = atan2f(dot3d(c, h) / hmag, dot3d(nvec, r)); // arg of latitude
        }
    } else {
        if (equatorial) {
            argp = atan2f(evec[1], evec[0]); // longitude of periapsis
            if (h[2] < 0) {
                argp = -argp;
            }
        } else {
            cross3d(nvec, evec, c);
            argp = atan2f(dot3d(c, h) / hmag, dot3d(nvec, evec));
        }
        cross3d(evec, r, c);
        nu = atan2f(dot3d(c, h) / hmag, dot3d(evec, r));
    }

    float E = 2.0f * atan2f(sqrtf(1.0f - e) * sinf(nu / 2.0f), sqrtf(1.0f + e) * cosf(nu / 2.0f));
    float M = E - e * sinf(E);

    const float twopi = (float)(2.0 * SGP4_PI);
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
 * \param[in]  dt       Time since element epoch (s), UTC-consistent clock. Take the
 *                      difference in double/integer time before casting to float.
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
int sgp4_propagate(const sgp4_sat_t *sat, float dt,
                   float *r_gcrf, float *v_gcrf, float *oe_osc, float *oe_mean,
                   float *r_teme, float *v_teme)
{
    memset(r_gcrf, 0, sizeof(float) * 3);
    memset(v_gcrf, 0, sizeof(float) * 3);
    memset(r_teme, 0, sizeof(float) * 3);
    memset(v_teme, 0, sizeof(float) * 3);
    memset(oe_osc, 0, sizeof(float) * 7);
    memset(oe_mean, 0, sizeof(float) * 6);

    if (sat->initErr != 0) {
        return sat->initErr;
    }

    const float twopi = (float)(2.0 * SGP4_PI);
    const float x2o3  = 2.0f / 3.0f;
    const float xke   = sat->xke;
    float tsince = dt / 60.0f; // minutes since epoch

    // secular gravity and atmospheric drag
    float xmdf   = sat->mo + sat->mdot * tsince;
    float argpdf = sat->argpo + sat->argpdot * tsince;
    float nodedf = sat->nodeo + sat->nodedot * tsince;
    float argpm  = argpdf;
    float mm     = xmdf;
    float t2     = tsince * tsince;
    float nodem  = nodedf + sat->nodecf * t2;
    float tempa  = 1.0f - sat->cc1 * tsince;
    float tempe  = sat->bstar * sat->cc4 * tsince;
    float templ  = sat->t2cof * t2;

    if (sat->isimp == 0) {
        float delomg   = sat->omgcof * tsince;
        float delmtemp = 1.0f + sat->eta * cosf(xmdf);
        float delm     = sat->xmcof * (delmtemp * delmtemp * delmtemp - sat->delmo);
        float temp     = delomg + delm;
        mm    = xmdf + temp;
        argpm = argpdf - temp;
        float t3 = t2 * tsince;
        float t4 = t3 * tsince;
        tempa = tempa - sat->d2 * t2 - sat->d3 * t3 - sat->d4 * t4;
        tempe = tempe + sat->bstar * sat->cc5 * (sinf(mm) - sat->sinmao);
        templ = templ + sat->t3cof * t3 + t4 * (sat->t4cof + tsince * sat->t5cof);
    }

    float nm    = sat->no_unkozai;
    float em    = sat->ecco;
    float inclm = sat->inclo;

    if (nm <= 0.0f) {
        return 2;
    }
    float am = powf(xke / nm, x2o3) * tempa * tempa;
    nm = xke / powf(am, 1.5f);
    em = em - tempe;

    if (em >= 1.0f || em < -0.001f) {
        return 1;
    }
    if (em < 1.0e-6f) {
        em = 1.0e-6f;
    }
    mm = mm + sat->no_unkozai * templ;
    float xlm = mm + argpm + nodem;

    nodem = fmodf(nodem, twopi);
    argpm = fmodf(argpm, twopi);
    xlm   = fmodf(xlm, twopi);
    mm    = fmodf(xlm - argpm - nodem, twopi);

    oe_mean[0] = am * sat->radiusearthkm;
    oe_mean[1] = em;
    oe_mean[2] = inclm;
    oe_mean[3] = mod_pos(nodem, twopi);
    oe_mean[4] = mod_pos(argpm, twopi);
    oe_mean[5] = mod_pos(mm, twopi);

    float sinip = sinf(inclm);
    float cosip = cosf(inclm);

    // long period periodics (near-Earth: no lunar-solar terms)
    float axnl = em * cosf(argpm);
    float temp = 1.0f / (am * (1.0f - em * em));
    float aynl = em * sinf(argpm) + temp * sat->aycof;
    float xl   = mm + argpm + nodem + temp * sat->xlcof * axnl;

    // solve Kepler's equation
    float u    = fmodf(xl - nodem, twopi);
    float eo1  = u;
    float tem5 = 9999.9f;
    int ktr = 1;
    float sineo1 = 0.0f, coseo1 = 1.0f;
    while (fabsf(tem5) >= 1.0e-6f && ktr <= 10) {
        sineo1 = sinf(eo1);
        coseo1 = cosf(eo1);
        tem5   = 1.0f - coseo1 * axnl - sineo1 * aynl;
        tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
        if (fabsf(tem5) >= 0.95f) {
            if (tem5 > 0.0f) {
                tem5 = 0.95f;
            } else {
                tem5 = -0.95f;
            }
        }
        eo1 = eo1 + tem5;
        ktr = ktr + 1;
    }

    // short period preliminary quantities
    float ecose = axnl * coseo1 + aynl * sineo1;
    float esine = axnl * sineo1 - aynl * coseo1;
    float el2   = axnl * axnl + aynl * aynl;
    float pl    = am * (1.0f - el2);
    if (pl < 0.0f) {
        return 4;
    }
    float rl     = am * (1.0f - ecose);
    float rdotl  = sqrtf(am) * esine / rl;
    float rvdotl = sqrtf(pl) / rl;
    float betal  = sqrtf(1.0f - el2);
    temp = esine / (1.0f + betal);
    float sinu  = am / rl * (sineo1 - aynl - axnl * temp);
    float cosu  = am / rl * (coseo1 - axnl + aynl * temp);
    float su    = atan2f(sinu, cosu);
    float sin2u = (cosu + cosu) * sinu;
    float cos2u = 1.0f - 2.0f * sinu * sinu;
    temp = 1.0f / pl;
    float temp1 = 0.5f * sat->j2 * temp;
    float temp2 = temp1 * temp;

    // update for short period periodics
    float mrt   = rl * (1.0f - 1.5f * temp2 * betal * sat->con41) + 0.5f * temp1 * sat->x1mth2 * cos2u;
    su = su - 0.25f * temp2 * sat->x7thm1 * sin2u;
    float xnode = nodem + 1.5f * temp2 * cosip * sin2u;
    float xinc  = inclm + 1.5f * temp2 * cosip * sinip * cos2u;
    float mvt   = rdotl - nm * temp1 * sat->x1mth2 * sin2u / xke;
    float rvdot = rvdotl + nm * temp1 * (sat->x1mth2 * cos2u + 1.5f * sat->con41) / xke;

    // orientation vectors
    float sinsu = sinf(su),    cossu = cosf(su);
    float snod  = sinf(xnode), cnod  = cosf(xnode);
    float sini  = sinf(xinc),  cosi  = cosf(xinc);
    float xmx   = -snod * cosi;
    float xmy   =  cnod * cosi;
    float uvec[3] = {xmx * sinsu + cnod * cossu, xmy * sinsu + snod * cossu, sini * sinsu};
    float vvec[3] = {xmx * cossu - cnod * sinsu, xmy * cossu - snod * sinsu, sini * cossu};

    for (int k = 0; k < 3; k++) {
        r_teme[k] = mrt * uvec[k] * sat->radiusearthkm;
        v_teme[k] = (mvt * uvec[k] + rvdot * vvec[k]) * sat->vkmpersec;
    }

    if (mrt < 1.0f) {
        return 6;
    }

    // TEME -> GCRF (held at the element epoch, see sgp4_init)
    matvec3(sat->R_teme2gcrf, r_teme, r_gcrf);
    matvec3(sat->R_teme2gcrf, v_teme, v_gcrf);

    rv2coe(r_gcrf, v_gcrf, 398600.4418f, oe_osc);
    return 0;
}
