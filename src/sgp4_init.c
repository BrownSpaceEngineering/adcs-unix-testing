#include "include/sgp4_init.h"
#include "math.h"

// ================= TUNABLE CONSTANTS =================================
// Gravity model: 72 = WGS-72 (use this for NORAD/Space-Track TLEs --
// they are generated with WGS-72), 84 = WGS-84.
#define GRAV_MODEL 72
// TAI-UTC leap seconds (37 s since 2017-01-01; update if IERS adds one)
#define DELTA_AT 37.0
// =====================================================================

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
 * \param[in]  epoch_jd  Julian date (UTC) of the element epoch. Only used for the
 *                       slowly varying TEME->GCRF rotation, so its precision is not
 *                       critical; the propagation itself only sees dt.
 * \param[out] sat       Precomputed constants for sgp4_propagate.
 *                       sat->initErr = 0 ok, 7 = deep-space orbit (period >= 225 min),
 *                       not supported.
 */
void sgp4_init(const double *oe_epoch, double bstar, double epoch_jd, sgp4_sat_t *sat)
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
    const double xke   = 60.0 / sqrt(radiusearthkm * radiusearthkm * radiusearthkm / mu); // sqrt(mu) in er^1.5/min
    const double j3oj2 = j3 / j2;
    const double twopi = 2.0 * SGP4_PI;
    const double x2o3  = 2.0 / 3.0;
    const double temp4 = 1.5e-12;

    double no_kozai = oe_epoch[0] * twopi / 1440.0; // rev/day -> rad/min
    double ecco  = oe_epoch[1];
    double inclo = oe_epoch[2];
    double nodeo = oe_epoch[3];
    double argpo = oe_epoch[4];
    double mo    = oe_epoch[5];

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

    // TT Julian centuries since J2000 at epoch (for TEME -> GCRF)
    double T0 = (epoch_jd + (DELTA_AT + 32.184) / 86400.0 - 2451545.0) / 36525.0;

    sat->initErr       = initErr;
    sat->isimp         = isimp;
    sat->radiusearthkm = radiusearthkm;
    sat->xke           = xke;
    sat->j2            = j2;
    sat->vkmpersec     = radiusearthkm * xke / 60.0;
    sat->bstar         = bstar;
    sat->ecco          = ecco;
    sat->inclo         = inclo;
    sat->nodeo         = nodeo;
    sat->argpo         = argpo;
    sat->mo            = mo;
    sat->no_unkozai    = no_unkozai;
    sat->mdot          = mdot;
    sat->argpdot       = argpdot;
    sat->nodedot       = nodedot;
    sat->nodecf        = nodecf;
    sat->cc1           = cc1;
    sat->cc4           = cc4;
    sat->cc5           = cc5;
    sat->t2cof         = t2cof;
    sat->omgcof        = omgcof;
    sat->xmcof         = xmcof;
    sat->eta           = eta;
    sat->delmo         = delmo;
    sat->sinmao        = sinmao;
    sat->d2            = d2;
    sat->d3            = d3;
    sat->d4            = d4;
    sat->t3cof         = t3cof;
    sat->t4cof         = t4cof;
    sat->t5cof         = t5cof;
    sat->con41         = con41;
    sat->x1mth2        = x1mth2;
    sat->x7thm1        = x7thm1;
    sat->xlcof         = xlcof;
    sat->aycof         = aycof;
    sat->T0            = T0;
}
