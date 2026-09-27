#ifndef SGP4_INIT_H
#define SGP4_INIT_H

#define SGP4_PI 3.14159265358979323846

/**
 * Precomputed SGP4 constants produced by sgp4_init and consumed by sgp4_propagate.
 * All values are double: SGP4 (and the Julian date epoch) needs more than float precision.
 */
typedef struct {
    int initErr;        // 0 ok, 7 = deep-space orbit (period >= 225 min), not supported
    int isimp;          // 1 if perigee < 220 km (higher-order drag terms dropped)
    double radiusearthkm;
    double xke;
    double j2;
    double vkmpersec;
    double bstar;
    double ecco;
    double inclo;
    double nodeo;
    double argpo;
    double mo;
    double no_unkozai;
    double mdot;
    double argpdot;
    double nodedot;
    double nodecf;
    double cc1;
    double cc4;
    double cc5;
    double t2cof;
    double omgcof;
    double xmcof;
    double eta;
    double delmo;
    double sinmao;
    double d2;
    double d3;
    double d4;
    double t3cof;
    double t4cof;
    double t5cof;
    double con41;
    double x1mth2;
    double x7thm1;
    double xlcof;
    double aycof;
    double T0;          // TT Julian centuries since J2000 at epoch
} sgp4_sat_t;

void sgp4_init(const double *oe_epoch, double bstar, double epoch_jd, sgp4_sat_t *sat);
#endif
