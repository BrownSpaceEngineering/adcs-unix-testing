#ifndef SGP4_INIT_H
#define SGP4_INIT_H

#define SGP4_PI 3.14159265358979323846

/**
 * Precomputed SGP4 constants produced by sgp4_init and consumed by sgp4_propagate.
 * sgp4_init computes these in double (it only runs once per uplink) and stores them as float,
 * so sgp4_propagate is all single precision. The float error grows ~130 m/day of dt, well
 * below SGP4's own ~1-3 km/day, as long as elements are uplinked daily.
 */
typedef struct {
    int initErr;        // 0 ok, 7 = deep-space orbit (period >= 225 min), not supported
    int isimp;          // 1 if perigee < 220 km (higher-order drag terms dropped)
    float radiusearthkm;
    float xke;
    float j2;
    float vkmpersec;
    float bstar;
    float ecco;
    float inclo;
    float nodeo;
    float argpo;
    float mo;
    float no_unkozai;
    float mdot;
    float argpdot;
    float nodedot;
    float nodecf;
    float cc1;
    float cc4;
    float cc5;
    float t2cof;
    float omgcof;
    float xmcof;
    float eta;
    float delmo;
    float sinmao;
    float d2;
    float d3;
    float d4;
    float t3cof;
    float t4cof;
    float t5cof;
    float con41;
    float x1mth2;
    float x7thm1;
    float xlcof;
    float aycof;
    float R_teme2gcrf[9]; // TEME -> GCRF rotation at epoch, row-major
} sgp4_sat_t;

void sgp4_init(const float *oe_epoch, float bstar, double epoch_jd, sgp4_sat_t *sat);
#endif
