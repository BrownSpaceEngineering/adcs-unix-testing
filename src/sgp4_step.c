#include "include/sgp4_step.h"
#include "include/sgp4_init.h"
#include "include/sgp4_propagate.h"
#include <stdbool.h>
#include <string.h>

/**
 * \fn sgp4_step
 *
 * \brief SGP4 entry point, called at the control rate (10 Hz).
 *
 * Hold the last uplinked elements on oe_epoch/bstar/epoch_jd and drive dt
 * with (current UTC - element epoch) in seconds. sgp4_init reruns only
 * when the held elements change (i.e. once per uplink); every other step
 * only runs sgp4_propagate.
 *
 * See sgp4_init / sgp4_propagate for units.
 *
 * \param[in]  oe_epoch  SGP4/TLE mean elements at epoch, 6 elements
 * \param[in]  bstar     B* drag term (1/earth radii)
 * \param[in]  epoch_jd  Julian date (UTC) of the element epoch
 * \param[in]  dt        Time since element epoch (s)
 * \param[out] r_gcrf    Position (km) in GCRF, 3 elements
 * \param[out] v_gcrf    Velocity (km/s) in GCRF, 3 elements
 * \param[out] oe_osc    Osculating Keplerian elements in GCRF, 7 elements
 *
 * \return errCode from sgp4_propagate (0 ok)
 */
int sgp4_step(const double *oe_epoch, double bstar, double epoch_jd, double dt,
              double *r_gcrf, double *v_gcrf, double *oe_osc)
{
    // Kept across calls so sgp4_init only reruns when new elements arrive
    static sgp4_sat_t sat;
    static double last_key[8];
    static bool initialized = false;

    double key[8];
    memcpy(key, oe_epoch, sizeof(double) * 6);
    key[6] = bstar;
    key[7] = epoch_jd;

    bool changed = !initialized;
    for (int k = 0; k < 8 && !changed; k++) {
        changed = (key[k] != last_key[k]);
    }

    if (changed) {
        memcpy(last_key, key, sizeof(key));
        sgp4_init(oe_epoch, bstar, epoch_jd, &sat);
        initialized = true;
    }

    double oe_mean[6], r_teme[3], v_teme[3];
    return sgp4_propagate(&sat, dt, r_gcrf, v_gcrf, oe_osc, oe_mean, r_teme, v_teme);
}
