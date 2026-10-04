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
 * \param[in]  dt        Time since element epoch (s); take the difference in double/integer
 *                       time before casting to float
 * \param[out] r_gcrf    Position (km) in GCRF, 3 elements
 * \param[out] v_gcrf    Velocity (km/s) in GCRF, 3 elements
 * \param[out] oe_osc    Osculating Keplerian elements in GCRF, 7 elements
 *
 * \return errCode from sgp4_propagate (0 ok)
 */
int sgp4_step(const float *oe_epoch, float bstar, double epoch_jd, float dt,
              float *r_gcrf, float *v_gcrf, float *oe_osc)
{
    // Kept across calls so sgp4_init only reruns when new elements arrive
    static sgp4_sat_t sat;
    static float last_elems[7];
    static double last_epoch_jd;
    static bool initialized = false;

    float elems[7];
    memcpy(elems, oe_epoch, sizeof(float) * 6);
    elems[6] = bstar;

    // bitwise compare: any change to the held elements means a new uplink
    bool changed = !initialized || memcmp(elems, last_elems, sizeof(elems)) != 0 ||
                   memcmp(&epoch_jd, &last_epoch_jd, sizeof(epoch_jd)) != 0;

    if (changed) {
        memcpy(last_elems, elems, sizeof(elems));
        last_epoch_jd = epoch_jd;
        sgp4_init(oe_epoch, bstar, epoch_jd, &sat);
        initialized = true;
    }

    float oe_mean[6], r_teme[3], v_teme[3];
    return sgp4_propagate(&sat, dt, r_gcrf, v_gcrf, oe_osc, oe_mean, r_teme, v_teme);
}
