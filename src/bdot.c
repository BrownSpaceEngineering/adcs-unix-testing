#include "include/bdot.h"
#include "arm_math.h"

/**
 * \fn Bdot
 *
 * \brief B-dot detumbling law, port of bDot.m: m = -k * dB/dt.
 *
 * The previous version returned +dB/dt (no gain, wrong sign), which pumps energy into the
 * tumble instead of removing it.
 *
 * \param[in] M_T Current magnetic field vector in body frame
 * \param[in] M_TMINUS1 Previous magnetic field vector in body frame
 * \param[in] k Positive detumbling gain
 * \param[in] dT Time step between measurements (s)
 * \param[out] moments Commanded magnetic dipole moment (body frame)
 */
void Bdot(const float32_t* M_T, const float32_t* M_TMINUS1, float32_t k, float32_t dT,
          float32_t* moments) {
    if (!(dT > 0.0f)) {
        moments[0] = 0.0f;
        moments[1] = 0.0f;
        moments[2] = 0.0f;
        return;
    }
    float32_t scale = -k / dT;
    for (int i = 0; i < 3; i++) {
        moments[i] = scale * (M_T[i] - M_TMINUS1[i]);
    }
}
