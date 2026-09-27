#include "include/down_quat.h"
#include "include/laextension.h"
#include "include/quat.h"
#include "math.h"
#include "arm_math.h"

#define EPS 1e-6f

/**
 * \fn normalize3
 *
 * \brief Normalizes a 3D vector
 *
 * \param[in,out] v The vector to normalize
 */
static void normalize3(float *v) {
    float n = l2_norm(v, 3);
    if (n > EPS) {
        v[0] /= n;
        v[1] /= n;
        v[2] /= n;
    }
}
/**
 * \fn pointing_error
 *
 * \brief Computes the error quaternion between the current attitude and the attitude
 *        where the body Z-axis points at Providence and the body Y-axis is normal to
 *        the plane of the target direction and the velocity.
 *
 * \param[in]  r_eci           ECI coordinates of satellite position (3 elements)
 * \param[in]  v_eci           ECI velocity of satellite (3 elements)
 * \param[in]  q_b2eci         Current orientation quaternion, body to ECI, ACTIVE convention, WXYZ (4 elements)
 * \param[in]  providence_eci  ECI coordinates of Providence, Rhode Island (3 elements)
 * \param[out] q_tgtb          Error quaternion, with axis in body coordinates, WXYZ (4 elements)
 * \param[out] z_want          Desired body Z-axis in ECI, unit vector (3 elements)
 */
void down_quat(float* r_eci, float* v_eci, float* q_b2eci, float* providence_eci,
                    float* q_tgtb, float* z_want)
{
    // Previous desired quaternion, kept across calls for sign continuity
    static float q_want_prev[4] = {1.0f, 0.0f, 0.0f, 0.0f};

    // z_want = (providence - r) / |providence - r|
    z_want[0] = providence_eci[0] - r_eci[0];
    z_want[1] = providence_eci[1] - r_eci[1];
    z_want[2] = providence_eci[2] - r_eci[2];
    normalize3(z_want);

    // y_want = z_want x v
    float y_want[3];
    cross(z_want, v_eci, y_want);
    normalize3(y_want);

    // x_want = y_want x z_want
    float x_want[3];
    cross(y_want, z_want, x_want);
    normalize3(x_want);

    // R_want_eci = [x_want, y_want, z_want] (columns), stored row-major
    float R_want_eci[9];
    R_want_eci[0] = x_want[0]; R_want_eci[1] = y_want[0]; R_want_eci[2] = z_want[0];
    R_want_eci[3] = x_want[1]; R_want_eci[4] = y_want[1]; R_want_eci[5] = z_want[1];
    R_want_eci[6] = x_want[2]; R_want_eci[7] = y_want[2]; R_want_eci[8] = z_want[2];

    float q_want_eci[4];
    rotm_to_quat(R_want_eci, q_want_eci);

    // Keep the same hemisphere as the previous timestep
    float d = q_want_eci[0] * q_want_prev[0] + q_want_eci[1] * q_want_prev[1]
            + q_want_eci[2] * q_want_prev[2] + q_want_eci[3] * q_want_prev[3];
    if (d < 0.0f) {
        q_want_eci[0] = -q_want_eci[0];
        q_want_eci[1] = -q_want_eci[1];
        q_want_eci[2] = -q_want_eci[2];
        q_want_eci[3] = -q_want_eci[3];
    }

    // Store for next timestep
    q_want_prev[0] = q_want_eci[0];
    q_want_prev[1] = q_want_eci[1];
    q_want_prev[2] = q_want_eci[2];
    q_want_prev[3] = q_want_eci[3];

    // q_tgtb = inv(q_b2eci) * q_want_eci
    float q_b2eci_inv[4];
    quat_inv(q_b2eci, q_b2eci_inv);
    quat_multiply(q_b2eci_inv, q_want_eci, q_tgtb);
}
