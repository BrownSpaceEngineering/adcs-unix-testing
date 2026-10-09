#include "include/down_quat.h"
#include "include/laextension.h"
#include "include/quat.h"
#include "arm_math.h"
#include "math.h"
#include <string.h>

#define EPS 1e-6f

static const float IDENTITY_Q[4] = {1.0f, 0.0f, 0.0f, 0.0f};

// Row-major rotation matrix with the given columns (body axes expressed in ECI)
static void rotm_from_columns(const float* x, const float* y, const float* z, float* R) {
    R[0] = x[0];
    R[1] = y[0];
    R[2] = z[0];
    R[3] = x[1];
    R[4] = y[1];
    R[5] = z[1];
    R[6] = x[2];
    R[7] = y[2];
    R[8] = z[2];
}

/**
 * \fn down_quat
 *
 * \brief Computes the body->ECI quaternion with the body Z axis pointing at the target.
 *
 * Fixes vs the previous version: the Z axis was from - to (pointing AWAY from the target), and
 * x = z cross y made a left-handed frame (det = -1), which is not a rotation at all.
 *
 * \param[in] from Current position vector of the satellite in ECI frame
 * \param[in] to Target position vector (e.g. Providence) in ECI frame
 * \param[in] q_body_to_eci Current orientation of the satellite (Body -> ECI)
 * \param[out] goal_q Desired orientation quaternion (Body -> ECI)
 */
void down_quat(const float* from, const float* to, const float* q_body_to_eci, float* goal_q) {
    float z[3] = {to[0] - from[0], to[1] - from[1], to[2] - from[2]};
    if (!normalize_vec(z, 3)) {
        memcpy(goal_q, IDENTITY_Q, sizeof(IDENTITY_Q));
        return;
    }

    // Current body Y-axis in ECI (2nd column of the body->ECI rotation matrix)
    float R[9];
    quat2rotm(q_body_to_eci, R);
    float y[3] = {R[1], R[4], R[7]};

    // Project y onto the plane perpendicular to z
    float proj = dot3(y, z);
    y[0] -= proj * z[0];
    y[1] -= proj * z[1];
    y[2] -= proj * z[2];

    if (!(l2_norm(y, 3) > EPS)) {
        // y nearly parallel to z: pick any perpendicular vector
        if (fabsf(z[0]) < 0.9f) {
            y[0] = 0.0f;
            y[1] = -z[2];
            y[2] = z[1];
        } else {
            y[0] = -z[1];
            y[1] = z[0];
            y[2] = 0.0f;
        }
    }
    normalize_vec(y, 3);

    // Right-handed: x = y cross z, then re-orthogonalize y = z cross x
    float x[3];
    cross(y, z, x);
    normalize_vec(x, 3);
    cross(z, x, y);

    float new_R[9];
    rotm_from_columns(x, y, z, new_R);
    rotm_to_quat(new_R, goal_q);
}

/**
 * \fn pointing_error
 *
 * \brief Port of pointing_error.m (see header).
 */
bool pointing_error(const float* r_eci, const float* v_eci, const float* q_b2eci,
                    const float* target_eci, float* q_want_prev, float* q_err, float* z_want) {
    float z[3] = {target_eci[0] - r_eci[0], target_eci[1] - r_eci[1], target_eci[2] - r_eci[2]};
    if (!normalize_vec(z, 3)) {
        memcpy(q_err, IDENTITY_Q, sizeof(IDENTITY_Q));
        return false;
    }
    if (z_want != NULL) {
        memcpy(z_want, z, sizeof(z));
    }

    float y[3];
    cross(z, v_eci, y);
    if (!(l2_norm(y, 3) > EPS * l2_norm(v_eci, 3)) || !normalize_vec(y, 3)) {
        // Velocity (nearly) along the line of sight; MATLAB would divide by zero here
        memcpy(q_err, IDENTITY_Q, sizeof(IDENTITY_Q));
        return false;
    }

    float x[3];
    cross(y, z, x);
    normalize_vec(x, 3);

    float R_want[9];
    rotm_from_columns(x, y, z, R_want);
    float q_want[4];
    rotm_to_quat(R_want, q_want);

    if (q_want_prev != NULL) {
        float d = q_want[0] * q_want_prev[0] + q_want[1] * q_want_prev[1]
                  + q_want[2] * q_want_prev[2] + q_want[3] * q_want_prev[3];
        if (d < 0.0f) {
            for (int i = 0; i < 4; i++) {
                q_want[i] = -q_want[i];
            }
        }
        memcpy(q_want_prev, q_want, sizeof(q_want));
    }

    float q_inv[4];
    quat_inv(q_b2eci, q_inv);
    quat_multiply(q_inv, q_want, q_err);
    return true;
}
