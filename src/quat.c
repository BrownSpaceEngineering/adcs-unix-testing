#include "include/quat.h"
#include "Include/arm_math_types.h"
#include "Include/dsp/quaternion_math_functions.h"
#include "arm_math.h"
#include "math.h"
#include <string.h>

/**
 * \fn quat_multiply
 * \brief Hamilton product q_left * q_right. Output may alias either input.
 */
void quat_multiply(const float32_t* q_left, const float32_t* q_right, float32_t* resulting_quat) {
    float32_t out[4];
    arm_quaternion_product_single_f32(q_left, q_right, out);
    memcpy(resulting_quat, out, sizeof(out));
}

/**
 * \fn quat_norm
 * \brief Normalizes a quaternion. A zero quaternion becomes the identity instead of NaN.
 */
void quat_norm(const float32_t* q, float32_t* resulting_quat) {
    float32_t n2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
    if (!(n2 > 0.0f)) {
        resulting_quat[0] = 1.0f;
        resulting_quat[1] = 0.0f;
        resulting_quat[2] = 0.0f;
        resulting_quat[3] = 0.0f;
        return;
    }
    float32_t out[4];
    arm_quaternion_normalize_f32(q, out, 1);
    memcpy(resulting_quat, out, sizeof(out));
}

/**
 * \fn quat_conj
 * \brief Conjugate [w, -x, -y, -z]. Equal to the inverse for unit quaternions.
 */
void quat_conj(const float32_t* q, float32_t* resulting_quat) {
    resulting_quat[0] = q[0];
    resulting_quat[1] = -q[1];
    resulting_quat[2] = -q[2];
    resulting_quat[3] = -q[3];
}

/**
 * \fn quat_inv
 * \brief Inverts a (not necessarily unit) quaternion.
 */
void quat_inv(const float32_t* q, float32_t* resulting_quat) {
    float32_t out[4];
    arm_quaternion_inverse_f32(q, out, 1);
    memcpy(resulting_quat, out, sizeof(out));
}

/**
 * \fn quat_apply
 *
 * \brief Rotates vec by q (active): q * [0, vec] * q^-1. Output may alias vec.
 */
void quat_apply(const float32_t* q, const float32_t* vec, float32_t* resulting_vec) {
    float32_t pure_vec[4] = {0, vec[0], vec[1], vec[2]};

    quat_t q_inv;
    arm_quaternion_inverse_f32(q, q_inv, 1);

    float32_t step_1[4];
    arm_quaternion_product_single_f32(q, pure_vec, step_1);

    float32_t step_2[4];
    arm_quaternion_product_single_f32(step_1, q_inv, step_2);

    resulting_vec[0] = step_2[1];
    resulting_vec[1] = step_2[2];
    resulting_vec[2] = step_2[3];
}

/**
 * \fn quat_diff
 *
 * \brief Calculates the quaternion needed to rotate from_q to get to to_q, i.e.
 * to_q * from_q^-1 (matches Quaternion.quat_diff in the MATLAB code).
 */
void quat_diff(const float32_t* from_q, const float32_t* to_q, float32_t* resulting_quat) {
    quat_t from_q_inv;
    arm_quaternion_inverse_f32(from_q, from_q_inv, 1);
    quat_multiply(to_q, from_q_inv, resulting_quat);
}

/**
 * \fn rotationvec2quat
 *
 * \brief Converts an angle-axis rotation vector (radians) to a unit rotation quaternion.
 */
void rotationvec2quat(const float32_t* vec, float32_t* resulting_quat) {
    float32_t angle = sqrtf(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    if (angle < 1e-6f) {
        // First-order expansion keeps tiny rotations instead of dropping them
        float32_t q[4] = {1.0f, 0.5f * vec[0], 0.5f * vec[1], 0.5f * vec[2]};
        quat_norm(q, resulting_quat);
        return;
    }
    float32_t half_angle = angle * 0.5f;
    float32_t s = sinf(half_angle) / angle;
    resulting_quat[0] = cosf(half_angle);
    resulting_quat[1] = s * vec[0];
    resulting_quat[2] = s * vec[1];
    resulting_quat[3] = s * vec[2];
}

/**
 * \fn quat2rotationvec
 *
 * \brief Converts a quaternion to the shortest angle-axis rotation vector (radians).
 *
 * Uses 2*atan2(|v|, w) rather than 2*acos(w): acos loses all precision near w = 1 in
 * float32 (anything under ~1e-3 rad rounds to zero), which is exactly the regime the
 * UKF's error quaternions live in. The input is not modified.
 */
void quat2rotationvec(const float32_t* q, float32_t* resulting_vec) {
    // q and -q are the same rotation; pick the one with w >= 0 (shortest path)
    float32_t sign = (q[0] < 0.0f) ? -1.0f : 1.0f;
    float32_t w = sign * q[0];
    float32_t x = sign * q[1];
    float32_t y = sign * q[2];
    float32_t z = sign * q[3];

    float32_t s = sqrtf(x * x + y * y + z * z);
    float32_t k;
    if (s < 1e-12f) {
        k = (w > 0.0f) ? 2.0f / w : 0.0f;
    } else {
        k = 2.0f * atan2f(s, w) / s;
    }
    resulting_vec[0] = k * x;
    resulting_vec[1] = k * y;
    resulting_vec[2] = k * z;
}

/**
 * \fn quat2rotm
 * \brief Converts a quaternion to a (row-major) rotation matrix R with R*v == quat_apply(q, v).
 */
void quat2rotm(const float32_t* q, float32_t* rotm) {
    float32_t qn[4];
    quat_norm(q, qn);
    arm_quaternion2rotation_f32(qn, rotm, 1);
}

/**
 * \fn rotm_to_quat
 * \brief Converts a (row-major, proper) rotation matrix to a unit quaternion with w >= 0.
 */
void rotm_to_quat(const float32_t* R, float32_t* q) {
    arm_rotation2quaternion_f32(R, q, 1);
    if (q[0] < 0.0f) {
        q[0] = -q[0];
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
    }
    quat_norm(q, q);
}
