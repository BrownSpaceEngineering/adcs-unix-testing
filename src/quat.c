#include "Include/arm_math_types.h"
#include "Include/dsp/quaternion_math_functions.h"
#include "arm_math.h"
#include "include/quat.h"
#include "Include/dsp/matrix_functions.h"
#include "Include/dsp/statistics_functions.h"
#include "Include/dsp/basic_math_functions.h"
#include "math.h"
#include <string.h>

/**
 * \fn quat_multiply
 * \brief Hamilton product of scalar-first quaternions.
 *
 * For q_left=[a,u] and q_right=[b,v]:
 *
 *   result = q_left (*) q_right
 *          = [a b - u dot v, a v + b u + u cross v]
 *
 * \param q_left Pointer to the left quaternion.
 * \param q_right Pointer to the right quaternion.
 * \param resulting_quat Pointer to the output quaternion.
 */
void quat_multiply(const float32_t q_left[4], const float32_t q_right[4],
                   float32_t resulting_quat[4]) {
    /*
     * For q_left=[a,u] and q_right=[b,v]:
     *
     *   q_left (*) q_right = [a b - u dot v,
     *                         a v + b u + u cross v]
     */
    arm_quaternion_product_single_f32(q_left, q_right, resulting_quat);
}

/**
 * \fn quat_norm
 *
 * \brief Normalizes a quaternion.
 *
 *   result = q / ||q||_2
 *
 * \param q Pointer to the input quaternion.
 * \param resulting_quat Pointer to the output normalized quaternion.
 */
void quat_norm(const float32_t q[4], float32_t resulting_quat[4]) {
    /* q_unit = q / ||q||_2. */
    arm_quaternion_normalize_f32(q, resulting_quat, 1);
}

/**
 * \fn quat_inv
 * \brief Inverts a quaternion.
 *
 * For q=[w,x,y,z]:
 *
 *   result = q^-1 = [w,-x,-y,-z] / ||q||_2^2
 *
 * \param q Pointer to the input quaternion.
 * \param resulting_quat Pointer to the output inverted quaternion.
 */
void quat_inv(const float32_t q[4], float32_t resulting_quat[4]) {
    /* For q=[w,x,y,z], q^-1 = [w,-x,-y,-z] / ||q||_2^2. */
    arm_quaternion_inverse_f32(q, resulting_quat, 1);
}

/**
 * \fn quat_apply
 * 
 * \brief Applies the quaternion q to a vector out-of-place.
 *
 *   pure_vec     = [0, vec]
 *   rotated_quat = q (*) pure_vec (*) q^-1
 *   result       = rotated_quat[1:4]
 *
 * \param q Pointer to the input quaternion.
 * \param vec Pointer to the input vector.
 * \param resulting_vec Pointer to the output vector.
 */
void quat_apply(const float32_t q[4], const float32_t vec[3], float32_t resulting_vec[3]) {
    /* pure_vec = [0, vec]. */
    float32_t pure_vec[4] = {0, vec[0], vec[1], vec[2]};

    /* q_inv = q^-1. */
    quat_t q_inv;
    arm_quaternion_inverse_f32(q, q_inv, 1);

    /* step_1 = q (*) [0, vec]. */
    float32_t step_1[4];
    arm_quaternion_product_single_f32(q, pure_vec, step_1);

    /* step_2 = step_1 (*) q^-1 = [0, v_rotated]. */
    float32_t step_2[4];
    arm_quaternion_product_single_f32(step_1, q_inv, step_2);

    resulting_vec[0] = step_2[1];
    resulting_vec[1] = step_2[2];
    resulting_vec[2] = step_2[3];
}

/**
 * \fn quat_diff
 *
 * \brief Calculates the quaternion needed to rotate from_q to get to to_q
 *
 *   result = to_q (*) from_q^-1
 *   to_q   = result (*) from_q
 *
 * \param from_q Pointer to the input from quaternion.
 * \param to_q Pointer to the input to quaternion.
 * \param resulting_quat Pointer to the output resulting quaternion.
 */
void quat_diff(const float32_t from_q[4], const float32_t to_q[4],
               float32_t resulting_quat[4]) {
    quat_t from_q_inv;

    /* from_q_inv = from_q^-1. */
    arm_quaternion_inverse_f32(from_q, from_q_inv, 1);

    /* difference = to_q (*) from_q_inv, so to_q = difference (*) from_q. */
    arm_quaternion_product_single_f32(to_q, from_q_inv, resulting_quat);
}

/**
 * \fn rotationvec2quat
 *
 * \brief Converts an angle-axis rotation vector to a rotation quaternion.
 *
 *   theta  = ||vec||_2
 *   scale  = sin(theta/2) / theta
 *   rotvec2quat(vec) = [cos(theta/2), scale vec]
 *
 * Near zero, scale = 1/2 - theta^2/48 + O(theta^4).
 *
 * \param vec Pointer to the input rotation vector.
 * \param resulting_quat Pointer to the output quaternion.
 */
void rotationvec2quat(const float32_t vec[3], float32_t resulting_quat[4]) {
    /* theta = ||rotation_vector||_2. */
    float32_t angle = sqrtf(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    /*
     * scale = sin(theta/2) / theta.
     * Near zero: scale = 1/2 - theta^2/48 + O(theta^4).
     */
    float32_t scale = angle < 1.0e-4f ? 0.5f - angle * angle / 48.0f
                                     : sinf(0.5f * angle) / angle;

    /* rotvec2quat(vec) = [cos(theta/2), vec sin(theta/2)/theta]. */
    resulting_quat[0] = cosf(0.5f * angle);
    for (int i = 0; i < 3; ++i) resulting_quat[i + 1] = scale * vec[i];
}

/**
 * \fn quat2rotationvec 
 * 
 * \brief Converts quaternion to an angle-axis rotation vector
 * 
 *   unit        = q / ||q||_2, with its sign chosen so unit.w >= 0
 *   vector_norm = ||unit.xyz||_2
 *   scale       = 2 atan2(vector_norm, unit.w) / vector_norm
 *   quat2rotvec(q) = scale unit.xyz
 *
 * At vector_norm = 0, scale uses its limiting value of 2.
 *
 * \param q the quaternion to convert 
 * \param resulting_vec pointer to store the rotation vector in 
 */
void quat2rotationvec(const float32_t q[4], float32_t resulting_vec[3]) {
    /* unit = q / ||q||_2. */
    quat_t unit;
    quat_norm(q, unit);

    /* q and -q encode the same rotation; choose the representative with w >= 0. */
    if (unit[0] < 0.0f) {
        for (int i = 0; i < 4; ++i) unit[i] = -unit[i];
    }
    /* vector_norm = ||unit[1:4]||_2. */
    float32_t vector_norm = sqrtf(unit[1] * unit[1] + unit[2] * unit[2]
                                + unit[3] * unit[3]);

    /* scale = 2 atan2(vector_norm, w) / vector_norm; its zero limit is 2. */
    float32_t scale = vector_norm < 1.0e-6f ? 2.0f
                    : 2.0f * atan2f(vector_norm, unit[0]) / vector_norm;

    /* quat2rotvec(q) = scale * unit[1:4]. */
    for (int i = 0; i < 3; ++i) resulting_vec[i] = scale * unit[i + 1];
}


/**
 * \fn quat2rotm
 * \brief Converts a quaternion to a rotation matrix.
 *
 * For a unit scalar-first quaternion q=[w,x,y,z]:
 *
 *   R = [1-2(y^2+z^2), 2(xy-wz),     2(xz+wy);
 *        2(xy+wz),     1-2(x^2+z^2), 2(yz-wx);
 *        2(xz-wy),     2(yz+wx),     1-2(x^2+y^2)]
 *
 * \param q Pointer to the input quaternion.
 * \param rotm Pointer to the output rotation matrix (row-major).
 */
void quat2rotm(const float32_t q[4], float32_t rotm[9]) {
    /*
     * R = [1-2(y^2+z^2), 2(xy-wz),     2(xz+wy);
     *      2(xy+wz),     1-2(x^2+z^2), 2(yz-wx);
     *      2(xz-wy),     2(yz+wx),     1-2(x^2+y^2)]
     */
    arm_quaternion2rotation_f32(q, rotm, 1);
}


/**
 * \fn rotm_to_quat
 * \brief Converts a rotation matrix to a quaternion.
 *
 * Return a unit scalar-first q=[w,x,y,z] satisfying:
 *
 *   R = [1-2(y^2+z^2), 2(xy-wz),     2(xz+wy);
 *        2(xy+wz),     1-2(x^2+z^2), 2(yz-wx);
 *        2(xz-wy),     2(yz+wx),     1-2(x^2+y^2)]
 *
 * \param R Pointer to the input rotation matrix (row-major).
 * \param q Pointer to the output quaternion.
 */
void rotm_to_quat(const float32_t R[9], float32_t q[4]) {
    /* Find unit q such that quat2rotm(q) = R. */
    arm_rotation2quaternion_f32(R, q, 1);
}
