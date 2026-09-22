#include "Include/arm_math_types.h"
#include "Include/dsp/quaternion_math_functions.h"
#include "arm_math.h"
#include "include/quat.h"
#include "Include/dsp/matrix_functions.h"
#include "Include/dsp/statistics_functions.h"
#include "Include/dsp/basic_math_functions.h"
#include "math.h"
#include <string.h>

/*Applies the Hamilton product*/
/**
 * \fn quat_multiply
 * \brief Multiplies two quaternions.
 * \param q_left Pointer to the left quaternion.
 * \param q_right Pointer to the right quaternion.
 * \param resulting_quat Pointer to the output quaternion.
 */
void quat_multiply(float32_t *q_left, float32_t *q_right, float32_t *resulting_quat) {
    arm_quaternion_product_single_f32(q_left, q_right, resulting_quat);
}

/**
 * \fn quat_norm
 *
 * \brief Normalizes a quaternion.
 *
 * \param q Pointer to the input quaternion.
 * \param resulting_quat Pointer to the output normalized quaternion.
 */
void quat_norm(float32_t *q, float32_t *resulting_quat) {
    arm_quaternion_normalize_f32(q, resulting_quat, 1);
}

/**
 * \fn quat_inv
 * \brief Inverts a quaternion.
 * \param q Pointer to the input quaternion.
 * \param resulting_quat Pointer to the output inverted quaternion.
 */
void quat_inv(float32_t *q, float32_t *resulting_quat) {
    arm_quaternion_inverse_f32(q, resulting_quat, 1);
}

/**
 * \fn quat_apply
 * 
 * \brief Applies the quaternion q to a vector out-of-place.
 *
 * \param q Pointer to the input quaternion.
 * \param vec Pointer to the input vector.
 * \param resulting_vec Pointer to the output vector.
 */
void quat_apply(quat_t q, float32_t *vec, float32_t *resulting_vec) {
    float32_t pure_vec[4] = {0, vec[0], vec[1], vec[2]}; // pure quaternion form of vec

    // q*pv*q^-1
    quat_t q_inv;
    arm_quaternion_inverse_f32(q, q_inv, 1);
    // quat_inv(q, q_inv);

    float32_t step_1[4];
    arm_quaternion_product_single_f32(q, pure_vec, step_1);
    // quat_multiply(q, pure_vec, step_1);

    float32_t step_2[4];
    arm_quaternion_product_single_f32(step_1, q_inv, step_2);
    // quat_multiply(step_1, q_inv, step_2);

    resulting_vec[0] = step_2[1];
    resulting_vec[1] = step_2[2];
    resulting_vec[2] = step_2[3];
}

/**
 * \fn quat_diff
 *
 * \brief Calculates the quaternion needed to rotate from_q to get to to_q
 *
 * \param from_q Pointer to the input from quaternion.
 * \param to_q Pointer to the input to quaternion.
 * \param resulting_quat Pointer to the output resulting quaternion.
 */
void quat_diff(quat_t from_q, quat_t to_q, quat_t resulting_quat) {
    
    quat_t from_q_inv;

    arm_quaternion_inverse_f32(from_q, from_q_inv, 1);
    arm_quaternion_product_single_f32(to_q, from_q_inv, resulting_quat);

    // quat_inv(from_q, from_q_inv);
    // quat_multiply(to_q, from_q_inv, resulting_quat);
}

/**
 * \fn rotationvec2quat
 *
 * \brief Converts an angle-axis rotation vector to a rotation quaternion.
 *
 * \param vec Pointer to the input rotation vector.
 * \param resulting_quat Pointer to the output quaternion.
 */
void rotationvec2quat(float32_t *vec, quat_t resulting_quat) {

    float32_t sum_squares;
    arm_power_f32(vec, 3, &sum_squares);
    float32_t angle;
    arm_sqrt_f32(sum_squares, &angle);
    // float angle = sqrtf(powf(vec[0], 2) + powf(vec[1], 2) + powf(vec[2], 2));
    if (angle < 1e-6f) {
        resulting_quat[0] = 1;
        resulting_quat[1] = 0;
        resulting_quat[2] = 0;
        resulting_quat[3] = 0;
    } else {
        // float axis[3] = {vec[0] / angle, vec[1] / angle, vec[2] / angle};
        float32_t axis[3] = {vec[0] / angle, vec[1] / angle, vec[2] / angle};
        float32_t half_angle = angle * 0.5f;

        float32_t working_array[3];

        resulting_quat[0] = arm_cos_f32(half_angle);
        working_array[0] = arm_sin_f32(half_angle);
        working_array[1] = working_array[0];
        working_array[2] = working_array[0];

        arm_mult_f32(working_array, axis, resulting_quat + 1, 3);

        // equivalent to:
        // resulting_quat[0] = cosf(angle / 2);
        // resulting_quat[1] = sinf(angle / 2) * axis[0];
        // resulting_quat[2] = sinf(angle / 2) * axis[1];
        // resulting_quat[3] = sinf(angle / 2) * axis[2];
    }
}

/**
 * \fn quat2rotationvec 
 * 
 * \brief Converts quaternion to an angle-axis rotation vector
 * 
 * \param q the quaternion to convert 
 * \param resulting_vec pointer to store the rotation vector in 
 */
void quat2rotationvec(quat_t q, float32_t* resulting_vec) {

    // NOTE: This is safe mathematically but not 
    // really in a programming sense, as 
    // -q represents the same rotation as q
    if (q[0] < 0) {
        q[0] = -q[0];
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
    }

    if (q[0] >= 1) {
        resulting_vec[0] = 0.0f;
        resulting_vec[1] = 0.0f;
        resulting_vec[2] = 0.0f;
        return;
    }
    
    float32_t theta = 2 * acosf(q[0]);
    if (fabsf(theta) < 1e-6f) {
        resulting_vec[0] = 0.0f;
        resulting_vec[1] = 0.0f;
        resulting_vec[2] = 0.0f;
        return;
    }

    float32_t sqrt_term;
    arm_sqrt_f32(1 - (q[0] * q[0]), &sqrt_term);
    float32_t denominators[3];
    denominators[0] = theta / sqrt_term;
    denominators[1] = denominators[0];
    denominators[2] = denominators[0];

    arm_mult_f32(q + 1, denominators, resulting_vec, 3);
    // should be equivalent to :
    // resulting_vec[0] = theta * q[1] / sqrtf(1 - pow(q[0], 2));
    // resulting_vec[1] = theta * q[2] / sqrtf(1 - pow(q[0], 2));
    // resulting_vec[2] = theta * q[3] / sqrtf(1 - pow(q[0], 2));
}


/**
 * \fn quat2rotm
 * \brief Converts a quaternion to a rotation matrix.
 * \param q Pointer to the input quaternion.
 * \param rotm Pointer to the output rotation matrix (row-major).
 */
void quat2rotm(float32_t *q, float32_t *rotm) {
    arm_quaternion2rotation_f32(q, rotm, 1);
}


/**
 * \fn rotm_to_quat
 * \brief Converts a rotation matrix to a quaternion.
 * \param R Pointer to the input rotation matrix (row-major).
 * \param q Pointer to the output quaternion.
 */
void rotm_to_quat(float32_t *R, float32_t *q) {
    arm_rotation2quaternion_f32(R, q, 1);
}