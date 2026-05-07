#include "Include/arm_math_types.h"
#include "Include/dsp/quaternion_math_functions.h"
#include "arm_math.h"
#include "include/quat.h"
#include "Include/dsp/matrix_functions.h"
#include "math.h"
#include <string.h>

/*Applies the Hamilton product*/
/**
 * \fn quat_multiply
 * \warning DEPRECATED: use arm_quaternion_product_f32 instead for better performance.
 * \brief Multiplies two quaternions.
 * \param q_left Pointer to the left quaternion.
 * \param q_right Pointer to the right quaternion.
 * \param resulting_quat Pointer to the output quaternion.
 */
// void quat_multiply(float* q_left, float* q_right, float* resulting_quat) {
    // float res[4];
    // res[0] = q_left[0] * q_right[0] - q_left[1] * q_right[1] - q_left[2] * q_right[2]
                        // - q_left[3] * q_right[3];
    // res[1] = q_left[0] * q_right[1] + q_left[1] * q_right[0] + q_left[2] * q_right[3]
                        // - q_left[3] * q_right[2];
    // res[2] = q_left[0] * q_right[2] - q_left[1] * q_right[3] + q_left[2] * q_right[0]
                        // + q_left[3] * q_right[1];
    // res[3] = q_left[0] * q_right[3] + q_left[1] * q_right[2] - q_left[2] * q_right[1]
                        // + q_left[3] * q_right[0];
    // memcpy(resulting_quat, res, sizeof(float) * 4);
// }

/**
 * \fn quat_norm
 *
 * \warning DEPRECATED: use arm_quaternion_normalize_f32 instead for better performance.
 *
 * \brief Normalizes a quaternion.
 *
 * \param q Pointer to the input quaternion.
 * \param resulting_quat Pointer to the output normalized quaternion.
 */
// void quat_norm(arm_matrix_instance_f32 *q, arm_matrix_instance_f32* resulting_quat) {
// 
    // float32_t mag = quat_mag(q->pData);
    // arm_mat_scale_f32(const arm_matrix_instance_f32 *pSrc, int scale, arm_matrix_instance_f32 *pDst)
    // resulting_quat[0] = q[0] / mag;
    // resulting_quat[1] = q[1] / mag;
    // resulting_quat[2] = q[2] / mag;
    // resulting_quat[3] = q[3] / mag;
// }

/**
 * \fn quat_mag
 * \warning DEPRECATED: use arm_quaternion_norm_f32 instead for better performance.
 * \brief Finds the magnitude of q.
 * \param q Pointer to the input quaternion.
 * \return The magnitude of the quaternion.
 */
// float quat_mag(float* q) { 
    // return sqrtf(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]); 
// }

/**
 * \fn quat_conj
 * \warning DEPRECATED: use arm_quaternion_conjugate_f32 instead for better performance.
 * \brief Computes the conjugate of a quaternion.
 * \param q Pointer to the input quaternion.
 * \param resulting_quat Pointer to the output conjugate quaternion.
 */
// void quat_conj(float* q, float* resulting_quat) {
    // resulting_quat[0] = q[0];
    // resulting_quat[1] = -q[1];
    // resulting_quat[2] = -q[2];
    // resulting_quat[3] = -q[3];
// }

/**
 * \fn quat_inv
 * \warning DEPRECATED: use arm_quaternion_inverse_f32 instead for better performance.
 * \brief Inverts a quaternion.
 * \param q Pointer to the input quaternion.
 * \param resulting_quat Pointer to the output inverted quaternion.
 */
// void quat_inv(float* q, float* resulting_quat) {
    // float res[4];
    // float mag = quat_mag(q);
    // quat_conj(q, res);
    // res[0] = res[0] / (mag * mag);
    // res[1] = res[1] / (mag * mag);
    // res[2] = res[2] / (mag * mag);
    // res[3] = res[3] / (mag * mag);
    // memcpy(resulting_quat, res, sizeof(float) * 4);
//}

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
    quat_t q_inv[4];
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
    
    quat_t from_q_inv[4];

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

    float32_t angle = arm_sqrt_f32(
        arm_mult_f32(vec[0], vec[0]) + arm_mult_f32(vec[1], vec[1]) + 
        arm_mult_f32(vec[2], vec[2])
    );
    // float angle = sqrtf(powf(vec[0], 2) + powf(vec[1], 2) + powf(vec[2], 2));
    if (angle < 1e-6f) {
        resulting_quat[0] = 1;
        resulting_quat[1] = 0;
        resulting_quat[2] = 0;
        resulting_quat[3] = 0;
    } else {
        // float axis[3] = {vec[0] / angle, vec[1] / angle, vec[2] / angle};
        float32_t axis[3] = {vec[0] / angle, vec[1] / angle, vec[2] / angle};  
        float32_t half_angle = arm_mult_f32(angle, 0.5f);
        
        float32_t working_array[3];
        
        resulting_quat_quat[0] = arm_cos_f32(half_angle);
        working_array[0] = arm_sin_f32(half_angle);
        working_array[1] = working_array[0];
        working_array[2] = working_array[0];

        arm_mult_f32(working_array[0], axis[0], resulting_quat[1], 3);
        
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
    if (f_eps_close_default(theta, 0.0f)) {
        resulting_vec[0] = 0.0f;
        resulting_vec[1] = 0.0f;
        resulting_vec[2] = 0.0f;
        return;
    }

    float32_t denominators[3]; 
    denominators[0] = theta / arm_sqrt_f32(1 - (q[0] * q[0])); 
    denominators[1] = denominators[0]; 
    denominators[2] = denominators[0]; 

    arm_mult_f32(q[1], denominators, resulting_vec, 3);
    // should be equivalent to : 
    // resulting_vec[0] = theta * q[1] / sqrtf(1 - pow(q[0], 2));
    // resulting_vec[1] = theta * q[2] / sqrtf(1 - pow(q[0], 2));
    // resulting_vec[2] = theta * q[3] / sqrtf(1 - pow(q[0], 2));
}


/**
 * \fn quat2rotm
 * \warning DEPRECATED: use arm_quaternion2rotation_f32 instead for better performance.
 * \brief Converts a quaternion to a rotation matrix.
 * \param q Pointer to the input quaternion.
 * \param rotm Pointer to the output rotation matrix.
 */
// void quat2rotm(float* q, float* rotm){
    // float qw = q[0];
    // float qx = q[1];
    // float qy = q[2];
    // float qz = q[3];
// 
    // //Precompute repeated terms
    // float qw2 = qw * qw;
    // float qx2 = qx * qx;
    // float qy2 = qy * qy;
    // float qz2 = qz * qz;
// 
    // float qxqy = qx * qy;
    // float qxqz = qx * qz;
    // float qyqz = qy * qz;
    // float qwqx = qw * qx;
    // float qwqy = qw * qy;
    // float qwqz = qw * qz;
// 
    // // Row-major 3x3 rotation matrix
    // rotm[0] = 1.0f - 2.0f * (qy2 + qz2);
    // rotm[1] = 2.0f * (qxqy - qwqz);
    // rotm[2] = 2.0f * (qxqz + qwqy);
// 
    // rotm[3] = 2.0f * (qxqy + qwqz);
    // rotm[4] = 1.0f - 2.0f * (qx2 + qz2);
    // rotm[5] = 2.0f * (qyqz - qwqx);
// 
    // rotm[6] = 2.0f * (qxqz - qwqy);
    // rotm[7] = 2.0f * (qyqz + qwqx);
    // rotm[8] = 1.0f - 2.0f * (qx2 + qy2);
// }


/**
 * \fn rotm_to_quat
 * \warning DEPRECATED: use arm_rotation2quaternion_f32 instead for better performance.
 * \brief Converts a rotation matrix to a quaternion.
 * \param R Pointer to the input rotation matrix.
 * \param q Pointer to the output quaternion.
 */
// void rotm_to_quat(float *R, float *q)
// {
    // float trace = R[0] + R[4] + R[8];
// 
    // if (trace > 0.0f)
    // {
        // float s = sqrtf(trace + 1.0f) * 2.0f; // s = 4*qw
        // q[0] = 0.25f * s;
        // q[1] = (R[7] - R[5]) / s;
        // q[2] = (R[2] - R[6]) / s;
        // q[3] = (R[3] - R[1]) / s;
    // }
    // else if ((R[0] > R[4]) && (R[0] > R[8]))
    // {
        // float s = sqrtf(1.0f + R[0] - R[4] - R[8]) * 2.0f; // s = 4*qx
        // q[0] = (R[7] - R[5]) / s;
        // q[1] = 0.25f * s;
        // q[2] = (R[1] + R[3]) / s;
        // q[3] = (R[2] + R[6]) / s;
    // }
    // else if (R[4] > R[8])
    // {
        // float s = sqrtf(1.0f + R[4] - R[0] - R[8]) * 2.0f; // s = 4*qy
        // q[0] = (R[2] - R[6]) / s;
        // q[1] = (R[1] + R[3]) / s;
        // q[2] = 0.25f * s;
        // q[3] = (R[5] + R[7]) / s;
    // }
    // else
    // {
        // float s = sqrtf(1.0f + R[8] - R[0] - R[4]) * 2.0f; // s = 4*qz
        // q[0] = (R[3] - R[1]) / s;
        // q[1] = (R[2] + R[6]) / s;
        // q[2] = (R[5] + R[7]) / s;
        // q[3] = 0.25f * s;
    // }
// 
    // // Optional: normalize to protect against numerical drift
    // float norm = sqrtf(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    // if (norm > 0.0f)
    // {
        // q[0] /= norm;
        // q[1] /= norm;
        // q[2] /= norm;
        // q[3] /= norm;
    // }
// }