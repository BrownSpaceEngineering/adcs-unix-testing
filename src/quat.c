#include "include/quat.h"
#include "math.h"
#include <string.h>

/*Applies the Hamilton product*/
void quat_multiply(float* q_left, float* q_right, float* resulting_quat) {
    float res[4];
    res[0] = q_left[0] * q_right[0] - q_left[1] * q_right[1] - q_left[2] * q_right[2]
                        - q_left[3] * q_right[3];
    res[1] = q_left[0] * q_right[1] + q_left[1] * q_right[0] + q_left[2] * q_right[3]
                        - q_left[3] * q_right[2];
    res[2] = q_left[0] * q_right[2] - q_left[1] * q_right[3] + q_left[2] * q_right[0]
                        + q_left[3] * q_right[1];
    res[3] = q_left[0] * q_right[3] + q_left[1] * q_right[2] - q_left[2] * q_right[1]
                        + q_left[3] * q_right[0];
    memcpy(resulting_quat, res, sizeof(float) * 4);
}
/*Out of place normalization of q*/
void quat_norm(float* q, float* resulting_quat) {
    float mag = quat_mag(q);
    resulting_quat[0] = q[0] / mag;
    resulting_quat[1] = q[1] / mag;
    resulting_quat[2] = q[2] / mag;
    resulting_quat[3] = q[3] / mag;
}
/*Finds the magnitude of q*/
float quat_mag(float* q) { return sqrtf(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]); }
/*Out-of-place conjugate of q*/
void quat_conj(float* q, float* resulting_quat) {
    resulting_quat[0] = q[0];
    resulting_quat[1] = -q[1];
    resulting_quat[2] = -q[2];
    resulting_quat[3] = -q[3];
}
/*Inverts q. For unit quaternions, this is the same as taking the conjugate*/
void quat_inv(float* q, float* resulting_quat) {
    float res[4];
    float mag = quat_mag(q);
    quat_conj(q, res);
    res[0] = res[0] / (mag * mag);
    res[1] = res[1] / (mag * mag);
    res[2] = res[2] / (mag * mag);
    res[3] = res[3] / (mag * mag);
    memcpy(resulting_quat, res, sizeof(float) * 4);
}
/*Applies the quaternion q to a vector out-of-place*/
void quat_apply(float* q, float* vec, float* resulting_vec) {
    float pure_vec[4] = {0, vec[0], vec[1], vec[2]}; // pure quaternion form of vec

    // q*pv*q^-1
    float q_inv[4];
    quat_inv(q, q_inv);

    float step_1[4];
    quat_multiply(q, pure_vec, step_1);

    float step_2[4];
    quat_multiply(step_1, q_inv, step_2);

    resulting_vec[0] = step_2[1];
    resulting_vec[1] = step_2[2];
    resulting_vec[2] = step_2[3];
}

/*Calculates the quaternion needed to rotate from_q to get to to_q */
void quat_diff(float* from_q, float* to_q, float* resulting_quat) {
    float from_q_inv[4];
    quat_inv(from_q, from_q_inv);
    quat_multiply(to_q, from_q_inv, resulting_quat);
}

/*Converts an angle-axis rotation vector to a rotation quaternion*/
void rotationvec2quat(float* vec, float* resulting_quat) {
    float angle = sqrtf(powf(vec[0], 2) + powf(vec[1], 2) + powf(vec[2], 2));
    if (angle < 1e-6f) {
        resulting_quat[0] = 1;
        resulting_quat[1] = 0;
        resulting_quat[2] = 0;
        resulting_quat[3] = 0;
    } else {
        float axis[3] = {vec[0] / angle, vec[1] / angle, vec[2] / angle};
        resulting_quat[0] = cosf(angle / 2);
        resulting_quat[1] = sinf(angle / 2) * axis[0];
        resulting_quat[2] = sinf(angle / 2) * axis[1];
        resulting_quat[3] = sinf(angle / 2) * axis[2];
    }
}
/*Converts q to an angle-axis rotation vector*/
void quat2rotationvec(float* q, float* resulting_vec) {
    if (q[0] < 0) {
        q[0] = -q[0];
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
    }

    if (q[0] >= 1) {
        resulting_vec[0] = 0;
        resulting_vec[1] = 0;
        resulting_vec[2] = 0;
        return;
    }
    
    float theta = 2 * acos(q[0]);
    if (theta == 0) {
        resulting_vec[0] = 0;
        resulting_vec[1] = 0;
        resulting_vec[2] = 0;
        return;
    }
    resulting_vec[0] = theta * q[1] / sqrtf(1 - pow(q[0], 2));
    resulting_vec[1] = theta * q[2] / sqrtf(1 - pow(q[0], 2));
    resulting_vec[2] = theta * q[3] / sqrtf(1 - pow(q[0], 2));
}