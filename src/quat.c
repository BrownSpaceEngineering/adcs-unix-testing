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
void quat2rotm(float* q, float* rotm){
    float qw = q[0];
    float qx = q[1];
    float qy = q[2];
    float qz = q[3];

    // Precompute repeated terms
    float qw2 = qw * qw;
    float qx2 = qx * qx;
    float qy2 = qy * qy;
    float qz2 = qz * qz;

    float qxqy = qx * qy;
    float qxqz = qx * qz;
    float qyqz = qy * qz;
    float qwqx = qw * qx;
    float qwqy = qw * qy;
    float qwqz = qw * qz;

    // Row-major 3x3 rotation matrix
    rotm[0] = 1.0f - 2.0f * (qy2 + qz2);
    rotm[1] = 2.0f * (qxqy - qwqz);
    rotm[2] = 2.0f * (qxqz + qwqy);

    rotm[3] = 2.0f * (qxqy + qwqz);
    rotm[4] = 1.0f - 2.0f * (qx2 + qz2);
    rotm[5] = 2.0f * (qyqz - qwqx);

    rotm[6] = 2.0f * (qxqz - qwqy);
    rotm[7] = 2.0f * (qyqz + qwqx);
    rotm[8] = 1.0f - 2.0f * (qx2 + qy2);
}
void rotm_to_quat(float *R, float *q)
{
    float trace = R[0] + R[4] + R[8];

    if (trace > 0.0f)
    {
        float s = sqrtf(trace + 1.0f) * 2.0f; // s = 4*qw
        q[0] = 0.25f * s;
        q[1] = (R[7] - R[5]) / s;
        q[2] = (R[2] - R[6]) / s;
        q[3] = (R[3] - R[1]) / s;
    }
    else if ((R[0] > R[4]) && (R[0] > R[8]))
    {
        float s = sqrtf(1.0f + R[0] - R[4] - R[8]) * 2.0f; // s = 4*qx
        q[0] = (R[7] - R[5]) / s;
        q[1] = 0.25f * s;
        q[2] = (R[1] + R[3]) / s;
        q[3] = (R[2] + R[6]) / s;
    }
    else if (R[4] > R[8])
    {
        float s = sqrtf(1.0f + R[4] - R[0] - R[8]) * 2.0f; // s = 4*qy
        q[0] = (R[2] - R[6]) / s;
        q[1] = (R[1] + R[3]) / s;
        q[2] = 0.25f * s;
        q[3] = (R[5] + R[7]) / s;
    }
    else
    {
        float s = sqrtf(1.0f + R[8] - R[0] - R[4]) * 2.0f; // s = 4*qz
        q[0] = (R[3] - R[1]) / s;
        q[1] = (R[2] + R[6]) / s;
        q[2] = (R[5] + R[7]) / s;
        q[3] = 0.25f * s;
    }

    // Optional: normalize to protect against numerical drift
    float norm = sqrtf(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    if (norm > 0.0f)
    {
        q[0] /= norm;
        q[1] /= norm;
        q[2] /= norm;
        q[3] /= norm;
    }
}