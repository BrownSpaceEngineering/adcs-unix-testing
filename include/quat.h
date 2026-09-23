#ifndef QUAT_H
#define QUAT_H
#include "arm_math.h"

typedef float32_t quat_t[4];
//All quaternions are assumed to be in WXYZ format
void quat_multiply(const float32_t q_left[4], const float32_t q_right[4],
                   float32_t resulting_quat[4]);
void quat_norm(const float32_t q[4], float32_t resulting_quat[4]);
// float quat_mag(float* q);
// void quat_conj(float* q, float* resulting_quat);
void quat_inv(const float32_t q[4], float32_t resulting_quat[4]);
void quat_apply(const float32_t q[4], const float32_t vec[3], float32_t resulting_vec[3]);
void quat_diff(const float32_t from_q[4], const float32_t to_q[4],
               float32_t resulting_quat[4]);
void rotationvec2quat(const float32_t vec[3], float32_t resulting_quat[4]);
void quat2rotationvec(const float32_t q[4], float32_t resulting_vec[3]);
void quat2rotm(const float32_t q[4], float32_t rotm[9]);
void rotm_to_quat(const float32_t R[9], float32_t q[4]);
#endif
