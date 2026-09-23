#ifndef QUAT_H
#define QUAT_H
#include "arm_math.h"

typedef float32_t quat_t[4];

// All quaternions are Hamilton, scalar-first [w, x, y, z], ACTIVE convention:
// quat_apply(q, v) = q * v * q^-1. A "body->ref" quaternion therefore maps a
// vector written in body coordinates to the same vector in ref coordinates.
void quat_multiply(const float* q_left, const float* q_right, float* resulting_quat);
void quat_norm(const float* q, float* resulting_quat);
void quat_conj(const float* q, float* resulting_quat);
void quat_inv(const float* q, float* resulting_quat);
void quat_apply(const float* q, const float* vec, float* resulting_vec);
// to_q * from_q^-1: the rotation (expressed in the outer/reference frame) that takes from_q to to_q
void quat_diff(const float* from_q, const float* to_q, float* resulting_quat);
void rotationvec2quat(const float* vec, float* resulting_quat);
void quat2rotationvec(const float* q, float* resulting_vec);
void quat2rotm(const float* q, float* rotm);
void rotm_to_quat(const float* R, float* q);
#endif
