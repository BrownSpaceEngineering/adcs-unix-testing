#ifndef QUAT_H
#define QUAT_H

//All quaternions are assumed to be in WXYZ format
void quat_multiply(float* q_left, float* q_right, float* resulting_quat);
void quat_norm(float* q, float* resulting_quat);
float quat_mag(float* q);
void quat_conj(float* q, float* resulting_quat);
void quat_inv(float* q, float* resulting_quat);
void quat_apply(float* q, float* vec, float* resulting_vec);
void quat_diff(float* from_q, float* to_q, float* resulting_quat);
void rotationvec2quat(float* vec, float* resulting_quat);
void quat2rotationvec(float* q, float* resulting_vec);

#endif