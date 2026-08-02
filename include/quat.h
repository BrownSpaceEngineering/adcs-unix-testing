#ifndef QUAT_H
#define QUAT_H

//All quaternions are assumed to be in WXYZ format
void quat_multiply(double* q_left, double* q_right, double* resulting_quat);
void quat_norm(double* q, double* resulting_quat);
double quat_mag(double* q);
void quat_conj(double* q, double* resulting_quat);
void quat_inv(double* q, double* resulting_quat);
void quat_apply(double* q, double* vec, double* resulting_vec);
void quat_diff(double* from_q, double* to_q, double* resulting_quat);
void rotationvec2quat(double* vec, double* resulting_quat);
void quat2rotationvec(double* q, double* resulting_vec);
void quat2rotm(double* q, double* rotm);
void rotm_to_quat(double*R, double*q);
#endif