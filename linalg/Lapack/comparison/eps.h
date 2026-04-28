#ifndef EPS_H
#define EPS_H

int f_eps_close(float a, float b, float epsilon);
int f_eps_close_default(float a, float b);
int f_eps_close_matrix(float* A, float* B, int row, int column, float epsilon);
int f_eps_close_matrix_default(float* A, float* B, int row, int column);
int dbl_eps_close(float a, float b, float epsilon);
int dbl_eps_close_default(float a, float b);
int dbl_eps_close_matrix(float* A, float* B, int row, int column, float epsilon);
int dbl_eps_close_matrix_default(float* A, float* B, int row, int column);

#endif // EPS_H