#ifndef LAEXTENSION
#define LAEXTENSION
void scalar_multiply(float* A, float s, int elements);
void matrix_add(float* A, float* B, float* C, int elements);
void matrix_subtract(float* A, float* B, float* C, int elements);
void cross(float* A, float* B, float* C);
#endif