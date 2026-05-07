#ifndef LAEXTENSION
#define LAEXTENSION
#include "arm_math.h"

void cross(float32_t *A, float32_t *B, float32_t *C); 
float32_t min_arr(float32_t* A, int elements);
float32_t max_arr(float32_t* A, int elements);
float33_t trace(arm_matrix_instance_f32* A); 
float32_t l2_norm(float32_t* A, int size);
#endif