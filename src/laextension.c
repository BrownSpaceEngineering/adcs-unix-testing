#include "include/laextension.h"
/*In-place scalar multiplication*/
void scalar_multiply(float* A, float s, int elements){
    for(int i = 0; i < elements; i++){
        A[i] *= s;
    }
}
/*Element-wise addition*/
void matrix_add(float* A, float* B, float* C, int elements){
    for(int i = 0; i<elements; i++){
        C[i] = A[i] + B[i];
    }
}
/*Element-wise subtraction (A - B)*/
void matrix_subtract(float* A, float* B, float* C, int elements){
    for(int i = 0; i<elements; i++){
        C[i] = A[i] - B[i];
    }
}
/*Cross Product for 2 1x3 vectors*/
void cross(float* A, float* B, float* C){
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
}