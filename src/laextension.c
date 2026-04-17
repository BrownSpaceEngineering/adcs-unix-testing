#include "include/laextension.h"
#include "math.h"
/*Cross Product for 2 1x3 vectors*/
void cross(float* A, float* B, float* C){
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
}
/*Max*/
float max_arr(float* A, int elements){
    if(elements == 0){
        return 0;
    }
    else{
        float max_found = A[0];
        for(int i = 0; i<elements; i++){
            max_found = fmaxf(max_found, A[i]);
        }
        return max_found;
    }
}
/*Mix*/
float min_arr(float* A, int elements){
    if(elements == 0){
        return 0;
    }
    else{
        float min_found = A[0];
        for(int i = 0; i<elements; i++){
            min_found = fminf(min_found, A[i]);
        }
        return min_found;
    }
}