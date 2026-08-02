#include "include/laextension.h"
#include "math.h"
/*Cross Product for 2 1x3 vectors*/
void cross(double* A, double* B, double* C){
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
}
/*Max*/
double max_arr(double* A, int elements){
    if(elements == 0){
        return 0;
    }
    else{
        double max_found = A[0];
        for(int i = 0; i<elements; i++){
            max_found = fmaxf(max_found, A[i]);
        }
        return max_found;
    }
}
/*Mix*/
double min_arr(double* A, int elements){
    if(elements == 0){
        return 0;
    }
    else{
        double min_found = A[0];
        for(int i = 0; i<elements; i++){
            min_found = fminf(min_found, A[i]);
        }
        return min_found;
    }
}
/*Trace*/
double trace(double* A, int sz){
    double n= 0;
    for(int i = 0; i<sz; i++){
        n += A[i * sz + i];
    }
    return n;
}
/*L2Norm*/
double l2_norm(double* A, int n){
    double x = 0;
    for (int i = 0; i<n; i++){
        x += pow(A[i], 2);
    }
    return sqrt(x);
}
