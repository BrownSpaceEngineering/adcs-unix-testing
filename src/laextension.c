#include "include/laextension.h"
#include "Include/arm_math_types.h"
#include "Include/dsp/matrix_functions.h"
#include "Include/dsp/statistics_functions.h"
#include "math.h"
#include "arm_math.h"

/**
 * \fn cross
 * 
 * \brief Compute the cross product of two 1x3 vectors
 *
 * \param[in] A: First vector
 * \param[in] B: Second vector
 * \param[out] C: Resultant vector
 * 
 * \return void
 */
void cross(float32_t *A, float32_t *B, float32_t *C){
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
}

/*Max*/
/**
 * \fn max_arr
 * 
 * \brief Find the maximum value in an array
 * 
 * \param[in] A: Array of floats
 * \param[in] elements: Number of elements in the array
 * 
 * \return float32_t, maximum value in the array
 */
inline float32_t max_arr(float32_t *A, int elements){
    
    float32_t max; 

    if(elements == 0){
        return 0;
    } else {
        arm_max_no_idx_f32(A, elements, &max);
    }
    return max;
}


/**
 * \fn min_arr 
 * 
 * \brief returns the minimum value in a flat array
 * 
 * \param[in] A: Array of floats
 * \param[in] elements: Number of elements in the array
 * 
 * \return float32_t, minimum value in the array 
 */
inline float32_t min_arr(float32_t *A, int elements){

    float32_t min;

    if(elements == 0){
        return 0;
    } else {
        arm_min_no_idx_f32(A, elements, &min);
    }
    return min;
}


/**
 * \fn trace
 * 
 * \brief Compute the trace of a square matrix
 * 
 * \param[in] A: Square matrix
 * \param[in] sz: Size of the matrix
 * 
 * \return float, trace of the matrix
 */
float32_t trace(arm_matrix_instance_f32* A){
    
    if (A->numRows != A->numCols){
        return 0;
    }
    
    int sz = A->numRows;
    float32_t *data = A->pData;
    float32_t trace = 0; 

    for (int i = 0; i < sz; i++) {
        trace += data[i * sz + i];
    }
    
    return trace;
}

/**
 * \fn l2_norm
 * 
 * \brief Compute the L2 norm of a vector
 * 
 * \param[in] A: Vector
 * \param[in] n: Number of elements in the vector
 * 
 * \return float, L2 norm of the vector
 */
inline float32_t l2_norm(float32_t *A, int n){
    
    float32_t sum_squares;
    float32_t norm;

    arm_power_f32(A, n, &sum_squares);
   
    arm_sqrt_f32(sum_squares, &norm);
    
    return norm;
}