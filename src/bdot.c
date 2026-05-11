#include "include/bdot.h"
#include "Include/arm_math_types.h"
#include "Include/dsp/matrix_functions.h"
#include "string.h"
#include "arm_math.h"

/**
 * \fn Bdot 
 * 
 * \brief Bdot algorithm for detumbling 
 * 
 * \param[in] M_T Current magnetic field vector in body frame
 * \param[in] M_TMINUS1 Previous magnetic field vector in body frame
 * \param[in] dT Time step between measurements
 * \param[out] moments Output magnetic moments to be applied to the reaction wheels
 */
void Bdot(float32_t* M_T, float32_t* M_TMINUS1, float32_t dT, float32_t* moments){

    float32_t k[3] = {1.0f / dT, 1.0f / dT, 1.0f / dT}; // Adjust these gains as needed
    float32_t working_array[3];

    arm_sub_f32(M_T, M_TMINUS1, working_array, 3); // working_array = M_T - M_TMINUS1

    arm_mult_f32(k, working_array, moments, 3); // moments = k * (M_T - M_TMINUS1)
}