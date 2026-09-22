#include "include/iterate.h"
#include "Include/arm_math_types.h"
#include "Include/dsp/matrix_functions.h"
#include "Include/dsp/basic_math_functions.h"
#include "include/laextension.h"
#include "math.h"
#include "include/quat.h"
#include "arm_math.h"
#include <assert.h>
#include <string.h>

/**
 * \fn ensure_psd
 * 
 * \brief Modifies a square matrix in-place to be positive semi-definite.
 * 
 * \param mat Pointer to a square matrix. 
 */
void ensure_psd(arm_matrix_instance_f32 *mat) {
    assert (mat->numRows == mat->numCols); // Must be square
    int size = mat->numRows; // Assuming square matrix, numRows == numCols

    float32_t *mat_data = mat->pData;
    float eps = 1e-6; // Adjust this if necessary (MATLAB used 1e-7)

    // 1. Force perfect symmetry: P = (P + P^T) / 2
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            // Average the mirrored elements
            float avg = (mat_data[i * size + j] + mat_data[j * size + i]) / 2.0f;
            
            // Assign the average back to both spots
            mat_data[i * size + j] = avg;
            mat_data[j * size + i] = avg;
        }
    }

    // 2. Add a tiny epsilon to the diagonal to ensure eigenvalues > 0
    // This physically guarantees Positive Definiteness as long as the 
    // variance hasn't completely blown up to negative infinity.
    for (int i = 0; i < size; i++) {
        mat_data[i * size + i] += eps;
        
        // Safety catch: if diagonal somehow became negative, hard-reset it
        if (mat_data[i * size + i] < 0.0f) {
            mat_data[i * size + i] = eps;
        }
    }
}

/**
 * \fn calculate_lambda
 * 
 * \brief Calculates the lambda parameter for sigma point generation.
 * 
 * \param n The number of state variables.
 * \param alpha The scaling parameter.
 *
 * \return The calculated lambda value.
 */
float32_t calculate_lambda(int n, float32_t alpha){
    int k = 3 - n;
    return alpha * alpha * (n + k) - n;
}

/**
 * \fn get_sigma_points
 * 
 * \brief Generates sigma points for the Unscented Kalman Filter.
 */
void get_sigma_points(float32_t lam, arm_matrix_instance_f32 *state, arm_matrix_instance_f32 *P, arm_matrix_instance_f32* sigmas){
    
    // instantiate a workiing copy of P so we don't mess with the original.
    float32_t P1_data [STATE_SIZE * STATE_SIZE];
    memcpy(P1_data, P->pData, sizeof(float32_t) * STATE_SIZE * STATE_SIZE);

    arm_matrix_instance_f32 P1 = {
        .numRows = STATE_SIZE,
        .numCols = STATE_SIZE,
        .pData = P1_data, 
    };   

    float32_t U_data[STATE_SIZE * STATE_SIZE];
    arm_matrix_instance_f32 U = {
        .numRows = STATE_SIZE,
        .numCols = STATE_SIZE,
        .pData = U_data, 
    }; 

    float32_t *sigmas_data = sigmas->pData;
    float32_t *state_data = state->pData;

    // TODO: is this safe to do in-place? I think so, but
    // we should check the implementation of arm_mat_scale_f32 to be sure.
    arm_mat_scale_f32(&P1, lam + STATE_SIZE, &P1);

    ensure_psd(&P1);
    
    // chol(P1, U, STATE_SIZE);

    arm_mat_cholesky_f32(&P1, &U);

    // TODO: figure out a way to do this with arm matrix functions.
    // populates the first  
    for(int i = 0; i<STATE_SIZE; i++){
        sigmas_data[i] = state_data[i];
    }
    for(int k = 0; k<STATE_SIZE; k++){
        int start = k * STATE_SIZE;
        int sigma_offset = STATE_SIZE;
        for(int i = 0; i<STATE_SIZE; i++){
            sigmas_data[sigma_offset + start + i] = state_data[i] + U_data[start + i];
        }
    }
    for(int k = 0; k<STATE_SIZE; k++){
        int start = k * STATE_SIZE;
        int sigma_offset = (STATE_SIZE + 1) * STATE_SIZE;
        for(int i = 0; i<STATE_SIZE; i++){
            sigmas_data[sigma_offset + start + i] = state_data[i] - U_data[start + i];
        }
    }
}

/**
 * \fn get_weights * 
 * 
 * \brief Calculates the mean and covariance weights for the Unscented Kalman Filter.
 * 
 * \param lambda The lambda parameter for sigma point generation.
 * \param n The number of state variables.
 * \param alpha The scaling parameter.
 * \param beta The distribution parameter.
 * \param cov_weights Output array for covariance weights (size NUM_SIGMAS).
 * \param mean_weights Output array for mean weights (size NUM_SIGMAS).
 * 
 */
void get_weights(float32_t lambda,
                 int n,
                 float32_t alpha,
                 float32_t beta,
                 float32_t *cov_weights,
                 float32_t *mean_weights) {

    float32_t c = (float32_t)0.5f / (n + lambda);

    for(int i = 0; i< NUM_SIGMAS; i++){
        cov_weights[i] = c;
        mean_weights[i] = c;
    }

    cov_weights[0] = lambda / (n + lambda) + (1 - alpha*alpha + beta);
    mean_weights[0] = lambda / (n + lambda);
}

/**
 * \fn error_sigmas_to_quat_sigmas
 *
 * \brief Converts error sigmas to quaternion sigmas.
 *
 * \param error_sigmas Pointer to the error sigmas array.
 * \param current_guess Pointer to the current quaternion guess.
 * \param quat_sigmas Pointer to the output quaternion sigmas array.
 */
void error_sigmas_to_quat_sigmas(float32_t *error_sigmas,
                                 float32_t *current_guess,
                                 float32_t *quat_sigmas){

    int n = NUM_SIGMAS;

    for(int i = 0; i<n; i++){
        float32_t error_rot[3] = {
            error_sigmas[i * STATE_SIZE],
            error_sigmas[i * STATE_SIZE + 1],
            error_sigmas[i * STATE_SIZE + 2]
        };
        float32_t error_quat[4];
        rotationvec2quat(error_rot, error_quat);
        float32_t new_rotation[4];
        quat_multiply(current_guess, error_quat, new_rotation);
        float32_t norm_rot[4];
        quat_norm(new_rotation, norm_rot);
        for(int j = 0; j<4; j++){
            quat_sigmas[i * (STATE_SIZE + 1) + j] = norm_rot[j];
        }
        for(int j = 4; j < STATE_SIZE + 1; j++){
            quat_sigmas[i * (STATE_SIZE + 1) + j] = error_sigmas[i * STATE_SIZE + (j - 1)];
        }
    }
}

/**
 * \fn propagate_sigmas
 *
 * \brief Propagates the quaternion sigmas through time.
 *
 * \param quat_sigmas Pointer to the quaternion sigmas array.
 * \param gyro Pointer to the gyroscope measurements array.
 * \param dt The time step.
 * \param propagated_sigmas Pointer to the output propagated sigmas array.
 */
void propagate_sigmas(float* quat_sigmas, float* gyro, float dt, float* propagated_sigmas){
    float p_sigmas[NUM_SIGMAS * (STATE_SIZE + 1)];
    for(int i = 0; i<NUM_SIGMAS; i++){
        int offset = (STATE_SIZE + 1) * i;
        float rotation[4];
        float bias[3];
        for(int j = 0; j<4; j++){
            rotation[j] = quat_sigmas[offset + j];
        }
        for(int j = 4; j<7; j++){
            bias[j - 4] = quat_sigmas[offset + j];
        }
        float omega[3];
        arm_scale_f32(bias, -1.0f, bias, 3);
        arm_add_f32(gyro, bias, omega, 3);
        arm_scale_f32(omega, dt, omega, 3);
        float delta_q[4];
        rotationvec2quat(omega, delta_q);

        float propagated_rot[4];
        quat_multiply(rotation, delta_q, propagated_rot);
        float norm_rot[4];
        quat_norm(propagated_rot, norm_rot);
        
        int p_offset = i * (STATE_SIZE + 1);
        for(int j = 0; j<4; j++){
            p_sigmas[j + p_offset] = norm_rot[j];
        }
        for(int j = 4; j<7; j++){
            p_sigmas[j + p_offset] = quat_sigmas[offset + j];
        }
    }
    memcpy(propagated_sigmas, p_sigmas, sizeof(float) * NUM_SIGMAS * (STATE_SIZE + 1));
}
void get_sigma_measurements(float* q_sigmas, float* ref, float* msmts){
    float m[NUM_SIGMAS * MSMT_SIZE];
    for(int i = 0; i< NUM_SIGMAS; i++){
        float sigma_rot[4];
        for(int j = 0; j<4; j++){
            sigma_rot[j] = q_sigmas[i * (STATE_SIZE + 1) + j];
        }
        float sigma_ref_to_body[4];
        quat_inv(sigma_rot, sigma_ref_to_body);
        float ref_1[3] = {ref[0], ref[1], ref[2]};
        float reading_1[3];
        quat_apply(sigma_ref_to_body, ref_1, reading_1);
        for(int j = 0; j<3; j++){
            m[i * MSMT_SIZE + j] = reading_1[j];
        }
        float ref_2[3] = {ref[3], ref[4], ref[5]};
        float reading_2[3];
        quat_apply(sigma_ref_to_body, ref_2, reading_2);
        for(int j = 0; j<3; j++){
            m[i * MSMT_SIZE + j + 3] = reading_2[j];
        }
    }
    memcpy(msmts, m, sizeof(float) * NUM_SIGMAS * MSMT_SIZE);
}
void get_error_vectors(int num_other_quats, float* other_quats, float* quat, float* error_vecs){
    float e[num_other_quats * 3];
    float x_inv[4];
    quat_inv(quat, x_inv);
    for(int i = 0; i<num_other_quats; i++){
        float other_quat[4];
        for(int j = 0; j<4; j++){
            other_quat[j] = other_quats[i * 4 + j];
        }
        float err_quat[4];
        quat_multiply(x_inv, other_quat, err_quat);
        float err_vec[3];
        quat2rotationvec(err_quat, err_vec);
        for(int j = 0; j<3; j++){
            e[i * 3 + j] = err_vec[j];
        }
    }
    memcpy(error_vecs, e, sizeof(float) * num_other_quats * 3);
}
void grad_descent(float* propagated_quats, float* current_quat, float* avg_quat, float* errors){
    float moving_avg[4];
    memcpy(moving_avg, current_quat, sizeof(float) * 4);
    
    for(int iter = 0; iter < 50; iter++){
        float error_vecs[3 * NUM_SIGMAS];
        get_error_vectors(NUM_SIGMAS, propagated_quats, moving_avg, error_vecs);
        memcpy(errors, error_vecs, sizeof(float) * 3 * NUM_SIGMAS);
        float avg_err[3] = {0,0,0};
        for(int i = 0; i<NUM_SIGMAS; i++){
            for(int j = 0; j < 3; j++){
                avg_err[j] += error_vecs[i * 3 + j] / NUM_SIGMAS;
            }
        }
        if(sqrtf(avg_err[0] * avg_err[0] + avg_err[1] * avg_err[1] + avg_err[2] * avg_err[2]) < 0.001){
            break;
        }
        float avg_error_quat[4];
        rotationvec2quat(avg_err, avg_error_quat);
        float new_moving_avg[4];
        memcpy(new_moving_avg, moving_avg, sizeof(float) * 4);
        quat_multiply(moving_avg, avg_error_quat, new_moving_avg);
        memcpy(moving_avg, new_moving_avg, sizeof(float) * 4);
        quat_norm(moving_avg, moving_avg);
    }
    memcpy(avg_quat, moving_avg, sizeof(float) * 4);
}


void iterate(float* error_state, float* quat_state, float* cov, float* body, float* ref, float* gyro, float* Q, float* R, float dt, float* new_err_state, float* new_quat_state, float* new_P){
    float P[STATE_SIZE * STATE_SIZE];
    memcpy(P, cov, sizeof(float) * STATE_SIZE * STATE_SIZE);

    arm_matrix_instance_f32 P_mat = {STATE_SIZE, STATE_SIZE, P};
    ensure_psd(&P_mat);

    float lambda = calculate_lambda(STATE_SIZE, ALPHA);
    float P_Q[STATE_SIZE * STATE_SIZE];
    arm_add_f32(P, Q, P_Q, STATE_SIZE * STATE_SIZE);
    float error_sigmas[STATE_SIZE * NUM_SIGMAS];
    arm_matrix_instance_f32 error_state_mat = {STATE_SIZE, 1, error_state};
    arm_matrix_instance_f32 P_Q_mat = {STATE_SIZE, STATE_SIZE, P_Q};
    arm_matrix_instance_f32 error_sigmas_mat = {STATE_SIZE, NUM_SIGMAS, error_sigmas};
    get_sigma_points(lambda, &error_state_mat, &P_Q_mat, &error_sigmas_mat);
    float quat_sigmas[(STATE_SIZE + 1) * NUM_SIGMAS];
    error_sigmas_to_quat_sigmas(error_sigmas, quat_state, quat_sigmas);
    float propagated_sigmas[(STATE_SIZE + 1) * NUM_SIGMAS];
    propagate_sigmas(quat_sigmas, gyro, dt, propagated_sigmas);
    

    float sigma_msmts[MSMT_SIZE * NUM_SIGMAS];
    get_sigma_measurements(propagated_sigmas, ref, sigma_msmts);

    float mean_weights[NUM_SIGMAS];
    float cov_weights[NUM_SIGMAS];
    get_weights(lambda, STATE_SIZE, ALPHA, BETA, cov_weights, mean_weights);

    float quaternions_in_sigmas[4 * NUM_SIGMAS];
    for(int i = 0; i< NUM_SIGMAS; i++){
        for(int j = 0; j < 4; j++){
            quaternions_in_sigmas[i * 4 + j] = propagated_sigmas[i * (STATE_SIZE + 1) + j];
        }
    }
    float avg_quat[4];
    float propagated_err_vecs[NUM_SIGMAS * 3];
    grad_descent(quaternions_in_sigmas, quat_state, avg_quat, propagated_err_vecs);

    float propagated_errors[NUM_SIGMAS * (STATE_SIZE)];
    for(int i = 0; i<NUM_SIGMAS; i++){
        for(int j = 0; j < STATE_SIZE; j++){
            if(j < 3){
                propagated_errors[i * STATE_SIZE + j] = propagated_err_vecs[i * 3 + j];
            }
            else{
                propagated_errors[i * STATE_SIZE + j] = propagated_sigmas[i * (STATE_SIZE + 1) + j + 1];
            }
        }
    }

    float mean_err[STATE_SIZE];
    for(int col = 0; col < STATE_SIZE; col++){
        float var_mean = 0;
        for(int row = 0; row < NUM_SIGMAS; row++){
            var_mean += mean_weights[row] * propagated_errors[row * STATE_SIZE + col];
        }
        mean_err[col] = var_mean;
    }

    float mean_msmt[MSMT_SIZE];
    for(int col = 0; col < MSMT_SIZE; col++){
        float var_mean = 0;
        for(int row = 0; row < NUM_SIGMAS; row++){
            var_mean += mean_weights[row] * sigma_msmts[row * MSMT_SIZE + col];
        }
        mean_msmt[col] = var_mean;
    }

    float P_hat[STATE_SIZE * STATE_SIZE] = {0};
    for(int i = 0; i< NUM_SIGMAS; i++){
        float err[STATE_SIZE] = {0};
        for(int j = 0; j< STATE_SIZE; j++){
            err[j] = propagated_errors[i * STATE_SIZE + j] - mean_err[j];
        }
        float errTerr[STATE_SIZE * STATE_SIZE];
        arm_matrix_instance_f32 err_col_mat = {STATE_SIZE, 1, err};
        arm_matrix_instance_f32 err_row_mat = {1, STATE_SIZE, err};
        arm_matrix_instance_f32 errTerr_mat = {STATE_SIZE, STATE_SIZE, errTerr};
        arm_mat_mult_f32(&err_col_mat, &err_row_mat, &errTerr_mat);
        arm_scale_f32(errTerr, cov_weights[i], errTerr, STATE_SIZE * STATE_SIZE);
        arm_add_f32(P_hat, errTerr, P_hat, STATE_SIZE * STATE_SIZE);
    }

    float P_xz[STATE_SIZE * MSMT_SIZE] = {0};
    for(int i = 0; i<NUM_SIGMAS; i++){
        float errT[STATE_SIZE] = {0};
        for(int j = 0; j< STATE_SIZE; j++){
            errT[j] = propagated_errors[i * STATE_SIZE + j] - mean_err[j];
        }

        float msmt_err[MSMT_SIZE] = {0};
        for(int j = 0; j<MSMT_SIZE; j++){
            msmt_err[j] = sigma_msmts[i * MSMT_SIZE + j] - mean_msmt[j];
        }

        float errT_msmterr[STATE_SIZE * MSMT_SIZE];
        arm_matrix_instance_f32 errT_col_mat = {STATE_SIZE, 1, errT};
        arm_matrix_instance_f32 msmt_err_row_mat = {1, MSMT_SIZE, msmt_err};
        arm_matrix_instance_f32 errT_msmterr_mat = {STATE_SIZE, MSMT_SIZE, errT_msmterr};
        arm_mat_mult_f32(&errT_col_mat, &msmt_err_row_mat, &errT_msmterr_mat);
        arm_scale_f32(errT_msmterr, cov_weights[i], errT_msmterr, STATE_SIZE * MSMT_SIZE);
        arm_add_f32(P_xz, errT_msmterr, P_xz, STATE_SIZE * MSMT_SIZE);
    }

    float P_zz[MSMT_SIZE * MSMT_SIZE] = {0};
    for(int i = 0; i< NUM_SIGMAS; i++){
        float msmt_err[MSMT_SIZE] = {0};
        for(int j = 0; j<MSMT_SIZE; j++){
            msmt_err[j] = sigma_msmts[i * MSMT_SIZE + j] - mean_msmt[j];
        }
        float msmtTmsmt[MSMT_SIZE * MSMT_SIZE];
        arm_matrix_instance_f32 msmt_col_mat = {MSMT_SIZE, 1, msmt_err};
        arm_matrix_instance_f32 msmt_row_mat = {1, MSMT_SIZE, msmt_err};
        arm_matrix_instance_f32 msmtTmsmt_mat = {MSMT_SIZE, MSMT_SIZE, msmtTmsmt};
        arm_mat_mult_f32(&msmt_col_mat, &msmt_row_mat, &msmtTmsmt_mat);
        arm_scale_f32(msmtTmsmt, cov_weights[i], msmtTmsmt, MSMT_SIZE * MSMT_SIZE);
        arm_add_f32(P_zz, msmtTmsmt, P_zz, MSMT_SIZE * MSMT_SIZE);
    }

    float P_vv[MSMT_SIZE * MSMT_SIZE];
    arm_add_f32(P_zz, R, P_vv, MSMT_SIZE * MSMT_SIZE);
    float P_vv_inv[MSMT_SIZE * MSMT_SIZE];
    arm_matrix_instance_f32 P_vv_mat = {MSMT_SIZE, MSMT_SIZE, P_vv};
    arm_matrix_instance_f32 P_vv_inv_mat = {MSMT_SIZE, MSMT_SIZE, P_vv_inv};
    arm_mat_inverse_f32(&P_vv_mat, &P_vv_inv_mat);

    float k[STATE_SIZE * MSMT_SIZE];
    arm_matrix_instance_f32 P_xz_mat = {STATE_SIZE, MSMT_SIZE, P_xz};
    arm_matrix_instance_f32 k_mat = {STATE_SIZE, MSMT_SIZE, k};
    arm_mat_mult_f32(&P_xz_mat, &P_vv_inv_mat, &k_mat);

    float innovation[MSMT_SIZE];
    for(int i = 0; i<MSMT_SIZE; i++){
        innovation[i] = body[i] - mean_msmt[i];
    }
    float k_innovation[STATE_SIZE];
    arm_matrix_instance_f32 innovation_mat = {MSMT_SIZE, 1, innovation};
    arm_matrix_instance_f32 k_innovation_mat = {STATE_SIZE, 1, k_innovation};
    arm_mat_mult_f32(&k_mat, &innovation_mat, &k_innovation_mat);

    float x_hat[STATE_SIZE];
    arm_add_f32(mean_err, k_innovation, x_hat, STATE_SIZE);

    // 1. Calculate K * P_vv
    float K_Pvv[STATE_SIZE * MSMT_SIZE];
    arm_matrix_instance_f32 K_Pvv_mat = {STATE_SIZE, MSMT_SIZE, K_Pvv};
    arm_mat_mult_f32(&k_mat, &P_vv_mat, &K_Pvv_mat);

    // 2. Transpose K
    float K_T[MSMT_SIZE * STATE_SIZE];
    arm_matrix_instance_f32 K_T_mat = {MSMT_SIZE, STATE_SIZE, K_T};
    arm_mat_trans_f32(&k_mat, &K_T_mat);

    // 3. Calculate (K * P_vv) * K^T
    float K_Pvv_KT[STATE_SIZE * STATE_SIZE];
    arm_matrix_instance_f32 K_Pvv_KT_mat = {STATE_SIZE, STATE_SIZE, K_Pvv_KT};
    arm_mat_mult_f32(&K_Pvv_mat, &K_T_mat, &K_Pvv_KT_mat);

    // 4. Update Covariance: P = P_hat - K_Pvv_KT
    arm_scale_f32(K_Pvv_KT, -1.0f, K_Pvv_KT, STATE_SIZE * STATE_SIZE); // Make it negative
    arm_add_f32(P_hat, K_Pvv_KT, P, STATE_SIZE * STATE_SIZE);

    ensure_psd(&P_mat);

    float x_hat_rot[4];
    float x_hat_vec[3] = {x_hat[0], x_hat[1], x_hat[2]};
    rotationvec2quat(x_hat_vec, x_hat_rot);
    float new_guess[4];
    quat_multiply(avg_quat, x_hat_rot, new_guess);
    quat_norm(new_guess, new_guess);
    
    memcpy(new_err_state, x_hat, sizeof(float) * STATE_SIZE);
    memcpy(new_quat_state, new_guess, sizeof(float) * 4);
    memcpy(new_P, P, sizeof(float) * STATE_SIZE * STATE_SIZE);
}

