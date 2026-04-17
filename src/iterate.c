#include "include/iterate.h"
#include "include/laextension.h"
#include "declareFunctions.h"
#include "math.h"
#include "include/quat.h"

void ensure_psd(float* mat, int size) {
    float eps = 1e-6; // Adjust this if necessary (MATLAB used 1e-7)

    // 1. Force perfect symmetry: P = (P + P^T) / 2
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            // Average the mirrored elements
            float avg = (mat[i * size + j] + mat[j * size + i]) / 2.0f;
            
            // Assign the average back to both spots
            mat[i * size + j] = avg;
            mat[j * size + i] = avg;
        }
    }

    // 2. Add a tiny epsilon to the diagonal to ensure eigenvalues > 0
    // This physically guarantees Positive Definiteness as long as the 
    // variance hasn't completely blown up to negative infinity.
    for (int i = 0; i < size; i++) {
        mat[i * size + i] += eps;
        
        // Safety catch: if diagonal somehow became negative, hard-reset it
        if (mat[i * size + i] < 0.0f) {
            mat[i * size + i] = eps;
        }
    }
}
float calculate_lambda(int n, float alpha){
    int k = 3 - n;
    return alpha * alpha * (n + k) - n;
}
void get_sigma_points(float lam, float* state, float* P, float* sigmas){
    float P1[STATE_SIZE * STATE_SIZE];
    memcpy(P1, P, sizeof(float) * STATE_SIZE * STATE_SIZE);
    scale(P1, lam + STATE_SIZE, STATE_SIZE, STATE_SIZE);
    ensure_psd(P1, STATE_SIZE);
    float U[STATE_SIZE * STATE_SIZE];
    chol(P1, U, STATE_SIZE);
    tran(U, STATE_SIZE, STATE_SIZE);

    for(int i = 0; i<STATE_SIZE; i++){
        sigmas[i] = state[i];
    }
    for(int k = 0; k<STATE_SIZE; k++){
        int start = k * STATE_SIZE;
        int sigma_offset = STATE_SIZE;
        for(int i = 0; i<STATE_SIZE; i++){
            sigmas[sigma_offset + start + i] = state[i] + U[start + i];
        }
    }
    for(int k = 0; k<STATE_SIZE; k++){
        int start = k * STATE_SIZE;
        int sigma_offset = (STATE_SIZE + 1) * STATE_SIZE;
        for(int i = 0; i<STATE_SIZE; i++){
            sigmas[sigma_offset + start + i] = state[i] - U[start + i];
        }
    }
}
void get_weights(float lambda, int n, float alpha, float beta, float* cov_weights, float* mean_weights){
    float c = 0.5f / (n + lambda);
    float c_weights[NUM_SIGMAS];
    float m_weights[NUM_SIGMAS];
    for(int i = 0; i< NUM_SIGMAS; i++){
        c_weights[i] = c;
        m_weights[i] = c;
    }
    c_weights[0] = lambda / (n + lambda) + (1 - alpha*alpha + beta);
    m_weights[0] = lambda / (n + lambda);
    memcpy(cov_weights, c_weights, sizeof(float) * NUM_SIGMAS);
    memcpy(mean_weights, m_weights, sizeof(float) * NUM_SIGMAS);
}
void error_sigmas_to_quat_sigmas(float* error_sigmas, float* current_guess, float* quat_sigmas){
    int n = NUM_SIGMAS;
    float q[n * (STATE_SIZE + 1)];
    for(int i = 0; i<n; i++){
        float error_rot[3] = {error_sigmas[i * STATE_SIZE],error_sigmas[i * STATE_SIZE + 1],error_sigmas[i * STATE_SIZE + 2]};
        float error_quat[4];
        rotationvec2quat(error_rot, error_quat);
        float new_rotation[4];
        quat_multiply(current_guess, error_quat, new_rotation);
        float norm_rot[4];
        quat_norm(new_rotation, norm_rot);
        for(int j = 0; j<4; j++){
            q[i * (STATE_SIZE + 1) + j] = norm_rot[j];
        }
        for(int j = 4; j < STATE_SIZE + 1; j++){
            q[i * (STATE_SIZE + 1) + j] = error_sigmas[i * STATE_SIZE + (j - 1)];
        }
    }
    memcpy(quat_sigmas, q, sizeof(float) * n * (STATE_SIZE + 1));
}
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
        scale(bias, -1, 1, 3);
        add(gyro, bias, omega, 1, 3, 3);
        scale(omega, dt, 1, 3);
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
    
    ensure_psd(P, STATE_SIZE);

    float lambda = calculate_lambda(STATE_SIZE, ALPHA);
    float P_Q[STATE_SIZE * STATE_SIZE];
    add(P, Q, P_Q, STATE_SIZE, STATE_SIZE, STATE_SIZE);
    float error_sigmas[STATE_SIZE * NUM_SIGMAS];
    get_sigma_points(lambda, error_state, P_Q, error_sigmas);
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
        float errT[STATE_SIZE];
        memcpy(errT, err, sizeof(float) * STATE_SIZE);
        float errTerr[STATE_SIZE * STATE_SIZE];
        mul(errT, err, false, errTerr, STATE_SIZE, 1, STATE_SIZE);
        scale(errTerr, cov_weights[i], STATE_SIZE, STATE_SIZE);
        float new_P_hat[STATE_SIZE * STATE_SIZE];
        add(P_hat, errTerr, new_P_hat, STATE_SIZE, STATE_SIZE, STATE_SIZE);
        memcpy(P_hat, new_P_hat, sizeof(float) * STATE_SIZE * STATE_SIZE);
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
        mul(errT, msmt_err, false, errT_msmterr, STATE_SIZE, 1, MSMT_SIZE);
        scale(errT_msmterr, cov_weights[i], STATE_SIZE, MSMT_SIZE);
        float new_P_xz[STATE_SIZE * MSMT_SIZE];
        add(P_xz, errT_msmterr, new_P_xz, STATE_SIZE, MSMT_SIZE, MSMT_SIZE);
        memcpy(P_xz, new_P_xz, sizeof(float) * STATE_SIZE * MSMT_SIZE);
    }

    float P_zz[MSMT_SIZE * MSMT_SIZE] = {0};
    for(int i = 0; i< NUM_SIGMAS; i++){
        float msmt_err[MSMT_SIZE] = {0};
        for(int j = 0; j<MSMT_SIZE; j++){
            msmt_err[j] = sigma_msmts[i * MSMT_SIZE + j] - mean_msmt[j];
        }
        float msmtT[MSMT_SIZE];
        memcpy(msmtT, msmt_err, sizeof(float) * MSMT_SIZE);
        float msmtTmsmt[MSMT_SIZE * MSMT_SIZE];
        mul(msmtT, msmt_err, false, msmtTmsmt, MSMT_SIZE, 1, MSMT_SIZE);
        scale(msmtTmsmt, cov_weights[i], MSMT_SIZE, MSMT_SIZE);
        float new_P_zz[MSMT_SIZE * MSMT_SIZE];
        add(P_zz, msmtTmsmt, new_P_zz, MSMT_SIZE, MSMT_SIZE, MSMT_SIZE);
        memcpy(P_zz, new_P_zz, sizeof(float) * MSMT_SIZE * MSMT_SIZE);
    }
    
    float P_vv[MSMT_SIZE * MSMT_SIZE];
    add(P_zz, R, P_vv, MSMT_SIZE, MSMT_SIZE, MSMT_SIZE);
    float P_vv_inv[MSMT_SIZE * MSMT_SIZE];
    memcpy(P_vv_inv, P_vv, sizeof(float) * MSMT_SIZE * MSMT_SIZE);
    pinv(P_vv_inv, MSMT_SIZE, MSMT_SIZE);

    float k[STATE_SIZE * MSMT_SIZE];
    mul(P_xz, P_vv_inv, false, k, STATE_SIZE, MSMT_SIZE, MSMT_SIZE);

    float innovation[MSMT_SIZE];
    for(int i = 0; i<MSMT_SIZE; i++){
        innovation[i] = body[i] - mean_msmt[i];
    }
    float k_innovation[STATE_SIZE];
    mul(k, innovation, false, k_innovation, STATE_SIZE, MSMT_SIZE, 1);

    float x_hat[STATE_SIZE];
    add(mean_err, k_innovation, x_hat, STATE_SIZE, 1, 1);

    // 1. Calculate K * P_vv
    float K_Pvv[STATE_SIZE * MSMT_SIZE];
    mul(k, P_vv, false, K_Pvv, STATE_SIZE, MSMT_SIZE, MSMT_SIZE);

    // 2. Transpose K
    float K_T[MSMT_SIZE * STATE_SIZE];
    memcpy(K_T, k, sizeof(float) * MSMT_SIZE * STATE_SIZE);
    tranf(K_T, STATE_SIZE, MSMT_SIZE);

    // 3. Calculate (K * P_vv) * K^T
    float K_Pvv_KT[STATE_SIZE * STATE_SIZE];
    mul(K_Pvv, K_T, false, K_Pvv_KT, STATE_SIZE, MSMT_SIZE, STATE_SIZE);

    // 4. Update Covariance: P = P_hat - K_Pvv_KT
    scale(K_Pvv_KT, -1.0f, STATE_SIZE, STATE_SIZE); // Make it negative
    add(P_hat, K_Pvv_KT, P, STATE_SIZE, STATE_SIZE, STATE_SIZE);
    
    ensure_psd(P, STATE_SIZE);

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

