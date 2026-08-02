#include "include/iterate.h"
#include "include/laextension.h"
#include "declareFunctions.h"
#include "math.h"
#include "include/quat.h"

// Gauss-Jordan in-place inverse for n×n double matrix.
// Replaces inv() which routes through LAPACK's single-precision linsolve.
static void double_inv(double* A, int n) {
    double aug[n * 2 * n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) aug[i * 2*n + j]   = A[i * n + j];
        for (int j = 0; j < n; j++) aug[i * 2*n + n+j] = (i == j) ? 1.0 : 0.0;
    }
    for (int col = 0; col < n; col++) {
        // partial pivot
        int pivot = col;
        for (int r = col+1; r < n; r++)
            if (fabs(aug[r * 2*n + col]) > fabs(aug[pivot * 2*n + col])) pivot = r;
        if (pivot != col)
            for (int j = 0; j < 2*n; j++) {
                double tmp = aug[col * 2*n + j];
                aug[col * 2*n + j] = aug[pivot * 2*n + j];
                aug[pivot * 2*n + j] = tmp;
            }
        double diag = aug[col * 2*n + col];
        if (fabs(diag) < 1e-15) continue;
        for (int j = 0; j < 2*n; j++) aug[col * 2*n + j] /= diag;
        for (int r = 0; r < n; r++) {
            if (r == col) continue;
            double factor = aug[r * 2*n + col];
            for (int j = 0; j < 2*n; j++) aug[r * 2*n + j] -= factor * aug[col * 2*n + j];
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) A[i * n + j] = aug[i * 2*n + n+j];
}

void ensure_psd(double* mat, int size) {
    double eps = 1e-9;

    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            double avg = (mat[i * size + j] + mat[j * size + i]) / 2.0;
            mat[i * size + j] = avg;
            mat[j * size + i] = avg;
        }
    }

    double L[6 * 6];
    memset(L, 0, sizeof(L));
    double min_residual = 1e300;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < i; j++) {
            double s = 0.0;
            for (int k = 0; k < j; k++) s += L[i * size + k] * L[j * size + k];
            double Ljj = L[j * size + j];
            L[i * size + j] = (Ljj > 0.0) ? (mat[i * size + j] - s) / Ljj : 0.0;
        }
        double s = 0.0;
        for (int k = 0; k < i; k++) s += L[i * size + k] * L[i * size + k];
        double residual = mat[i * size + i] - s;
        if (residual < min_residual) min_residual = residual;
        L[i * size + i] = (residual > 0.0) ? sqrt(residual) : 0.0;
    }

    if (min_residual < eps) {
        double add = eps - min_residual + eps;
        for (int i = 0; i < size; i++) mat[i * size + i] += add;
    }
}

double calculate_lambda(int n, double alpha){
    int k = 9 - n;
    return alpha * alpha * (n + k) - n;
}

void get_sigma_points(double lam, double* state, double* P, double* sigmas){
    double P1[STATE_SIZE * STATE_SIZE];
    memcpy(P1, P, sizeof(double) * STATE_SIZE * STATE_SIZE);
    scale(P1, lam + STATE_SIZE, STATE_SIZE, STATE_SIZE);
    ensure_psd(P1, STATE_SIZE);
    double U[STATE_SIZE * STATE_SIZE];
    chol(P1, U, STATE_SIZE);

    for(int i = 0; i<STATE_SIZE; i++){
        sigmas[i] = state[i];
    }
    for(int k = 0; k<STATE_SIZE; k++){
        int start = k * STATE_SIZE;
        int sigma_offset = STATE_SIZE;
        for(int i = 0; i<STATE_SIZE; i++){
            sigmas[sigma_offset + start + i] = state[i] + U[i * STATE_SIZE + k];
        }
    }
    for(int k = 0; k<STATE_SIZE; k++){
        int start = k * STATE_SIZE;
        int sigma_offset = (STATE_SIZE + 1) * STATE_SIZE;
        for(int i = 0; i<STATE_SIZE; i++){
            sigmas[sigma_offset + start + i] = state[i] - U[i * STATE_SIZE + k];
        }
    }
}

void get_weights(double lambda, int n, double alpha, double beta, double* cov_weights, double* mean_weights){
    double c = 0.5 / (n + lambda);
    double c_weights[NUM_SIGMAS];
    double m_weights[NUM_SIGMAS];
    for(int i = 0; i< NUM_SIGMAS; i++){
        c_weights[i] = c;
        m_weights[i] = c;
    }
    c_weights[0] = lambda / (n + lambda) + (1 - alpha*alpha + beta);
    m_weights[0] = lambda / (n + lambda);
    memcpy(cov_weights, c_weights, sizeof(double) * NUM_SIGMAS);
    memcpy(mean_weights, m_weights, sizeof(double) * NUM_SIGMAS);
}

void error_sigmas_to_quat_sigmas(double* error_sigmas, double* current_guess, double* quat_sigmas){
    int n = NUM_SIGMAS;
    double q[n * (STATE_SIZE + 1)];
    for(int i = 0; i<n; i++){
        double error_rot[3] = {error_sigmas[i * STATE_SIZE],error_sigmas[i * STATE_SIZE + 1],error_sigmas[i * STATE_SIZE + 2]};
        double error_quat[4];
        rotationvec2quat(error_rot, error_quat);
        double new_rotation[4];
        quat_multiply(error_quat, current_guess, new_rotation);
        double norm_rot[4];
        quat_norm(new_rotation, norm_rot);
        for(int j = 0; j<4; j++){
            q[i * (STATE_SIZE + 1) + j] = norm_rot[j];
        }
        for(int j = 4; j < STATE_SIZE + 1; j++){
            q[i * (STATE_SIZE + 1) + j] = error_sigmas[i * STATE_SIZE + (j - 1)];
        }
    }
    memcpy(quat_sigmas, q, sizeof(double) * n * (STATE_SIZE + 1));
}

void propagate_sigmas(double* quat_sigmas, double* gyro, double dt, double* propagated_sigmas){
    double p_sigmas[NUM_SIGMAS * (STATE_SIZE + 1)];
    for(int i = 0; i<NUM_SIGMAS; i++){
        int offset = (STATE_SIZE + 1) * i;
        double rotation[4];
        double bias[3];
        for(int j = 0; j<4; j++){
            rotation[j] = quat_sigmas[offset + j];
        }
        for(int j = 4; j<7; j++){
            bias[j - 4] = quat_sigmas[offset + j];
        }
        double omega[3];
        scale(bias, -1, 1, 3);
        add(gyro, bias, omega, 1, 3, 3);
        scale(omega, dt, 1, 3);
        double delta_q[4];
        rotationvec2quat(omega, delta_q);

        double propagated_rot[4];
        quat_multiply(rotation, delta_q, propagated_rot);
        double norm_rot[4];
        quat_norm(propagated_rot, norm_rot);

        int p_offset = i * (STATE_SIZE + 1);
        for(int j = 0; j<4; j++){
            p_sigmas[j + p_offset] = norm_rot[j];
        }
        for(int j = 4; j<7; j++){
            p_sigmas[j + p_offset] = quat_sigmas[offset + j];
        }
    }
    memcpy(propagated_sigmas, p_sigmas, sizeof(double) * NUM_SIGMAS * (STATE_SIZE + 1));
}

void get_sigma_measurements(double* q_sigmas, double* ref, int num_vecs, double* msmts){
    int msmt_size = 3 * num_vecs;
    double m[NUM_SIGMAS * MAX_MSMT_SIZE];
    for(int i = 0; i < NUM_SIGMAS; i++){
        double sigma_rot[4];
        for(int j = 0; j < 4; j++){
            sigma_rot[j] = q_sigmas[i * (STATE_SIZE + 1) + j];
        }
        double sigma_ref_to_body[4];
        quat_inv(sigma_rot, sigma_ref_to_body);
        for(int v = 0; v < num_vecs; v++){
            double ref_v[3] = {ref[v*3], ref[v*3+1], ref[v*3+2]};
            double reading_v[3];
            quat_apply(sigma_ref_to_body, ref_v, reading_v);
            for(int j = 0; j < 3; j++){
                m[i * msmt_size + v*3 + j] = reading_v[j];
            }
        }
    }
    memcpy(msmts, m, sizeof(double) * NUM_SIGMAS * msmt_size);
}

void get_error_vectors(int num_other_quats, double* other_quats, double* quat, double* error_vecs){
    double e[num_other_quats * 3];
    double x_inv[4];
    quat_inv(quat, x_inv);
    for(int i = 0; i<num_other_quats; i++){
        double other_quat[4];
        for(int j = 0; j<4; j++){
            other_quat[j] = other_quats[i * 4 + j];
        }
        double err_quat[4];
        quat_multiply(other_quat, x_inv, err_quat);
        double err_vec[3];
        quat2rotationvec(err_quat, err_vec);
        for(int j = 0; j<3; j++){
            e[i * 3 + j] = err_vec[j];
        }
    }
    memcpy(error_vecs, e, sizeof(double) * num_other_quats * 3);
}

void grad_descent(double* propagated_quats, double* current_quat, double* mean_weights, double* avg_quat, double* errors){
    double safe_w[NUM_SIGMAS];
    double wsum = 0.0;
    for(int i = 0; i < NUM_SIGMAS; i++){
        safe_w[i] = (mean_weights[i] > 0.0) ? mean_weights[i] : 0.0;
        wsum += safe_w[i];
    }
    if(wsum > 0.0){
        for(int i = 0; i < NUM_SIGMAS; i++) safe_w[i] /= wsum;
    } else {
        for(int i = 0; i < NUM_SIGMAS; i++) safe_w[i] = 1.0 / NUM_SIGMAS;
    }

    double moving_avg[4];
    memcpy(moving_avg, current_quat, sizeof(double) * 4);

    for(int iter = 0; iter < 100; iter++){
        double error_vecs[3 * NUM_SIGMAS];
        get_error_vectors(NUM_SIGMAS, propagated_quats, moving_avg, error_vecs);
        memcpy(errors, error_vecs, sizeof(double) * 3 * NUM_SIGMAS);
        double avg_err[3] = {0,0,0};
        for(int i = 0; i<NUM_SIGMAS; i++){
            for(int j = 0; j < 3; j++){
                avg_err[j] += safe_w[i] * error_vecs[i * 3 + j];
            }
        }
        if(sqrt(avg_err[0]*avg_err[0] + avg_err[1]*avg_err[1] + avg_err[2]*avg_err[2]) < 1e-6){
            break;
        }
        double avg_error_quat[4];
        rotationvec2quat(avg_err, avg_error_quat);
        double new_moving_avg[4];
        quat_multiply(avg_error_quat, moving_avg, new_moving_avg);
        memcpy(moving_avg, new_moving_avg, sizeof(double) * 4);
        quat_norm(moving_avg, moving_avg);
    }
    memcpy(avg_quat, moving_avg, sizeof(double) * 4);
}


void iterate(double* error_state, double* quat_state, double* cov,
             double* body, double* ref, int num_vecs,
             double* gyro, double* Q, double* R, double dt,
             double* new_err_state, double* new_quat_state, double* new_P){
    int msmt_size = 3 * num_vecs;

    double P[STATE_SIZE * STATE_SIZE];
    memcpy(P, cov, sizeof(double) * STATE_SIZE * STATE_SIZE);

    ensure_psd(P, STATE_SIZE);

    double lambda = calculate_lambda(STATE_SIZE, ALPHA);
    double error_sigmas[STATE_SIZE * NUM_SIGMAS];
    get_sigma_points(lambda, error_state, P, error_sigmas);
    double quat_sigmas[(STATE_SIZE + 1) * NUM_SIGMAS];
    error_sigmas_to_quat_sigmas(error_sigmas, quat_state, quat_sigmas);
    double propagated_sigmas[(STATE_SIZE + 1) * NUM_SIGMAS];
    propagate_sigmas(quat_sigmas, gyro, dt, propagated_sigmas);

    double mean_weights[NUM_SIGMAS];
    double cov_weights[NUM_SIGMAS];
    get_weights(lambda, STATE_SIZE, ALPHA, BETA, cov_weights, mean_weights);

    double quaternions_in_sigmas[4 * NUM_SIGMAS];
    for(int i = 0; i< NUM_SIGMAS; i++){
        for(int j = 0; j < 4; j++){
            quaternions_in_sigmas[i * 4 + j] = propagated_sigmas[i * (STATE_SIZE + 1) + j];
        }
    }
    double avg_quat[4];
    double propagated_err_vecs[NUM_SIGMAS * 3];
    grad_descent(quaternions_in_sigmas, quat_state, mean_weights, avg_quat, propagated_err_vecs);

    double propagated_errors[NUM_SIGMAS * STATE_SIZE];
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

    double mean_err[STATE_SIZE];
    for(int col = 0; col < STATE_SIZE; col++){
        double var_mean = 0;
        for(int row = 0; row < NUM_SIGMAS; row++){
            var_mean += mean_weights[row] * propagated_errors[row * STATE_SIZE + col];
        }
        mean_err[col] = var_mean;
    }

    double P_hat[STATE_SIZE * STATE_SIZE] = {0};
    for(int i = 0; i< NUM_SIGMAS; i++){
        double err[STATE_SIZE] = {0};
        for(int j = 0; j< STATE_SIZE; j++){
            err[j] = propagated_errors[i * STATE_SIZE + j] - mean_err[j];
        }
        double errTerr[STATE_SIZE * STATE_SIZE];
        mul(err, err, false, errTerr, STATE_SIZE, 1, STATE_SIZE);
        scale(errTerr, cov_weights[i], STATE_SIZE, STATE_SIZE);
        double new_P_hat[STATE_SIZE * STATE_SIZE];
        add(P_hat, errTerr, new_P_hat, STATE_SIZE, STATE_SIZE, STATE_SIZE);
        memcpy(P_hat, new_P_hat, sizeof(double) * STATE_SIZE * STATE_SIZE);
    }
    double P_hat_Q[STATE_SIZE * STATE_SIZE];
    add(P_hat, Q, P_hat_Q, STATE_SIZE, STATE_SIZE, STATE_SIZE);
    memcpy(P_hat, P_hat_Q, sizeof(double) * STATE_SIZE * STATE_SIZE);

    if(num_vecs == 0 || body == NULL || ref == NULL){
        ensure_psd(P_hat, STATE_SIZE);
        memcpy(new_err_state, mean_err, sizeof(double) * STATE_SIZE);
        memcpy(new_quat_state, avg_quat, sizeof(double) * 4);
        memcpy(new_P, P_hat, sizeof(double) * STATE_SIZE * STATE_SIZE);
        return;
    }

    double sigma_msmts[MAX_MSMT_SIZE * NUM_SIGMAS];
    get_sigma_measurements(propagated_sigmas, ref, num_vecs, sigma_msmts);

    double mean_msmt[MAX_MSMT_SIZE];
    for(int col = 0; col < msmt_size; col++){
        double var_mean = 0;
        for(int row = 0; row < NUM_SIGMAS; row++){
            var_mean += mean_weights[row] * sigma_msmts[row * msmt_size + col];
        }
        mean_msmt[col] = var_mean;
    }

    double P_xz[STATE_SIZE * MAX_MSMT_SIZE] = {0};
    for(int i = 0; i<NUM_SIGMAS; i++){
        double errT[STATE_SIZE] = {0};
        for(int j = 0; j< STATE_SIZE; j++){
            errT[j] = propagated_errors[i * STATE_SIZE + j] - mean_err[j];
        }
        double msmt_err[MAX_MSMT_SIZE] = {0};
        for(int j = 0; j<msmt_size; j++){
            msmt_err[j] = sigma_msmts[i * msmt_size + j] - mean_msmt[j];
        }
        double errT_msmterr[STATE_SIZE * MAX_MSMT_SIZE];
        mul(errT, msmt_err, false, errT_msmterr, STATE_SIZE, 1, msmt_size);
        scale(errT_msmterr, cov_weights[i], STATE_SIZE, msmt_size);
        double new_P_xz[STATE_SIZE * MAX_MSMT_SIZE];
        add(P_xz, errT_msmterr, new_P_xz, STATE_SIZE, msmt_size, msmt_size);
        memcpy(P_xz, new_P_xz, sizeof(double) * STATE_SIZE * msmt_size);
    }

    double P_zz[MAX_MSMT_SIZE * MAX_MSMT_SIZE] = {0};
    for(int i = 0; i< NUM_SIGMAS; i++){
        double msmt_err[MAX_MSMT_SIZE] = {0};
        for(int j = 0; j<msmt_size; j++){
            msmt_err[j] = sigma_msmts[i * msmt_size + j] - mean_msmt[j];
        }
        double msmtTmsmt[MAX_MSMT_SIZE * MAX_MSMT_SIZE];
        mul(msmt_err, msmt_err, false, msmtTmsmt, msmt_size, 1, msmt_size);
        scale(msmtTmsmt, cov_weights[i], msmt_size, msmt_size);
        double new_P_zz[MAX_MSMT_SIZE * MAX_MSMT_SIZE];
        add(P_zz, msmtTmsmt, new_P_zz, msmt_size, msmt_size, msmt_size);
        memcpy(P_zz, new_P_zz, sizeof(double) * msmt_size * msmt_size);
    }

    double P_vv[MAX_MSMT_SIZE * MAX_MSMT_SIZE];
    add(P_zz, R, P_vv, msmt_size, msmt_size, msmt_size);
    double P_vv_inv[MAX_MSMT_SIZE * MAX_MSMT_SIZE];
    memcpy(P_vv_inv, P_vv, sizeof(double) * msmt_size * msmt_size);
    double_inv(P_vv_inv, msmt_size);

    double k[STATE_SIZE * MAX_MSMT_SIZE];
    mul(P_xz, P_vv_inv, false, k, STATE_SIZE, msmt_size, msmt_size);

    double innovation[MAX_MSMT_SIZE];
    for(int i = 0; i<msmt_size; i++){
        innovation[i] = body[i] - mean_msmt[i];
    }
    double k_innovation[STATE_SIZE];
    mul(k, innovation, false, k_innovation, STATE_SIZE, msmt_size, 1);

    double x_hat[STATE_SIZE];
    add(mean_err, k_innovation, x_hat, STATE_SIZE, 1, 1);

    double K_Pvv[STATE_SIZE * MAX_MSMT_SIZE];
    mul(k, P_vv, false, K_Pvv, STATE_SIZE, msmt_size, msmt_size);

    double K_T[MAX_MSMT_SIZE * STATE_SIZE];
    memcpy(K_T, k, sizeof(double) * STATE_SIZE * msmt_size);
    tranf(K_T, STATE_SIZE, msmt_size);

    double K_Pvv_KT[STATE_SIZE * STATE_SIZE];
    mul(K_Pvv, K_T, false, K_Pvv_KT, STATE_SIZE, msmt_size, STATE_SIZE);

    scale(K_Pvv_KT, -1.0, STATE_SIZE, STATE_SIZE);
    double P_new[STATE_SIZE * STATE_SIZE];
    add(P_hat, K_Pvv_KT, P_new, STATE_SIZE, STATE_SIZE, STATE_SIZE);
    ensure_psd(P_new, STATE_SIZE);

    double x_hat_rot[4];
    double x_hat_vec[3] = {x_hat[0], x_hat[1], x_hat[2]};
    rotationvec2quat(x_hat_vec, x_hat_rot);
    double new_guess[4];
    quat_multiply(x_hat_rot, avg_quat, new_guess);
    quat_norm(new_guess, new_guess);

    memcpy(new_err_state, x_hat, sizeof(double) * STATE_SIZE);
    memcpy(new_quat_state, new_guess, sizeof(double) * 4);
    memcpy(new_P, P_new, sizeof(double) * STATE_SIZE * STATE_SIZE);
}
