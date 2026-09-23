#include "include/iterate.h"
#include "Include/arm_math_types.h"
#include "Include/dsp/matrix_functions.h"
#include "include/laextension.h"
#include "include/quat.h"
#include "arm_math.h"
#include "math.h"
#include <string.h>

#define QUAT_SIGMA_SIZE (STATE_SIZE + 1) // [q (4), bias (3)]
#define GRAD_DESCENT_MAX_ITERS 20
#define GRAD_DESCENT_TOL 1e-6f

/**
 * \fn cholesky_lower
 *
 * \brief L * L^T = A for a symmetric positive definite n x n A. Unlike arm_mat_cholesky_f32
 * alone, this also zeroes the strict upper triangle of L (CMSIS leaves it untouched).
 */
static bool cholesky_lower(const float32_t* A, float32_t* L, int n) {
    memset(L, 0, sizeof(float32_t) * n * n);
    arm_matrix_instance_f32 A_mat = {(uint16_t)n, (uint16_t)n, (float32_t*)A};
    arm_matrix_instance_f32 L_mat = {(uint16_t)n, (uint16_t)n, L};
    if (arm_mat_cholesky_f32(&A_mat, &L_mat) != ARM_MATH_SUCCESS) {
        return false;
    }
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            L[i * n + j] = 0.0f;
        }
    }
    return all_finite(L, n * n);
}

/**
 * \fn ensure_psd
 *
 * \brief Makes a square matrix (n <= STATE_SIZE) symmetric positive definite in place.
 *
 * Port of ensure_positive_definite() from simulink.m. MATLAB uses eig(); CMSIS has no symmetric
 * eigen-solver, so instead we test with a Cholesky factorization and add growing diagonal jitter
 * until it succeeds. Unlike the old version, a healthy matrix is left unchanged (no epsilon is
 * added every call).
 *
 * \return false if the matrix contains NaN/inf or could not be repaired
 */
bool ensure_psd(float32_t* P, int n) {
    if (n <= 0 || n > STATE_SIZE || !all_finite(P, n * n)) {
        return false;
    }

    // 1. Force perfect symmetry: P = (P + P^T) / 2
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            float avg = 0.5f * (P[i * n + j] + P[j * n + i]);
            P[i * n + j] = avg;
            P[j * n + i] = avg;
        }
    }

    float32_t L[STATE_SIZE * STATE_SIZE];
    if (cholesky_lower(P, L, n)) {
        return true;
    }

    // 2. Add jitter relative to the matrix scale until the factorization succeeds
    float scale = 0.0f;
    for (int i = 0; i < n; i++) {
        scale = fmaxf(scale, fabsf(P[i * n + i]));
    }
    if (!(scale > 0.0f)) {
        scale = 1.0f;
    }
    float jitter = 1e-7f * scale;
    for (int attempt = 0; attempt < 12; attempt++) {
        for (int i = 0; i < n; i++) {
            P[i * n + i] += jitter;
        }
        if (cholesky_lower(P, L, n)) {
            return true;
        }
        jitter *= 10.0f;
    }
    return false;
}

/**
 * \fn calculate_lambda
 *
 * \brief lambda = alpha^2 (n + kappa) - n
 */
float32_t calculate_lambda(int n, float32_t alpha, float32_t kappa) {
    return alpha * alpha * ((float32_t)n + kappa) - (float32_t)n;
}

/**
 * \fn get_sigma_points
 *
 * \brief Generates the 2n+1 sigma points x, x + col_k(L), x - col_k(L) with
 * L L^T = (n + lambda) P.
 *
 * The previous version used rows of the (lower-triangular) CMSIS Cholesky factor together
 * with its uninitialized upper triangle. Rows of L only reproduce P when P is diagonal, so
 * every step after the first produced the wrong spread. Columns of L (= rows of MATLAB's
 * upper-triangular chol(P)) are the correct offsets.
 */
bool get_sigma_points(float32_t lam, const float32_t* state, const float32_t* P, float32_t* sigmas) {
    float32_t P1[STATE_SIZE * STATE_SIZE];
    for (int i = 0; i < STATE_SIZE * STATE_SIZE; i++) {
        P1[i] = P[i] * (lam + STATE_SIZE);
    }
    if (!ensure_psd(P1, STATE_SIZE)) {
        return false;
    }

    float32_t L[STATE_SIZE * STATE_SIZE];
    if (!cholesky_lower(P1, L, STATE_SIZE)) {
        return false;
    }

    memcpy(sigmas, state, sizeof(float32_t) * STATE_SIZE);
    for (int k = 0; k < STATE_SIZE; k++) {
        float32_t* plus = &sigmas[(1 + k) * STATE_SIZE];
        float32_t* minus = &sigmas[(1 + STATE_SIZE + k) * STATE_SIZE];
        for (int i = 0; i < STATE_SIZE; i++) {
            float32_t offset = L[i * STATE_SIZE + k]; // column k
            plus[i] = state[i] + offset;
            minus[i] = state[i] - offset;
        }
    }
    return true;
}

/**
 * \fn get_weights
 *
 * \brief Calculates the mean and covariance weights for the Unscented Kalman Filter.
 */
void get_weights(float32_t lambda, int n, float32_t alpha, float32_t beta, float32_t* cov_weights,
                 float32_t* mean_weights) {
    float32_t c = 0.5f / ((float32_t)n + lambda);

    for (int i = 0; i < 2 * n + 1; i++) {
        cov_weights[i] = c;
        mean_weights[i] = c;
    }

    mean_weights[0] = lambda / ((float32_t)n + lambda);
    cov_weights[0] = mean_weights[0] + (1 - alpha * alpha + beta);
}

/**
 * \fn error_sigmas_to_quat_sigmas
 *
 * \brief Applies each sigma's (local) attitude error to the current estimate:
 * q_i = q_est * dq(err_i). Output rows are [q (4), bias (3)].
 */
static void error_sigmas_to_quat_sigmas(const float32_t* error_sigmas, const float32_t* current_guess,
                                        float32_t* quat_sigmas) {
    for (int i = 0; i < NUM_SIGMAS; i++) {
        const float32_t* err = &error_sigmas[i * STATE_SIZE];
        float32_t* out = &quat_sigmas[i * QUAT_SIGMA_SIZE];
        float32_t error_quat[4];
        rotationvec2quat(err, error_quat);
        quat_multiply(current_guess, error_quat, out);
        quat_norm(out, out);
        out[4] = err[3];
        out[5] = err[4];
        out[6] = err[5];
    }
}

/**
 * \fn propagate_sigmas
 *
 * \brief Propagates each quaternion sigma by its own bias-corrected gyro rate:
 * q <- q * dq((gyro - bias) dt).
 */
static void propagate_sigmas(float32_t* quat_sigmas, const float32_t* gyro, float dt) {
    for (int i = 0; i < NUM_SIGMAS; i++) {
        float32_t* s = &quat_sigmas[i * QUAT_SIGMA_SIZE];
        float32_t omega_dt[3] = {(gyro[0] - s[4]) * dt, (gyro[1] - s[5]) * dt, (gyro[2] - s[6]) * dt};
        float32_t delta_q[4];
        rotationvec2quat(omega_dt, delta_q);
        quat_multiply(s, delta_q, s);
        quat_norm(s, s);
    }
}

/**
 * \fn get_sigma_measurements
 *
 * \brief Predicted body-frame measurements: ref vectors rotated by each sigma's ref->body.
 */
static void get_sigma_measurements(const float32_t* quat_sigmas, const float32_t* ref, int num_vecs,
                                   float32_t* msmts) {
    int m = 3 * num_vecs;
    for (int i = 0; i < NUM_SIGMAS; i++) {
        float32_t ref_to_body[4];
        quat_conj(&quat_sigmas[i * QUAT_SIGMA_SIZE], ref_to_body);
        for (int v = 0; v < num_vecs; v++) {
            quat_apply(ref_to_body, &ref[3 * v], &msmts[i * m + 3 * v]);
        }
    }
}

/**
 * \fn get_error_vectors
 *
 * \brief Local error rotation vectors of each sigma quaternion about x: rv(x^-1 * q_i).
 */
static void get_error_vectors(const float32_t* quat_sigmas, const float32_t* x, float32_t* error_vecs) {
    float32_t x_inv[4];
    quat_conj(x, x_inv);
    for (int i = 0; i < NUM_SIGMAS; i++) {
        float32_t err_quat[4];
        quat_multiply(x_inv, &quat_sigmas[i * QUAT_SIGMA_SIZE], err_quat);
        quat2rotationvec(err_quat, &error_vecs[i * 3]);
    }
}

/**
 * \fn grad_descent
 *
 * \brief Weighted quaternion mean of the sigma points (as in simulink.m's gradient_descent,
 * but in the local-error convention). error_vecs always correspond to the returned average.
 */
static void grad_descent(const float32_t* quat_sigmas, const float32_t* initial, const float32_t* weights,
                         float32_t* avg_quat, float32_t* error_vecs) {
    float32_t w_sum = 0.0f;
    for (int i = 0; i < NUM_SIGMAS; i++) {
        w_sum += weights[i];
    }
    if (fabsf(w_sum) < 1e-12f) {
        w_sum = 1.0f;
    }

    float32_t moving_avg[4];
    quat_norm(initial, moving_avg);
    for (int iter = 0; iter < GRAD_DESCENT_MAX_ITERS; iter++) {
        get_error_vectors(quat_sigmas, moving_avg, error_vecs);
        float32_t avg_err[3] = {0, 0, 0};
        for (int i = 0; i < NUM_SIGMAS; i++) {
            for (int j = 0; j < 3; j++) {
                avg_err[j] += weights[i] * error_vecs[i * 3 + j] / w_sum;
            }
        }
        if (l2_norm(avg_err, 3) < GRAD_DESCENT_TOL) {
            memcpy(avg_quat, moving_avg, sizeof(float32_t) * 4);
            return;
        }
        float32_t avg_error_quat[4];
        rotationvec2quat(avg_err, avg_error_quat);
        quat_multiply(moving_avg, avg_error_quat, moving_avg);
        quat_norm(moving_avg, moving_avg);
    }
    get_error_vectors(quat_sigmas, moving_avg, error_vecs);
    memcpy(avg_quat, moving_avg, sizeof(float32_t) * 4);
}

ukf_status_t iterate(const float* error_state, const float* quat_state, const float* cov,
                     const float* body, const float* ref, int num_vecs, const float* gyro,
                     const float* Q, const float* R, float dt, float* new_err_state,
                     float* new_quat_state, float* new_P) {
    // Copy inputs first so outputs may alias them, and so failures can return the inputs
    float x_in[STATE_SIZE];
    float q_in[4];
    float P_in[STATE_SIZE * STATE_SIZE];
    memcpy(x_in, error_state, sizeof(x_in));
    memcpy(q_in, quat_state, sizeof(q_in));
    memcpy(P_in, cov, sizeof(P_in));

#define UKF_FAIL(code)                                                                             \
    do {                                                                                           \
        memcpy(new_err_state, x_in, sizeof(x_in));                                                 \
        memcpy(new_quat_state, q_in, sizeof(q_in));                                                \
        memcpy(new_P, P_in, sizeof(P_in));                                                         \
        return (code);                                                                             \
    } while (0)

    if (num_vecs < 0 || num_vecs > MAX_MSMT_VECS || !(dt >= 0.0f) || !all_finite(x_in, STATE_SIZE)
        || !all_finite(q_in, 4) || !all_finite(gyro, 3) || !all_finite(Q, STATE_SIZE * STATE_SIZE)) {
        UKF_FAIL(UKF_ERR_ARGS);
    }
    const int m = 3 * num_vecs;
    if (num_vecs > 0 && (!all_finite(body, m) || !all_finite(ref, m) || !all_finite(R, m * m))) {
        UKF_FAIL(UKF_ERR_ARGS);
    }

    // --- Predict ---
    float P_Q[STATE_SIZE * STATE_SIZE];
    for (int i = 0; i < STATE_SIZE * STATE_SIZE; i++) {
        P_Q[i] = P_in[i] + Q[i];
    }
    if (!ensure_psd(P_Q, STATE_SIZE)) {
        UKF_FAIL(UKF_ERR_CHOLESKY);
    }

    float lambda = calculate_lambda(STATE_SIZE, UKF_ALPHA, UKF_KAPPA);
    float mean_weights[NUM_SIGMAS];
    float cov_weights[NUM_SIGMAS];
    get_weights(lambda, STATE_SIZE, UKF_ALPHA, UKF_BETA, cov_weights, mean_weights);

    float error_sigmas[NUM_SIGMAS * STATE_SIZE];
    if (!get_sigma_points(lambda, x_in, P_Q, error_sigmas)) {
        UKF_FAIL(UKF_ERR_CHOLESKY);
    }

    float quat_sigmas[NUM_SIGMAS * QUAT_SIGMA_SIZE];
    error_sigmas_to_quat_sigmas(error_sigmas, q_in, quat_sigmas);
    propagate_sigmas(quat_sigmas, gyro, dt);

    // Mean attitude, starting from the propagated central sigma point
    float avg_quat[4];
    float err_vecs[NUM_SIGMAS * 3];
    grad_descent(quat_sigmas, &quat_sigmas[0], mean_weights, avg_quat, err_vecs);

    float propagated_errors[NUM_SIGMAS * STATE_SIZE];
    for (int i = 0; i < NUM_SIGMAS; i++) {
        for (int j = 0; j < 3; j++) {
            propagated_errors[i * STATE_SIZE + j] = err_vecs[i * 3 + j];
            propagated_errors[i * STATE_SIZE + 3 + j] = quat_sigmas[i * QUAT_SIGMA_SIZE + 4 + j];
        }
    }

    float mean_err[STATE_SIZE] = {0};
    for (int i = 0; i < NUM_SIGMAS; i++) {
        for (int j = 0; j < STATE_SIZE; j++) {
            mean_err[j] += mean_weights[i] * propagated_errors[i * STATE_SIZE + j];
        }
    }

    float dx[NUM_SIGMAS * STATE_SIZE];
    for (int i = 0; i < NUM_SIGMAS; i++) {
        for (int j = 0; j < STATE_SIZE; j++) {
            dx[i * STATE_SIZE + j] = propagated_errors[i * STATE_SIZE + j] - mean_err[j];
        }
    }

    float P_hat[STATE_SIZE * STATE_SIZE] = {0};
    for (int i = 0; i < NUM_SIGMAS; i++) {
        for (int r = 0; r < STATE_SIZE; r++) {
            for (int c = 0; c < STATE_SIZE; c++) {
                P_hat[r * STATE_SIZE + c] += cov_weights[i] * dx[i * STATE_SIZE + r] * dx[i * STATE_SIZE + c];
            }
        }
    }

    float x_hat[STATE_SIZE];
    float P[STATE_SIZE * STATE_SIZE];
    memcpy(x_hat, mean_err, sizeof(x_hat));
    memcpy(P, P_hat, sizeof(P));

    // --- Update ---
    if (num_vecs > 0) {
        float sigma_msmts[NUM_SIGMAS * MAX_MSMT_SIZE];
        get_sigma_measurements(quat_sigmas, ref, num_vecs, sigma_msmts);

        float mean_msmt[MAX_MSMT_SIZE] = {0};
        for (int i = 0; i < NUM_SIGMAS; i++) {
            for (int j = 0; j < m; j++) {
                mean_msmt[j] += mean_weights[i] * sigma_msmts[i * m + j];
            }
        }

        float dz[NUM_SIGMAS * MAX_MSMT_SIZE];
        for (int i = 0; i < NUM_SIGMAS; i++) {
            for (int j = 0; j < m; j++) {
                dz[i * m + j] = sigma_msmts[i * m + j] - mean_msmt[j];
            }
        }

        float P_xz[STATE_SIZE * MAX_MSMT_SIZE] = {0};
        float P_vv[MAX_MSMT_SIZE * MAX_MSMT_SIZE];
        memcpy(P_vv, R, sizeof(float) * m * m);
        for (int i = 0; i < NUM_SIGMAS; i++) {
            for (int r = 0; r < STATE_SIZE; r++) {
                for (int c = 0; c < m; c++) {
                    P_xz[r * m + c] += cov_weights[i] * dx[i * STATE_SIZE + r] * dz[i * m + c];
                }
            }
            for (int r = 0; r < m; r++) {
                for (int c = 0; c < m; c++) {
                    P_vv[r * m + c] += cov_weights[i] * dz[i * m + r] * dz[i * m + c];
                }
            }
        }

        // arm_mat_inverse_f32 destroys its input, so invert a scratch copy. (The old code
        // inverted P_vv in place and then reused the destroyed P_vv in K P_vv K^T.)
        float P_vv_scratch[MAX_MSMT_SIZE * MAX_MSMT_SIZE];
        float P_vv_inv[MAX_MSMT_SIZE * MAX_MSMT_SIZE];
        memcpy(P_vv_scratch, P_vv, sizeof(float) * m * m);
        arm_matrix_instance_f32 P_vv_scratch_mat = {(uint16_t)m, (uint16_t)m, P_vv_scratch};
        arm_matrix_instance_f32 P_vv_inv_mat = {(uint16_t)m, (uint16_t)m, P_vv_inv};
        if (arm_mat_inverse_f32(&P_vv_scratch_mat, &P_vv_inv_mat) != ARM_MATH_SUCCESS
            || !all_finite(P_vv_inv, m * m)) {
            UKF_FAIL(UKF_ERR_SINGULAR);
        }

        // K = P_xz P_vv^-1
        float K[STATE_SIZE * MAX_MSMT_SIZE];
        arm_matrix_instance_f32 P_xz_mat = {STATE_SIZE, (uint16_t)m, P_xz};
        arm_matrix_instance_f32 K_mat = {STATE_SIZE, (uint16_t)m, K};
        arm_mat_mult_f32(&P_xz_mat, &P_vv_inv_mat, &K_mat);

        for (int r = 0; r < STATE_SIZE; r++) {
            float corr = 0.0f;
            for (int c = 0; c < m; c++) {
                corr += K[r * m + c] * (body[c] - mean_msmt[c]);
            }
            x_hat[r] += corr;
        }

        // P = P_hat - K P_vv K^T  (== P_hat - K P_xz^T)
        for (int r = 0; r < STATE_SIZE; r++) {
            for (int c = 0; c < STATE_SIZE; c++) {
                float acc = 0.0f;
                for (int j = 0; j < m; j++) {
                    acc += K[r * m + j] * P_xz[c * m + j];
                }
                P[r * STATE_SIZE + c] -= acc;
            }
        }
    }

    if (!ensure_psd(P, STATE_SIZE)) {
        UKF_FAIL(UKF_ERR_CHOLESKY);
    }

    // Fold the attitude correction into the quaternion (local convention)
    float x_hat_rot[4];
    rotationvec2quat(x_hat, x_hat_rot);
    float new_guess[4];
    quat_multiply(avg_quat, x_hat_rot, new_guess);
    quat_norm(new_guess, new_guess);
    x_hat[0] = 0.0f;
    x_hat[1] = 0.0f;
    x_hat[2] = 0.0f;

    if (!all_finite(x_hat, STATE_SIZE) || !all_finite(new_guess, 4) || !all_finite(P, STATE_SIZE * STATE_SIZE)) {
        UKF_FAIL(UKF_ERR_NONFINITE);
    }

    memcpy(new_err_state, x_hat, sizeof(x_hat));
    memcpy(new_quat_state, new_guess, sizeof(new_guess));
    memcpy(new_P, P, sizeof(P));
    return UKF_OK;
#undef UKF_FAIL
}
