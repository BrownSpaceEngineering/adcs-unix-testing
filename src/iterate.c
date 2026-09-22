#include "include/iterate.h"
#include "include/quat.h"
#include "arm_math.h"

#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <string.h>

enum { QUAT_SIGMA_SIZE = STATE_SIZE + 1, MAX_VECTOR_COUNT = 2 };

/*
 * UKF notation used in this file
 * ==============================
 *
 *   n       = STATE_SIZE = 6
 *   x       = [dtheta_x, dtheta_y, dtheta_z, bias_x, bias_y, bias_z]
 *   P       = covariance of x
 *   Q, R    = process-noise and measurement-noise covariance
 *   q       = [w, x, y, z], rotating body coordinates into reference coordinates
 *   a (*) b = Hamilton quaternion product, implemented by quat_multiply(a, b)
 *   q^-1    = quat_inv(q)
 *   rotvec2quat(r) = rotationvec2quat(r)
 *   quat2rotvec(q) = quat2rotationvec(q)
 *   Wm, Wc  = unscented-transform mean and covariance weights
 *   X[i, :] = Euclidean sigma point i
 *   z       = stacked body-frame vector observations
 *   r       = stacked reference-frame vectors
 *
 * A quaternion sigma point is stored as:
 *
 *   [q_w, q_x, q_y, q_z, bias_x, bias_y, bias_z]
 *
 * Matrices are flat row-major arrays, so A[row, column] is given by:
 *
 *   A[row * column_count + column]
 */
static const float32_t UKF_ALPHA = 1.0f;
static const float32_t UKF_BETA = 2.0f;
static const float32_t UKF_KAPPA = 0.0f;
static const float32_t COV_EPSILON = 1.0e-9f;

/*
 * Compute the scaled unscented-transform parameter:
 *
 *   lambda = alpha^2 (n + kappa) - n
 */
static float32_t sigma_lambda(void) {
    /* lambda = alpha^2 (n + kappa) - n. */
    return UKF_ALPHA * UKF_ALPHA * (STATE_SIZE + UKF_KAPPA) - STATE_SIZE;
}

/*
 * Return true exactly when every input is finite:
 *
 *   result = AND(i = 0 ... count-1, isfinite(values[i]))
 */
static bool all_finite(const float32_t *values, int count) {
    if (values == NULL) return false;
    for (int i = 0; i < count; ++i) {
        if (!isfinite(values[i])) return false;
    }
    return true;
}

/*
 * Normalize a three-dimensional vector:
 *
 *   norm = ||source||_2 = sqrt(sum(i = 0 ... 2, source[i]^2))
 *   unit = source / norm
 */
static arm_status normalize_vector3(const float32_t *source, float32_t *unit) {
    if (unit == NULL) return ARM_MATH_ARGUMENT_ERROR;
    if (!all_finite(source, 3)) return ARM_MATH_NANINF;

    /* norm = ||source||_2 = sqrt(sum_i source[i]^2). */
    float32_t norm = sqrtf(source[0] * source[0]
                           + source[1] * source[1]
                           + source[2] * source[2]);
    if (!isfinite(norm)) return ARM_MATH_NANINF;
    if (norm <= 0.0f) return ARM_MATH_ARGUMENT_ERROR;

    /* unit = source / norm. */
    for (int axis = 0; axis < 3; ++axis) {
        unit[axis] = source[axis] / norm;
    }
    return all_finite(unit, 3) ? ARM_MATH_SUCCESS : ARM_MATH_NANINF;
}

/*
 * Prepare a covariance for Cholesky decomposition:
 *
 *   P_sym = (P + P^T) / 2
 *   jitter_0 = max(epsilon, 8 FLT_EPSILON max(|diag(P_sym)|))
 *   P_try = P_sym                            on the first attempt
 *   P_try <- P_try + jitter I                after a failed attempt
 *   jitter <- 10 jitter                      after each addition
 *
 * Success means that P_try = L L^T. Exhausting the bounded attempts reports
 * an indefinite covariance instead of applying an unbounded repair.
 */
static arm_status ensure_positive_definite(float32_t *matrix) {
    if (!all_finite(matrix, STATE_SIZE * STATE_SIZE)) return ARM_MATH_NANINF;

    /* P_sym = (P + P^T) / 2. */
    for (int r = 0; r < STATE_SIZE; ++r) {
        for (int c = r + 1; c < STATE_SIZE; ++c) {
            float32_t average = 0.5f * (matrix[r * STATE_SIZE + c]
                                    + matrix[c * STATE_SIZE + r]);
            matrix[r * STATE_SIZE + c] = matrix[c * STATE_SIZE + r] = average;
        }
    }
    if (!all_finite(matrix, STATE_SIZE * STATE_SIZE)) return ARM_MATH_NANINF;
    /*
     * Try P_try <- P_try + jitter I until P_try = L L^T exists. If the bounded
     * retries cannot repair P, report the decomposition failure.
     */
    float32_t scale = 0.0f;
    for (int r = 0; r < STATE_SIZE; ++r) {
        if (fabsf(matrix[r * STATE_SIZE + r]) > scale)
            scale = fabsf(matrix[r * STATE_SIZE + r]);
    }
    float32_t jitter = fmaxf(COV_EPSILON, 8.0f * FLT_EPSILON * scale);
    for (int attempt = 0; attempt < 6; ++attempt) {
        float32_t input_data[STATE_SIZE * STATE_SIZE], lower[STATE_SIZE * STATE_SIZE] = {0};
        memcpy(input_data, matrix, sizeof(input_data));
        arm_matrix_instance_f32 input = {STATE_SIZE, STATE_SIZE, input_data};
        arm_matrix_instance_f32 factor = {STATE_SIZE, STATE_SIZE, lower};
        if (arm_mat_cholesky_f32(&input, &factor) == ARM_MATH_SUCCESS
            && all_finite(lower, STATE_SIZE * STATE_SIZE)) return ARM_MATH_SUCCESS;
        for (int r = 0; r < STATE_SIZE; ++r) matrix[r * STATE_SIZE + r] += jitter;
        jitter *= 10.0f;
    }
    return ARM_MATH_DECOMPOSITION_FAILURE;
}

/*
 * Compute scaled unscented-transform weights for 2n + 1 sigma points:
 *
 *   lambda = alpha^2 (n + kappa) - n
 *   Wm[0]  = lambda / (n + lambda)
 *   Wc[0]  = Wm[0] + (1 - alpha^2 + beta)
 *   Wm[i]  = Wc[i] = 1 / (2 (n + lambda)),  i = 1 ... 2n
 */
static void sigma_weights(float32_t *mean_weights, float32_t *cov_weights) {
    /* lambda = alpha^2 (n + kappa) - n. */
    const float32_t lambda = sigma_lambda();

    /* Wm[i] = Wc[i] = 1 / (2 (n + lambda)), i = 1 ... 2n. */
    const float32_t side_weight = 0.5f / (STATE_SIZE + lambda);
    for (int i = 0; i < NUM_SIGMAS; ++i) {
        mean_weights[i] = side_weight;
        cov_weights[i] = side_weight;
    }
    /* Wm[0] = lambda / (n + lambda). */
    mean_weights[0] = lambda / (STATE_SIZE + lambda);

    /* Wc[0] = Wm[0] + (1 - alpha^2 + beta). */
    cov_weights[0] = mean_weights[0] + 1.0f - UKF_ALPHA * UKF_ALPHA + UKF_BETA;
}

/*
 * Generate Euclidean sigma points from state x and covariance P:
 *
 *   P_sym     = (P + P^T) / 2
 *   L L^T    = (n + lambda) P_sym
 *   X[0]     = x
 *   X[k+1]   = x + L[:, k]
 *   X[k+1+n] = x - L[:, k],  k = 0 ... n-1
 */
static arm_status sigma_points(const float32_t *state, const float32_t *covariance,
                               float32_t *sigmas) {
    float32_t scaled[STATE_SIZE * STATE_SIZE];
    float32_t lower[STATE_SIZE * STATE_SIZE] = {0};
    memcpy(scaled, covariance, sizeof(scaled));

    /* scaled = P_sym = (P + P^T) / 2. */
    for (int r = 0; r < STATE_SIZE; ++r) {
        for (int c = r + 1; c < STATE_SIZE; ++c) {
            float32_t average = 0.5f * (scaled[r * STATE_SIZE + c]
                                    + scaled[c * STATE_SIZE + r]);
            scaled[r * STATE_SIZE + c] = scaled[c * STATE_SIZE + r] = average;
        }
    }
    /* scaled = (n + lambda) P_sym. */
    const float32_t lambda = sigma_lambda();
    for (int i = 0; i < STATE_SIZE * STATE_SIZE; ++i) scaled[i] *= STATE_SIZE + lambda;

    /* L L^T = scaled = (n + lambda) P_sym. */
    arm_matrix_instance_f32 input = {STATE_SIZE, STATE_SIZE, scaled};
    arm_matrix_instance_f32 factor = {STATE_SIZE, STATE_SIZE, lower};
    arm_status status = arm_mat_cholesky_f32(&input, &factor);
    if (status != ARM_MATH_SUCCESS) return status;
    if (!all_finite(lower, STATE_SIZE * STATE_SIZE)) return ARM_MATH_NANINF;

    /* X[0] = x. */
    memcpy(sigmas, state, sizeof(float32_t) * STATE_SIZE);

    /*
     * X[k+1]   = x + L[:, k]
     * X[k+1+n] = x - L[:, k],  k = 0 ... n-1
     */
    for (int k = 0; k < STATE_SIZE; ++k) {
        for (int j = 0; j < STATE_SIZE; ++j) {
            float32_t upper_row_value = lower[j * STATE_SIZE + k];
            sigmas[(k + 1) * STATE_SIZE + j] = state[j] + upper_row_value;
            sigmas[(k + 1 + STATE_SIZE) * STATE_SIZE + j] = state[j] - upper_row_value;
        }
    }
    return all_finite(sigmas, NUM_SIGMAS * STATE_SIZE)
           ? ARM_MATH_SUCCESS : ARM_MATH_NANINF;
}

/*
 * Lift Euclidean error sigma points onto SO(3) with a left error:
 *
 *   error_q[i] = rotvec2quat(dtheta[i])
 *   q[i]       = normalize(error_q[i] (*) q_guess)
 *   bias[i]    = X[i, 3:6]
 */
static arm_status quaternion_sigmas(const float32_t *error_sigmas,
                                    const float32_t *guess, float32_t *result) {
    for (int i = 0; i < NUM_SIGMAS; ++i) {
        float32_t error_q[4], q[4];
        if (!all_finite(&error_sigmas[i * STATE_SIZE], 3)) return ARM_MATH_NANINF;

        /* error_q = rotvec2quat(dtheta_i). */
        rotationvec2quat(&error_sigmas[i * STATE_SIZE], error_q);
        if (!all_finite(error_q, 4)) return ARM_MATH_NANINF;

        /* q_i = normalize(rotvec2quat(dtheta_i) (*) q_guess). */
        quat_multiply(error_q, guess, q); /* left-side attitude error */
        if (!all_finite(q, 4)) return ARM_MATH_NANINF;
        quat_norm(q, &result[i * QUAT_SIGMA_SIZE]);
        if (!all_finite(&result[i * QUAT_SIGMA_SIZE], 4)) return ARM_MATH_NANINF;
        memcpy(&result[i * QUAT_SIGMA_SIZE + 4], &error_sigmas[i * STATE_SIZE + 3],
               3 * sizeof(float32_t));
    }
    return ARM_MATH_SUCCESS;
}

/*
 * Propagate every quaternion sigma with a body-frame gyro increment:
 *
 *   omega[i] = gyro - bias[i]
 *   dq[i]    = rotvec2quat(omega[i] dt)
 *   q[i]+    = normalize(q[i] (*) dq[i])
 *   bias[i]+ = bias[i]
 */
static arm_status propagate_sigmas(const float32_t *input, const float32_t *gyro,
                                   float32_t dt, float32_t *result) {
    for (int i = 0; i < NUM_SIGMAS; ++i) {
        int offset = i * QUAT_SIGMA_SIZE;
        float32_t rotation[3], delta[4], propagated[4];

        /* rotation_i = omega_i dt = (gyro - bias_i) dt. */
        for (int j = 0; j < 3; ++j) rotation[j] = (gyro[j] - input[offset + 4 + j]) * dt;
        if (!all_finite(rotation, 3)) return ARM_MATH_NANINF;

        /* dq_i = rotvec2quat(rotation_i). */
        rotationvec2quat(rotation, delta);
        if (!all_finite(delta, 4)) return ARM_MATH_NANINF;

        /* q_i+ = normalize(q_i (*) dq_i); body-frame dq_i is on the right. */
        quat_multiply(&input[offset], delta, propagated); /* body-frame gyro */
        if (!all_finite(propagated, 4)) return ARM_MATH_NANINF;
        quat_norm(propagated, &result[offset]);
        if (!all_finite(&result[offset], 4)) return ARM_MATH_NANINF;

        /* copy the bias over */
        memcpy(&result[offset + 4], &input[offset + 4], 3 * sizeof(float32_t));
    }
    return ARM_MATH_SUCCESS;
}

/*
 * Express quaternion sigma points as left rotation-vector errors about q_mean:
 *
 *   relative[i] = q[i] (*) q_mean^-1
 *   dtheta[i]   = quat2rotvec(relative[i])
 */
static arm_status error_vectors(const float32_t *sigmas, const float32_t *mean,
                                float32_t *errors) {
    /* inverse = q_mean^-1. */
    float32_t inverse[4];
    quat_inv(mean, inverse);
    if (!all_finite(inverse, 4)) return ARM_MATH_NANINF;
    for (int i = 0; i < NUM_SIGMAS; ++i) {
        float32_t relative[4];
        /* relative_i = q_i (*) q_mean^-1. */
        quat_multiply(&sigmas[i * QUAT_SIGMA_SIZE], inverse, relative);
        if (!all_finite(relative, 4)) return ARM_MATH_NANINF;

        /* dtheta_i = quat2rotvec(relative_i). */
        quat2rotationvec(relative, &errors[i * 3]);
        if (!all_finite(&errors[i * 3], 3)) return ARM_MATH_NANINF;
    }
    return ARM_MATH_SUCCESS;
}

/*
 * Compute a weighted quaternion mean in its tangent space. Iterate:
 *
 *   e[i]   = quat2rotvec(q[i] (*) q_mean^-1)
 *   e_bar  = sum(i = 0 ... 2n, Wm[i] e[i])
 *   q_mean <- normalize(rotvec2quat(e_bar) (*) q_mean)
 *
 * Stop when ||e_bar||_2 < 1e-9, then recompute every e[i] about the final mean.
 */
static arm_status quaternion_mean(const float32_t *sigmas, const float32_t *guess,
                                  const float32_t *weights, float32_t *mean,
                                  float32_t *errors) {
    memcpy(mean, guess, 4 * sizeof(float32_t));
    if (!all_finite(mean, 4)) return ARM_MATH_NANINF;
    quat_norm(mean, mean);
    if (!all_finite(mean, 4)) return ARM_MATH_NANINF;
    arm_status status;
    for (int iteration = 0; iteration < 100; ++iteration) {
        /* e_i = quat2rotvec(q_i (*) q_mean^-1). */
        status = error_vectors(sigmas, mean, errors);
        if (status != ARM_MATH_SUCCESS) return status;

        /* average_err = sum_i Wm[i] e_i. */
        float32_t average[3] = {0};
        for (int i = 0; i < NUM_SIGMAS; ++i) {
            for (int j = 0; j < 3; ++j) average[j] += weights[i] * errors[i * 3 + j];
        }
        /* Stop when ||average_err||_2 < 1e-9. */
        float32_t length = sqrtf(average[0] * average[0] + average[1] * average[1]
                             + average[2] * average[2]);
        if (!isfinite(length)) return ARM_MATH_NANINF;
        if (length < 1.0e-9f) break;
        float32_t correction[4], updated[4];

        /* q_mean <- normalize(rotvec2quat(average_err) (*) q_mean). */
        rotationvec2quat(average, correction);
        if (!all_finite(correction, 4)) return ARM_MATH_NANINF;
        quat_multiply(correction, mean, updated);
        if (!all_finite(updated, 4)) return ARM_MATH_NANINF;
        quat_norm(updated, mean);
        if (!all_finite(mean, 4)) return ARM_MATH_NANINF;
    }
    /* Errors used by P_hat must be relative to the final mean. */
    return error_vectors(sigmas, mean, errors);
}

/*
 * Predict each reference vector in body coordinates:
 *
 *   z[i,v] = q[i]^-1 (*) [0, reference[v]] (*) q[i]
 *
 * Here q[i] maps body coordinates to reference coordinates.
 */
static arm_status sigma_measurements(const float32_t *sigmas,
                                     const float32_t *references,
                                     int vector_count, float32_t *measurements) {
    int size = 3 * vector_count;
    for (int i = 0; i < NUM_SIGMAS; ++i) {
        float32_t inverse[4];

        /* q_i maps body -> reference, so reference -> body is q_i^-1. */
        quat_inv(&sigmas[i * QUAT_SIGMA_SIZE], inverse);
        if (!all_finite(inverse, 4)) return ARM_MATH_NANINF;
        for (int v = 0; v < vector_count; ++v) {
            /* z_i,v = q_i^-1 (*) [0, reference_v] (*) q_i. */
            quat_apply(inverse, &references[3 * v], &measurements[i * size + 3 * v]);
        }
    }
    return all_finite(measurements, NUM_SIGMAS * size)
           ? ARM_MATH_SUCCESS : ARM_MATH_NANINF;
}

/*
 * Compute the predicted Euclidean error state and covariance:
 *
 *   Xerr[i] = [quat2rotvec(q[i] (*) q_mean^-1), bias[i]]
 *   x_hat   = sum(i = 0 ... 2n, Wm[i] Xerr[i])
 *   dx[i]   = Xerr[i] - x_hat
 *   P_hat   = sum(i = 0 ... 2n, Wc[i] dx[i] dx[i]^T) + Q
 */
static arm_status predicted_error_statistics(const float32_t *propagated_sigmas,
                                             const float32_t *attitude_errors,
                                             const float32_t *mean_weights,
                                             const float32_t *covariance_weights,
                                             const float32_t *process_noise,
                                             float32_t *sigma_errors,
                                             float32_t *predicted_state,
                                             float32_t *predicted_covariance) {
    memset(predicted_state, 0, STATE_SIZE * sizeof(float32_t));
    memset(predicted_covariance, 0, STATE_SIZE * STATE_SIZE * sizeof(float32_t));

    for (int sigma = 0; sigma < NUM_SIGMAS; ++sigma) {
        /* Xerr[i] = [quat2rotvec(q_i (*) q_mean^-1), bias_i]. */
        for (int axis = 0; axis < 3; ++axis) {
            sigma_errors[sigma * STATE_SIZE + axis] = attitude_errors[sigma * 3 + axis];
            sigma_errors[sigma * STATE_SIZE + axis + 3]
                = propagated_sigmas[sigma * QUAT_SIGMA_SIZE + axis + 4];
        }
        /* x_hat = sum_i Wm[i] Xerr[i]. */
        for (int axis = 0; axis < STATE_SIZE; ++axis) {
            predicted_state[axis] += mean_weights[sigma]
                                     * sigma_errors[sigma * STATE_SIZE + axis];
        }
    }

    /* P_hat = sum_i Wc[i] (Xerr[i] - x_hat)(Xerr[i] - x_hat)^T. */
    for (int sigma = 0; sigma < NUM_SIGMAS; ++sigma) {
        for (int row = 0; row < STATE_SIZE; ++row) {
            float32_t row_error = sigma_errors[sigma * STATE_SIZE + row] - predicted_state[row];
            for (int column = 0; column < STATE_SIZE; ++column) {
                float32_t column_error = sigma_errors[sigma * STATE_SIZE + column]
                                         - predicted_state[column];
                predicted_covariance[row * STATE_SIZE + column]
                    += covariance_weights[sigma] * row_error * column_error;
            }
        }
    }
    /* P_hat <- P_hat + Q, exactly once per filter call. */
    for (int i = 0; i < STATE_SIZE * STATE_SIZE; ++i) {
        predicted_covariance[i] += process_noise[i];
    }
    return all_finite(predicted_covariance, STATE_SIZE * STATE_SIZE)
           ? ARM_MATH_SUCCESS : ARM_MATH_NANINF;
}

/*
 * Correct the predicted state with one or two observed body vectors:
 *
 *   z_hat = sum(i = 0 ... 2n, Wm[i] z[i])
 *   dx[i] = Xerr[i] - x_hat
 *   dz[i] = z[i] - z_hat
 *   P_xz  = sum(i = 0 ... 2n, Wc[i] dx[i] dz[i]^T)
 *   S     = sum(i = 0 ... 2n, Wc[i] dz[i] dz[i]^T) + R
 *   K     = P_xz S^-1
 *   x_new = x_hat + K (z - z_hat)
 *   P_new = P_hat - K S K^T
 */
static arm_status correct_with_measurements(const float32_t *propagated_sigmas,
                                            const float32_t *sigma_errors,
                                            const float32_t *predicted_state,
                                            const float32_t *mean_weights,
                                            const float32_t *covariance_weights,
                                            const float32_t *body_vectors,
                                            const float32_t *reference_vectors,
                                            const float32_t *measurement_noise,
                                            int vector_count,
                                            float32_t *corrected_state,
                                            float32_t *corrected_covariance) {
    const int measurement_size = 3 * vector_count;
    float32_t measurements[NUM_SIGMAS * MSMT_SIZE];
    float32_t measurement_mean[MSMT_SIZE] = {0};
    float32_t state_measurement_covariance[STATE_SIZE * MSMT_SIZE] = {0};
    float32_t innovation_covariance[MSMT_SIZE * MSMT_SIZE] = {0};

    arm_status status = sigma_measurements(propagated_sigmas, reference_vectors,
                                           vector_count, measurements);
    if (status != ARM_MATH_SUCCESS) return status;

    /* z_hat = sum_i Wm[i] z_i. */
    for (int sigma = 0; sigma < NUM_SIGMAS; ++sigma) {
        for (int component = 0; component < measurement_size; ++component) {
            measurement_mean[component] += mean_weights[sigma]
                                            * measurements[sigma * measurement_size + component];
        }
    }

    /*
     * dx_i = Xerr[i] - x_hat
     * dz_i = z_i - z_hat
     * P_xz = sum_i Wc[i] dx_i dz_i^T
     * P_zz = sum_i Wc[i] dz_i dz_i^T
     */
    for (int sigma = 0; sigma < NUM_SIGMAS; ++sigma) {
        for (int row = 0; row < STATE_SIZE; ++row) {
            float32_t state_error = sigma_errors[sigma * STATE_SIZE + row]
                                    - predicted_state[row];
            for (int column = 0; column < measurement_size; ++column) {
                float32_t measurement_error = measurements[sigma * measurement_size + column]
                                              - measurement_mean[column];
                state_measurement_covariance[row * measurement_size + column]
                    += covariance_weights[sigma] * state_error * measurement_error;
            }
        }
        for (int row = 0; row < measurement_size; ++row) {
            float32_t row_error = measurements[sigma * measurement_size + row]
                                  - measurement_mean[row];
            for (int column = 0; column < measurement_size; ++column) {
                float32_t column_error = measurements[sigma * measurement_size + column]
                                         - measurement_mean[column];
                innovation_covariance[row * measurement_size + column]
                    += covariance_weights[sigma] * row_error * column_error;
            }
        }
    }
    /* S = P_zz + R. */
    for (int i = 0; i < measurement_size * measurement_size; ++i) {
        innovation_covariance[i] += measurement_noise[i];
    }
    if (!all_finite(innovation_covariance, measurement_size * measurement_size)) {
        return ARM_MATH_NANINF;
    }

    /* S_inverse = S^-1. CMSIS overwrites its input, so invert a copy. */
    float32_t inverse_input[MSMT_SIZE * MSMT_SIZE];
    float32_t inverse[MSMT_SIZE * MSMT_SIZE];
    memcpy(inverse_input, innovation_covariance,
           (size_t)measurement_size * measurement_size * sizeof(float32_t));
    arm_matrix_instance_f32 matrix = {(uint16_t)measurement_size,
                                      (uint16_t)measurement_size, inverse_input};
    arm_matrix_instance_f32 inverse_matrix = {(uint16_t)measurement_size,
                                              (uint16_t)measurement_size, inverse};
    status = arm_mat_inverse_f32(&matrix, &inverse_matrix);
    if (status != ARM_MATH_SUCCESS) return status;
    if (!all_finite(inverse, measurement_size * measurement_size)) return ARM_MATH_NANINF;

    /* K = P_xz S^-1; x_new = x_hat + K(z - z_hat). */
    float32_t kalman_gain[STATE_SIZE * MSMT_SIZE] = {0};
    for (int row = 0; row < STATE_SIZE; ++row) {
        for (int column = 0; column < measurement_size; ++column) {
            for (int k = 0; k < measurement_size; ++k) {
                kalman_gain[row * measurement_size + column]
                    += state_measurement_covariance[row * measurement_size + k]
                       * inverse[k * measurement_size + column];
            }
            corrected_state[row] += kalman_gain[row * measurement_size + column]
                                    * (body_vectors[column] - measurement_mean[column]);
        }
    }

    /* P_new = P_hat - K S K^T. */
    for (int row = 0; row < STATE_SIZE; ++row) {
        for (int column = 0; column < STATE_SIZE; ++column) {
            float32_t reduction = 0.0f;
            for (int j = 0; j < measurement_size; ++j) {
                for (int k = 0; k < measurement_size; ++k) {
                    reduction += kalman_gain[row * measurement_size + j]
                                 * innovation_covariance[j * measurement_size + k]
                                 * kalman_gain[column * measurement_size + k];
                }
            }
            corrected_covariance[row * STATE_SIZE + column] -= reduction;
        }
    }
    return ARM_MATH_SUCCESS;
}

/*
 * Perform one complete UKF step.
 *
 * Inputs at step k:
 *
 *   x_k = error_state = [dtheta_k, bias_k]
 *   q_k = attitude_quaternion
 *   P_k = covariance
 *
 * 1. Construct Euclidean sigma points:
 *
 *   lambda = alpha^2 (n + kappa) - n
 *   Wm[0]  = lambda / (n + lambda)
 *   Wc[0]  = Wm[0] + (1 - alpha^2 + beta)
 *   Wm[i]  = Wc[i] = 1 / (2 (n + lambda)),  i = 1 ... 2n
 *   L L^T  = (n + lambda) P_guarded
 *   X[0,:] = x_k
 *   X[j+1,:]   = x_k + L[:, j]
 *   X[j+1+n,:] = x_k - L[:, j],  j = 0 ... n-1
 *
 * 2. Lift the attitude errors into quaternions and propagate each sigma:
 *
 *   q[i]    = normalize(rotvec2quat(X[i,0:3]) (*) q_k)
 *   b[i]    = X[i,3:6]
 *   q[i]+   = normalize(q[i] (*) rotvec2quat((gyro - b[i]) dt))
 *   b[i]+   = b[i]
 *
 * The attitude error is multiplied on the left. The body-frame gyro
 * increment is multiplied on the right.
 *
 * 3. Recover the predicted mean and covariance:
 *
 *   q_bar solves sum(i, Wm[i] quat2rotvec(q[i]+ (*) q_bar^-1)) = 0
 *   Xerr[i,:] = [quat2rotvec(q[i]+ (*) q_bar^-1), b[i]+]
 *   x_hat     = sum(i, Wm[i] Xerr[i,:])
 *   dx[i,:]   = Xerr[i,:] - x_hat
 *   P_hat     = sum(i, Wc[i] dx[i,:]^T dx[i,:]) + Q
 *
 * 4. If do_update is true, apply the vector-measurement correction:
 *
 *   z[v]   <- z[v] / ||z[v]||_2
 *   r[v]   <- r[v] / ||r[v]||_2
 *   Z[i,v] = q[i]+^-1 (*) [0, r[v]] (*) q[i]+
 *   z_hat  = sum(i, Wm[i] Z[i,:])
 *   dz[i,:] = Z[i,:] - z_hat
 *   P_xz   = sum(i, Wc[i] dx[i,:]^T dz[i,:])
 *   S      = sum(i, Wc[i] dz[i,:]^T dz[i,:]) + R
 *   K      = P_xz S^-1
 *   x_new  = x_hat + K (z - z_hat)
 *   P_new  = P_hat - K S K^T
 *
 * For a predict-only call, x_new = x_hat and P_new = P_hat.
 *
 * 5. Inject the corrected attitude error into the nominal quaternion:
 *
 *   q_new = normalize(rotvec2quat(x_new[0:3]) (*) q_bar)
 *
 * The returned error state still contains x_new[0:3] so callers can inspect
 * the applied correction. Before the next call, the caller sets those three
 * components to zero and retains x_new[3:6] as the gyro-bias estimate.
 *
 * Code-name mapping:
 *
 *   X             = error_sigma_points
 *   [q[i], b[i]]  = quaternion_sigma_points
 *   [q[i]+, b[i]+]= propagated_sigma_points
 *   q_bar         = mean_quaternion
 *   Xerr          = sigma_errors
 *   x_hat, P_hat  = predicted_state, predicted_covariance
 *   x_new          = corrected_state
 */
arm_status iterate(const float32_t *error_state, const float32_t *attitude_quaternion,
                   const float32_t *covariance, const float32_t *body_vectors,
                   const float32_t *reference_vectors, const float32_t *gyro_measurement,
                   const float32_t *process_noise, const float32_t *measurement_noise, float32_t dt,
                   int vector_count, bool do_update,
                   float32_t *new_error_state, float32_t *new_attitude_quaternion,
                   float32_t *new_covariance) {
    /* Validate first: outputs remain untouched if any stage fails. */
    if (vector_count < 1 || vector_count > MAX_VECTOR_COUNT || dt <= 0.0f
        || error_state == NULL || attitude_quaternion == NULL || covariance == NULL
        || gyro_measurement == NULL || process_noise == NULL
        || new_error_state == NULL || new_attitude_quaternion == NULL
        || new_covariance == NULL) {
        return ARM_MATH_ARGUMENT_ERROR;
    }
    if (!isfinite(dt) || !all_finite(error_state, STATE_SIZE)
        || !all_finite(attitude_quaternion, 4) || !all_finite(covariance, STATE_SIZE * STATE_SIZE)
        || !all_finite(gyro_measurement, 3) || !all_finite(process_noise, STATE_SIZE * STATE_SIZE))
        return ARM_MATH_NANINF;

    /* Require ||q_k||_2^2 > 0. */
    float32_t quaternion_norm_squared = 0.0f;
    for (int axis = 0; axis < 4; ++axis) {
        quaternion_norm_squared += attitude_quaternion[axis] * attitude_quaternion[axis];
    }
    if (!isfinite(quaternion_norm_squared)) return ARM_MATH_NANINF;
    if (quaternion_norm_squared <= 0.0f) return ARM_MATH_ARGUMENT_ERROR;

    const int measurement_size = 3 * vector_count;
    float32_t normalized_body_vectors[MSMT_SIZE];
    float32_t normalized_reference_vectors[MSMT_SIZE];
    if (do_update) {
        if (body_vectors == NULL || reference_vectors == NULL || measurement_noise == NULL) {
            return ARM_MATH_ARGUMENT_ERROR;
        }
        if (!all_finite(body_vectors, measurement_size)
            || !all_finite(reference_vectors, measurement_size)
            || !all_finite(measurement_noise, measurement_size * measurement_size)) {
            return ARM_MATH_NANINF;
        }
        /* Step 0: z[v] <- z[v]/||z[v]||_2 and r[v] <- r[v]/||r[v]||_2. */
        for (int vector = 0; vector < vector_count; ++vector) {
            arm_status vector_status = normalize_vector3(&body_vectors[3 * vector],
                                                         &normalized_body_vectors[3 * vector]);
            if (vector_status != ARM_MATH_SUCCESS) return vector_status;
            vector_status = normalize_vector3(&reference_vectors[3 * vector],
                                              &normalized_reference_vectors[3 * vector]);
            if (vector_status != ARM_MATH_SUCCESS) return vector_status;
        }
    }

    /*
     * Step 1a: Basic degenerate eigenvalue numerical guard
     */
    float32_t guarded_covariance[STATE_SIZE * STATE_SIZE];
    memcpy(guarded_covariance, covariance, sizeof(guarded_covariance));
    arm_status status = ensure_positive_definite(guarded_covariance);
    if (status != ARM_MATH_SUCCESS) return status;

    /*
     * Step 1b: Generate error sigma points
     */
    float32_t error_sigma_points[NUM_SIGMAS * STATE_SIZE];
    status = sigma_points(error_state, guarded_covariance, error_sigma_points);
    if (status != ARM_MATH_SUCCESS) return status;

    float32_t quaternion_sigma_points[NUM_SIGMAS * QUAT_SIGMA_SIZE];
    float32_t propagated_sigma_points[NUM_SIGMAS * QUAT_SIGMA_SIZE];

    /*
     * Step 2a: Convert error sigma points to quaternion sigma points
     */
    status = quaternion_sigmas(error_sigma_points, attitude_quaternion, quaternion_sigma_points);
    if (status != ARM_MATH_SUCCESS) return status;

    /*
     * Step 2b: Propagate quaternion sigma points using gyro measurements
     */
    status = propagate_sigmas(quaternion_sigma_points, gyro_measurement, dt, propagated_sigma_points);
    if (status != ARM_MATH_SUCCESS) return status;

    /* Step 3a: calculate W_m[i] (mean weights) and W_c[i] (covariance weights) from alpha, beta, kappa, and n. */
    float32_t mean_weights[NUM_SIGMAS], covariance_weights[NUM_SIGMAS];
    sigma_weights(mean_weights, covariance_weights);
    float32_t mean_quaternion[4], attitude_errors[NUM_SIGMAS * 3];

    /*
     * Step 3b: Find the mean propagated quaternion using gradient descent. The attitude errors between all propagated quaternions and the mean is stored in attitude_errors
     */
    status = quaternion_mean(propagated_sigma_points, attitude_quaternion, mean_weights,
                             mean_quaternion, attitude_errors);
    if (status != ARM_MATH_SUCCESS) return status;

    float32_t sigma_errors[NUM_SIGMAS * STATE_SIZE];
    float32_t predicted_state[STATE_SIZE];
    float32_t predicted_covariance[STATE_SIZE * STATE_SIZE];

    /*
     * Step 3c: Calculate
     *
     *   Xerr[i,:] = [attitude_errors[i,:], b[i]+]
     *   x_hat = sum(i, Wm[i] Xerr[i,:])
     *   P_hat = sum(i, Wc[i] dx[i,:]^T dx[i,:]) + Q
     */
    status = predicted_error_statistics(propagated_sigma_points, attitude_errors,
                                        mean_weights, covariance_weights, process_noise,
                                        sigma_errors, predicted_state,
                                        predicted_covariance);
    if (status != ARM_MATH_SUCCESS) return status;

    /* Predict-only default: x_new = x_hat and P_new = P_hat.*/
    float32_t corrected_state[STATE_SIZE];
    memcpy(corrected_state, predicted_state, sizeof(corrected_state));

    /* If there is an update, do the Kalman update step */
    if (do_update) {
        /*
         *   Z[i,v] = q[i]+^-1 (*) [0,r[v]] (*) q[i]+
         *   z_hat  = sum(i, Wm[i] Z[i,:])
         *   K      = P_xz (P_zz + R)^-1
         *   x_new  = x_hat + K(z - z_hat)
         *   P_new  = P_hat - K(P_zz + R)K^T
         */
        status = correct_with_measurements(propagated_sigma_points, sigma_errors,
                                           predicted_state, mean_weights, covariance_weights,
                                           normalized_body_vectors, normalized_reference_vectors,
                                           measurement_noise, vector_count, corrected_state,
                                           predicted_covariance);
        if (status != ARM_MATH_SUCCESS) return status;
    }
    if (!all_finite(corrected_state, STATE_SIZE)) return ARM_MATH_NANINF;

    /* Reguard the outgoing predicted covariance */
    status = ensure_positive_definite(predicted_covariance);
    if (status != ARM_MATH_SUCCESS) return status;
    float32_t correction[4], result_quaternion[4];

    /* Step 5a: Turn the correction error vector into a quaternion */
    rotationvec2quat(corrected_state, correction);
    if (!all_finite(correction, 4)) return ARM_MATH_NANINF;

    /* Step 5b: result_quaternion = correction (*) q_bar. */
    quat_multiply(correction, mean_quaternion, result_quaternion);
    if (!all_finite(result_quaternion, 4)) return ARM_MATH_NANINF;

    /* Step 5c: Renormalize the resulting quaternion */
    quat_norm(result_quaternion, result_quaternion);
    if (!all_finite(result_quaternion, 4)) return ARM_MATH_NANINF;

    memcpy(new_error_state, corrected_state, sizeof(corrected_state));
    memcpy(new_attitude_quaternion, result_quaternion, sizeof(result_quaternion));
    memcpy(new_covariance, predicted_covariance, sizeof(predicted_covariance));
    return ARM_MATH_SUCCESS;
}
