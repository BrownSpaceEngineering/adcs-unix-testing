/*
 * iterate.c — Multiplicative Unscented Kalman Filter (MUKF)
 *
 * Implements a 6-state attitude + gyro-bias UKF using the multiplicative
 * quaternion error parameterisation described in:
 *   Kraft, E. (2003). "A Quaternion-based Unscented Kalman Filter for
 *   Orientation Tracking." UPenn KodLab technical report.
 *
 * State vector:  x = [dq_x, dq_y, dq_z,  b_x, b_y, b_z]
 *   dq_xyz  — Rodrigues (rotation-vector) attitude error relative to the
 *              current quaternion estimate.  Reset to zero after each update.
 *   b_xyz   — gyroscope bias estimate [rad/s].
 *
 * Quaternion convention: body-to-ECI, scalar-first [w, x, y, z].
 * Left-compose convention: sigma_q = error_q ⊗ q_mean  (matches the Python
 * reference implementation MUKF_FINAL.py).
 */

#include "include/iterate.h"
#include "include/laextension.h"
#include "declareFunctions.h"
#include "math.h"
#include "include/quat.h"

/* ── double_inv ────────────────────────────────────────────────────────────────
 * Gauss-Jordan in-place matrix inverse for an n×n double-precision matrix A.
 * Builds the augmented matrix [A | I], applies full row-reduction with partial
 * pivoting, then extracts the right half as A^{-1}.
 *
 * This replaces the linalg inv() call, which internally routes through a
 * single-precision LAPACK linsolve and would corrupt results on double arrays.
 */
static void double_inv(double* A, int n) {
    double aug[n * 2 * n];
    /* Build augmented matrix [A | I] */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) aug[i * 2*n + j]   = A[i * n + j];
        for (int j = 0; j < n; j++) aug[i * 2*n + n+j] = (i == j) ? 1.0 : 0.0;
    }
    for (int col = 0; col < n; col++) {
        /* Partial pivoting: swap the row with the largest pivot element
         * into the current position to improve numerical stability. */
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
        if (fabs(diag) < 1e-15) continue; /* singular column — output is undefined; caller must ensure invertibility */
        /* Normalise pivot row so the diagonal becomes 1 */
        for (int j = 0; j < 2*n; j++) aug[col * 2*n + j] /= diag;
        /* Eliminate all other rows in this column */
        for (int r = 0; r < n; r++) {
            if (r == col) continue;
            double factor = aug[r * 2*n + col];
            for (int j = 0; j < 2*n; j++) aug[r * 2*n + j] -= factor * aug[col * 2*n + j];
        }
    }
    /* Extract right half of augmented matrix — that is A^{-1} */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) A[i * n + j] = aug[i * 2*n + n+j];
}

/* ── ensure_psd ────────────────────────────────────────────────────────────────
 * Enforce positive semi-definiteness on a symmetric matrix in-place.
 *
 * Step 1 — Symmetrise: P ← (P + Pᵀ)/2.  Floating-point arithmetic can break
 *           exact symmetry; this restores it before the Cholesky test.
 *
 * Step 2 — Modified-Cholesky test: partially factor P and track the minimum
 *           diagonal residual r_i = P[i][i] − Σ_{k<i} L[i][k]².  For a PSD
 *           matrix every r_i > 0; for a non-PSD matrix at least one r_i ≤ 0.
 *           If min(r_i) < ε, add (ε − min_r + ε)·I to the diagonal — the
 *           smallest diagonal shift that guarantees the next Cholesky succeeds.
 *           This mirrors Python's ensure_positive_definite:
 *             P += (epsilon − np.min(eigvals) + epsilon) * I
 *           but without a full eigendecomposition.
 */
void ensure_psd(double* mat, int size) {
    double eps = 1e-9;

    /* Step 1: symmetrise */
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            double avg = (mat[i * size + j] + mat[j * size + i]) / 2.0;
            mat[i * size + j] = avg;
            mat[j * size + i] = avg;
        }
    }

    /* Step 2: modified-Cholesky residual test */
    double L[size * size];   /* lower-triangular work array */
    memset(L, 0, sizeof(L));
    double min_residual = 1e300;
    for (int i = 0; i < size; i++) {
        /* Off-diagonal entries of row i */
        for (int j = 0; j < i; j++) {
            double s = 0.0;
            for (int k = 0; k < j; k++) s += L[i * size + k] * L[j * size + k];
            double Ljj = L[j * size + j];
            L[i * size + j] = (Ljj > 0.0) ? (mat[i * size + j] - s) / Ljj : 0.0;
        }
        /* Diagonal residual for row i */
        double s = 0.0;
        for (int k = 0; k < i; k++) s += L[i * size + k] * L[i * size + k];
        double residual = mat[i * size + i] - s;
        if (residual < min_residual) min_residual = residual;
        L[i * size + i] = (residual > 0.0) ? sqrt(residual) : 0.0;
    }

    /* Add the minimum necessary diagonal shift to restore PSD */
    if (min_residual < eps) {
        double add = eps - min_residual + eps;
        for (int i = 0; i < size; i++) mat[i * size + i] += add;
    }
}

/* ── calculate_lambda ─────────────────────────────────────────────────────────
 * Computes the UKF scaling parameter λ = α²(n + κ) − n.
 *
 * κ = 9 − n  (tuned so that n + κ = 9, matching the Python reference which
 *              uses kappa = 9 − n for a 6-state filter).
 * α = 0.1    (controls sigma-point spread; set in iterate.h).
 *
 * For n=6: λ ≈ 0.01·9 − 6 = −5.91.  The negative value is intentional; it
 * places sigma points closer together relative to the covariance spread.
 */
double calculate_lambda(int n, double alpha){
    int k = 9 - n;
    return alpha * alpha * (n + k) - n;
}

/* ── get_sigma_points ─────────────────────────────────────────────────────────
 * Draws the 2n+1 deterministic UKF sigma points from the error-state mean
 * (zero vector) and covariance P.
 *
 * Sigma-point layout (each is a STATE_SIZE-element error vector):
 *   sigma_0          = state              (centre point)
 *   sigma_{1..n}     = state + col_k(L)   (positive spread)
 *   sigma_{n+1..2n}  = state − col_k(L)   (negative spread)
 * where L is the lower-Cholesky factor of (λ+n)·P, and col_k(L) is its k-th
 * column.  Using columns (not rows) of L ensures the sigma-point set correctly
 * reproduces the full covariance P including off-diagonal terms.
 */
void get_sigma_points(double lam, double* state, double* P, double* sigmas){
    /* Scale the covariance: P1 = (λ+n)·P */
    double P1[STATE_SIZE * STATE_SIZE];
    memcpy(P1, P, sizeof(double) * STATE_SIZE * STATE_SIZE);
    scale(P1, lam + STATE_SIZE, STATE_SIZE, STATE_SIZE);
    ensure_psd(P1, STATE_SIZE);

    /* Cholesky decomposition: P1 = L·Lᵀ  →  lower-triangular factor stored in L */
    double L[STATE_SIZE * STATE_SIZE];
    chol(P1, L, STATE_SIZE);

    /* Centre sigma point */
    for(int i = 0; i<STATE_SIZE; i++){
        sigmas[i] = state[i];
    }
    /* Positive sigma points: state + column k of L */
    for(int k = 0; k<STATE_SIZE; k++){
        int start = k * STATE_SIZE;
        int sigma_offset = STATE_SIZE;
        for(int i = 0; i<STATE_SIZE; i++){
            sigmas[sigma_offset + start + i] = state[i] + L[i * STATE_SIZE + k];
        }
    }
    /* Negative sigma points: state − column k of L */
    for(int k = 0; k<STATE_SIZE; k++){
        int start = k * STATE_SIZE;
        int sigma_offset = (STATE_SIZE + 1) * STATE_SIZE;
        for(int i = 0; i<STATE_SIZE; i++){
            sigmas[sigma_offset + start + i] = state[i] - L[i * STATE_SIZE + k];
        }
    }
}

/* ── get_weights ──────────────────────────────────────────────────────────────
 * Computes the 2n+1 mean and covariance weights for the UKF.
 *
 * All off-centre weights are equal: w_i = 1 / (2(n+λ))  for i ≥ 1.
 *
 * Centre weights differ between mean and covariance:
 *   mean_weights[0]  = λ / (n+λ)
 *   cov_weights[0]   = λ / (n+λ) + (1 − α² + β)
 * where β=2 is optimal for Gaussian distributions (adds fourth-order
 * moment information to the covariance weight).
 *
 * Note: with λ ≈ −5.91, mean_weights[0] ≈ −65.7.  This large negative centre
 * weight is intentional and is handled by gradient_descent which clips
 * negative weights before computing the weighted quaternion mean.
 */
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

/* ── error_sigmas_to_quat_sigmas ──────────────────────────────────────────────
 * Lifts each 6-element error-state sigma point to a 7-element quaternion sigma
 * point by composing the attitude error with the current quaternion estimate.
 *
 * For each sigma point i with error vector e = sigma[0:3]:
 *   1. Convert e to a unit quaternion δq via the rotation-vector map.
 *   2. Left-compose with the current estimate: q_sigma = δq ⊗ q_current.
 *      (Left-compose convention matches the Python reference implementation.)
 *   3. Normalise q_sigma.
 *   4. Pass the bias components sigma[3:6] through unchanged.
 *
 * Output layout per sigma: [q_w, q_x, q_y, q_z, b_x, b_y, b_z]  (7 elements).
 */
void error_sigmas_to_quat_sigmas(double* error_sigmas, double* current_guess, double* quat_sigmas){
    int n = NUM_SIGMAS;
    double q[n * (STATE_SIZE + 1)];
    for(int i = 0; i<n; i++){
        /* Extract 3-element rotation-vector attitude error */
        double error_rot[3] = {error_sigmas[i * STATE_SIZE],error_sigmas[i * STATE_SIZE + 1],error_sigmas[i * STATE_SIZE + 2]};
        /* Map rotation vector → unit quaternion δq */
        double error_quat[4];
        rotationvec2quat(error_rot, error_quat);
        /* Left-compose: q_sigma = δq ⊗ q_current */
        double new_rotation[4];
        quat_multiply(error_quat, current_guess, new_rotation);
        double norm_rot[4];
        quat_norm(new_rotation, norm_rot);
        /* Store quaternion part */
        for(int j = 0; j<4; j++){
            q[i * (STATE_SIZE + 1) + j] = norm_rot[j];
        }
        /* Pass bias components through unchanged */
        for(int j = 4; j < STATE_SIZE + 1; j++){
            q[i * (STATE_SIZE + 1) + j] = error_sigmas[i * STATE_SIZE + (j - 1)];
        }
    }
    memcpy(quat_sigmas, q, sizeof(double) * n * (STATE_SIZE + 1));
}

/* ── propagate_sigmas ─────────────────────────────────────────────────────────
 * Propagates each quaternion sigma point forward by one time step using
 * the gyroscope kinematic model (gyro-integration propagation):
 *
 *   ω = gyro − bias_sigma          (remove per-sigma bias estimate)
 *   δq = rotationvec2quat(ω · dt)  (small-angle rotation increment)
 *   q_new = q_sigma ⊗ δq           (right-compose body-frame increment)
 *
 * Bias components are propagated as a random walk (held constant here;
 * process noise Q is added to P_hat after propagation in iterate()).
 */
void propagate_sigmas(double* quat_sigmas, double* gyro, double dt, double* propagated_sigmas){
    double p_sigmas[NUM_SIGMAS * (STATE_SIZE + 1)];
    for(int i = 0; i<NUM_SIGMAS; i++){
        int offset = (STATE_SIZE + 1) * i;
        /* Extract quaternion and bias from this sigma point */
        double rotation[4];
        double bias[3];
        for(int j = 0; j<4; j++){
            rotation[j] = quat_sigmas[offset + j];
        }
        for(int j = 4; j<7; j++){
            bias[j - 4] = quat_sigmas[offset + j];
        }
        /* ω = (gyro − bias) · dt */
        double omega[3];
        scale(bias, -1, 1, 3);          /* bias ← −bias (in-place) */
        add(gyro, bias, omega, 1, 3, 3); /* omega = gyro + (−bias) */
        scale(omega, dt, 1, 3);          /* omega ← omega · dt */
        /* Incremental rotation quaternion from angular velocity */
        double delta_q[4];
        rotationvec2quat(omega, delta_q);
        /* Integrate: q_new = q ⊗ δq */
        double propagated_rot[4];
        quat_multiply(rotation, delta_q, propagated_rot);
        double norm_rot[4];
        quat_norm(propagated_rot, norm_rot);

        int p_offset = i * (STATE_SIZE + 1);
        for(int j = 0; j<4; j++){
            p_sigmas[j + p_offset] = norm_rot[j];
        }
        /* Bias held constant (random-walk; Q added to P_hat separately) */
        for(int j = 4; j<7; j++){
            p_sigmas[j + p_offset] = quat_sigmas[offset + j];
        }
    }
    memcpy(propagated_sigmas, p_sigmas, sizeof(double) * NUM_SIGMAS * (STATE_SIZE + 1));
}

/* ── get_sigma_measurements ───────────────────────────────────────────────────
 * Passes each propagated quaternion sigma point through the measurement model:
 * for each reference vector ref_v (ECI frame), rotate it into the body frame
 * using the sigma-point quaternion to get the expected body-frame reading.
 *
 * Supports num_vecs = 1 (magnetometer only, eclipse) or
 *          num_vecs = 2 (magnetometer + sun/reference, sunlit).
 *
 * Output: msmts[i * msmt_size + v*3 + j] — expected body-frame measurement
 *         for sigma point i, reference vector v, component j.
 */
void get_sigma_measurements(double* q_sigmas, double* ref, int num_vecs, double* msmts){
    int msmt_size = 3 * num_vecs;
    double m[NUM_SIGMAS * MAX_MSMT_SIZE];
    for(int i = 0; i < NUM_SIGMAS; i++){
        /* Extract the quaternion for this sigma point */
        double sigma_rot[4];
        for(int j = 0; j < 4; j++){
            sigma_rot[j] = q_sigmas[i * (STATE_SIZE + 1) + j];
        }
        /* q^{-1} rotates ECI vectors into body frame */
        double sigma_ref_to_body[4];
        quat_inv(sigma_rot, sigma_ref_to_body);
        /* Rotate each reference vector into the body frame */
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

/* ── get_error_vectors ────────────────────────────────────────────────────────
 * Computes the Lie-algebra (rotation-vector) attitude errors between a set of
 * quaternions and a reference (mean) quaternion.
 *
 * Left-decomposition convention (matches Python reference):
 *   error_i = log(q_sigma_i ⊗ q_mean^{-1})
 * so that q_sigma_i = error_i ⊗ q_mean for each sigma point i.
 */
void get_error_vectors(int num_other_quats, double* other_quats, double* quat, double* error_vecs){
    double e[num_other_quats * 3];
    double x_inv[4];
    quat_inv(quat, x_inv);   /* x_inv = q_mean^{-1} */
    for(int i = 0; i<num_other_quats; i++){
        double other_quat[4];
        for(int j = 0; j<4; j++){
            other_quat[j] = other_quats[i * 4 + j];
        }
        /* error_q = q_sigma ⊗ q_mean^{-1} */
        double err_quat[4];
        quat_multiply(other_quat, x_inv, err_quat);
        /* Convert to rotation vector (Rodrigues parameterisation) */
        double err_vec[3];
        quat2rotationvec(err_quat, err_vec);
        for(int j = 0; j<3; j++){
            e[i * 3 + j] = err_vec[j];
        }
    }
    memcpy(error_vecs, e, sizeof(double) * num_other_quats * 3);
}

/* ── grad_descent ─────────────────────────────────────────────────────────────
 * Iterative gradient-descent algorithm to compute the weighted Fréchet mean
 * of a set of unit quaternions on SO(3).
 *
 * Algorithm (Markley et al. / Kraft):
 *   1. Initialise the running mean estimate as the prior quaternion.
 *   2. Compute per-sigma attitude error vectors relative to the current mean.
 *   3. Form the weighted average error vector ē.
 *   4. Convert ē → δq and left-compose: mean ← δq ⊗ mean.
 *   5. Repeat until |ē| < 1e-6 or 100 iterations.
 *
 * Weight clipping: mean_weights[0] ≈ −65.7 (large negative due to λ<0).
 * The negative centre weight is clipped to zero and the remaining weights
 * renormalised before computing ē, matching Python's safe_w treatment.
 * This prevents the centre sigma point from destabilising the mean estimate.
 *
 * Outputs: avg_quat — the converged weighted mean quaternion.
 *          errors   — final per-sigma rotation-vector errors relative to mean
 *                     (used by iterate() to build P_hat and P_xz).
 */
void grad_descent(double* propagated_quats, double* current_quat, double* mean_weights, double* avg_quat, double* errors){
    /* Clip negative weights and renormalise (safe_w in Python reference) */
    double safe_w[NUM_SIGMAS];
    double wsum = 0.0;
    for(int i = 0; i < NUM_SIGMAS; i++){
        safe_w[i] = (mean_weights[i] > 0.0) ? mean_weights[i] : 0.0;
        wsum += safe_w[i];
    }
    if(wsum > 0.0){
        for(int i = 0; i < NUM_SIGMAS; i++) safe_w[i] /= wsum;
    } else {
        /* Fallback: equal weights if all weights are non-positive */
        for(int i = 0; i < NUM_SIGMAS; i++) safe_w[i] = 1.0 / NUM_SIGMAS;
    }

    double moving_avg[4];
    memcpy(moving_avg, current_quat, sizeof(double) * 4);

    for(int iter = 0; iter < 100; iter++){
        /* Compute rotation-vector errors from current mean estimate */
        double error_vecs[3 * NUM_SIGMAS];
        get_error_vectors(NUM_SIGMAS, propagated_quats, moving_avg, error_vecs);
        memcpy(errors, error_vecs, sizeof(double) * 3 * NUM_SIGMAS);
        /* Weighted average error vector */
        double avg_err[3] = {0,0,0};
        for(int i = 0; i<NUM_SIGMAS; i++){
            for(int j = 0; j < 3; j++){
                avg_err[j] += safe_w[i] * error_vecs[i * 3 + j];
            }
        }
        /* Convergence check: |ē| < 1 μrad */
        if(sqrt(avg_err[0]*avg_err[0] + avg_err[1]*avg_err[1] + avg_err[2]*avg_err[2]) < 1e-6){
            break;
        }
        /* Left-compose mean update: mean ← δq(ē) ⊗ mean */
        double avg_error_quat[4];
        rotationvec2quat(avg_err, avg_error_quat);
        double new_moving_avg[4];
        quat_multiply(avg_error_quat, moving_avg, new_moving_avg);
        memcpy(moving_avg, new_moving_avg, sizeof(double) * 4);
        quat_norm(moving_avg, moving_avg);
    }
    memcpy(avg_quat, moving_avg, sizeof(double) * 4);
}


/* ── iterate ──────────────────────────────────────────────────────────────────
 * Main MUKF step: propagates the filter by one time step and optionally applies
 * a measurement update.
 *
 * Propagation (always):
 *   1. Draw 2n+1 sigma points from the current error covariance P.
 *   2. Lift each error sigma to a quaternion sigma (error_sigmas_to_quat_sigmas).
 *   3. Propagate each quaternion sigma through the gyro kinematic model.
 *   4. Recover the predicted quaternion mean via gradient descent on SO(3).
 *   5. Build the predicted error covariance P_hat from the sigma spread, then
 *      add process noise Q (Q is added *after* propagation, not before sigma
 *      generation, matching Python timing).
 *
 * Measurement update (when num_vecs > 0):
 *   6. Project each propagated sigma through the measurement model h(σ) to get
 *      expected body-frame reference-vector readings.
 *   7. Compute predicted measurement mean ẑ.
 *   8. Compute cross-covariance P_xz and innovation covariance P_zz.
 *   9. Kalman gain:  K = P_xz · (P_zz + R)^{-1}
 *  10. State update: x̂ = mean_err + K·(z − ẑ)
 *  11. Covariance update: P_new = P_hat − K·(P_zz+R)·Kᵀ
 *  12. Quaternion update (left-compose):  q_new = rot(x̂_att) ⊗ q_mean
 *
 * Parameters:
 *   error_state   — current 6-element error state [dq_xyz, bias_xyz]
 *   quat_state    — current quaternion estimate [w, x, y, z]
 *   cov           — current 6×6 error covariance P
 *   body          — stacked body-frame measurements (3·num_vecs elements)
 *   ref           — stacked ECI reference vectors  (3·num_vecs elements)
 *   num_vecs      — 0 = propagation only, 1 = single-vector, 2 = two-vector
 *   gyro          — gyroscope measurement [rad/s]
 *   Q             — 6×6 process noise covariance
 *   R             — (3·num_vecs)×(3·num_vecs) measurement noise covariance
 *   dt            — time step [s]
 *   new_err_state — output: updated error state
 *   new_quat_state— output: updated quaternion
 *   new_P         — output: updated covariance
 */
void iterate(double* error_state, double* quat_state, double* cov,
             double* body, double* ref, int num_vecs,
             double* gyro, const double* Q, const double* R, double dt,
             double* new_err_state, double* new_quat_state, double* new_P){
    int msmt_size = 3 * num_vecs;

    double P[STATE_SIZE * STATE_SIZE];
    memcpy(P, cov, sizeof(double) * STATE_SIZE * STATE_SIZE);
    ensure_psd(P, STATE_SIZE); /* guarantee Cholesky succeeds before sigma generation */

    /* ── Sigma-point generation ────────────────────────────────────────────── */
    double lambda = calculate_lambda(STATE_SIZE, ALPHA);
    double error_sigmas[STATE_SIZE * NUM_SIGMAS];
    get_sigma_points(lambda, error_state, P, error_sigmas);

    /* Lift error sigma points to full quaternion + bias sigma points */
    double quat_sigmas[(STATE_SIZE + 1) * NUM_SIGMAS];
    error_sigmas_to_quat_sigmas(error_sigmas, quat_state, quat_sigmas);

    /* ── Propagation ───────────────────────────────────────────────────────── */
    double propagated_sigmas[(STATE_SIZE + 1) * NUM_SIGMAS];
    propagate_sigmas(quat_sigmas, gyro, dt, propagated_sigmas);

    /* Weights for mean and covariance reconstruction */
    double mean_weights[NUM_SIGMAS];
    double cov_weights[NUM_SIGMAS];
    get_weights(lambda, STATE_SIZE, ALPHA, BETA, cov_weights, mean_weights);

    /* Extract quaternion components from propagated sigmas for mean computation */
    double quaternions_in_sigmas[4 * NUM_SIGMAS];
    for(int i = 0; i< NUM_SIGMAS; i++){
        for(int j = 0; j < 4; j++){
            quaternions_in_sigmas[i * 4 + j] = propagated_sigmas[i * (STATE_SIZE + 1) + j];
        }
    }

    /* Compute weighted Fréchet mean quaternion via gradient descent on SO(3).
     * Returns the mean quaternion and per-sigma attitude error vectors. */
    double avg_quat[4];
    double propagated_err_vecs[NUM_SIGMAS * 3];
    grad_descent(quaternions_in_sigmas, quat_state, mean_weights, avg_quat, propagated_err_vecs);

    /* Assemble full propagated error matrix [att_errors | bias_components] */
    double propagated_errors[NUM_SIGMAS * STATE_SIZE];
    for(int i = 0; i<NUM_SIGMAS; i++){
        for(int j = 0; j < STATE_SIZE; j++){
            if(j < 3){
                /* Attitude error from gradient descent */
                propagated_errors[i * STATE_SIZE + j] = propagated_err_vecs[i * 3 + j];
            }
            else{
                /* Bias components passed through from propagated sigma */
                propagated_errors[i * STATE_SIZE + j] = propagated_sigmas[i * (STATE_SIZE + 1) + j + 1];
            }
        }
    }

    /* Weighted mean of error state (using unclipped mean_weights) */
    double mean_err[STATE_SIZE];
    for(int col = 0; col < STATE_SIZE; col++){
        double var_mean = 0;
        for(int row = 0; row < NUM_SIGMAS; row++){
            var_mean += mean_weights[row] * propagated_errors[row * STATE_SIZE + col];
        }
        mean_err[col] = var_mean;
    }

    /* Predicted covariance P_hat = Σ w_c_i · (eᵢ − ē)(eᵢ − ē)ᵀ + Q
     * Q is added *after* propagation (Python timing), not before sigma
     * generation, so it does not inflate the sigma-point spread. */
    double P_hat[STATE_SIZE * STATE_SIZE] = {0};
    for(int i = 0; i< NUM_SIGMAS; i++){
        double err[STATE_SIZE] = {0};
        for(int j = 0; j< STATE_SIZE; j++){
            err[j] = propagated_errors[i * STATE_SIZE + j] - mean_err[j];
        }
        /* Outer product: err · errᵀ */
        double errTerr[STATE_SIZE * STATE_SIZE];
        mul(err, err, false, errTerr, STATE_SIZE, 1, STATE_SIZE);
        scale(errTerr, cov_weights[i], STATE_SIZE, STATE_SIZE);
        double new_P_hat[STATE_SIZE * STATE_SIZE];
        add(P_hat, errTerr, new_P_hat, STATE_SIZE, STATE_SIZE, STATE_SIZE);
        memcpy(P_hat, new_P_hat, sizeof(double) * STATE_SIZE * STATE_SIZE);
    }
    /* Add process noise Q */
    double P_hat_Q[STATE_SIZE * STATE_SIZE];
    add(P_hat, (double*)Q, P_hat_Q, STATE_SIZE, STATE_SIZE, STATE_SIZE);
    memcpy(P_hat, P_hat_Q, sizeof(double) * STATE_SIZE * STATE_SIZE);

    /* ── Propagation-only exit (no measurement) ────────────────────────────── */
    if(num_vecs == 0 || body == NULL || ref == NULL){
        ensure_psd(P_hat, STATE_SIZE);
        memcpy(new_err_state, mean_err, sizeof(double) * STATE_SIZE);
        memcpy(new_quat_state, avg_quat, sizeof(double) * 4);
        memcpy(new_P, P_hat, sizeof(double) * STATE_SIZE * STATE_SIZE);
        return;
    }

    /* ── Measurement update ────────────────────────────────────────────────── */

    /* Predicted measurements for each sigma point */
    double sigma_msmts[MAX_MSMT_SIZE * NUM_SIGMAS];
    get_sigma_measurements(propagated_sigmas, ref, num_vecs, sigma_msmts);

    /* Predicted measurement mean ẑ = Σ w_m_i · h(σᵢ) */
    double mean_msmt[MAX_MSMT_SIZE];
    for(int col = 0; col < msmt_size; col++){
        double var_mean = 0;
        for(int row = 0; row < NUM_SIGMAS; row++){
            var_mean += mean_weights[row] * sigma_msmts[row * msmt_size + col];
        }
        mean_msmt[col] = var_mean;
    }

    /* Cross-covariance P_xz = Σ w_c_i · (eᵢ − ē)(h(σᵢ) − ẑ)ᵀ */
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

    /* Innovation covariance P_zz = Σ w_c_i · (h(σᵢ) − ẑ)(h(σᵢ) − ẑ)ᵀ */
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

    /* Total measurement noise covariance: P_vv = P_zz + R */
    double P_vv[MAX_MSMT_SIZE * MAX_MSMT_SIZE];
    add(P_zz, (double*)R, P_vv, msmt_size, msmt_size, msmt_size);

    /* Invert P_vv with the double-precision Gauss-Jordan inverse.
     * (LAPACK's inv() routes through single-precision linsolve internally,
     *  which would corrupt these double arrays.) */
    double P_vv_inv[MAX_MSMT_SIZE * MAX_MSMT_SIZE];
    memcpy(P_vv_inv, P_vv, sizeof(double) * msmt_size * msmt_size);
    double_inv(P_vv_inv, msmt_size);

    /* Kalman gain: K = P_xz · P_vv^{-1} */
    double k[STATE_SIZE * MAX_MSMT_SIZE];
    mul(P_xz, P_vv_inv, false, k, STATE_SIZE, msmt_size, msmt_size);

    /* Innovation (measurement residual): ν = z − ẑ */
    double innovation[MAX_MSMT_SIZE];
    for(int i = 0; i<msmt_size; i++){
        innovation[i] = body[i] - mean_msmt[i];
    }

    /* State correction: K·ν */
    double k_innovation[STATE_SIZE];
    mul(k, innovation, false, k_innovation, STATE_SIZE, msmt_size, 1);

    /* Updated error state: x̂ = ē + K·ν */
    double x_hat[STATE_SIZE];
    add(mean_err, k_innovation, x_hat, STATE_SIZE, 1, 1);

    /* Updated covariance: P_new = P_hat − K·P_vv·Kᵀ */
    double K_Pvv[STATE_SIZE * MAX_MSMT_SIZE];
    mul(k, P_vv, false, K_Pvv, STATE_SIZE, msmt_size, msmt_size);
    double K_T[MAX_MSMT_SIZE * STATE_SIZE];
    memcpy(K_T, k, sizeof(double) * STATE_SIZE * msmt_size);
    tranf(K_T, STATE_SIZE, msmt_size);          /* in-place transpose → Kᵀ */
    double K_Pvv_KT[STATE_SIZE * STATE_SIZE];
    mul(K_Pvv, K_T, false, K_Pvv_KT, STATE_SIZE, msmt_size, STATE_SIZE);
    scale(K_Pvv_KT, -1.0, STATE_SIZE, STATE_SIZE);
    double P_new[STATE_SIZE * STATE_SIZE];
    add(P_hat, K_Pvv_KT, P_new, STATE_SIZE, STATE_SIZE, STATE_SIZE);
    ensure_psd(P_new, STATE_SIZE);

    /* Quaternion update (left-compose, matching Python reference):
     * q_new = rot(x̂_att) ⊗ q_mean
     * The attitude error x̂[0:3] is converted to a rotation quaternion and
     * left-composed onto the propagated mean; the caller then resets
     * x̂[0:3] = 0 so the error state stays near zero. */
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
