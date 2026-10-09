#include "include/quest.h"
#include "include/laextension.h"
#include "include/quat.h"
#include "arm_math.h"
#include <math.h>
#include <string.h>

#define QUEST_NEWTON_ITERS 20

static float newton_raphson(float proposed_eigen, float a, float b, float c, float d, float sigma) {
    float l2 = proposed_eigen * proposed_eigen;
    float f = l2 * l2 - (a + b) * l2 - c * proposed_eigen + (a * b + c * sigma - d);
    float df = 4.0f * l2 * proposed_eigen - 2.0f * (a + b) * proposed_eigen - c;
    if (fabsf(df) < 1e-20f) {
        return proposed_eigen;
    }
    return proposed_eigen - f / df;
}

// Core QUEST solve. b/r are unit vectors, w sums to 1. q_out is the UNNORMALIZED body->ref
// quaternion [gamma, X] (scalar first); its magnitude shrinks to 0 near a 180 degree rotation.
static void quest_core(const float* b, const float* r, const float* w, int n, float* q_out) {
    float32_t B[3 * 3] = {0};
    float32_t Z[3] = {0};
    for (int i = 0; i < n; i++) {
        const float* bi = &b[3 * i];
        const float* ri = &r[3 * i];
        // Attitude profile matrix B = sum w b r^T
        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) {
                B[row * 3 + col] += w[i] * bi[row] * ri[col];
            }
        }
        float32_t c_prod[3];
        cross(bi, ri, c_prod);
        Z[0] += w[i] * c_prod[0];
        Z[1] += w[i] * c_prod[1];
        Z[2] += w[i] * c_prod[2];
    }

    float32_t S[3 * 3];
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            S[row * 3 + col] = B[row * 3 + col] + B[col * 3 + row];
        }
    }

    // 3x3 determinant, expanded along the first row.
    float32_t delta = S[0] * (S[4] * S[8] - S[5] * S[7]) - S[1] * (S[3] * S[8] - S[5] * S[6])
                      + S[2] * (S[3] * S[7] - S[4] * S[6]);
    // kappa = tr(adj(S)); this form remains valid even if S is singular.
    float32_t kappa = S[0] * S[4] + S[4] * S[8] + S[0] * S[8] - S[1] * S[3] - S[2] * S[6] - S[5] * S[7];
    float32_t sigma = 0.5f * (S[0] + S[4] + S[8]);

    float32_t SZ[3], SSZ[3];
    for (int row = 0; row < 3; row++) {
        SZ[row] = S[row * 3 + 0] * Z[0] + S[row * 3 + 1] * Z[1] + S[row * 3 + 2] * Z[2];
    }
    for (int row = 0; row < 3; row++) {
        SSZ[row] = S[row * 3 + 0] * SZ[0] + S[row * 3 + 1] * SZ[1] + S[row * 3 + 2] * SZ[2];
    }

    float32_t d = dot3(Z, SSZ);
    float32_t c = delta + dot3(Z, SZ);
    float32_t bb = sigma * sigma + dot3(Z, Z);
    float32_t a = sigma * sigma - kappa;

    // Largest eigenvalue is ~ sum of weights (= 1) for consistent measurements
    float32_t proposed_eigen = 1.0f;
    for (int i = 0; i < QUEST_NEWTON_ITERS; i++) {
        proposed_eigen = newton_raphson(proposed_eigen, a, bb, c, d, sigma);
    }
    float32_t alpha = proposed_eigen * proposed_eigen - sigma * sigma + kappa;
    float32_t beta = proposed_eigen - sigma;
    float32_t gamma = (proposed_eigen + sigma) * alpha - delta;

    // X = (alpha I + beta S + S^2) Z
    float32_t X[3];
    for (int row = 0; row < 3; row++) {
        X[row] = alpha * Z[row] + beta * SZ[row] + SSZ[row];
    }

    q_out[0] = gamma;
    q_out[1] = X[0];
    q_out[2] = X[1];
    q_out[3] = X[2];
}

bool quest_weighted(const float* body, const float* ref, const float* weights, int msmt_ct,
                    float* result) {
    static const float identity[4] = {1, 0, 0, 0};
    float b[3 * QUEST_MAX_MSMTS];
    float r[3 * QUEST_MAX_MSMTS];
    float w[QUEST_MAX_MSMTS];
    int n = 0;
    float w_sum = 0.0f;

    if (msmt_ct > QUEST_MAX_MSMTS) {
        msmt_ct = QUEST_MAX_MSMTS;
    }
    for (int i = 0; i < msmt_ct; i++) {
        memcpy(&b[3 * n], &body[3 * i], sizeof(float) * 3);
        memcpy(&r[3 * n], &ref[3 * i], sizeof(float) * 3);
        float wi = (weights == NULL) ? 1.0f : weights[i];
        // QUEST assumes unit vectors; skip anything we can't normalize
        if (!(wi > 0.0f) || !normalize_vec(&b[3 * n], 3) || !normalize_vec(&r[3 * n], 3)) {
            continue;
        }
        w[n] = wi;
        w_sum += wi;
        n++;
    }

    // Need two non-parallel observations for a unique attitude
    bool observable = false;
    for (int i = 0; i < n && !observable; i++) {
        for (int j = i + 1; j < n; j++) {
            float c[3];
            cross(&r[3 * i], &r[3 * j], c);
            if (l2_norm(c, 3) > 1e-4f) {
                observable = true;
                break;
            }
        }
    }
    if (!observable) {
        memcpy(result, identity, sizeof(identity));
        return false;
    }
    for (int i = 0; i < n; i++) {
        w[i] /= w_sum;
    }

    // QUEST is singular for 180 degree rotations (the unnormalized solution shrinks to float
    // noise). Shuster's method of sequential rotations: also solve with the reference frame
    // rotated 180 deg about each coordinate axis and keep the best-conditioned solution
    // (largest unnormalized scalar part). At least one of the four frames has
    // |cos(theta'/2)| >= 1/2 for any attitude. Each solve is a handful of 3x3 operations.
    float q[4] = {1, 0, 0, 0};
    float best_gamma = -1.0f;
    for (int axis = 0; axis <= 3; axis++) {
        float r_rot[3 * QUEST_MAX_MSMTS];
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < 3; k++) {
                // 180 deg about coordinate axis `axis`: keep that component, negate the others
                r_rot[3 * i + k] = (axis == 0 || k == axis - 1) ? r[3 * i + k] : -r[3 * i + k];
            }
        }
        float q_rot_raw[4];
        quest_core(b, r_rot, w, n, q_rot_raw);
        if (fabsf(q_rot_raw[0]) > best_gamma) {
            best_gamma = fabsf(q_rot_raw[0]);
            float q_rot[4];
            quat_norm(q_rot_raw, q_rot);
            if (axis == 0) {
                memcpy(q, q_rot, sizeof(q));
            } else {
                // q_rot = q_axis * q, so q = q_axis^-1 * q_rot
                float q_axis_inv[4] = {0, 0, 0, 0};
                q_axis_inv[axis] = -1.0f;
                quat_multiply(q_axis_inv, q_rot, q);
            }
        }
    }

    if (q[0] < 0.0f) {
        q[0] = -q[0];
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
    }
    if (!all_finite(q, 4)) {
        memcpy(result, identity, sizeof(identity));
        return false;
    }
    quat_norm(q, result);
    return true;
}

bool quest(const float* body, const float* ref, int msmt_ct, float* result) {
    return quest_weighted(body, ref, NULL, msmt_ct, result);
}
