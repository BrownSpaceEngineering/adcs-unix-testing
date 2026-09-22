#include "include/quest.h"
#include "include/laextension.h"
#include "Include/dsp/matrix_functions.h"
#include "Include/dsp/basic_math_functions.h"
#include "arm_math.h"
#include <math.h>
#include <string.h>

float newton_raphson(float proposed_eigen, float a, float b, float c, float d, float sigma){
    return proposed_eigen - (powf(proposed_eigen, 4) - (a + b) * (proposed_eigen * proposed_eigen) - c * proposed_eigen + (a * b + c * sigma - d)) / (4 * pow(proposed_eigen, 3) - 2 * (a+b) * proposed_eigen - c);
}

void quest(float* body, float* ref, int msmt_ct, float* result){
    int iters = 20;
    float32_t weights[msmt_ct];
    for(int i = 0; i<msmt_ct; i++){
        weights[i] = 1.0f/(float32_t)msmt_ct;
    }
    float32_t B[3 * 3] = {0};
    float32_t Z[3] = {0};
    for(int i = 0; i<msmt_ct; i++){
        float32_t body_msmt[3] = {body[i * 3], body[i * 3 + 1], body[i * 3 + 2]};
        float32_t ref_msmt[3] = {ref[i * 3], ref[i * 3 + 1], ref[i * 3 + 2]};

        float32_t outer[3 * 3];
        arm_matrix_instance_f32 body_msmt_mat = {3, 1, body_msmt};
        arm_matrix_instance_f32 ref_msmt_mat = {1, 3, ref_msmt};
        arm_matrix_instance_f32 outer_mat = {3, 3, outer};
        arm_mat_mult_f32(&body_msmt_mat, &ref_msmt_mat, &outer_mat);
        arm_scale_f32(outer, weights[i], outer, 9);

        float32_t c_prod[3];
        cross(body_msmt, ref_msmt, c_prod);
        arm_scale_f32(c_prod, weights[i], c_prod, 3);

        arm_add_f32(B, outer, B, 9);
        arm_add_f32(Z, c_prod, Z, 3);
    }

    float32_t B_T[3*3];
    arm_matrix_instance_f32 B_mat = {3, 3, B};
    arm_matrix_instance_f32 B_T_mat = {3, 3, B_T};
    arm_mat_trans_f32(&B_mat, &B_T_mat);

    float32_t S[3*3];
    arm_add_f32(B, B_T, S, 9);

    // 3x3 determinant, expanded along the first row.
    float32_t delta = S[0] * (S[4] * S[8] - S[5] * S[7])
        - S[1] * (S[3] * S[8] - S[5] * S[6])
        + S[2] * (S[3] * S[7] - S[4] * S[6]);
    // kappa = tr(adj(S)); this form remains valid even if S is singular.
    float32_t kappa = S[0] * S[4] + S[4] * S[8] + S[0] * S[8]
        - S[1] * S[3] - S[2] * S[6] - S[5] * S[7];
    arm_matrix_instance_f32 S_mat = {3, 3, S};
    float32_t sigma = 0.5f * trace(&S_mat);

    float32_t SZ[3];
    arm_matrix_instance_f32 Z_mat = {3, 1, Z};
    arm_matrix_instance_f32 SZ_mat = {3, 1, SZ};
    arm_mat_mult_f32(&S_mat, &Z_mat, &SZ_mat);

    float32_t SSZ[3];
    arm_matrix_instance_f32 SSZ_mat = {3, 1, SSZ};
    arm_mat_mult_f32(&S_mat, &SZ_mat, &SSZ_mat);

    float32_t d;
    arm_dot_prod_f32(Z, SSZ, 3, &d);

    float32_t zsz;
    arm_dot_prod_f32(Z, SZ, 3, &zsz);
    float32_t c = delta + zsz;

    float32_t zz;
    arm_dot_prod_f32(Z, Z, 3, &zz);
    float32_t b = sigma * sigma + zz;

    float32_t a = sigma * sigma - kappa;

    float32_t proposed_eigen = 1;
    for(int i = 0; i<iters; i++){
        proposed_eigen = newton_raphson(proposed_eigen, a, b, c, d, sigma);
    }
    float32_t alpha = proposed_eigen * proposed_eigen - sigma * sigma + kappa;
    float32_t beta = proposed_eigen - sigma;
    float32_t gamma = (proposed_eigen + sigma) * alpha - delta;

    float32_t S_S[3 * 3];
    arm_matrix_instance_f32 S_S_mat = {3, 3, S_S};
    arm_mat_mult_f32(&S_mat, &S_mat, &S_S_mat);

    float32_t b_S[3*3];
    arm_scale_f32(S, beta, b_S, 9);

    float32_t a_I[3*3];
    eye(a_I, 3);
    arm_scale_f32(a_I, alpha, a_I, 9);

    float32_t X_1[3*3];
    arm_add_f32(a_I, b_S, X_1, 9);
    float32_t X_2[3*3];
    arm_add_f32(X_1, S_S, X_2, 9);

    float32_t X[3];
    arm_matrix_instance_f32 X_2_mat = {3, 3, X_2};
    arm_matrix_instance_f32 X_mat = {3, 1, X};
    arm_mat_mult_f32(&X_2_mat, &Z_mat, &X_mat);
    float32_t norm_X = l2_norm(X, 3);

    float32_t res[4] = {gamma, X[0], X[1], X[2]};
    float32_t factor = 1.0f / sqrtf(gamma * gamma + norm_X * norm_X);
    arm_scale_f32(res, factor, res, 4);
    memcpy(result, res, sizeof(float32_t) * 4);
}
