#include "include/quest.h"
#include "declareFunctions.h"
#include "stdlib.h"
#include "include/laextension.h"
float newton_raphson(float proposed_eigen, float a, float b, float c, float d, float sigma){
    return proposed_eigen - (powf(proposed_eigen, 4) - (a + b) * (proposed_eigen * proposed_eigen) - c * proposed_eigen + (a * b + c * sigma - d)) / (4 * pow(proposed_eigen, 3) - 2 * (a+b) * proposed_eigen - c);
}

void quest(float* body, float* ref, int msmt_ct, float* result){
    int iters = 20;
    float weights[msmt_ct];
    for(int i = 0; i<msmt_ct; i++){
        weights[i] = 1.0f/(float)msmt_ct;
    }
    float B[3 * 3] = {0};
    float Z[3] = {0};
    for(int i = 0; i<msmt_ct; i++){
        float body_msmt[3] = {body[i * 3], body[i * 3 + 1], body[i * 3 + 2]};
        float ref_msmt[3] = {ref[i * 3], ref[i * 3 + 1], ref[i * 3 + 2]};
        float outer[3 * 3];
        mul(body_msmt, ref_msmt, false, outer, 3, 1, 3);
        scale(outer, weights[i], 3, 3);
        float c_prod[3];
        cross(body_msmt, ref_msmt, c_prod);
        scale(c_prod, weights[i], 1, 3);
        float new_B[3 * 3];
        float new_Z[3];
        add(B, outer, new_B, 3, 3, 3);
        add(Z, c_prod, new_Z, 1, 3, 3);
        memcpy(B, new_B, sizeof(float) * 3 * 3);
        memcpy(Z, new_Z, sizeof(float) * 3);
    }

    float B_T[3*3];
    memcpy(B_T, B, sizeof(float) * 3 * 3);
    tranf(B_T, 3, 3);

    float S[3*3];
    add(B, B_T, S, 3, 3, 3);

    float delta = det(S, 3);
    // kappa = tr(adj(S)); this form remains valid even if S is singular.
    float kappa = S[0] * S[4] + S[4] * S[8] + S[0] * S[8]
        - S[1] * S[3] - S[2] * S[6] - S[5] * S[7];
    float sigma = 0.5 * trace(S, 3);

    float Z_T[1 * 3];
    memcpy(Z_T, Z, sizeof(float) * 3);
    
    float SZT[3 * 1];
    mul(S, Z_T, false, SZT, 3, 3, 1);
    float SSZT[3 * 1];
    mul(S, SZT, false, SSZT, 3, 3, 1);
    float ZSSZT[1 * 1];
    mul(Z, SSZT, false, ZSSZT, 1, 3, 1);
    float d = ZSSZT[0];

    float ZSZT[1 * 1];
    mul(Z, SZT, false, ZSZT, 1, 3, 1);
    float c = delta + ZSZT[0];
    
    float ZZT[1*1];
    mul(Z, Z_T, false, ZZT, 1, 3, 1);
    float b = sigma * sigma + ZZT[0];

    float a = sigma * sigma - kappa;

    float proposed_eigen = 1;
    for(int i = 0; i<iters; i++){
        proposed_eigen = newton_raphson(proposed_eigen, a, b, c, d, sigma);
    }
    float alpha = proposed_eigen * proposed_eigen - sigma * sigma + kappa;
    float beta = proposed_eigen - sigma;
    float gamma = (proposed_eigen + sigma) * alpha - delta;
    
    float S_S[3 * 3];
    mul(S, S, false, S_S, 3, 3, 3);

    float b_S[3*3];
    memcpy(b_S, S, sizeof(float) * 3 * 3);
    scale(b_S, beta, 3 ,3);

    float a_I[3*3];
    eye(a_I, 3, 3);
    scale(a_I, alpha, 3, 3);
    
    float X_1[3*3];
    add(a_I, b_S, X_1, 3, 3, 3);
    float X_2[3*3];
    add(X_1, S_S, X_2, 3, 3, 3);

    float X[3];
    mul(X_2, Z, false, X, 3, 3, 1);
    float norm_X = l2_norm(X, 3);

    float res[4] = {gamma, X[0], X[1], X[2]};
    float factor = 1 / sqrtf(gamma * gamma + norm_X * norm_X);
    scale(res, factor, 1, 4);
    memcpy(result, res, sizeof(float) * 4);
}