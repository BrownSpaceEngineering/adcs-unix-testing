#include "include/quest.h"
#include "declareFunctions.h"
#include "stdlib.h"
#include "include/laextension.h"
double newton_raphson(double proposed_eigen, double a, double b, double c, double d, double sigma){
    return proposed_eigen - (pow(proposed_eigen, 4) - (a + b) * (proposed_eigen * proposed_eigen) - c * proposed_eigen + (a * b + c * sigma - d)) / (4 * pow(proposed_eigen, 3) - 2 * (a+b) * proposed_eigen - c);
}

void quest(double* body, double* ref, int msmt_ct, double* result){
    int iters = 20;
    double weights[msmt_ct];
    for(int i = 0; i<msmt_ct; i++){
        weights[i] = 1.0/(double)msmt_ct;
    }
    double B[3 * 3] = {0};
    double Z[3] = {0};
    for(int i = 0; i<msmt_ct; i++){
        double body_msmt[3] = {body[i * 3], body[i * 3 + 1], body[i * 3 + 2]};
        double ref_msmt[3] = {ref[i * 3], ref[i * 3 + 1], ref[i * 3 + 2]};
        double outer[3 * 3];
        mul(body_msmt, ref_msmt, false, outer, 3, 1, 3);
        scale(outer, weights[i], 3, 3);
        double c_prod[3];
        cross(body_msmt, ref_msmt, c_prod);
        scale(c_prod, weights[i], 1, 3);
        double new_B[3 * 3];
        double new_Z[3];
        add(B, outer, new_B, 3, 3, 3);
        add(Z, c_prod, new_Z, 1, 3, 3);
        memcpy(B, new_B, sizeof(double) * 3 * 3);
        memcpy(Z, new_Z, sizeof(double) * 3);
    }

    double B_T[3*3];
    memcpy(B_T, B, sizeof(double) * 3 * 3);
    tranf(B_T, 3, 3);

    double S[3*3];
    add(B, B_T, S, 3, 3, 3);

    double delta = det(S, 3);
    // kappa = tr(adj(S)); this form remains valid even if S is singular.
    double kappa = S[0] * S[4] + S[4] * S[8] + S[0] * S[8]
        - S[1] * S[3] - S[2] * S[6] - S[5] * S[7];
    double sigma = 0.5 * trace(S, 3);

    double Z_T[1 * 3];
    memcpy(Z_T, Z, sizeof(double) * 3);
    
    double SZT[3 * 1];
    mul(S, Z_T, false, SZT, 3, 3, 1);
    double SSZT[3 * 1];
    mul(S, SZT, false, SSZT, 3, 3, 1);
    double ZSSZT[1 * 1];
    mul(Z, SSZT, false, ZSSZT, 1, 3, 1);
    double d = ZSSZT[0];

    double ZSZT[1 * 1];
    mul(Z, SZT, false, ZSZT, 1, 3, 1);
    double c = delta + ZSZT[0];
    
    double ZZT[1*1];
    mul(Z, Z_T, false, ZZT, 1, 3, 1);
    double b = sigma * sigma + ZZT[0];

    double a = sigma * sigma - kappa;

    double proposed_eigen = 1;
    for(int i = 0; i<iters; i++){
        proposed_eigen = newton_raphson(proposed_eigen, a, b, c, d, sigma);
    }
    double alpha = proposed_eigen * proposed_eigen - sigma * sigma + kappa;
    double beta = proposed_eigen - sigma;
    double gamma = (proposed_eigen + sigma) * alpha - delta;
    
    double S_S[3 * 3];
    mul(S, S, false, S_S, 3, 3, 3);

    double b_S[3*3];
    memcpy(b_S, S, sizeof(double) * 3 * 3);
    scale(b_S, beta, 3 ,3);

    double a_I[3*3];
    eye(a_I, 3, 3);
    scale(a_I, alpha, 3, 3);
    
    double X_1[3*3];
    add(a_I, b_S, X_1, 3, 3, 3);
    double X_2[3*3];
    add(X_1, S_S, X_2, 3, 3, 3);

    double X[3];
    mul(X_2, Z, false, X, 3, 3, 1);
    double norm_X = l2_norm(X, 3);

    double res[4] = {gamma, X[0], X[1], X[2]};
    double factor = 1.0 / sqrt(gamma * gamma + norm_X * norm_X);
    scale(res, factor, 1, 4);
    memcpy(result, res, sizeof(double) * 4);
}
