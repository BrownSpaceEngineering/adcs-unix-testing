#include "include/test.h"
#include "declareFunctions.h"
#include "include/quat.h"
#include "stdlib.h"
#include "include/iterate.h"
// put test function definitions here

void test_run_all(void) { 
    test_iteration();
}

void test_matrix_product(void) {

    printf("----- testing matrix product -----\n");

    // test case for 2*2 matrix product
    float A[4] = {1., 2., 3., 4.};
    float B[4] = {5., 6., 7., 8.};
    float C[4] = {0.};
    float C_expected[4] = {19., 22., 43., 50.};

    mul(A, B, false, C, 2, 2, 2);
    debug_matrix(C, 2, 2);

    if (dbl_eps_close_matrix(C, C_expected, 2, 2, DBL_EPSILON)) {
        printf("2 * 2 matrix product test passed!\n");
    } else {
        printf("2 * 2 matrix product test failed!\n");
    }

    // 4*4 identity matrix test
    float identity[4 * 4] = {0.};
    float large_result[4 * 4] = {0.};

    eye(identity, 4, 4);

    mul(identity, identity, false, large_result, 4, 4, 4);

    if (dbl_eps_close_matrix(identity, large_result, 4, 4, DBL_EPSILON)) {
        printf("Identity matrix squared test passed!\n");
    } else {
        printf("Identity matrix squared test failed!\n");
    }

    float A_large[4 * 4] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.};

    mul(A_large, identity, false, large_result, 4, 4, 4);

    if (dbl_eps_close_matrix(A_large, large_result, 4, 4, DBL_EPSILON)) {
        printf("Identity post-multiplication test passed!\n");
    } else {
        printf("Identity post-multiplication test failed!\n");
    }

    mul(identity, A_large, false, large_result, 4, 4, 4);

    if (dbl_eps_close_matrix(A_large, large_result, 4, 4, DBL_EPSILON)) {
        printf("Identity pre-multiplication test passed!\n");
    } else {
        printf("Identity pre-multiplication test failed!\n");
    }

    float B_large[4 * 4] = {5.24829, 6.21496, 3.27374, 3.49223, 1.52040, 3.70849, 7.21884, 0.41667,
                            7.77438, 8.24807, 8.63347, 2.01096, 8.29170, 1.46735, 8.53606, 5.14221};

    float C_large[4 * 4] = {1.92304, 0.14043, 1.64762, 4.97396, 0.68077, 4.99275, 7.04041, 2.44857,
                            8.22049, 1.66745, 7.94150, 0.56302, 2.68638, 7.59450, 1.43236, 3.59834};

    float large_multiplication_expected[4 * 4]
        = {50.6168, 63.7473, 83.4036, 55.732,  65.9102, 33.9305, 86.5396, 22.2066,
           96.939,  71.9404, 142.322, 70.9624, 100.929, 61.7765, 99.1469, 68.1449};

    mul(B_large, C_large, false, large_result, 4, 4, 4);

    if (dbl_eps_close_matrix(large_result, large_multiplication_expected, 4, 4, DBL_EPSILON)) {
        printf("Large matrix product test passed!\n");
    } else {
        printf("Large matrix product test failed!\n");
    }
}

void test_quaternion(void) {
    // Multiplication Tests
    float res_quat[4];
    float res_vec[3];
    float test_qd[4];

    // Identity
    float id_quat[4] = {1, 0, 0, 0};
    quat_multiply(id_quat, id_quat, res_quat);
    float expected_res_id[4] = {1, 0, 0, 0};
    if (f_eps_close_matrix(res_quat, expected_res_id, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Multiplication Identity Test Passed\n");
    } else {
        printf("Quaternion Multiplication Identity Test Failed\n");
    }

    // Random Pair
    float rp_quat_1[4] = {0.59513523, 0.53820488, 0.00657071, 0.5967465};
    float rp_quat_2[4] = {0.78805728, 0.07242287, 0.48362778, 0.37393157};
    quat_multiply(rp_quat_1, rp_quat_2, res_quat);
    float expected_res_rp[4] = {0.203702175001, 0.181091486099, 0.134968324367, 0.952625236179};
    if (f_eps_close_matrix(res_quat, expected_res_rp, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Multiplication Random Pair Test Passed\n");
    } else {
        printf("Quaternion Multiplication Random Pair Test Failed\n");
    }

    // Reverse Random Pair
    quat_multiply(rp_quat_2, rp_quat_1, res_quat);
    float expected_res_rp_inv[4] = {0.203702175001, 0.753383864322, 0.451035727503, 0.432995312932};
    if (f_eps_close_matrix(res_quat, expected_res_rp_inv, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Multiplication Reverse Random Pair Test Passed\n");
    } else {
        printf("Quaternion Multiplication Reverse Random Pair Test Failed\n");
    }

    // Test normalization
    float unnorm_quat[4] = {0.44245429, 0.97492992, 0.24490942, 0.74845996};
    float expected_norm_quat[4] = {0.33290518, 0.73354294, 0.18427127, 0.56314562};
    quat_norm(unnorm_quat, res_quat);
    if (f_eps_close_matrix(res_quat, expected_norm_quat, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Normalization Test Passed\n");
    } else {
        printf("Quaternion Normalization Test Failed\n");
    }

    // Test inversion
    float quat_to_invert[4] = {7, 4, 5, 9};
    float expected_inverted_quat[4] = {0.0409357, -0.0233918, -0.0292398, -0.0526316};
    quat_inv(quat_to_invert, res_quat);
    if (f_eps_close_matrix(res_quat, expected_inverted_quat, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Inversion Test Passed\n");
    } else {
        printf("Quaternion Inversion Test Failed\n");
    }

    // Test quat2rotvec using one of the random pair quaternions
    float expected_rp_1_rotvec[3] = {1.2502, 0.0153, 1.3862};
    quat2rotationvec(rp_quat_1, res_vec);
    if (f_eps_close_matrix(res_vec, expected_rp_1_rotvec, 1, 3, 1e-4)) {
        printf("Quaternion2RotVec Test Passed\n");
    } else {
        printf("Quaternion2RotVec Test Failed\n");
    }

    // Test quat_diff using the random pair quats from before
    float expected_quat_diff[4] = {0.7343, -0.66718, 0.12461, 0.012084};
    quat_diff(rp_quat_1, rp_quat_2, res_quat);
    if (f_eps_close_matrix(res_quat, expected_quat_diff, 1, 4, 1e-4)) {
        printf("Quaternion Diff Test Passed\n");
    } else {
        printf("Quaternion Diff Test Failed\n");
    }

    // Test the robustness of quatdiff
    float robustness_test[4];
    quat_multiply(res_quat, rp_quat_1, robustness_test);
    if (f_eps_close_matrix(robustness_test, rp_quat_2, 1, 4, 1e-4)) {
        printf("Quaternion Robustness Test Passed\n");
    } else {
        printf("Quaternion Robustness Test Failed\n");
    }

    // Test Rotation Vec --> Quat
    float random_vec[3] = {0.96667295, 0.7002543, 0.61082435};
    float expected_quat_conv[4] = {0.78355, 0.44793, 0.32448, 0.28304};
    rotationvec2quat(random_vec, res_quat);
    if (f_eps_close_matrix(res_quat, expected_quat_conv, 1, 4, 1e-4)) {
        printf("RotVec2Quat Test Passed\n");
    } else {
        printf("RotVec2Quat Test Failed\n");
    }

    // Test quat apply
    float expected_rotated_vec[3] = {0.1828, 0.1028, 1.3244};
    quat_apply(rp_quat_1, random_vec, res_vec);
    if (f_eps_close_matrix(res_vec, expected_rotated_vec, 1, 3, 1e-4)) {
        printf("Quat Apply Test Passed\n");
    } else {
        printf("Quat Apply Test Failed\n");
        printf("Got: %f, %f, %f\n", res_vec[0], res_vec[1], res_vec[2]);
    }
}

void test_iteration(void){
    float dt = 0.5;
    float true_body_to_ref[4] = {0.2812, -0.6998, 0.5497, -0.3592};
    float true_ref_to_body[4];
    quat_norm(true_body_to_ref, true_body_to_ref);
    quat_inv(true_body_to_ref, true_ref_to_body);

    float current_guess[4] = {0.2812, -0.6998, 0.5497, -0.3592};
    quat_norm(current_guess, current_guess);

    float true_omega[3] = {2, 0.5, -1};
    float true_bias[3] = {1.0, -0.5, 0.01};
    scale(true_omega, M_PI/180.0f, 1, 3);
    scale(true_bias, M_PI/180.0f, 1, 3);

    float gyro_noise = 0.001;
    float msmt_noise = 0.03;

    float state[6] = {0,0,0,0,0,0};
    float cov[36];
    float Q[36];
    float R[36];
    eye(cov, 6, 6);
    scale(cov, 0.1, 6, 6);
    eye(R, 6, 6);
    scale(R, 0.1, 6, 6);
    eye(Q, 6, 6);
    scale(Q, 0.01, 6, 6);

    float ref[6] = {40, 0, 0, 0, 40, 0};
    float ref_1[3] = {40,0,0};
    float ref_2[3] = {0,40,0};

    for(int iter = 0; iter < 1000; iter++){
        float simulated_gyro_measurement[3];
        for(int i = 0; i<3; i++){
            float uniform_noise = (2.0f * ((float) rand() / (float) RAND_MAX)) - 1.0f;
            simulated_gyro_measurement[i] = true_omega[i] + true_bias[i] + uniform_noise * gyro_noise;
        }
        
        float delta_vec[3];
        memcpy(delta_vec, true_omega, sizeof(float) * 3);
        scale(delta_vec, dt, 1, 3);

        float delta_q[4];
        rotationvec2quat(delta_vec, delta_q);

        float new_true_body_to_ref[4];
        quat_multiply(true_body_to_ref, delta_q, new_true_body_to_ref);
        quat_norm(new_true_body_to_ref, true_body_to_ref);
        float new_ref_to_body[4];
        quat_inv(true_body_to_ref, new_ref_to_body);

        
        float body_1[3];
        float body_2[3];
        quat_apply(new_ref_to_body, ref_1, body_1);
        quat_apply(new_ref_to_body, ref_2, body_2);
        float body[6] = {body_1[0], body_1[1], body_1[2], body_2[0], body_2[1], body_2[2]};
        for(int i = 0; i<6; i++){
            float uniform_noise = (2.0f * ((float) rand() / (float) RAND_MAX)) - 1.0f;
            body[i] = body[i] + uniform_noise * msmt_noise;
        }

        float new_error_state[6];
        float new_quat[4];
        float new_P[36];
        iterate(state, current_guess, cov, body, ref, simulated_gyro_measurement, Q, R, dt, new_error_state, new_quat, new_P);

        new_error_state[0] = 0;
        new_error_state[1] = 0;
        new_error_state[2] = 0;

        memcpy(state, new_error_state, sizeof(float) * 6);
        memcpy(cov, new_P, sizeof(float) * 36);
        memcpy(current_guess, new_quat, sizeof(float) * 4);
        quat_norm(current_guess, current_guess);
        if(iter % 50 == 0){
            float current_q_err[4];
            quat_diff(true_body_to_ref, current_guess, current_q_err);
            float current_err[3];
            quat2rotationvec(current_q_err, current_err);
            printf("Msmt Error: %f\n", (sqrtf(current_err[0]*current_err[0] + current_err[1] * current_err[1] + current_err[2] * current_err[2])));
            printf("Estimated Bias: ");
            float bias[3] = {state[3], state[4], state[5]};
            scale(bias, 180.0f/M_PI, 1, 3);
            print(bias, 1, 3);
            printf("\n");
        }
    }
}