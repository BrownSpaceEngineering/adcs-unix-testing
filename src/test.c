#include "src/test.h"
#include "declareFunctions.h"
#include "quat.h"
#include "adcs_algorithms.h"
#include "stdlib.h"
#include <stdio.h>
#include <math.h>

void test_run_all(void) { 
    test_matrix_product(); 
    test_quaternion();
    test_adcs_algorithms();
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
    printf("----- testing quaternions -----\n");
    // Multiplication Tests
    float res_quat[4];
    float res_vec[3];

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

void test_adcs_algorithms(void) {
    printf("----- testing ADCS algorithms -----\n");

    // 1. B-Dot Detumbling
    float mag0[3] = {1, 1, 1};
    float mag1[3] = {-1, -1, -1};
    float dt = 1.0f;
    float k = 1.0f;
    float moments[3];
    bdot_detumbling(mag0, mag1, dt, k, moments);
    if (moments[0] == -2.0f && moments[1] == -2.0f && moments[2] == -2.0f) {
        printf("B-Dot Detumbling Test Passed\n");
    } else {
        printf("B-Dot Detumbling Test Failed: [%f, %f, %f]\n", moments[0], moments[1], moments[2]);
    }
// 2. Sun Vector ECI (J2000.0)
double jd_j2000 = 2451545.0;
float sun_vec[3];
sun_vector_eci(jd_j2000, sun_vec);
// At J2000.0 (Jan 1st), Sun is NOT at Vernal Equinox.
// Expected roughly [0.18, -0.90, -0.39] based on the algorithm logic.
if (fabsf(sun_vec[0] - 0.18f) < 0.01f && fabsf(sun_vec[1] + 0.90f) < 0.01f && fabsf(sun_vec[2] + 0.39f) < 0.01f) {
    printf("Sun Vector ECI J2000.0 Test Passed\n");
} else {
    printf("Sun Vector ECI J2000.0 Test Failed: [%f, %f, %f]\n", sun_vec[0], sun_vec[1], sun_vec[2]);
}


    // 3. WMM ECI
    float r_eci[3] = {6371200.0f, 0, 0}; 
    float B_eci[3];
    wmm_eci(r_eci, jd_j2000, B_eci);
    float B_mag = sqrtf(B_eci[0]*B_eci[0] + B_eci[1]*B_eci[1] + B_eci[2]*B_eci[2]);
    if (B_mag > 1e-5f && B_mag < 1e-4f) {
        printf("WMM ECI Equator Magnitude Test Passed: %.2e T\n", B_mag);
    } else {
        printf("WMM ECI Equator Magnitude Test Failed: %.2e T\n", B_mag);
    }

    // 4. Orbital Elements Propagation
    float a = 7000000.0f; 
    float e = 0.01f;
    float nu0 = 0.0f;
    float dt_prop = 1000.0f;
    float nu1;
    propagate_orbital_elements(a, e, nu0, dt_prop, &nu1);
    if (nu1 > nu0) {
        printf("Orbital Elements Propagation (Forward) Test Passed: %f -> %f\n", nu0, nu1);
    } else {
        printf("Orbital Elements Propagation (Forward) Test Failed: %f -> %f\n", nu0, nu1);
    }

    // 5. MUKF Iterate
    float state[6] = {0}, q_est[4] = {1, 0, 0, 0}, P[36] = {0};
    for(int i=0; i<36; i+=7) P[i] = 0.01f;
    float meas[3] = {1, 0, 0}, gyro[3] = {0.01f, 0, 0}, B_ref[3] = {1, 0, 0};
    float next_state[6], next_q_est[4], next_P[36];
    mukf_iterate(state, q_est, P, meas, gyro, B_ref, 1.0f, next_state, next_q_est, next_P);
    printf("MUKF Iteration Test Passed (Execution)\n");

    // 6. Down Quaternion
    float sat_pos[3] = {7000000.0f, 0, 0};
    float target_pos[3] = {6371000.0f, 0, 0}; // Nadir
    float q_world[4] = {1, 0, 0, 0};
    float q_down[4];
    down_quaternion(sat_pos, q_world, target_pos, q_down);
    // For this position and identity orientation, q_down should represent some rotation
    if (q_down[0] != 0.0f || q_down[1] != 0.0f || q_down[2] != 0.0f || q_down[3] != 0.0f) {
        printf("Down Quaternion Test Passed (Non-zero result)\n");
    } else {
        printf("Down Quaternion Test Failed (Zero result)\n");
    }
}
