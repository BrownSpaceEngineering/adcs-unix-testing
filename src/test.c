#include "src/test.h"
#include "declareFunctions.h"
#include "quat.h"
#include "magnetosphere.h"
#include "quat.h"
#include "stdlib.h"
// put test function definitions here


void test_matrix_product(void) {

    printf("----- testing matrix product -----\n");

    // test case for 2*2 matrix product
    double A[4] = {1., 2., 3., 4.};
    double B[4] = {5., 6., 7., 8.};
    double C[4] = {0.};
    double C_expected[4] = {19., 22., 43., 50.};

    mul(A, B, false, C, 2, 2, 2);
    debug_matrix(C, 2, 2);

    if (dbl_eps_close_matrix(C, C_expected, 2, 2, DBL_EPSILON)) {
        printf("2 * 2 matrix product test passed!\n");
    } else {
        printf("2 * 2 matrix product test failed!\n");
    }

    // 4*4 identity matrix test
    double identity[4 * 4] = {0.};
    double large_result[4 * 4] = {0.};

    eye(identity, 4, 4);

    mul(identity, identity, false, large_result, 4, 4, 4);

    if (dbl_eps_close_matrix(identity, large_result, 4, 4, DBL_EPSILON)) {
        printf("Identity matrix squared test passed!\n");
    } else {
        printf("Identity matrix squared test failed!\n");
    }

    double A_large[4 * 4] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.};

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

    double B_large[4 * 4]
        = {5.24829, 6.21496, 3.27374, 3.49223, 1.52040, 3.70849, 7.21884, 0.41667,
           7.77438, 8.24807, 8.63347, 2.01096, 8.29170, 1.46735, 8.53606, 5.14221};

    double C_large[4 * 4]
        = {1.92304, 0.14043, 1.64762, 4.97396, 0.68077, 4.99275, 7.04041, 2.44857,
           8.22049, 1.66745, 7.94150, 0.56302, 2.68638, 7.59450, 1.43236, 3.59834};

    double large_multiplication_expected[4 * 4]
        = {50.6168, 63.7473, 83.4036, 55.732,  65.9102, 33.9305, 86.5396, 22.2066,
           96.939,  71.9404, 142.322, 70.9624, 100.929, 61.7765, 99.1469, 68.1449};

    mul(B_large, C_large, false, large_result, 4, 4, 4);

    if (dbl_eps_close_matrix(large_result, large_multiplication_expected, 4, 4, DBL_EPSILON)) {
        printf("Large matrix product test passed!\n");
    } else {
        printf("Large matrix product test failed!\n");
    }
}

void test_quaternion(void){
    //Multiplication Tests
    float res_quat[4];
    float res_vec[3];
    float test_qd[4];

    //Identity
    float id_quat[4] = {1, 0, 0, 0};
    quat_multiply(id_quat, id_quat, res_quat);
    float expected_res_id[4] = {1, 0, 0, 0};
    if (f_eps_close_matrix(res_quat, expected_res_id, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Multiplication Identity Test Passed\n");
    }
    else{
        printf("Quaternion Multiplication Identity Test Failed\n");
    }

    //Random Pair
    float rp_quat_1[4] = {0.59513523, 0.53820488, 0.00657071, 0.5967465};
    float rp_quat_2[4] = {0.78805728, 0.07242287, 0.48362778, 0.37393157};
    quat_multiply(rp_quat_1, rp_quat_2, res_quat);
    float expected_res_rp[4] = {0.203702175001, 0.181091486099, 0.134968324367, 0.952625236179};
    if (f_eps_close_matrix(res_quat, expected_res_rp, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Multiplication Random Pair Test Passed\n");
    }
    else{
        printf("Quaternion Multiplication Random Pair Test Failed\n");
    }

    //Reverse Random Pair
    quat_multiply(rp_quat_2, rp_quat_1, res_quat);
    float expected_res_rp_inv[4] = {0.203702175001, 0.753383864322, 0.451035727503, 0.432995312932};
    if (f_eps_close_matrix(res_quat, expected_res_rp_inv, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Multiplication Reverse Random Pair Test Passed\n");
    }
    else{
        printf("Quaternion Multiplication Reverse Random Pair Test Failed\n");
    }
    
    //Test normalization
    float unnorm_quat[4] = {0.44245429, 0.97492992, 0.24490942, 0.74845996};
    float expected_norm_quat[4] = {0.33290518, 0.73354294, 0.18427127, 0.56314562};
    quat_norm(unnorm_quat, res_quat);
    if (f_eps_close_matrix(res_quat, expected_norm_quat, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Normalization Test Passed\n");
    }
    else{
        printf("Quaternion Normalization Test Failed\n");
    }
    
    //Test inversion
    float quat_to_invert[4] = {7,4,5,9};
    float expected_inverted_quat[4] = {0.0409357,-0.0233918,-0.0292398,-0.0526316};
    quat_inv(quat_to_invert, res_quat);
    if (f_eps_close_matrix(res_quat, expected_inverted_quat, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Inversion Test Passed\n");
    }
    else{
        printf("Quaternion Inversion Test Failed\n");
    }

    //Test quat2rotvec using one of the random pair quaternions
    float expected_rp_1_rotvec[3] = {1.2502, 0.0153, 1.3862};
    quat_to_rotation_vec(rp_quat_1, res_vec);
    if (f_eps_close_matrix(res_vec, expected_rp_1_rotvec, 1, 3, 1e-4)) {
        printf("Quaternion2RotVec Test Passed\n");
    }
    else{
        printf("Quaternion2RotVec Test Failed\n");
    }

    //Test quat_diff using the random pair quats from before
    float expected_quat_diff[4] = {0.7343,  -0.66718,  0.12461, 0.012084};
    quat_diff(rp_quat_1, rp_quat_2, res_quat);
    if (f_eps_close_matrix(res_quat, expected_quat_diff, 1, 4, 1e-4)) {
        printf("Quaternion Diff Test Passed\n");
    }
    else{
        printf("Quaternion Diff Test Failed\n");
    }

    //Test the robustness of quatdiff
    float robustness_test[4];
    quat_multiply(res_quat, rp_quat_1, robustness_test);
    if (f_eps_close_matrix(robustness_test, rp_quat_2, 1, 4, 1e-4)) {
        printf("Quaternion Robustness Test Passed\n");
    }
    else{
        printf("Quaternion Robustness Test Failed\n");
    }

    //Test Rotation Vec --> Quat
    float random_vec[3] = {0.96667295, 0.7002543,  0.61082435};
    float expected_quat_conv[4] = {0.78355, 0.44793, 0.32448, 0.28304};
    rotation_vec_to_quat(random_vec, res_quat);
    if (f_eps_close_matrix(res_quat, expected_quat_conv, 1, 4, 1e-4)) {
        printf("RotVec2Quat Test Passed\n");
    }
    else{
        printf("RotVec2Quat Test Failed\n");
    }

    //Test quat apply
    float expected_rotated_vec[3] = {0.1828, 0.1028, 1.3244};
    quat_apply(rp_quat_1, random_vec, res_vec);
    if (f_eps_close_matrix(res_vec, expected_rotated_vec, 1, 3, 1e-4)) {
        printf("Quat Apply Test Passed\n");
    }
    else{
        printf("Quat Apply Test Failed\n");
        printf("Got: %f, %f, %f\n", res_vec[0], res_vec[1], res_vec[2]);
    }

}


/* =========================================================================
 * TEST CASES (Main Function)
 * ========================================================================= */
typedef struct {
    const char *name;
    float r[3];
    float expected_ned[3];
} magnetosphere_test_t;

void test_magnetosphere(void)
{
    /* High-precision JD used by both builds */
    double jd_d = (double)2461116 + 0.223974;

    /* Array of the 10 orbital test cases with expected Matlab NED vectors */
    magnetosphere_test_t tests[] = {
        {"SAA Center",  { 4.5e6f, -3.7e6f, -3.4e6f}, { 2.0598e+04f,  6.0145e+03f, -1.8521e+04f}},
        {"ISS Max Lat", { 4.2e6f,  0.0f,    5.3e6f}, { 1.2509e+04f, -2.9673e+03f,  4.3991e+04f}},
        {"Mag N Pole",  { 0.5e6f,  0.0f,    6.75e6f}, { 1.2450e+03f, -1.3110e+03f,  4.7655e+04f}},
        {"Mag S Pole",  {-2.0e6f,  2.0e6f, -6.1e6f}, { 6.2109e+03f, -1.3507e+04f, -3.8299e+04f}},
        {"Eq 90 deg E", { 0.0f,    6.77e6f, 0.0f},   { 2.3820e+04f, -2.0337e+02f, -1.1299e+04f}},
        {"Bermuda",     { 2.5e6f, -5.3e6f,  3.4e6f}, { 2.0728e+04f,  4.2250e+03f,  2.6296e+04f}},
        {"Pacific 180", {-6.77e6f, 0.0f,    0.0f},   { 3.3273e+04f, -9.4567e+01f, -9.5107e+03f}},
        {"45N 45E",     { 3.4e6f,  3.4e6f,  4.8e6f}, { 1.8653e+04f, -3.0303e+03f,  3.3657e+04f}},
        {"SSO Polar",   { 1.1e6f,  0.0f,   -6.6e6f}, { 1.3022e+04f,  8.2517e+03f, -3.9577e+04f}},
        {"Vandenberg",  {-2.8e6f, -4.8e6f,  3.8e6f}, { 2.3925e+04f,  2.3729e+02f,  2.5251e+04f}}
    };

    int num_tests = sizeof(tests) / sizeof(tests[0]);

    /* Coefficient arrays — type matches the active implementation */
#ifdef MAG_SPHERE_ORIGINAL
    double g[14][14];
    double h[14][14];
    load_wmm2025(jd_d, g, h);
#elif defined(MAG_SPHERE_UPDATE)
    float g[14][14];
    float h[14][14];
    load_wmm2025(jd_d, g, h);
#endif

    printf("Test Name       | Input Mag    | Output B_ECI [Bx, By, Bz] (nT)\n");
    printf("----------------------------------------------------------------\n");

    for (int i = 0; i < num_tests; ++i) {

        /* Per-iteration working scalars — typed per implementation */
#ifdef MAG_SPHERE_ORIGINAL
        double r_check[3] = { tests[i].r[0], tests[i].r[1], tests[i].r[2] };
        double B_vec[3];
        wmm_eci_embedded(r_check, jd_d, B_vec);

        double theta  = gmst_from_jd(jd_d);
        double cos_t  = cos(theta);
        double sin_t  = sin(theta);
        double r_ecef[3];
        r_ecef[0] =  cos_t * r_check[0] + sin_t * r_check[1];
        r_ecef[1] = -sin_t * r_check[0] + cos_t * r_check[1];
        r_ecef[2] =  r_check[2];
        double lat_rad, lon_rad, alt_m;
        ecef_to_geodetic(r_ecef, &lat_rad, &lon_rad, &alt_m);
        double B_N, B_E, B_D;
        synthesize_mag_field(lat_rad, lon_rad, alt_m, g, h, &B_N, &B_E, &B_D);
#elif defined(MAG_SPHERE_UPDATE)
        float r_check[3] = { tests[i].r[0], tests[i].r[1], tests[i].r[2] };
        float B_vec[3];
        wmm_eci_embedded_v2(r_check, 2461116, 0.223974f, B_vec);

        float theta  = gmst_from_jd(jd_d);
        float cos_t  = cosf(theta);
        float sin_t  = sinf(theta);
        float r_ecef[3];
        r_ecef[0] =  cos_t * r_check[0] + sin_t * r_check[1];
        r_ecef[1] = -sin_t * r_check[0] + cos_t * r_check[1];
        r_ecef[2] =  r_check[2];
        float lat_rad, lon_rad, alt_m;
        ecef_to_geodetic(r_ecef, &lat_rad, &lon_rad, &alt_m);
        float B_N, B_E, B_D;
        synthesize_mag_field(lat_rad, lon_rad, alt_m, g, h, &B_N, &B_E, &B_D);
#endif

        /* Accuracy metric — computed and printed in float for both builds */
        float mat_N = tests[i].expected_ned[0];
        float mat_E = tests[i].expected_ned[1];
        float mat_D = tests[i].expected_ned[2];

        float dot_prod = (float)(mat_N * B_N + mat_E * B_E + mat_D * B_D);
        float norm_mat = sqrtf(mat_N*mat_N + mat_E*mat_E + mat_D*mat_D);
        float norm_our = (float)sqrt(B_N*B_N + B_E*B_E + B_D*B_D);

        float accuracy = dot_prod / (norm_mat * norm_our);
        if (accuracy > 1.0f) accuracy = 1.0f;
        if (accuracy < -1.0f) accuracy = -1.0f;
        float angle_error_deg = acosf(accuracy) * 180.0f / (float)M_PI;

        float r_mag = sqrtf(tests[i].r[0]*tests[i].r[0] + tests[i].r[1]*tests[i].r[1] + tests[i].r[2]*tests[i].r[2]) / 1000.0f;
        float bx_nt = (float)(B_vec[0] * 1.0e9);
        float by_nt = (float)(B_vec[1] * 1.0e9);
        float bz_nt = (float)(B_vec[2] * 1.0e9);
        float b_mag = sqrtf(bx_nt*bx_nt + by_nt*by_nt + bz_nt*bz_nt);

        printf("Matlab NED Vector (Correct magnetic field): [%.4e, %.4e, %.4e]\n", mat_N, mat_E, mat_D);
        printf("Our NED Vector: [%.4e, %.4e, %.4e]\n", (float)B_N, (float)B_E, (float)B_D);
        printf("Our accuracy: %f\n", accuracy);
        printf("Our angular error (degrees): %f\n", angle_error_deg);
        printf("Our Starting ECI position: [%.4e, %.4e, %.4e]\n", tests[i].r[0], tests[i].r[1], tests[i].r[2]);
        printf("%-15s | %6.0f km   | [%8.1f, %8.1f, %8.1f] (Total: %.1f nT)\n",
               tests[i].name, r_mag, bx_nt, by_nt, bz_nt, b_mag);
        printf("----------------------------------------------------------------\n");
    }
}

void test_run_all(void) {
    test_matrix_product(); 
    test_quaternion();
    test_magnetosphere();
}