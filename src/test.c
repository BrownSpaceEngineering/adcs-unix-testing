#include "include/test.h"
#include "declareFunctions.h"
#include "include/quat.h"
#include "stdlib.h"
#include "include/iterate.h"
#include "include/quest.h"
#include "include/laextension.h"
#include "include/magnetosphere.h"

static double uniform_01(void) {
    return ((double) rand() + 1.0) / ((double) RAND_MAX + 2.0);
}

static double normal_sample(double mean, double stddev) {
    double u1 = uniform_01();
    double u2 = uniform_01();
    double mag = sqrt(-2.0 * log(u1));
    double z0 = mag * cos(2.0 * (double) M_PI * u2);
    return mean + stddev * z0;
}

// put test function definitions here

void test_run_all(void) { 
    //test_quest();
    test_long();
}

void test_quest(void){
    float true_ref_to_body[4] = {0.78355, 0.44793, 0.32448, 0.28304};
    float true_body_to_ref[4];
    quat_inv(true_ref_to_body, true_body_to_ref);
    float ref_1[3] = {1, 2, 3};
    float ref_2[3] = {-4,3, -6};
    float ref[6] = {1,2,3,-4,3,-6};

    float body_1[3];
    quat_apply(true_ref_to_body, ref_1, body_1);
    float body_2[3];
    quat_apply(true_ref_to_body, ref_2, body_2);
    float body[6] = {body_1[0], body_1[1], body_1[2], body_2[0], body_2[1], body_2[2]};

    float guess[4];
    quest(body, ref, 2, guess);

    if (f_eps_close_matrix(guess, true_body_to_ref, 1, 4, 1e-4)) {
        printf("QuEST Test Passed\n");
        print(guess, 1, 4);
        print(true_body_to_ref, 1, 4);
    } else {
        printf("QuEST Test Failed\n");
        print(guess, 1, 4);
        print(true_body_to_ref, 1, 4);
    }

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

void test_iteration_2vec(void){ /* unused – replaced by test_long */ }

void test_iteration_1vec(void){ /* unused - replaced by test_long */ }

void test_long(void){
    // ── simulation parameters (mirrors MUKF_FINAL.py) ─────────────────────────
    double dt = 0.1;
    int update_every = 100;     // measurement update every 100 steps = every 10 s (0.1 Hz)
    int switch_every = 27000;   // toggle 1-vec / 2-vec every 45 min  (45*60/0.1)
    int simulation_time = 108000; // 3 hours of 10 Hz  (180*60/0.1)

    double sigma_gyro = 0.0005;
    double sigma_mag  = 0.001;
    double sigma_sun  = 0.01;
    double true_bias[3] = {0.002, -0.001, 0.0015};

    double Q_filt[36] = {0};
    {
        double att_var = sigma_gyro * dt * sigma_gyro * dt;
        Q_filt[0]  = att_var; Q_filt[7]  = att_var; Q_filt[14] = att_var;
        Q_filt[21] = 1e-10;   Q_filt[28] = 1e-10;   Q_filt[35] = 1e-10;
    }

    double R_2vec[36] = {0};
    R_2vec[0]  = 2e-5; R_2vec[7]  = 2e-5; R_2vec[14] = 2e-5;
    R_2vec[21] = 1e-4; R_2vec[28] = 1e-4; R_2vec[35] = 1e-4;

    double R_1vec[9] = {0};
    R_1vec[0] = 2e-5; R_1vec[4] = 2e-5; R_1vec[8] = 2e-5;

    // ── orbital parameters ──────────────────────────────────────────────────────
    double orbital_rate = 2.0 * (double)M_PI / (92.0 * 60.0);
    double incl_rad = 51.6 * (double)M_PI / 180.0;
    double s_incl = sin(incl_rad), c_incl = cos(incl_rad);
    double r_ISS_m  = 6787000.0;
    int32_t jd_int  = 2460676;
    double  jd_frac = 0.5;

    // Reference vector 2 — fixed in ECI (matches Python's ref_vec_2 = [1,0,0])
    double ref_vec2[3] = {1.0, 0.0, 0.0};

    // ── true initial attitude ────────────────────────────────────────────────────
    double nadir_quat_0[4];
    {
        double bX[3]={0.0,-s_incl,c_incl}, bY[3]={0.0,c_incl,s_incl}, bZ[3]={-1.0,0.0,0.0};
        double Rm[9]={bX[0],bY[0],bZ[0],bX[1],bY[1],bZ[1],bX[2],bY[2],bZ[2]};
        double tr0=Rm[0]+Rm[4]+Rm[8];
        double sv=0.5/sqrt(tr0+1.0);
        nadir_quat_0[0]=0.25/sv; nadir_quat_0[1]=(Rm[7]-Rm[5])*sv;
        nadir_quat_0[2]=(Rm[2]-Rm[6])*sv; nadir_quat_0[3]=(Rm[3]-Rm[1])*sv;
        quat_norm(nadir_quat_0, nadir_quat_0);
    }

    double true_quat[4];
    memcpy(true_quat, nadir_quat_0, 4*sizeof(double));

    srand(42);
    double perturb_vec[3] = {
        normal_sample(0.0, 0.05),
        normal_sample(0.0, 0.05),
        normal_sample(0.0, 0.05)
    };
    double perturb_q[4];
    rotationvec2quat(perturb_vec, perturb_q);
    double current_quat[4];
    quat_multiply(perturb_q, true_quat, current_quat);
    quat_norm(current_quat, current_quat);

    double state[6]  = {0};
    double cov[36];
    memset(cov, 0, sizeof(double) * 36);
    // Initialise full covariance to 0.01 * I (matches Python's P = np.eye(6) * 0.01)
    cov[0] = 0.01; cov[7] = 0.01; cov[14] = 0.01;
    cov[21] = 0.01; cov[28] = 0.01; cov[35] = 0.01;

    int vector_mode = 1;

    for(int i = 0; i < simulation_time; i++){
        if(i > 0 && i % switch_every == 0){
            vector_mode = (vector_mode == 1) ? 2 : 1;
        }

        double theta = (double)i * orbital_rate * dt;
        double ct = cos(theta), st = sin(theta);

        double bX[3] = {0.0,  -s_incl,      c_incl     };
        double bY[3] = {-st,   ct*c_incl,   ct*s_incl  };
        double bZ[3] = {-ct,  -st*c_incl,  -st*s_incl  };

        double Rmat[9] = {bX[0],bY[0],bZ[0], bX[1],bY[1],bZ[1], bX[2],bY[2],bZ[2]};
        double trv = Rmat[0]+Rmat[4]+Rmat[8];
        double nadir_q[4];
        if (trv > 0.0) {
            double sv=0.5/sqrt(trv+1.0);
            nadir_q[0]=0.25/sv; nadir_q[1]=(Rmat[7]-Rmat[5])*sv;
            nadir_q[2]=(Rmat[2]-Rmat[6])*sv; nadir_q[3]=(Rmat[3]-Rmat[1])*sv;
        } else if (Rmat[0]>Rmat[4] && Rmat[0]>Rmat[8]) {
            double sv=2.0*sqrt(1.0+Rmat[0]-Rmat[4]-Rmat[8]);
            nadir_q[0]=(Rmat[7]-Rmat[5])/sv; nadir_q[1]=0.25*sv;
            nadir_q[2]=(Rmat[1]+Rmat[3])/sv; nadir_q[3]=(Rmat[2]+Rmat[6])/sv;
        } else if (Rmat[4] > Rmat[8]) {
            double sv=2.0*sqrt(1.0+Rmat[4]-Rmat[0]-Rmat[8]);
            nadir_q[0]=(Rmat[2]-Rmat[6])/sv; nadir_q[1]=(Rmat[1]+Rmat[3])/sv;
            nadir_q[2]=0.25*sv; nadir_q[3]=(Rmat[5]+Rmat[7])/sv;
        } else {
            double sv=2.0*sqrt(1.0+Rmat[8]-Rmat[0]-Rmat[4]);
            nadir_q[0]=(Rmat[3]-Rmat[1])/sv; nadir_q[1]=(Rmat[2]+Rmat[6])/sv;
            nadir_q[2]=(Rmat[5]+Rmat[7])/sv; nadir_q[3]=0.25*sv;
        }
        quat_norm(nadir_q, nadir_q);

        double q_inv_prev[4], delta_q_step[4], omega_vec[3];
        quat_inv(true_quat, q_inv_prev);
        quat_multiply(q_inv_prev, nadir_q, delta_q_step);
        quat2rotationvec(delta_q_step, omega_vec);
        double true_omega[3] = {omega_vec[0]/dt, omega_vec[1]/dt, omega_vec[2]/dt};
        memcpy(true_quat, nadir_q, 4*sizeof(double));

        double r_ECI_m[3] = {r_ISS_m*ct, r_ISS_m*st*c_incl, r_ISS_m*st*s_incl};
        double B_ECI[3];
        wmm_eci_embedded_v2(r_ECI_m, jd_int, jd_frac, B_ECI);
        double B_mag = l2_norm(B_ECI, 3);
        double ref_vec1[3] = {B_ECI[0]/B_mag, B_ECI[1]/B_mag, B_ECI[2]/B_mag};

        double gyro[3];
        for(int j = 0; j < 3; j++){
            gyro[j] = true_omega[j] + true_bias[j] + normal_sample(0.0, sigma_gyro);
        }

        double new_state[6];
        double new_quat[4];
        double new_cov[36];

        if(i % update_every == 0){
            double ref_to_body[4];
            quat_inv(true_quat, ref_to_body);

            if(vector_mode == 2){
                double b1[3], b2[3];
                quat_apply(ref_to_body, ref_vec1, b1);
                quat_apply(ref_to_body, ref_vec2, b2);
                for(int j = 0; j < 3; j++){
                    b1[j] += normal_sample(0.0, sigma_mag);
                    b2[j] += normal_sample(0.0, sigma_sun);
                }
                double n1 = l2_norm(b1, 3), n2 = l2_norm(b2, 3);
                for(int j = 0; j < 3; j++){ b1[j] /= n1; b2[j] /= n2; }

                double body[6] = {b1[0],b1[1],b1[2], b2[0],b2[1],b2[2]};
                double ref[6]  = {ref_vec1[0],ref_vec1[1],ref_vec1[2],
                                 ref_vec2[0],ref_vec2[1],ref_vec2[2]};
                iterate(state, current_quat, cov, body, ref, 2,
                        gyro, Q_filt, R_2vec, dt, new_state, new_quat, new_cov);
            } else {
                // 1-vector update: magnetometer only (eclipse mode)
                double b1[3];
                quat_apply(ref_to_body, ref_vec1, b1);
                for(int j = 0; j < 3; j++) b1[j] += normal_sample(0.0, sigma_mag);
                double n1 = l2_norm(b1, 3);
                for(int j = 0; j < 3; j++) b1[j] /= n1;

                iterate(state, current_quat, cov, b1, ref_vec1, 1,
                        gyro, Q_filt, R_1vec, dt, new_state, new_quat, new_cov);
            }
        } else {
            iterate(state, current_quat, cov, NULL, NULL, 0,
                    gyro, Q_filt, NULL, dt, new_state, new_quat, new_cov);
        }

        new_state[0] = 0.0;
        new_state[1] = 0.0;
        new_state[2] = 0.0;

        memcpy(state, new_state, sizeof(double) * 6);
        memcpy(cov, new_cov, sizeof(double) * 36);
        memcpy(current_quat, new_quat, sizeof(double) * 4);
        quat_norm(current_quat, current_quat);

        if(i % 1000 == 0){
            double err_q[4];
            quat_diff(true_quat, current_quat, err_q);
            double err_vec[3];
            quat2rotationvec(err_q, err_vec);
            double err_deg = (180.0 / (double)M_PI) * l2_norm(err_vec, 3);
            double time_s = i * dt;
            double be[3] = {state[3]-true_bias[0], state[4]-true_bias[1], state[5]-true_bias[2]};
            printf("%5.0f s | Mode: %dv | att: %6.3f deg"
                   " | bias est [%+7.5f %+7.5f %+7.5f]"
                   " | bias err [%+7.5f %+7.5f %+7.5f] rad/s\n",
                   time_s, vector_mode, err_deg,
                   state[3], state[4], state[5],
                   be[0], be[1], be[2]);
        }
    }
}
