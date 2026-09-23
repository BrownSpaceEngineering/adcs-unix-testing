#include "include/test.h"
#include "Include/dsp/basic_math_functions.h"
#include "Include/dsp/matrix_functions.h"
#include "arm_math.h"
#include "include/iterate.h"
#include "include/laextension.h"
#include "include/quat.h"
#include "include/quest.h"
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI PI
#endif

/*
 * Draw a value strictly inside (0, 1):
 *
 *   u = (rand() + 1) / (RAND_MAX + 2)
 */
static float32_t uniform_01(void) {
    /* u = (rand() + 1) / (RAND_MAX + 2). */
    return ((float32_t)rand() + 1.0f) / ((float32_t)RAND_MAX + 2.0f);
}

/*
 * Draw one normal sample with the Box-Muller transform:
 *
 *   u1, u2 ~ Uniform(0, 1)
 *   z0     = sqrt(-2 log(u1)) cos(2 pi u2)
 *   sample = mean + stddev z0
 */
static float32_t normal_sample(float32_t mean, float32_t stddev) {
    float32_t u1 = uniform_01();
    float32_t u2 = uniform_01();

    /* z0 = sqrt(-2 log(u1)) cos(2 pi u2). */
    float32_t mag = sqrtf(-2.0f * logf(u1));
    float32_t z0 = mag * cosf(2.0f * (float32_t)M_PI * u2);

    /* sample = mean + stddev z0. */
    return mean + stddev * z0;
}

static void print_matrix(float32_t* A, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%f ", A[i * cols + j]);
        }
        printf("\n");
    }
}

static bool eps_close_matrix(float32_t* A, float32_t* B, int rows, int cols, float32_t eps) {
    for (int i = 0; i < rows * cols; i++) {
        if (fabsf(A[i] - B[i]) > eps) {
            return false;
        }
    }
    return true;
}

// put test function definitions here

void test_run_all(void) {
    test_iteration();
}

void test_quest(void) {
    float true_ref_to_body[4] = {0.78355, 0.44793, 0.32448, 0.28304};
    float true_body_to_ref[4];
    quat_inv(true_ref_to_body, true_body_to_ref);
    float ref_1[3] = {1, 2, 3};
    float ref_2[3] = {-4, 3, -6};
    float ref[6] = {1, 2, 3, -4, 3, -6};

    float body_1[3];
    quat_apply(true_ref_to_body, ref_1, body_1);
    float body_2[3];
    quat_apply(true_ref_to_body, ref_2, body_2);
    float body[6] = {body_1[0], body_1[1], body_1[2], body_2[0], body_2[1], body_2[2]};

    float guess[4];
    quest(body, ref, 2, guess);

    if (eps_close_matrix(guess, true_body_to_ref, 1, 4, 1e-4)) {
        printf("QuEST Test Passed\n");
        print_matrix(guess, 1, 4);
        print_matrix(true_body_to_ref, 1, 4);
    } else {
        printf("QuEST Test Failed\n");
        print_matrix(guess, 1, 4);
        print_matrix(true_body_to_ref, 1, 4);
    }
}
void test_matrix_product(void) {

    printf("----- testing matrix product -----\n");

    // test case for 2*2 matrix product
    float A[4] = {1., 2., 3., 4.};
    float B[4] = {5., 6., 7., 8.};
    float C[4] = {0.};
    float C_expected[4] = {19., 22., 43., 50.};

    arm_matrix_instance_f32 A_mat = {2, 2, A};
    arm_matrix_instance_f32 B_mat = {2, 2, B};
    arm_matrix_instance_f32 C_mat = {2, 2, C};
    arm_mat_mult_f32(&A_mat, &B_mat, &C_mat);
    print_matrix(C, 2, 2);

    if (eps_close_matrix(C, C_expected, 2, 2, FLT_EPSILON)) {
        printf("2 * 2 matrix product test passed!\n");
    } else {
        printf("2 * 2 matrix product test failed!\n");
    }

    // 4*4 identity matrix test
    float identity[4 * 4] = {0.};
    float large_result[4 * 4] = {0.};

    eye(identity, 4);

    arm_matrix_instance_f32 identity_mat = {4, 4, identity};
    arm_matrix_instance_f32 large_result_mat = {4, 4, large_result};
    arm_mat_mult_f32(&identity_mat, &identity_mat, &large_result_mat);

    if (eps_close_matrix(identity, large_result, 4, 4, FLT_EPSILON)) {
        printf("Identity matrix squared test passed!\n");
    } else {
        printf("Identity matrix squared test failed!\n");
    }

    float A_large[4 * 4] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.};
    arm_matrix_instance_f32 A_large_mat = {4, 4, A_large};

    arm_mat_mult_f32(&A_large_mat, &identity_mat, &large_result_mat);

    if (eps_close_matrix(A_large, large_result, 4, 4, FLT_EPSILON)) {
        printf("Identity post-multiplication test passed!\n");
    } else {
        printf("Identity post-multiplication test failed!\n");
    }

    arm_mat_mult_f32(&identity_mat, &A_large_mat, &large_result_mat);

    if (eps_close_matrix(A_large, large_result, 4, 4, FLT_EPSILON)) {
        printf("Identity pre-multiplication test passed!\n");
    } else {
        printf("Identity pre-multiplication test failed!\n");
    }

    float B_large[4 * 4] = {5.24829, 6.21496, 3.27374, 3.49223, 1.52040, 3.70849, 7.21884, 0.41667,
                            7.77438, 8.24807, 8.63347, 2.01096, 8.29170, 1.46735, 8.53606, 5.14221};
    arm_matrix_instance_f32 B_large_mat = {4, 4, B_large};

    float C_large[4 * 4] = {1.92304, 0.14043, 1.64762, 4.97396, 0.68077, 4.99275, 7.04041, 2.44857,
                            8.22049, 1.66745, 7.94150, 0.56302, 2.68638, 7.59450, 1.43236, 3.59834};

    float large_multiplication_expected[4 * 4]
        = {50.6168, 63.7473, 83.4036, 55.732,  65.9102, 33.9305, 86.5396, 22.2066,
           96.939,  71.9404, 142.322, 70.9624, 100.929, 61.7765, 99.1469, 68.1449};

    arm_matrix_instance_f32 C_large_mat = {4, 4, C_large};
    arm_mat_mult_f32(&B_large_mat, &C_large_mat, &large_result_mat);

    if (eps_close_matrix(large_result, large_multiplication_expected, 4, 4, FLT_EPSILON)) {
        printf("Large matrix product test passed!\n");
    } else {
        printf("Large matrix product test failed!\n");
    }
}

void test_quaternion(void) {
    // Multiplication Tests
    float res_quat[4];
    float res_vec[3];

    // Identity
    float id_quat[4] = {1, 0, 0, 0};
    quat_multiply(id_quat, id_quat, res_quat);
    float expected_res_id[4] = {1, 0, 0, 0};
    if (eps_close_matrix(res_quat, expected_res_id, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Multiplication Identity Test Passed\n");
    } else {
        printf("Quaternion Multiplication Identity Test Failed\n");
    }

    // Random Pair
    float rp_quat_1[4] = {0.59513523, 0.53820488, 0.00657071, 0.5967465};
    float rp_quat_2[4] = {0.78805728, 0.07242287, 0.48362778, 0.37393157};
    quat_multiply(rp_quat_1, rp_quat_2, res_quat);
    float expected_res_rp[4] = {0.203702175001, 0.181091486099, 0.134968324367, 0.952625236179};
    if (eps_close_matrix(res_quat, expected_res_rp, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Multiplication Random Pair Test Passed\n");
    } else {
        printf("Quaternion Multiplication Random Pair Test Failed\n");
    }

    // Reverse Random Pair
    quat_multiply(rp_quat_2, rp_quat_1, res_quat);
    float expected_res_rp_inv[4] = {0.203702175001, 0.753383864322, 0.451035727503, 0.432995312932};
    if (eps_close_matrix(res_quat, expected_res_rp_inv, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Multiplication Reverse Random Pair Test Passed\n");
    } else {
        printf("Quaternion Multiplication Reverse Random Pair Test Failed\n");
    }

    // Test normalization
    float unnorm_quat[4] = {0.44245429, 0.97492992, 0.24490942, 0.74845996};
    float expected_norm_quat[4] = {0.33290518, 0.73354294, 0.18427127, 0.56314562};
    quat_norm(unnorm_quat, res_quat);
    if (eps_close_matrix(res_quat, expected_norm_quat, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Normalization Test Passed\n");
    } else {
        printf("Quaternion Normalization Test Failed\n");
    }

    // Test inversion
    float quat_to_invert[4] = {7, 4, 5, 9};
    float expected_inverted_quat[4] = {0.0409357, -0.0233918, -0.0292398, -0.0526316};
    quat_inv(quat_to_invert, res_quat);
    if (eps_close_matrix(res_quat, expected_inverted_quat, 1, 4, FLT_EPSILON)) {
        printf("Quaternion Inversion Test Passed\n");
    } else {
        printf("Quaternion Inversion Test Failed\n");
    }

    // Test quat2rotvec using one of the random pair quaternions
    float expected_rp_1_rotvec[3] = {1.2502, 0.0153, 1.3862};
    quat2rotationvec(rp_quat_1, res_vec);
    if (eps_close_matrix(res_vec, expected_rp_1_rotvec, 1, 3, 1e-4)) {
        printf("Quaternion2RotVec Test Passed\n");
    } else {
        printf("Quaternion2RotVec Test Failed\n");
    }

    // Test quat_diff using the random pair quats from before
    float expected_quat_diff[4] = {0.7343, -0.66718, 0.12461, 0.012084};
    quat_diff(rp_quat_1, rp_quat_2, res_quat);
    if (eps_close_matrix(res_quat, expected_quat_diff, 1, 4, 1e-4)) {
        printf("Quaternion Diff Test Passed\n");
    } else {
        printf("Quaternion Diff Test Failed\n");
    }

    // Test the robustness of quatdiff
    float robustness_test[4];
    quat_multiply(res_quat, rp_quat_1, robustness_test);
    if (eps_close_matrix(robustness_test, rp_quat_2, 1, 4, 1e-4)) {
        printf("Quaternion Robustness Test Passed\n");
    } else {
        printf("Quaternion Robustness Test Failed\n");
    }

    // Test Rotation Vec --> Quat
    float random_vec[3] = {0.96667295, 0.7002543, 0.61082435};
    float expected_quat_conv[4] = {0.78355, 0.44793, 0.32448, 0.28304};
    rotationvec2quat(random_vec, res_quat);
    if (eps_close_matrix(res_quat, expected_quat_conv, 1, 4, 1e-4)) {
        printf("RotVec2Quat Test Passed\n");
    } else {
        printf("RotVec2Quat Test Failed\n");
    }

    // Test quat apply
    float expected_rotated_vec[3] = {0.1828, 0.1028, 1.3244};
    quat_apply(rp_quat_1, random_vec, res_vec);
    if (eps_close_matrix(res_vec, expected_rotated_vec, 1, 3, 1e-4)) {
        printf("Quat Apply Test Passed\n");
    } else {
        printf("Quat Apply Test Failed\n");
        printf("Got: %f, %f, %f\n", res_vec[0], res_vec[1], res_vec[2]);
    }
}


#define UKF_CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: %s failed\n", __func__, __LINE__, #condition); \
        exit(EXIT_FAILURE); \
    } \
} while (0)

/*
 * Initialize the simulation covariances:
 *
 *   P0 = diag([0.01, 0.01, 0.01, 2.5e-5, 2.5e-5, 2.5e-5])
 *   Q  = diag([4e-6, 4e-6, 4e-6, 1e-10, 1e-10, 1e-10])
 *   R1 = 2e-5 I3
 *   R2 = diag([2e-5, 2e-5, 2e-5, 1e-4, 1e-4, 1e-4])
 */
static void initialize_filter_matrices(float32_t *covariance, float32_t *process_noise,
                                       float32_t *two_vector_noise,
                                       float32_t *one_vector_noise) {
    memset(covariance, 0, 36 * sizeof(float32_t));
    memset(process_noise, 0, 36 * sizeof(float32_t));
    memset(two_vector_noise, 0, 36 * sizeof(float32_t));
    memset(one_vector_noise, 0, 9 * sizeof(float32_t));
    /*
     * P0 = diag([0.01, 0.01, 0.01, 2.5e-5, 2.5e-5, 2.5e-5])
     * Q  = diag([4e-6, 4e-6, 4e-6, 1e-10, 1e-10, 1e-10])
     * R1 = 2e-5 I3
     * R2 = diag([2e-5, 2e-5, 2e-5, 1e-4, 1e-4, 1e-4])
     */
    for (int i = 0; i < 3; ++i) {
        covariance[i * 6 + i] = 0.01f;
        covariance[(i + 3) * 6 + i + 3] = 2.5e-5f;
        process_noise[i * 6 + i] = 4e-6f;
        process_noise[(i + 3) * 6 + i + 3] = 1e-10f;
        two_vector_noise[i * 6 + i] = one_vector_noise[i * 3 + i] = 2e-5f;
        two_vector_noise[(i + 3) * 6 + i + 3] = 1e-4f;
    }
}

/*
 * Normalize a three-dimensional vector:
 *
 *   output = input / ||input||_2
 */
static void unit_vector3(const float32_t *input, float32_t *output) {
    /* output = input / ||input||_2. */
    float32_t length = sqrtf(input[0] * input[0] + input[1] * input[1]
                         + input[2] * input[2]);
    UKF_CHECK(isfinite(length) && length > 0.0f);
    for (int i = 0; i < 3; ++i) output[i] = input[i] / length;
}

/*
 * Simulate noisy unit-vector observations in body coordinates:
 *
 *   q_reference_to_body = q_body_to_reference^-1
 *   body_mag = unit(q^-1 (*) magnetic_field (*) q + N(0, 0.01^2 I))
 *   body_sun = unit(q^-1 (*) sun_field      (*) q + N(0, 0.01^2 I))
 *
 * The sun equation is evaluated only in two-vector mode.
 */
static void simulate_body_vectors(const float32_t *truth_quaternion,
                                  const float32_t *magnetic_field,
                                  const float32_t *sun_field,
                                  bool two_vectors,
                                  float32_t *body_vectors) {
    float32_t reference_to_body_quaternion[4], noisy_vector[3];

    /* q_reference_to_body = q_body_to_reference^-1. */
    quat_inv(truth_quaternion, reference_to_body_quaternion);

    /* body_mag = unit(q^-1 (*) magnetic_field (*) q + N(0, 0.01^2 I)). */
    quat_apply(reference_to_body_quaternion, magnetic_field, noisy_vector);
    for (int axis = 0; axis < 3; ++axis) {
        noisy_vector[axis] += normal_sample(0.0f, 0.01f);
    }
    unit_vector3(noisy_vector, body_vectors);

    if (two_vectors) {
        /* body_sun = unit(q^-1 (*) sun_field (*) q + N(0, 0.01^2 I)). */
        quat_apply(reference_to_body_quaternion, sun_field, noisy_vector);
        for (int axis = 0; axis < 3; ++axis) {
            noisy_vector[axis] += normal_sample(0.0f, 0.01f);
        }
        unit_vector3(noisy_vector, body_vectors + 3);
    }
}

/*
 * Return the shortest angular distance between two attitudes:
 *
 *   difference = normalize(left (*) right^-1)
 *   angle      = 2 atan2(||difference.xyz||_2, |difference.w|)
 */
static float32_t quaternion_angle(const float32_t *left, const float32_t *right) {
    float32_t inverse[4], difference[4];

    /* difference = left (*) right^-1. */
    quat_inv(right, inverse);
    quat_multiply(left, inverse, difference);
    quat_norm(difference, difference);
    /* angle = 2 atan2(||difference.xyz||_2, |difference.w|). */
    float32_t vector_norm = sqrtf(difference[1] * difference[1]
                              + difference[2] * difference[2]
                              + difference[3] * difference[3]);
    return 2.0f * atan2f(vector_norm, fabsf(difference[0]));
}

static int compare_float32(const void *left, const void *right) {
    float32_t left_value = *(const float32_t *)left;
    float32_t right_value = *(const float32_t *)right;
    return (left_value > right_value) - (left_value < right_value);
}

/* Return a linearly interpolated percentile from a sorted array. */
static float32_t percentile(const float32_t *sorted_values, int count, float32_t fraction) {
    /* position = fraction * (count - 1). */
    float32_t position = fraction * (count - 1);
    int lower = (int)floorf(position);
    int upper = lower + 1 < count ? lower + 1 : lower;
    float32_t weight = position - lower;

    /* result = lower_value + weight * (upper_value - lower_value). */
    return sorted_values[lower]
           + weight * (sorted_values[upper] - sorted_values[lower]);
}

/*
 * Run the fixed-seed UKF simulation. At each step:
 *
 *   truth_rotation = rotvec2quat(truth_rate dt)
 *   q_truth+       = normalize(q_truth (*) truth_rotation)
 *   gyro_measured  = truth_rate + truth_bias + N(0, 0.001^2 I)
 *
 * Every 100 steps, the active body-vector measurements are:
 *
 *   body[v] = unit(q_truth^-1 (*) reference[v] (*) q_truth
 *                  + N(0, 0.01^2 I))
 *
 * After the first 1000 warmup steps, store every post-step attitude error in
 * degrees, then report each mode's:
 *
 *   average = sum(k, e[k]) / count
 *   p50     = percentile(e, 0.50)
 *   p99.5   = percentile(e, 0.995)
 *   maximum = max(k, e[k])
 */
void test_iteration(void) {
    enum {
        MODE_STEPS = 27000,
        MEASUREMENT_PERIOD_STEPS = 100,
        CYCLES = 10
    };
    const float32_t dt = 0.1f;
    const float32_t radians_to_degrees = 180.0f / (float32_t)M_PI;
    const int samples_per_mode = MODE_STEPS * CYCLES;
    float32_t covariance[36], process_noise[36], two_vector_noise[36], one_vector_noise[9];
    float32_t body_vectors[6], reference_vectors[6];
    float32_t next_error_state[6], next_quaternion[4], next_covariance[36];
    float32_t *one_vector_errors = malloc((size_t)samples_per_mode * sizeof(float32_t));
    float32_t *two_vector_errors = malloc((size_t)samples_per_mode * sizeof(float32_t));
    UKF_CHECK(one_vector_errors != NULL && two_vector_errors != NULL);

    /* MATLAB's truth bias and sensor-noise model, with a fixed seed. */
    const float32_t truth_bias[3] = {0.002f, -0.001f, 0.0015f};
    const float32_t truth_rate[3] = {0.018f, -0.011f, 0.014f};
    const float32_t magnetic_field[3] = {22.0f, -12.0f, 34.0f}; /* microtesla */
    const float32_t sun_field[3] = {0.2f, 0.8f, 0.4f};
    unit_vector3(magnetic_field, reference_vectors);
    unit_vector3(sun_field, reference_vectors + 3);
    float32_t truth_quaternion[4] = {1, 0, 0, 0};
    float32_t estimated_quaternion[4] = {1, 0, 0, 0};
    float32_t error_state[6] = {0};
    int two_vector_samples = 0, one_vector_samples = 0;
    float32_t two_vector_sum = 0.0f, one_vector_sum = 0.0f;
    float32_t two_vector_maximum = 0.0f, one_vector_maximum = 0.0f;
    initialize_filter_matrices(covariance, process_noise, two_vector_noise, one_vector_noise);
    srand(42);
    for (int step = 0; step < 2 * MODE_STEPS * CYCLES; ++step) {
        bool two_vectors = (step % (2 * MODE_STEPS)) < MODE_STEPS;
        const int vector_count = two_vectors ? 2 : 1;

        float32_t truth_rotation_step[3], truth_rotation_quaternion[4];
        float32_t next_truth_quaternion[4], measured_gyro[3];

        /*
         * truth_rotation_step = truth_rate * dt
         * measured_gyro = truth_rate + truth_bias + N(0, 0.001^2 I)
         */
        for (int axis = 0; axis < 3; ++axis) {
            truth_rotation_step[axis] = truth_rate[axis] * dt;
            measured_gyro[axis] = truth_rate[axis] + truth_bias[axis]
                                  + normal_sample(0.0f, 0.001f);
        }
        /* truth_rotation_quaternion = rotvec2quat(truth_rate * dt). */
        rotationvec2quat(truth_rotation_step, truth_rotation_quaternion);

        /* next_truth_quaternion = q_truth (*) truth_rotation_quaternion. */
        quat_multiply(truth_quaternion, truth_rotation_quaternion, next_truth_quaternion);

        /* q_truth_next = next_truth_quaternion / ||next_truth_quaternion||_2. */
        quat_norm(next_truth_quaternion, truth_quaternion);

        bool update = step % MEASUREMENT_PERIOD_STEPS == 0;
        if (update) {
            simulate_body_vectors(truth_quaternion, magnetic_field, sun_field,
                                  two_vectors, body_vectors);
        }
        const float32_t *active_noise = two_vectors ? two_vector_noise : one_vector_noise;
        const float32_t *body_input = update ? body_vectors : NULL;
        const float32_t *reference_input = update ? reference_vectors : NULL;
        const float32_t *noise_input = update ? active_noise : NULL;
        arm_status status = iterate(error_state, estimated_quaternion, covariance,
                                    body_input, reference_input, measured_gyro,
                                    process_noise, noise_input, dt, vector_count, update,
                                    next_error_state, next_quaternion, next_covariance);
        if (status != ARM_MATH_SUCCESS) {
            fprintf(stderr, "UKF failed: step=%d mode=%s update=%d status=%d\n",
                    step, two_vectors ? "2vec" : "1vec", update, status);
            exit(EXIT_FAILURE);
        }
        memcpy(error_state, next_error_state, sizeof(error_state));

        /* dtheta was injected into next_quaternion, so retain only the bias. */
        memset(error_state, 0, 3 * sizeof(float32_t));
        memcpy(estimated_quaternion, next_quaternion, sizeof(estimated_quaternion));
        memcpy(covariance, next_covariance, sizeof(covariance));
        /* Include prediction steps as well as measurement updates in each mode. */
        float32_t error_degrees = quaternion_angle(estimated_quaternion, truth_quaternion)
                                  * radians_to_degrees;
        UKF_CHECK(isfinite(error_degrees));
        if (two_vectors) {
            UKF_CHECK(two_vector_samples < samples_per_mode);
            two_vector_errors[two_vector_samples] = error_degrees;
            two_vector_sum += error_degrees;
            if (error_degrees > two_vector_maximum) two_vector_maximum = error_degrees;
            ++two_vector_samples;
        } else {
            UKF_CHECK(one_vector_samples < samples_per_mode);
            one_vector_errors[one_vector_samples] = error_degrees;
            one_vector_sum += error_degrees;
            if (error_degrees > one_vector_maximum) one_vector_maximum = error_degrees;
            ++one_vector_samples;
        }
    }

    UKF_CHECK(one_vector_samples == samples_per_mode);
    qsort(one_vector_errors, (size_t)one_vector_samples, sizeof(float32_t), compare_float32);
    qsort(two_vector_errors, (size_t)two_vector_samples, sizeof(float32_t), compare_float32);
    float32_t one_vector_p50 = percentile(one_vector_errors, one_vector_samples, 0.50f);
    float32_t one_vector_p995 = percentile(one_vector_errors, one_vector_samples, 0.995f);
    float32_t two_vector_p50 = percentile(two_vector_errors, two_vector_samples, 0.50f);
    float32_t two_vector_p995 = percentile(two_vector_errors, two_vector_samples, 0.995f);

    printf("seeded simulation attitude error (deg):\n");
    printf("  1vec: average %.6f, 50th percentile %.6f, 99.5th percentile %.6f, max %.6f\n",
           one_vector_sum / one_vector_samples, one_vector_p50,
           one_vector_p995, one_vector_maximum);
    printf("  2vec: average %.6f, 50th percentile %.6f, 99.5th percentile %.6f, max %.6f\n",
           two_vector_sum / two_vector_samples, two_vector_p50,
           two_vector_p995, two_vector_maximum);

    free(one_vector_errors);
    free(two_vector_errors);
}
