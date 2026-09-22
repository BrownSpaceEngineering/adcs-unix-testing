#ifndef ITER
#define ITER
#include "arm_math.h"
#include <stdbool.h>
#define STATE_SIZE 6
#define MSMT_SIZE 6
#define NUM_SIGMAS (STATE_SIZE * 2 + 1)
/* One UKF step. All arrays contain float32_t values.
 * error_state: [attitude rotation error (3), gyro bias (3)]
 * attitude_quaternion: scalar-first, body coordinates to reference coordinates
 * covariance and process_noise: row-major 6x6 matrices
 * body_vectors and reference_vectors: one or two stacked 3D vectors
 * measurement_noise: row-major 3x3 or 6x6 matrix
 *
 * Every call propagates. A call with do_update=false may pass NULL for the
 * vectors and measurement_noise. The caller resets the first three output
 * error components before the next call because they were injected into the
 * output quaternion; the bias components persist. On failure, no output array
 * is modified. */
arm_status iterate(const float32_t *error_state, const float32_t *attitude_quaternion,
                   const float32_t *covariance, const float32_t *body_vectors,
                   const float32_t *reference_vectors, const float32_t *gyro_measurement,
                   const float32_t *process_noise, const float32_t *measurement_noise,
                   float32_t dt, int vector_count, bool do_update,
                   float32_t *new_error_state, float32_t *new_attitude_quaternion,
                   float32_t *new_covariance);
#endif
