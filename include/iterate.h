#ifndef ITER
#define ITER
static const int STATE_SIZE = 6;
static const int MSMT_SIZE = 6;
static const int MAX_MSMT_SIZE = 6;
static const int NUM_SIGMAS = STATE_SIZE * 2 + 1;
static const double ALPHA = 0.1;
static const double BETA = 2.0;
// num_vecs: 0 = propagation only, 1 = single-vector update, 2 = two-vector update
// body and ref are 3*num_vecs doubles; R is (3*num_vecs)x(3*num_vecs); pass NULL when num_vecs=0
void iterate(double* error_state, double* quat_state, double* cov,
             double* body, double* ref, int num_vecs,
             double* gyro, const double* Q, const double* R, double dt,
             double* new_err_state, double* new_quat_state, double* new_P);
#endif
