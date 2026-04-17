#ifndef ITER
#define ITER
static const int STATE_SIZE = 6;
static const int MSMT_SIZE = 6;
static const int NUM_SIGMAS = STATE_SIZE * 2 + 1;
static const float ALPHA = 0.01;
static const float BETA = 2;
void iterate(float* error_state, float* quat_state, float* cov, float* body, float* ref, float* gyro, float* Q, float* R, float dt, float* new_err_state, float* new_quat_state, float* new_P);
#endif
