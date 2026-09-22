#ifndef FILTERS
#define FILTERS
#define GAUSSIAN_RADIUS 5
void gaussian_smooth_1d(float* src, float* dst, int n, float sigma);
void exponential_filter(float* src, float* dst, int n, float alpha);
void holtz_double_exp_filter(float* src, float* dst, int n, float alpha, float beta);
#endif