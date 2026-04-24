#include "include/filters.h"
#include "math.h"
void gaussian_smooth_1d(float* src, float* dst, int n, float sigma) {
    int radius = (int)(3.0 * sigma + 0.5);
    if (radius > GAUSSIAN_RADIUS) radius = GAUSSIAN_RADIUS;
    int ksize = 2 * radius + 1;
    float kernel[2 * GAUSSIAN_RADIUS + 1];
    float sum = 0.0;

    for (int i = 0; i < ksize; i++) {
        float x = i - radius;
        kernel[i] = expf(-0.5 * (x / sigma) * (x / sigma));
        sum += kernel[i];
    }
    for (int i = 0; i < ksize; i++)
        kernel[i] /= sum;

    for (int i = 0; i < n; i++) {
        float val = 0.0;
        for (int t = 0; t < ksize; t++) {
            int j = i + t - radius;
            if (j < 0)  j = -j;             // reflect left
            if (j >= n) j = 2*n - 2 - j;    // reflect right
            val += src[j] * kernel[t];
        }
        dst[i] = val;
    }
}
void holtz_double_exp_filter(float* src, float* dst, int n, float alpha, float beta){
    double level = src[0];
    double trend = (n > 1) ? src[1] - src[0] : 0.0;

    dst[0] = level;

    for (int i = 1; i < n; i++) {
        double prev_level = level;
        level = alpha * src[i] + (1.0 - alpha) * (level + trend);
        trend = beta  * (level - prev_level) + (1.0 - beta) * trend;
        dst[i] = level;
    }
}
void exponential_filter(float* src, float* dst, int n, float alpha){
     dst[0] = src[0];

    for (int i = 1; i < n; i++){
        dst[i] = alpha * src[i] + (1.0 - alpha) * dst[i - 1];
    }
}