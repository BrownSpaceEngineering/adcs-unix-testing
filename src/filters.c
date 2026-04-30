#include "include/filters.h"
#include "Include/dsp/controller_functions.h"
#include "Include/dsp/matrix_functions.h"
#include "math.h"

/**
 * \fn gaussian_smooth_1d
 * 
 * \brief Applies a 1D Gaussian smoothing filter to the input data.
 * 
 * \param src Pointer to the input data array.
 * \param dst Pointer to the output data array (must be pre-allocated).
 * \param n The number of elements in the input and output arrays.
 * \param sigma The standard deviation of the Gaussian kernel.
 */
void gaussian_smooth_1d(arm_matrix_instance_f32* src, arm_matrix_instance_f32* dst, int n, float32_t sigma) {
    int radius = (int)(3.0 * sigma + 0.5);
    if (radius > GAUSSIAN_RADIUS) { radius = GAUSSIAN_RADIUS; }
    int ksize = 2 * radius + 1;
    float32_t kernel[2 * GAUSSIAN_RADIUS + 1];
    float32_t sum = 0.0;

    float32_t *src_data = src->pData;
    float32_t *dst_data = dst->pData;

    float32_t x; 

    // generate the gaussian kernel 
    for (int i = 0; i < ksize; i++) {
        x = i - radius;
        kernel[i] = expf(-0.5 * (x / sigma) * (x / sigma));
        sum += kernel[i];
    }

    // and now normalise it.
    for (int i = 0; i < ksize; i++) {
        kernel[i] /= sum;
    }

    // apply the kernel to the data with reflection at the boundaries
    for (int i = 0; i < n; i++) {
        float32_t val = 0.0;
        for (int t = 0; t < ksize; t++) {
            int j = i + t - radius;
            if (j < 0)  j = -j;             // reflect left
            if (j >= n) j = 2*n - 2 - j;    // reflect right
            val += src_data[j] * kernel[t];
        }
        dst_data[i] = val;
    }
}

/**
 * \fn holtz_double_exp_filter
 * 
 * \brief Applies a Holt-Winters double exponential smoothing filter to the input data.
 * 
 * \param src Pointer to the input data array.
 * \param dst Pointer to the output data array (must be pre-allocated).
 * \param n The number of elements in the input and output arrays.
 * \param alpha The smoothing parameter for the level.
 * \param beta The smoothing parameter for the trend.
 */
// TODO: Why does this use doubles? Should we make a float version?
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

/** 
 * \fn exponential_filter
 * 
 * \brief Applies an exponential smoothing filter to the input data.
 * 
 * \param src Pointer to the input data array.
 * \param dst Pointer to the output data array (must be pre-allocated).
 * \param n The number of elements in the input and output arrays.
 * \param alpha The smoothing parameter.
 */
void exponential_filter(float* src, float* dst, int n, float alpha){
     dst[0] = src[0];

    for (int i = 1; i < n; i++){
        dst[i] = alpha * src[i] + (1.0 - alpha) * dst[i - 1];
    }
}