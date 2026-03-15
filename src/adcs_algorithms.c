#include "adcs_algorithms.h"
#include "linalg/LinearAlgebra/declareFunctions.h"
#include "quat.h"
#include "wmm_coeffs.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923
#endif

// Helper for degree-based trig
static float deg2rad(float deg) { return deg * (M_PI / 180.0f); }
static float rad2deg(float rad) { return rad * (180.0f / M_PI); }
static float sind(float deg) { return sinf(deg2rad(deg)); }
static float cosd(float deg) { return cosf(deg2rad(deg)); }

void bdot_detumbling(const float* mag_t0, const float* mag_t1, float dt, float k, float* moments) {
    if (dt <= 0.0f) {
        moments[0] = moments[1] = moments[2] = 0.0f;
        return;
    }
    for (int i = 0; i < 3; i++) {
        float bdot = (mag_t0[i] - mag_t1[i]) / dt;
        moments[i] = -k * bdot;
    }
}

// Helper: cross product
static void cross_product(const float* a, const float* b, float* out) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

// Helper: rotm2quat (WXYZ)
static void rotm2quat(const float* R, float* q) {
    // R is 3x3 in row-major
    float trace = R[0] + R[4] + R[8];
    if (trace > 0) {
        float s = 0.5f / sqrtf(trace + 1.0f);
        q[0] = 0.25f / s;
        q[1] = (R[7] - R[5]) * s;
        q[2] = (R[2] - R[6]) * s;
        q[3] = (R[3] - R[1]) * s;
    } else {
        if (R[0] > R[4] && R[0] > R[8]) {
            float s = 2.0f * sqrtf(1.0f + R[0] - R[4] - R[8]);
            q[0] = (R[7] - R[5]) / s;
            q[1] = 0.25f * s;
            q[2] = (R[1] + R[3]) / s;
            q[3] = (R[2] + R[6]) / s;
        } else if (R[4] > R[8]) {
            float s = 2.0f * sqrtf(1.0f + R[4] - R[0] - R[8]);
            q[0] = (R[2] - R[6]) / s;
            q[1] = (R[1] + R[3]) / s;
            q[2] = 0.25f * s;
            q[3] = (R[5] + R[7]) / s;
        } else {
            float s = 2.0f * sqrtf(1.0f + R[8] - R[0] - R[4]);
            q[0] = (R[3] - R[1]) / s;
            q[1] = (R[2] + R[6]) / s;
            q[2] = (R[5] + R[7]) / s;
            q[3] = 0.25f * s;
        }
    }
    // Normalize
    float mag = sqrtf(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    if (mag > 1e-12f) {
        q[0] /= mag;
        q[1] /= mag;
        q[2] /= mag;
        q[3] /= mag;
    } else {
        q[0] = 1.0f;
        q[1] = q[2] = q[3] = 0.0f;
    }
}

void down_quaternion(const float* sat_pos_eci, const float* q_world, const float* target_pos_eci,
                     float* q_down) {
    float nadir[3];
    for (int i = 0; i < 3; i++) {
        nadir[i] = sat_pos_eci[i] - target_pos_eci[i];
    }
    float n_mag = sqrtf(nadir[0] * nadir[0] + nadir[1] * nadir[1] + nadir[2] * nadir[2]);
    if (n_mag < 1e-6f) {
        q_down[0] = 1.0f;
        q_down[1] = q_down[2] = q_down[3] = 0.0f;
        return;
    }
    for (int i = 0; i < 3; i++)
        nadir[i] /= n_mag;

    float body_y_local[3] = {0.0f, 1.0f, 0.0f};
    float body_y_eci[3];
    quat_apply((float*)q_world, body_y_local, body_y_eci);

    float y_dot_n = body_y_eci[0] * nadir[0] + body_y_eci[1] * nadir[1] + body_y_eci[2] * nadir[2];
    float new_y[3];
    for (int i = 0; i < 3; i++)
        new_y[i] = body_y_eci[i] - y_dot_n * nadir[i];

    float ny_mag = sqrtf(new_y[0] * new_y[0] + new_y[1] * new_y[1] + new_y[2] * new_y[2]);
    if (ny_mag < 1e-6f) {
        float fallback_y[3] = {1.0f, 0.0f, 0.0f};
        if (fabsf(nadir[0]) > 0.9f) {
            fallback_y[0] = 0.0f;
            fallback_y[1] = 1.0f;
        }
        float f_dot_n
            = fallback_y[0] * nadir[0] + fallback_y[1] * nadir[1] + fallback_y[2] * nadir[2];
        for (int i = 0; i < 3; i++)
            new_y[i] = fallback_y[i] - f_dot_n * nadir[i];
        ny_mag = sqrtf(new_y[0] * new_y[0] + new_y[1] * new_y[1] + new_y[2] * new_y[2]);
    }
    for (int i = 0; i < 3; i++)
        new_y[i] /= ny_mag;

    float x_axis[3];
    cross_product(nadir, new_y, x_axis);

    float R[9] = {x_axis[0], new_y[0],  nadir[0], x_axis[1], new_y[1],
                  nadir[1],  x_axis[2], new_y[2], nadir[2]};

    rotm2quat(R, q_down);
}

void sun_vector_eci(double jd, float* sun_vec) {
    double julianOffset = jd - 2451545.0;
    double julianC = julianOffset / 36525.0;

    float meanAnomaly = 357.529f + 35999.050f * (float)julianC;
    float meanLongitude = 280.459f + 36000.770f * (float)julianC;

    meanAnomaly = fmodf(meanAnomaly, 360.0f);
    meanLongitude = fmodf(meanLongitude, 360.0f);

    float sunCenter
        = (1.914602f - 0.004817f * (float)julianC - 0.000014f * (float)(julianC * julianC))
              * sind(meanAnomaly)
          + (0.019993f - 0.000101f * (float)julianC) * sind(2.0f * meanAnomaly)
          + 0.000289f * sind(3.0f * meanAnomaly);

    float eclipticLongitude = meanLongitude + sunCenter;
    eclipticLongitude = fmodf(eclipticLongitude, 360.0f);

    float obliquityEcliptic = 23.0f + 26.0f / 60.0f + 21.448f / 3600.0f
                              - (46.8150f * (float)julianC + 0.00059f * (float)(julianC * julianC)
                                 - 0.001813f * (float)(julianC * julianC * julianC))
                                    / 3600.0f;

    sun_vec[0] = cosd(eclipticLongitude);
    sun_vec[1] = cosd(obliquityEcliptic) * sind(eclipticLongitude);
    sun_vec[2] = sind(obliquityEcliptic) * sind(eclipticLongitude);
}

void orbital_to_eci(float smAxis, float eccentricity, float inclination, float aNodeLongitude,
                    float periapsisArg, float trueAnomaly, float* pos, float* vel) {
    float r0
        = smAxis * (1.0f - eccentricity * eccentricity) / (1.0f + eccentricity * cosf(trueAnomaly));
    float r1[3] = {r0 * cosf(trueAnomaly), r0 * sinf(trueAnomaly), 0.0f};

    float gravConst = 398600.4418f;
    float h = sqrtf(gravConst * smAxis * (1.0f - eccentricity * eccentricity));

    float v1[3] = {-gravConst / h * sinf(trueAnomaly),
                   gravConst / h * (eccentricity + cosf(trueAnomaly)), 0.0f};

    float qPeri[4] = {cosf(periapsisArg / 2.0f), 0.0f, 0.0f, sinf(periapsisArg / 2.0f)};
    float qInc[4] = {cosf(inclination / 2.0f), sinf(inclination / 2.0f), 0.0f, 0.0f};
    float qLong[4] = {cosf(aNodeLongitude / 2.0f), 0.0f, 0.0f, sinf(aNodeLongitude / 2.0f)};

    float qTmp[4], qComb[4];
    quat_multiply(qInc, qPeri, qTmp);
    quat_multiply(qLong, qTmp, qComb);

    quat_apply(qComb, r1, pos);
    quat_apply(qComb, v1, vel);
}

void pd_control_loop(const float* q_error_aa, const float* omega_aa, const float* B,
                     const float* Kp, const float* Kd, float* m) {
    float tau[3];
    for (int i = 0; i < 3; i++) {
        tau[i] = Kp[i] * q_error_aa[i] + Kd[i] * omega_aa[i];
    }

    float B_mag2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];
    if (B_mag2 < 1e-12f) {
        m[0] = m[1] = m[2] = 0.0f;
        return;
    }
    float cross_B_tau[3];
    cross_product(B, tau, cross_B_tau);
    for (int i = 0; i < 3; i++) {
        m[i] = cross_B_tau[i] / B_mag2;
    }
}

// WMM Helper: GMST from Julian Date
static float gmst_from_jd(double JD) {
    double T = (JD - 2451545.0) / 36525.0;
    double GMST_sec = 67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * T + 0.093104 * T * T
                      - 6.2e-6 * T * T * T;
    double GMST_deg = fmod(GMST_sec / 240.0, 360.0);
    if (GMST_deg < 0)
        GMST_deg += 360.0;
    return deg2rad((float)GMST_deg);
}

// WMM Helper: ECEF to Geodetic (WGS-84)
static void ecef2geodetic(const float* r_ecef, float* lat, float* lon, float* alt) {
    float a = 6378137.0f;
    float f = 1.0f / 298.257223563f;
    float b = a * (1.0f - f);
    float e2 = f * (2.0f - f);
    float ep2 = e2 / (1.0f - e2);

    float x = r_ecef[0];
    float y = r_ecef[1];
    float z = r_ecef[2];

    *lon = atan2f(y, x);
    float p = sqrtf(x * x + y * y);

    if (p < 1e-12f) {
        *lon = 0.0f;
        if (z >= 0) {
            *lat = M_PI_2;
            *alt = z - b;
        } else {
            *lat = -M_PI_2;
            *alt = -z - b;
        }
        return;
    }

    float theta = atan2f(z * a, p * b);
    for (int iter = 0; iter < 5; iter++) {
        float st = sinf(theta), ct = cosf(theta);
        *lat = atan2f(z + ep2 * b * st * st * st, p - e2 * a * ct * ct * ct);
        float sl = sinf(*lat);
        float N = a / sqrtf(1.0f - e2 * sl * sl);
        *alt = p / cosf(*lat) - N;
        float theta_new = atan2f(z * (1.0f - e2 * N / (N + *alt)), p);
        if (fabsf(theta_new - theta) < 1e-12f)
            break;
        theta = theta_new;
    }
}

static void synthesize_mag_field(float lat, float lon, float alt, double jd, float* B_ned) {
    float a = 6371200.0f;
    float Re = 6378137.0f;
    float f = 1.0f / 298.257223563f;

    float sin_lat = sinf(lat), cos_lat = cosf(lat);
    float rc = Re / sqrtf(1.0f - (2.0f * f - f * f) * sin_lat * sin_lat);
    float p = (rc + alt) * cos_lat;
    float z = (rc * (1.0f - f) * (1.0f - f) + alt) * sin_lat;
    float r = sqrtf(p * p + z * z);
    float phi_prime = asinf(z / r);
    float psi = lat - phi_prime;
    float theta = M_PI_2 - phi_prime;
    float cos_theta = cosf(theta), sin_theta = sinf(theta);

    float decimalYear = (float)(2000.0 + (jd - 2451545.0) / 365.25);
    float dt = decimalYear - 2025.0f;

    static float P[WMM_MAX_N + 1][WMM_MAX_N + 1];
    static float dP[WMM_MAX_N + 1][WMM_MAX_N + 1];
    memset(P, 0, sizeof(P));
    memset(dP, 0, sizeof(dP));

    P[0][0] = 1.0f;
    P[1][0] = cos_theta;
    P[1][1] = sin_theta;
    dP[0][0] = 0.0f;
    dP[1][0] = -sin_theta;
    dP[1][1] = cos_theta;

    for (int n = 2; n <= WMM_MAX_N; n++) {
        for (int m = 0; m <= n; m++) {
            if (m == n) {
                float factor = sqrtf(1.0f - 1.0f / (2.0f * n));
                P[n][n] = sin_theta * P[n - 1][n - 1] * factor;
                dP[n][n] = (sin_theta * dP[n - 1][n - 1] + cos_theta * P[n - 1][n - 1]) * factor;
            } else {
                float denom = sqrtf((float)(n * n - m * m));
                float K = (float)(2 * n - 1) / denom;
                float M = sqrtf((float)((n - 1) * (n - 1) - m * m)) / denom;
                P[n][m] = K * cos_theta * P[n - 1][m] - M * P[n - 2][m];
                dP[n][m]
                    = K * (cos_theta * dP[n - 1][m] - sin_theta * P[n - 1][m]) - M * dP[n - 2][m];
            }
        }
    }

    float Br = 0, Bt = 0, Bphi = 0;
    for (int i = 0; i < sizeof(wmm2025_data) / sizeof(wmm_coeff_t); i++) {
        int n = wmm2025_data[i].n;
        int m = wmm2025_data[i].m;
        float g = wmm2025_data[i].g + dt * wmm2025_data[i].gdot;
        float h = wmm2025_data[i].h + dt * wmm2025_data[i].hdot;
        float ratio = powf(a / r, (float)(n + 2));
        float cos_mlon = cosf(m * lon), sin_mlon = sinf(m * lon);
        float xy_comp = g * cos_mlon + h * sin_mlon;
        Br += (float)(n + 1) * ratio * xy_comp * P[n][m];
        Bt += ratio * xy_comp * dP[n][m];
        if (m > 0 && fabsf(sin_theta) > 1e-10f) {
            float z_comp = g * sin_mlon - h * cos_mlon;
            Bphi += ratio * (float)m * z_comp * P[n][m] / sin_theta;
        }
    }

    float B_X_geo = -Bt, B_Y_geo = Bphi, B_Z_geo = -Br;
    B_ned[0] = -B_X_geo * cosf(psi) - B_Z_geo * sinf(psi);
    B_ned[1] = B_Y_geo;
    B_ned[2] = B_X_geo * sinf(psi) + B_Z_geo * cosf(psi);
}

void wmm_eci(const float* r_eci, double jd, float* B_eci) {
    float theta = gmst_from_jd(jd);
    float r_ecef[3] = {r_eci[0] * cosf(theta) + r_eci[1] * sinf(theta),
                       -r_eci[0] * sinf(theta) + r_eci[1] * cosf(theta), r_eci[2]};
    float lat, lon, alt;
    ecef2geodetic(r_ecef, &lat, &lon, &alt);
    float B_ned[3];
    synthesize_mag_field(lat, lon, alt, jd, B_ned);
    float B_ecef[3] = {
        -sinf(lat) * cosf(lon) * B_ned[0] - sinf(lon) * B_ned[1] - cosf(lat) * cosf(lon) * B_ned[2],
        -sinf(lat) * sinf(lon) * B_ned[0] + cosf(lon) * B_ned[1] - cosf(lat) * sinf(lon) * B_ned[2],
        cosf(lat) * B_ned[0] - sinf(lat) * B_ned[2]};
    B_eci[0] = B_ecef[0] * cosf(theta) - B_ecef[1] * sinf(theta);
    B_eci[1] = B_ecef[0] * sinf(theta) + B_ecef[1] * cosf(theta);
    B_eci[2] = B_ecef[2];
    for (int i = 0; i < 3; i++)
        B_eci[i] *= 1e-9f;
}

void mukf_iterate(float* state, float* q_est, float* P, const float* measurement,
                  const float* gyro_meas, const float* B_ref, float dt, float* next_state,
                  float* next_q_est, float* next_P) {
    const int n = 6;
    const int num_sigmas = 13;
    float alpha = 1e-3f, beta = 2.0f;
    float kappa = 3.0f - (float)n;
    float lam = alpha * alpha * ((float)n + kappa) - (float)n;

    float Q[36] = {0}; // Process noise
    for (int i = 0; i < 6; i++)
        Q[i * 6 + i] = 1e-5f;
    float R[9] = {0}; // Measurement noise
    for (int i = 0; i < 3; i++)
        R[i * 3 + i] = 1e-3f;

    float P_plus_Q[36];
    for (int i = 0; i < 36; i++)
        P_plus_Q[i] = P[i] + Q[i];
    float L[36];
    chol(P_plus_Q, L, n);

    float sigmas[13][6];
    memcpy(sigmas[0], state, 6 * sizeof(float));
    float factor = sqrtf(lam + (float)n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sigmas[i + 1][j] = state[j] + factor * L[j * n + i];
            sigmas[i + 1 + n][j] = state[j] - factor * L[j * n + i];
        }
    }

    float prop_sigmas_q[13][4];
    float prop_sigmas_bias[13][3];
    for (int i = 0; i < num_sigmas; i++) {
        float err_q[4], bias[3];
        rotationvec2quat(&sigmas[i][0], err_q);
        quat_multiply(err_q, q_est, prop_sigmas_q[i]);
        memcpy(bias, &sigmas[i][3], 3 * sizeof(float));
        float omega[3];
        for (int j = 0; j < 3; j++)
            omega[j] = (gyro_meas[j] - bias[j]) * dt;
        float rot_q[4], q_tmp[4];
        rotationvec2quat(omega, rot_q);
        quat_multiply(prop_sigmas_q[i], rot_q, q_tmp);
        memcpy(prop_sigmas_q[i], q_tmp, 4 * sizeof(float));
        memcpy(prop_sigmas_bias[i], bias, 3 * sizeof(float));
    }

    float w_m[13], w_c[13];
    w_m[0] = lam / (lam + (float)n);
    w_c[0] = w_m[0] + (1.0f - alpha * alpha + beta);
    for (int i = 1; i < num_sigmas; i++)
        w_m[i] = w_c[i] = 0.5f / (lam + (float)n);

    float q_avg[4];
    memcpy(q_avg, prop_sigmas_q[0], 4 * sizeof(float));
    for (int iter = 0; iter < 10; iter++) {
        float err_avg[3] = {0, 0, 0};
        for (int i = 0; i < num_sigmas; i++) {
            float dq[4], ev[3];
            quat_diff(q_avg, prop_sigmas_q[i], dq);
            quat2rotationvec(dq, ev);
            for (int j = 0; j < 3; j++)
                err_avg[j] += w_m[i] * ev[j];
        }
        float dq_avg[4], q_tmp[4];
        rotationvec2quat(err_avg, dq_avg);
        quat_multiply(dq_avg, q_avg, q_tmp);
        memcpy(q_avg, q_tmp, 4 * sizeof(float));
        if (err_avg[0] * err_avg[0] + err_avg[1] * err_avg[1] + err_avg[2] * err_avg[2] < 1e-8f)
            break;
    }
    memcpy(next_q_est, q_avg, 4 * sizeof(float));

    float bias_avg[3] = {0, 0, 0};
    for (int i = 0; i < num_sigmas; i++) {
        for (int j = 0; j < 3; j++)
            bias_avg[j] += w_m[i] * prop_sigmas_bias[i][j];
    }

    float prop_errors[13][6];
    for (int i = 0; i < num_sigmas; i++) {
        float dq[4], ev[3];
        quat_diff(q_avg, prop_sigmas_q[i], dq);
        quat2rotationvec(dq, ev);
        memcpy(&prop_errors[i][0], ev, 3 * sizeof(float));
        for (int j = 0; j < 3; j++)
            prop_errors[i][3 + j] = prop_sigmas_bias[i][j] - bias_avg[j];
    }

    float pred_meas[13][3];
    float meas_avg[3] = {0, 0, 0};
    for (int i = 0; i < num_sigmas; i++) {
        float q_inv[4];
        quat_inv(prop_sigmas_q[i], q_inv);
        quat_apply(q_inv, (float*)B_ref, pred_meas[i]);
        float m_mag = sqrtf(pred_meas[i][0] * pred_meas[i][0] + pred_meas[i][1] * pred_meas[i][1]
                            + pred_meas[i][2] * pred_meas[i][2]);
        if (m_mag > 1e-12f)
            for (int j = 0; j < 3; j++)
                pred_meas[i][j] /= m_mag;
        for (int j = 0; j < 3; j++)
            meas_avg[j] += w_m[i] * pred_meas[i][j];
    }

    float P_hat[36] = {0}, P_zz[9] = {0}, P_xz[18] = {0};
    for (int i = 0; i < num_sigmas; i++) {
        for (int row = 0; row < 6; row++)
            for (int col = 0; col < 6; col++)
                P_hat[row * 6 + col] += w_c[i] * prop_errors[i][row] * prop_errors[i][col];
        float dz[3];
        for (int j = 0; j < 3; j++)
            dz[j] = pred_meas[i][j] - meas_avg[j];
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++)
                P_zz[row * 3 + col] += w_c[i] * dz[row] * dz[col];
        for (int row = 0; row < 6; row++)
            for (int col = 0; col < 3; col++)
                P_xz[row * 3 + col] += w_c[i] * prop_errors[i][row] * dz[col];
    }
    for (int i = 0; i < 9; i++)
        P_zz[i] += R[i];

    float P_zz_inv[9];
    memcpy(P_zz_inv, P_zz, 9 * sizeof(float));
    inv(P_zz_inv, 3);

    float K[18];
    mul(P_xz, P_zz_inv, false, K, 6, 3, 3);

    float innovation[3];
    for (int i = 0; i < 3; i++)
        innovation[i] = measurement[i] - meas_avg[i];
    float dx[6];
    memset(dx, 0, 6 * sizeof(float));
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            dx[i] += K[i * 3 + j] * innovation[j];

    for (int i = 0; i < 6; i++)
        next_state[i] = dx[i]; // Reset state to innovation for next iter?
    // In many UKF for quat, the state is reset to 0 after update, and q_est is updated.
    float dq_final[4], q_final[4];
    rotationvec2quat(dx, dq_final);
    quat_multiply(dq_final, q_avg, q_final);
    memcpy(next_q_est, q_final, 4 * sizeof(float));
    for (int i = 0; i < 3; i++)
        next_state[3 + i] = bias_avg[i] + dx[3 + i];
    memset(next_state, 0, 3 * sizeof(float)); // Reset rotation error

    float K_Pzz[18];
    mul(K, P_zz, false, K_Pzz, 6, 3, 3);
    float K_Pzz_KT[36];
    float KT[18];
    for (int r = 0; r < 6; r++)
        for (int c = 0; c < 3; c++)
            KT[c * 6 + r] = K[r * 3 + c];
    mul(K_Pzz, KT, false, K_Pzz_KT, 6, 3, 6);
    for (int i = 0; i < 36; i++)
        next_P[i] = P_hat[i] - K_Pzz_KT[i];
}

void moment_to_current(const float* m, float Imax, float* I) {
    // Basic scaling for current. G is assumed 1.0 for now.
    float n[3] = {1.0f, 1.0f, 1.0f};
    float A[3] = {1.0f, 1.0f, 1.0f};
    for (int i = 0; i < 3; i++) {
        I[i] = m[i] / (n[i] * A[i]);
        if (I[i] > Imax)
            I[i] = Imax;
        if (I[i] < -Imax)
            I[i] = -Imax;
    }
}

void ecef_to_eci(const float* r_ecef, double unix_seconds, float* r_eci) {
    double JD = unix_seconds / 86400.0 + 2440587.5;
    float theta = gmst_from_jd(JD);
    float cos_t = cosf(theta), sin_t = sinf(theta);
    // R_z(theta)
    r_eci[0] = cos_t * r_ecef[0] - sin_t * r_ecef[1];
    r_eci[1] = sin_t * r_ecef[0] + cos_t * r_ecef[1];
    r_eci[2] = r_ecef[2];
}

void propagate_orbital_elements(float semi_major_axis, float eccentricity, float true_anomaly_deg,
                                float time_delta, float* new_true_anomaly_deg) {
    float u = 3.986004418e14f;
    float a = semi_major_axis;
    float e = eccentricity;
    float f = deg2rad(true_anomaly_deg);

    float E = 2.0f * atan2f(sqrtf(1.0f - e) * sinf(f / 2.0f), sqrtf(1.0f + e) * cosf(f / 2.0f));
    float M = E - e * sinf(E);
    float n = sqrtf(u / (a * a * a));
    float M_new = M + n * time_delta;

    float E_curr = M_new;
    for (int i = 0; i < 10; i++) {
        float E_change = (E_curr - e * sinf(E_curr) - M_new) / (1.0f - e * cosf(E_curr));
        E_curr -= E_change;
        if (fabsf(E_change) < 1e-10f)
            break;
    }

    float f_new
        = 2.0f
          * atan2f(sqrtf(1.0f + e) * sinf(E_curr / 2.0f), sqrtf(1.0f - e) * cosf(E_curr / 2.0f));
    *new_true_anomaly_deg = rad2deg(f_new);
}
