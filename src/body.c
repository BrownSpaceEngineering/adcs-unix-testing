
#include "string.h"
#include "math.h"
#include "body.h"
#include "laextension.h"
#include "quest.h"
#include "bdot.h"
#include "moments2currents.h"
#include "photodiode_determination.h"
#include "magnetosphere.h"
#include "orbital2eci.h"
#include "sunvec.h"
#include "kepler.h"
#include "iterate.h"
#include "quat.h"
#include "pd.h"
#include "ecef2eci.h"
#include "down_quat.h"
#include "torque2moments.h"
#include "declareFunctions.h"

float PVD_ECEF[3] = {6.3761e6, -0.1387e6, 0.0807e6};
const float DETUMBLING_GYRO_THRESHOLD = 0.17f;
const float DETUMBLING_MAG_THRESHOLD = 0.01f;
float Imax[3] = {1, 1, 1};//TEMP
const float INIT_P[6 * 6] = {0.01, 0, 0, 0, 0, 0,
                                0, 0.01, 0, 0, 0, 0,
                                0, 0, 0.01, 0, 0, 0,
                                0, 0, 0, 0.01, 0, 0,
                                0, 0, 0, 0, 0.01, 0,
                                0, 0, 0, 0, 0, 0.01};
const float Q[6 * 6] =  {0.01, 0, 0, 0, 0, 0,
                                0, 0.01, 0, 0, 0, 0,
                                0, 0, 0.01, 0, 0, 0,
                                0, 0, 0, 0.01, 0, 0,
                                0, 0, 0, 0, 0.01, 0,
                                0, 0, 0, 0, 0, 0.01};
const float R_IN_SUN[6 * 6] =  {0.01, 0, 0, 0, 0, 0,
                                        0, 0.01, 0, 0, 0, 0,
                                        0, 0, 0.01, 0, 0, 0,
                                        0, 0, 0, 0.01, 0, 0,
                                        0, 0, 0, 0, 0.01, 0,
                                        0, 0, 0, 0, 0, 0.01};
const float R_IN_SHADOW[6 * 6] = {0.01, 0, 0, 0, 0, 0,
                                        0, 0.01, 0, 0, 0, 0,
                                        0, 0, 0.01, 0, 0, 0,
                                        0, 0, 0, 1.0, 0, 0,
                                        0, 0, 0, 0, 1.0, 0,
                                        0, 0, 0, 0, 0, 1.0};

bool filter_failure(float* new_quat, float* gyro, float* new_P, float dt) {
    
    //Assess properties of new quat
    float diff[4];
    
    quat_diff(estimated_quat, new_quat, diff);
    
    float diff_rot[3];
    
    quat2rotationvec(diff, diff_rot);
    scale(diff_rot, 1 / dt, 1, 3);
    
    if(l2_norm(diff_rot, 3) / l2_norm(gyro, 3) > 10){
        return true;
    }

    //Assess properties of cov
    for(int i = 0; i < 6 * 6; i++) {
        if (isnan(error_quat_cov[i])) {
            return true;
        }
    }

    return false;
}

bool ready_to_point(float* photodiode_measurements, float* gyro, float* last_magnetometer_measurements, float* magnetometer_measurements, float dt){
    if(gyro == NULLPTR || last_magnetometer_measurements == NULLPTR || magnetometer_measurements == NULLPTR){
        return false;
    }
    float dM[3] = {magnetometer_measurements[0] - last_magnetometer_measurements[0],
                magnetometer_measurements[1] - last_magnetometer_measurements[2],
                magnetometer_measurements[1] - last_magnetometer_measurements[2]};
    scale(dM, dt, 1, 3);
    float photodiode_sun_vec[3];
    bool in_sun = get_vec_from_photodiode_readings(photodiode_measurements, photodiode_sun_vec);
    
    return (l2_norm(dM, 3) < DETUMBLING_MAG_THRESHOLD || l2_norm(gyro, 3) < DETUMBLING_GYRO_THRESHOLD) && in_sun;
}
void body(float* last_magnetometer_measurements, //1x3
    float* magnetometer_measurements,  //1x3
    float* gyro_measurements,//1x3, radians
    float* photodiode_measurements, //1xNUM_DIODES, raw readings
    float* posn_update,//1x6, may be NULLPTR
    int unix,//unix time
    int jd_scalar,//scalar JD
    float jd_frac,//fractional JD
    float dt,//time since last update
    float* output_currents)
{
    //Update position
    if(posn_update != NULLPTR){
        memcpy(kepler_posn, posn_update, sizeof(float) * 6);
        posn_initialized = true;
    }else{
        if(posn_initialized){
            float new_posn[6];
            propogateOrbitalElements(kepler_posn[0], kepler_posn[1], kepler_posn[2], kepler_posn[3], kepler_posn[4], kepler_posn[5], dt, new_posn);
            memcpy(kepler_posn, new_posn, sizeof(float) * 6);
        }
    }
    
    //Make sure position exists before doing anything else
    if(!posn_initialized){
        return;
    }

    //During detumbling...
    if(!pointing){

        //Case 1: Measurements tell us to stop detumbling
        if(ready_to_point(photodiode_measurements, gyro_measurements, last_magnetometer_measurements, magnetometer_measurements, dt)){
            //Assumptions: We only detumble if A) we just got punted out or B) our covariance matrix has failed and we need to reset the filter. In either case, we must reset our estimate for both quat and bias before detumbling, so we can't consider bias here
            float photodiode_sun_vec[3];
            get_vec_from_photodiode_readings(photodiode_measurements, photodiode_sun_vec);
    
            //Gets expected magnetometer reading
            float expected_mag[3];
            float eci[3];
            orbitalToECI(kepler_posn, eci);
            wmm_eci_embedded_v2(eci, jd_scalar, jd_frac, expected_mag);
            
            //Gets expected sun vector reading
            float expected_sun_vec[3];
            sun_vec(unix, expected_sun_vec);

            //Assembles ref & body measurements
            float body[6] = {magnetometer_measurements[0], magnetometer_measurements[1], magnetometer_measurements[2], photodiode_sun_vec[0], photodiode_sun_vec[1], photodiode_sun_vec[2]};
            float ref[6] = {expected_mag[0], expected_mag[1], expected_mag[2], expected_sun_vec[0], expected_sun_vec[1], expected_sun_vec[2]};

            //Uses QUEST to get initial attitude estimate
            float estimated_q[4];
            quest(body, ref, 2, estimated_q);
            attitude_initialized = true;

            //Prepares filter
            memcpy(estimated_quat, estimated_q, sizeof(float) * 4);
            memset(estimated_gyro_bias, 0, sizeof(float) * 3);
            memset(error_quat_state, 0, sizeof(float) * 6);
            memcpy(error_quat_cov, INIT_P, sizeof(float) * 6 * 6);

            //Enables pointing
            pointing = true;

            //Sets currents to 0 to prepare
            memset(output_currents, 0, sizeof(float) * 3);
            return;
        }
        else{//rely on BDot to keep detumbling or to hold our attitude steady
            float moments[3];
            bdot(magnetometer_measurements, last_magnetometer_measurements, dt, moments);
            moment2current3axis(moments, Imax, output_currents);
            return;
        }
    }
    //During pointing
    else{
        float photodiode_sun_vec[3];
        bool in_sun = get_vec_from_photodiode_readings(photodiode_measurements, photodiode_sun_vec);
       
        //Gets expected magnetometer reading
        float expected_mag[3];
        float eci[3];
        orbitalToECI(kepler_posn, eci);
        wmm_eci_embedded_v2(eci, jd_scalar, jd_frac, expected_mag);

        //Filter during the sun
        if(in_sun){
            //Gets expected sun vector reading
            float expected_sun_vec[3];
            sun_vec(unix, expected_sun_vec);

            //Assembles ref & body measurements
            float body[6] = {magnetometer_measurements[0], magnetometer_measurements[1], magnetometer_measurements[2], photodiode_sun_vec[0], photodiode_sun_vec[1], photodiode_sun_vec[2]};
            float ref[6] = {expected_mag[0], expected_mag[1], expected_mag[2], expected_sun_vec[0], expected_sun_vec[1], expected_sun_vec[2]};
            
            //Run the filter
            float new_error_state[6];
            float new_estimated_quat[4];
            float new_P[6 * 6];
            iterate(error_quat_state, estimated_quat, error_quat_cov, body, ref, gyro_measurements, Q, R_IN_SUN, dt, new_error_state, new_estimated_quat, new_P);

            //TODO: add a better failure check here
            if(filter_failure(new_estimated_quat, gyro_measurements, new_P, dt)){
                //Reset everything
                memset(estimated_quat, 0, sizeof(float) * 4);
                memset(estimated_gyro_bias, 0, sizeof(float) * 3);
                memset(error_quat_state, 0, sizeof(float) * 6);
                memset(error_quat_cov, 0, sizeof(float) * 6 * 6);

                //Return to detumbling just in case
                attitude_initialized = false;
                pointing = false;

                //Sets currents to 0 to prepare
                memset(output_currents, 0, sizeof(float) * 3);
                return;
            }

            //Prepare for next run
            error_quat_state[0] = 0;
            error_quat_state[1] = 0;
            error_quat_state[2] = 0;
            error_quat_state[3] = new_error_state[3];
            error_quat_state[4] = new_error_state[4];
            error_quat_state[5] = new_error_state[5];

            memcpy(error_quat_cov, new_P, sizeof(float) * 6 * 6);
            memcpy(estimated_quat, new_estimated_quat, sizeof(float) * 4);

            //Track the bias
            float bias[3] = {new_error_state[3], new_error_state[4], new_error_state[5]};
            memcpy(estimated_gyro_bias, bias, sizeof(float) * 3);
        }
        //Filter in the shade
        else if(last_ref_magnetosphere_initialized){
            //Simulates a reading using the last magnetosphere initialization measurement
            float finite_difference_reference[3];
            sub(expected_mag, last_ref_magnetosphere, finite_difference_reference, 1, 3, 3);
            scale(finite_difference_reference, dt, 1, 3);

            float finite_difference_body[3];
            sub(magnetometer_measurements, last_magnetometer_measurements, finite_difference_body, 1, 3, 3);
            scale(finite_difference_body, dt, 1, 3);

            float cross_mag_body[3];
            float unbiased_gyro[3];
            sub(gyro_measurements, estimated_gyro_bias, unbiased_gyro, 1, 3, 3);
            cross(unbiased_gyro, magnetometer_measurements, cross_mag_body);
            float finite_difference_body_updated[3];
            add(finite_difference_body, cross_mag_body, finite_difference_body_updated, 1, 3, 3);

            //Assembles ref & body measurements
            float body[6] = {magnetometer_measurements[0], magnetometer_measurements[1], magnetometer_measurements[2], finite_difference_body_updated[0], finite_difference_body_updated[1], finite_difference_body_updated[2]};

            float ref[6] = {expected_mag[0], expected_mag[1], expected_mag[2], finite_difference_reference[0], finite_difference_reference[1], finite_difference_reference[2]};

            //Run the filter with UNBIASED gyro
            float new_error_state[6];
            float new_estimated_quat[4];
            float new_P[6 * 6];
            iterate(error_quat_state, estimated_quat, error_quat_cov, body, ref, unbiased_gyro, Q, R_IN_SUN, dt, new_error_state, new_estimated_quat, new_P);

            //TODO: add a better failure check here
            if(filter_failure(new_estimated_quat, gyro_measurements, new_P, dt)){
                //Reset everything
                memset(estimated_quat, 0, sizeof(float) * 4);
                memset(estimated_gyro_bias, 0, sizeof(float) * 3);
                memset(error_quat_state, 0, sizeof(float) * 6);
                memset(error_quat_cov, 0, sizeof(float) * 6 * 6);

                //Return to detumbling just in case
                attitude_initialized = false;
                pointing = false;

                //Sets currents to 0 to prepare
                memset(output_currents, 0, sizeof(float) * 3);
                return;
            }

            //Prepare for next run w/out tracking bias(in other words, we're relying on previous gyro bias tracking during our time in the sun)
            error_quat_state[0] = 0;
            error_quat_state[1] = 0;
            error_quat_state[2] = 0;
            error_quat_state[3] = 0;
            error_quat_state[4] = 0;
            error_quat_state[5] = 0;

            memcpy(error_quat_cov, new_P, sizeof(float) * 6 * 6);
            memcpy(estimated_quat, new_estimated_quat, sizeof(float) * 4);
        }

        //Sets the last reference magnetosphere model
        last_ref_magnetosphere_initialized = true;
        memcpy(last_ref_magnetosphere, expected_mag, sizeof(float) * 3);

        //Uses estimated quat to get goal q
        float pvd_eci[3];
        ecef_2_eci(PVD_ECEF, pvd_eci, unix);
        float goal_q[4];
        down_quat(eci, pvd_eci, estimated_quat, goal_q);

        //Find the error quat
        float q_diff[4];
        quat_diff(estimated_quat, goal_q, q_diff);

        //Find the rotation vec
        float rotvec[3];
        quat2rotationvec(q_diff, rotvec);

        //Use PD to find torques
        float t[3];
        PD_loop(rotvec, gyro_measurements, t);

        //Use torques to find moments
        float m[3];
        torque_2_moments(magnetometer_measurements, t, m);

        //Use moments to find currents
        moment2current3axis(m, Imax, output_currents);
        return;
    }

    
}