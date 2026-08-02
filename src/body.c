
#include "string.h"
#include "math.h"
#include "include/body.h"
#include "include/laextension.h"
#include "include/quest.h"
#include "include/bdot.h"
#include "include/moments2currents.h"
#include "include/photodiode_determination.h"
#include "include/magnetosphere.h"
#include "include/orbital2eci.h"
#include "include/sunvec.h"
#include "include/kepler.h"
#include "include/iterate.h"
#include "include/quat.h"
#include "include/pd.h"
#include "include/ecef2eci.h"
#include "include/down_quat.h"
#include "include/torque2moments.h"
#include "declareFunctions.h"

double PVD_ECEF[3] = {6.3761e6, -0.1387e6, 0.0807e6};
const double DETUMBLING_GYRO_THRESHOLD = 0.17;
const double DETUMBLING_MAG_THRESHOLD = 0.01;
// iterate() is called on every body() invocation for gyro propagation (mode 0),
// and every FILTER_UPDATE_EVERY calls a full measurement update (mode 1 or 2) is added.
// At 10 Hz this gives a 0.1 Hz measurement update rate, matching test_long.c.
static const int FILTER_UPDATE_EVERY = 100;
double Imax[3] = {1, 1, 1};//TEMP
const double INIT_P[6 * 6] = {0.01, 0, 0, 0, 0, 0,
                                0, 0.01, 0, 0, 0, 0,
                                0, 0, 0.01, 0, 0, 0,
                                0, 0, 0, 0.01, 0, 0,
                                0, 0, 0, 0, 0.01, 0,
                                0, 0, 0, 0, 0, 0.01};
const double Q[6 * 6] =  {0.01, 0, 0, 0, 0, 0,
                                0, 0.01, 0, 0, 0, 0,
                                0, 0, 0.01, 0, 0, 0,
                                0, 0, 0, 0.01, 0, 0,
                                0, 0, 0, 0, 0.01, 0,
                                0, 0, 0, 0, 0, 0.01};
const double R_IN_SUN[6 * 6] =  {0.01, 0, 0, 0, 0, 0,
                                        0, 0.01, 0, 0, 0, 0,
                                        0, 0, 0.01, 0, 0, 0,
                                        0, 0, 0, 0.01, 0, 0,
                                        0, 0, 0, 0, 0.01, 0,
                                        0, 0, 0, 0, 0, 0.01};
const double R_IN_SHADOW[3 * 3] = {0.01, 0, 0,
                                         0, 0.01, 0,
                                         0, 0, 0.01};
bool filter_failure(double* new_quat, double* gyro, double* new_P, double dt){
    //Fail if quaternion jump implies angular rate inconsistent with gyro measurement
    // estimated_quat and new_quat are already double [body→ECI]; use directly
    double diff[4]; // inter-step rotation: new_quat ⊗ estimated_quat^-1
    quat_diff(estimated_quat, new_quat, diff);
    double diff_rot[3];
    quat2rotationvec(diff, diff_rot);
    scale(diff_rot, 1.0 / dt, 1, 3);
    double gyro_norm = l2_norm(gyro, 3);
    if(gyro_norm > 0.0 && l2_norm(diff_rot, 3) / gyro_norm > 10){
        return true;
    }

    //Fail if new covariance contains NaN
    for(int i = 0; i < 6 * 6; i++){
        if(isnan(new_P[i])){
            return true;
        }
    }
    return false;
}
bool ready_to_point(double* photodiode_measurements, double* gyro, double* last_magnetometer_measurements, double* magnetometer_measurements, double dt){
    if(gyro == NULLPTR || last_magnetometer_measurements == NULLPTR || magnetometer_measurements == NULLPTR || photodiode_measurements == NULLPTR){
        return false;
    }
    // Component-wise magnetometer delta, then divide by dt to get rate of change [T/s]
    double dM[3] = {magnetometer_measurements[0] - last_magnetometer_measurements[0],
                    magnetometer_measurements[1] - last_magnetometer_measurements[1],
                    magnetometer_measurements[2] - last_magnetometer_measurements[2]};
    scale(dM, 1.0 / dt, 1, 3);
    double sun_d[3];
    bool in_sun = get_vec_from_photodiode_readings(photodiode_measurements, sun_d);
    return (l2_norm(dM, 3) < DETUMBLING_MAG_THRESHOLD && l2_norm(gyro, 3) < DETUMBLING_GYRO_THRESHOLD) && in_sun;
}
void body(double* last_magnetometer_measurements, //1x3
    double* magnetometer_measurements,  //1x3
    double* gyro_measurements,//1x3, radians
    double* photodiode_measurements, //1xNUM_DIODES, raw readings
    double* posn_update,//1x6, may be NULLPTR
    int unix,//unix time
    int jd_scalar,//scalar JD
    double jd_frac,//fractional JD
    double dt,//time since last update
    double* output_currents)
{
    //Update position
    if(posn_update != NULLPTR){
        memcpy(kepler_posn, posn_update, sizeof(double) * 6);
        posn_initialized = true;
    }else{
        if(posn_initialized){
            double new_posn[6];
            propogateOrbitalElements(kepler_posn[0], kepler_posn[1], kepler_posn[2], kepler_posn[3], kepler_posn[4], kepler_posn[5], dt, new_posn);
            memcpy(kepler_posn, new_posn, sizeof(double) * 6);
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
            double photodiode_sun_vec[3];
            get_vec_from_photodiode_readings(photodiode_measurements, photodiode_sun_vec);
    
            //Gets expected magnetometer reading
            double eci_qi[3];
            orbitalToECI(kepler_posn, eci_qi);
            double expected_mag_qi[3];
            wmm_eci_embedded_v2(eci_qi, jd_scalar, jd_frac, expected_mag_qi);
            
            //Gets expected sun vector reading
            double expected_sun_vec[3];
            sun_vec(unix, expected_sun_vec);

            //Assembles ref & body measurements
            double body_d[6] = {magnetometer_measurements[0], magnetometer_measurements[1], magnetometer_measurements[2], photodiode_sun_vec[0], photodiode_sun_vec[1], photodiode_sun_vec[2]};
            double ref_d[6] = {expected_mag_qi[0], expected_mag_qi[1], expected_mag_qi[2], expected_sun_vec[0], expected_sun_vec[1], expected_sun_vec[2]};

            //Uses QUEST to get initial attitude estimate
            double estimated_q[4]; // body→ECI  (QUEST convention: body→ref)
            quest(body_d, ref_d, 2, estimated_q);

            //Prepares filter
            memcpy(estimated_quat, estimated_q, sizeof(double) * 4);
            memset(error_quat_state, 0, sizeof(double) * 6);
            memcpy(error_quat_cov, INIT_P, sizeof(double) * 6 * 6);

            //Enables pointing
            filter_step = 0;
            pointing = true;

            //Sets currents to 0 to prepare
            memset(output_currents, 0, sizeof(double) * 3);
            return;
        }
        else{// not yet ready to point; apply B-dot to continue reducing angular velocity
            double moments_bd[3];
            Bdot(magnetometer_measurements, last_magnetometer_measurements, dt, moments_bd);
            moment2current3axis(moments_bd, Imax, output_currents);
            return;
        }
    }
    //During pointing: iterate() runs on every call (10 Hz gyro propagation).
    //Every FILTER_UPDATE_EVERY calls a measurement update is also applied.
    else{
        double photodiode_sun_vec[3];
        bool in_sun = get_vec_from_photodiode_readings(photodiode_measurements, photodiode_sun_vec);
       
        //Gets expected magnetometer reading
        double eci[3];
        orbitalToECI(kepler_posn, eci);
        double expected_mag[3];
        wmm_eci_embedded_v2(eci, jd_scalar, jd_frac, expected_mag);

        if(filter_step % FILTER_UPDATE_EVERY == 0){
            //Measurement update step (every ~10 s at 10 Hz)

            //Filter during the sun (2-vector: magnetometer + photodiode sun vector)
            if(in_sun){
                //Gets expected sun vector reading
                double expected_sun_vec[3];
                sun_vec(unix, expected_sun_vec);

                //Assembles ref & body measurements
                double body[6] = {magnetometer_measurements[0], magnetometer_measurements[1], magnetometer_measurements[2], photodiode_sun_vec[0], photodiode_sun_vec[1], photodiode_sun_vec[2]};
                double ref[6] = {expected_mag[0], expected_mag[1], expected_mag[2], expected_sun_vec[0], expected_sun_vec[1], expected_sun_vec[2]};
                
                //Run the filter
                double new_error_state[6];
                double new_estimated_quat[4]; // updated estimate [body→ECI]
                double new_P[6 * 6];
                iterate(error_quat_state, estimated_quat, error_quat_cov, body, ref, 2, gyro_measurements, Q, R_IN_SUN, dt, new_error_state, new_estimated_quat, new_P);

                //TODO: add a better failure check here
                if(filter_failure(new_estimated_quat, gyro_measurements, new_P, dt)){
                    //Reset everything
                    estimated_quat[0] = 1.0; estimated_quat[1] = 0.0;
                    estimated_quat[2] = 0.0; estimated_quat[3] = 0.0;
                    memset(error_quat_state, 0, sizeof(double) * 6);
                    memset(error_quat_cov, 0, sizeof(double) * 6 * 6);

                    //Return to detumbling just in case
                    pointing = false;

                    //Sets currents to 0 to prepare
                    memset(output_currents, 0, sizeof(double) * 3);
                    return;
                }

                //Prepare for next run
                error_quat_state[0] = 0;
                error_quat_state[1] = 0;
                error_quat_state[2] = 0;
                error_quat_state[3] = new_error_state[3];
                error_quat_state[4] = new_error_state[4];
                error_quat_state[5] = new_error_state[5];

                memcpy(error_quat_cov, new_P, sizeof(double) * 6 * 6);
                memcpy(estimated_quat, new_estimated_quat, sizeof(double) * 4);
            }
            //Filter in eclipse: single-vector magnetometer only (no sun reference available).
            else{
                double body[3] = {magnetometer_measurements[0], magnetometer_measurements[1], magnetometer_measurements[2]};
                double ref[3]  = {expected_mag[0], expected_mag[1], expected_mag[2]};
                double n_body = sqrt(body[0]*body[0] + body[1]*body[1] + body[2]*body[2]);
                double n_ref  = sqrt(ref[0]*ref[0] + ref[1]*ref[1] + ref[2]*ref[2]);
                for(int i = 0; i < 3; i++){ body[i] /= n_body; ref[i] /= n_ref; }

                double new_error_state[6];
                double new_estimated_quat[4]; // updated estimate [body→ECI]
                double new_P[6 * 6];
                iterate(error_quat_state, estimated_quat, error_quat_cov, body, ref, 1,
                        gyro_measurements, Q, R_IN_SHADOW, dt, new_error_state, new_estimated_quat, new_P);

                //TODO: add a better failure check here
                if(filter_failure(new_estimated_quat, gyro_measurements, new_P, dt)){
                    //Reset everything
                    estimated_quat[0] = 1.0; estimated_quat[1] = 0.0;
                    estimated_quat[2] = 0.0; estimated_quat[3] = 0.0;
                    memset(error_quat_state, 0, sizeof(double) * 6);
                    memset(error_quat_cov, 0, sizeof(double) * 6 * 6);

                    //Return to detumbling just in case
                    pointing = false;

                    //Sets currents to 0 to prepare
                    memset(output_currents, 0, sizeof(double) * 3);
                    return;
                }

                //Reset attitude error; preserve and update bias estimate (mirrors test_long.c behaviour)
                error_quat_state[0] = 0;
                error_quat_state[1] = 0;
                error_quat_state[2] = 0;
                error_quat_state[3] = new_error_state[3];
                error_quat_state[4] = new_error_state[4];
                error_quat_state[5] = new_error_state[5];

                memcpy(error_quat_cov, new_P, sizeof(double) * 6 * 6);
                memcpy(estimated_quat, new_estimated_quat, sizeof(double) * 4);
            }
        }
        else{
            //Propagation-only step: integrate gyro dynamics, advance covariance, skip sensor correction
            double new_error_state[6];
            double new_estimated_quat[4]; // updated estimate [body→ECI]
            double new_P[6 * 6];
            iterate(error_quat_state, estimated_quat, error_quat_cov, NULL, NULL, 0,
                    gyro_measurements, Q, NULL, dt, new_error_state, new_estimated_quat, new_P);

            error_quat_state[0] = 0;
            error_quat_state[1] = 0;
            error_quat_state[2] = 0;
            error_quat_state[3] = new_error_state[3];
            error_quat_state[4] = new_error_state[4];
            error_quat_state[5] = new_error_state[5];
            memcpy(error_quat_cov, new_P, sizeof(double) * 6 * 6);
            memcpy(estimated_quat, new_estimated_quat, sizeof(double) * 4);
        }

        filter_step++;

        double pvd_eci[3];
        ecef_2_eci(PVD_ECEF, pvd_eci, unix);
        double goal_q[4]; // desired attitude [body→ECI]
        // nadir = pvd_eci - eci: body z-axis points toward pvd (nadir/target direction)
        down_quat(eci, pvd_eci, estimated_quat, goal_q);

        //Find the body-frame attitude error: estimated_quat^-1 ⊗ goal_q
        // Body-frame error keeps rotvec in body frame, matching gyro_measurements
        double est_q_inv[4]; // ECI→body
        quat_inv(estimated_quat, est_q_inv);
        double q_diff_d[4]; // attitude error in body frame [body→ECI delta]
        quat_multiply(est_q_inv, goal_q, q_diff_d);

        //Find the rotation vec (body frame)
        double rotvec[3];
        quat2rotationvec(q_diff_d, rotvec);

        //Use PD to find torques
        double t[3];
        PD_loop(rotvec, gyro_measurements, t);

        //Use torques to find moments
        double m[3];
        torque_2_moments(magnetometer_measurements, t, m);

        //Use moments to find currents
        moment2current3axis(m, Imax, output_currents);
        return;
    }

    
}
