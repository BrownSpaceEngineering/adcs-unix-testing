#include "include/photodiode_determination.h"
#include "declareFunctions.h"
#include "math.h"

const int NUM_PAIRS = 9;
const int NUM_SELECTED_DIODES = 5;

const double MAX_READING = 1.7;
const double MIN_READING = 0.17;
const double PHOTODIODES[18][3] = {
    {  0.8660254,  -0.27968387, -0.41445981 },
    { -0.8660254,   0.27968387,  0.41445981 },

    {  0.8660254,  -0.16748457,  0.47111455 },
    { -0.8660254,   0.16748457, -0.47111455 },

    {  0.8660254,   0.4999767,   0.00482655 },
    { -0.8660254,  -0.4999767,  -0.00482655 },

    { -0.07387608,  0.8660254,   0.49451221 },
    {  0.07387608, -0.8660254,  -0.49451221 },

    { -0.48736734,  0.8660254,  -0.11168291 },
    {  0.48736734, -0.8660254,   0.11168291 },

    {  0.32616791,  0.8660254,  -0.37896503 },
    { -0.32616791, -0.8660254,   0.37896503 },

    { -0.36819159, -0.33828235, -0.8660254  },
    {  0.36819159,  0.33828235,  0.8660254  },

    {  0.48601695,  0.1174203,  -0.8660254  },
    { -0.48601695, -0.1174203,   0.8660254  },

    { -0.23581003,  0.44090093, -0.8660254  },
    {  0.23581003, -0.44090093,  0.8660254  }
};

// Assumes photodiode_readings are passed in the same order as PHOTODIODES.
// Pairs are 0-1, 2-3, ..., 16-17.
bool get_vec_from_photodiode_readings(double* photodiode_readings,
                                      double* estimated_sun_vector) {
    int valid_readings = 0;
    int selected_indices[NUM_SELECTED_DIODES];
    double selected_readings[NUM_SELECTED_DIODES];

    for (int i = 0; i < NUM_SELECTED_DIODES; i++) {
        selected_indices[i] = -1;
        selected_readings[i] = -1.0;
    }

    // For each pair, choose the brighter diode.
    // Keep only the NUM_SELECTED_DIODES brightest pair-winners
    for (int pair = 0; pair < NUM_PAIRS; pair++) {
        int i0 = 2 * pair;
        int i1 = 2 * pair + 1;

        int brighter_index;
        double brighter_reading;

        if (photodiode_readings[i0] >= photodiode_readings[i1]) {
            brighter_index = i0;
            brighter_reading = photodiode_readings[i0];
        } else {
            brighter_index = i1;
            brighter_reading = photodiode_readings[i1];
        }

        if(brighter_reading > MIN_READING){
            valid_readings += 1;
        }

        int min_pos = 0;
        for (int i = 1; i < NUM_SELECTED_DIODES; i++) {
            if (selected_readings[i] < selected_readings[min_pos]) {
                min_pos = i;
            }
        }

        if (brighter_reading > selected_readings[min_pos]) {
            selected_readings[min_pos] = brighter_reading;
            selected_indices[min_pos] = brighter_index;
        }
    }

    if(valid_readings < 3){
        return false;
    }

    // Construct system of equations
    double A[NUM_SELECTED_DIODES * 3];
    double b[NUM_SELECTED_DIODES];

    for (int i = 0; i < NUM_SELECTED_DIODES; i++) {
        int diode_index = selected_indices[i];

        A[i * 3 + 0] = PHOTODIODES[diode_index][0];
        A[i * 3 + 1] = PHOTODIODES[diode_index][1];
        A[i * 3 + 2] = PHOTODIODES[diode_index][2];

        b[i] = selected_readings[i] / MAX_READING;
    }

    // Least-squares inverse
    pinv(A, NUM_SELECTED_DIODES, 3);

    // Solving using inverse
    mul(A, b, false, estimated_sun_vector, 3, NUM_SELECTED_DIODES, 1);

    // Normalize output vector
    double mag = sqrt(
        estimated_sun_vector[0] * estimated_sun_vector[0] +
        estimated_sun_vector[1] * estimated_sun_vector[1] +
        estimated_sun_vector[2] * estimated_sun_vector[2]
    );

    estimated_sun_vector[0] /= mag;
    estimated_sun_vector[1] /= mag;
    estimated_sun_vector[2] /= mag;

    return true;

}
