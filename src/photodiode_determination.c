#include "include/photodiode_determination.h"
#include "declareFunctions.h"
#include "math.h"

const int NUM_PAIRS = 9;
const int NUM_SELECTED_DIODES = 5;


const float MAX_READING = 1.7f;
const float MIN_READING = 0.17f;
const float PHOTODIODES[18][3] = {
    {  0.8660254f,  -0.27968387f, -0.41445981f },
    { -0.8660254f,   0.27968387f,  0.41445981f },

    {  0.8660254f,  -0.16748457f,  0.47111455f },
    { -0.8660254f,   0.16748457f, -0.47111455f },

    {  0.8660254f,   0.4999767f,   0.00482655f },
    { -0.8660254f,  -0.4999767f,  -0.00482655f },

    { -0.07387608f,  0.8660254f,   0.49451221f },
    {  0.07387608f, -0.8660254f,  -0.49451221f },

    { -0.48736734f,  0.8660254f,  -0.11168291f },
    {  0.48736734f, -0.8660254f,   0.11168291f },

    {  0.32616791f,  0.8660254f,  -0.37896503f },
    { -0.32616791f, -0.8660254f,   0.37896503f },

    { -0.36819159f, -0.33828235f, -0.8660254f  },
    {  0.36819159f,  0.33828235f,  0.8660254f  },

    {  0.48601695f,  0.1174203f,  -0.8660254f  },
    { -0.48601695f, -0.1174203f,   0.8660254f  },

    { -0.23581003f,  0.44090093f, -0.8660254f  },
    {  0.23581003f, -0.44090093f,  0.8660254f  }
};

// Assumes photodiode_readings are passed in the same order as PHOTODIODES.
// Pairs are 0-1, 2-3, ..., 16-17.
bool get_vec_from_photodiode_readings(float* photodiode_readings,
                                      float* estimated_sun_vector) {
    int valid_readings = 0;
    int selected_indices[NUM_SELECTED_DIODES];
    float selected_readings[NUM_SELECTED_DIODES];

    for (int i = 0; i < NUM_SELECTED_DIODES; i++) {
        selected_indices[i] = -1;
        selected_readings[i] = -1.0f;
    }

    // For each pair, choose the brighter diode.
    // Keep only the NUM_SELECTED_DIODES brightest pair-winners
    for (int pair = 0; pair < NUM_PAIRS; pair++) {
        int i0 = 2 * pair;
        int i1 = 2 * pair + 1;

        int brighter_index;
        float brighter_reading;

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

    float A[NUM_SELECTED_DIODES * 3];
    float b[NUM_SELECTED_DIODES];

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
    mul(A, b, false,
        estimated_sun_vector,
        3,
        NUM_SELECTED_DIODES,
        1);

    // Normalize output vector
    float mag = sqrtf(
        estimated_sun_vector[0] * estimated_sun_vector[0] +
        estimated_sun_vector[1] * estimated_sun_vector[1] +
        estimated_sun_vector[2] * estimated_sun_vector[2]
    );

    estimated_sun_vector[0] /= mag;
    estimated_sun_vector[1] /= mag;
    estimated_sun_vector[2] /= mag;

    return true;
}