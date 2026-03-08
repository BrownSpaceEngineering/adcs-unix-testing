/*
 * inv.c
 *
 *  Created on: 12 feb. 2019
 *      Author: Daniel Mårtensson
 */

#include "declareFunctions.h"

/*
 * Turn A into A^-1. A have the size row x row
 */

void inv(double* A, int row) {

    // Create identity matrix
    double I[row * row];
    eye(I, row, row);

    // Do inverse of A
    linsolve(A, A, I, row, row);
}

void invf(float* A, int row) {

    // Create identity matrix
    double I[row * row];
    eye(I, row, row);

    // Do inverse of A
    linsolvef(A, A, I, row, row);
}
