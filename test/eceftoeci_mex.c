#include "matrix.h"
#include "mex.h"
#include "include/ecef2eci.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nlhs > 1) {
        mexErrMsgIdAndTxt("ecef2eci nlhs", "number of output error");
    }
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("ecef2eci nrhs", "number of input error");
    }
    if (!mxIsSingle(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfElements(prhs[0]) != 3)
        mexErrMsgIdAndTxt("ecef2eci:ecef", "ecef must be a real single 3-vector.");
    if (!mxIsInt32(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1)
        mexErrMsgIdAndTxt("ecef2eci:time", "unix_time must be an int32 scalar.");
    float* ecef = mxGetSingles(prhs[0]);
    int* unix_time = mxGetInt32s(prhs[1]);

    plhs[0] = mxCreateNumericMatrix(3, 1, mxSINGLE_CLASS, mxREAL);
    float* eci = mxGetSingles(plhs[0]);
    ecef_2_eci(ecef, eci, *unix_time);
}


    

