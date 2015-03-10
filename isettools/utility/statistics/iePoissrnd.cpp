/*
 * iePoissrnd.c - Poisson random number generator
 *
 *   Inputs:
 *     lambda  - Poisson mean
 *     nSamp   - number of samples
 *
 * This is a MEX-file for MATLAB
 * (HJ) April, 2014
*/

#include "mex.h"
#include "matrix.h"
#include <random>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 1)
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                          "Two inputs required.");
    double *lambda = mxGetPr(prhs[0]);;
    double *outputPr;
    int nSamp = 1;
    std::mt19937 mt; // Mersenne Twister generator
    
    if (mxGetNumberOfElements(prhs[0]) == 1) {
        // Single lambda, multiple samples
        // Output is the same size as samples
        std::poisson_distribution<> poisson(*lambda);
        if (nrhs < 2){
            plhs[0] = mxCreateDoubleScalar(poisson(mt));
        }
        else{
            int ndim = mxGetNumberOfElements(prhs[1]);
            int totProd = 1;
            mwSize *dimSize = (mwSize*)mxCalloc(ndim, sizeof(mwSize));
            double *dims = (double*)mxGetData(prhs[1]);
            for (int ii = 0; ii < ndim; ii++){
                dimSize[ii] = (mwSize)dims[ii];
                totProd *= dims[ii];
            }
            
            plhs[0] = mxCreateNumericArray(ndim, dimSize, mxDOUBLE_CLASS, mxREAL);
            outputPr = mxGetPr(plhs[0]);
            
            for (int ii = 0; ii < totProd; ii++){
                outputPr[ii] = poisson(mt);
            }
            mxFree(dimSize);
        }
    }
    else{
        // Lambda matrix, output is the same size as lambda
        const mwSize* dims = mxGetDimensions(prhs[0]);
        int ndim = mxGetNumberOfDimensions(prhs[0]);
        plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
        outputPr = mxGetPr(plhs[0]);
        for (int ii=0; ii<mxGetNumberOfElements(prhs[0]); ii++) {
            std::poisson_distribution<> poisson(lambda[ii]);
            outputPr[ii] = poisson(mt);
        }
    }
}
