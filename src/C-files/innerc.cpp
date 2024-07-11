#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the proper number of input and output arguments
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:inner:nrhs", "Three inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:inner:nlhs", "One output required.");
    }

    // Get the scalar nelems
    double nelems = mxGetScalar(prhs[0]);

    // Create a pointer for the output
    plhs[0] = mxCreateDoubleScalar(0);
    double *val = mxGetPr(plhs[0]);

    // Loop through the elements
    for (int i = 0; i < (int)nelems; i++) {
        mxArray *uCell = mxGetCell(prhs[1], i);
        mxArray *vCell = mxGetCell(prhs[2], i);
        
        double *u = mxGetPr(uCell);
        double *v = mxGetPr(vCell);

        mwSize numel = mxGetNumberOfElements(uCell);

        for (int j = 0; j < numel; j++) {
            *val += u[j] * v[j];
        }
    }
}
