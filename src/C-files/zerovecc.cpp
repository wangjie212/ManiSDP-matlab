#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the proper number of input and output arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:zerovec:nrhs", "Two inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:zerovec:nlhs", "One output required.");
    }

    // Get the scalar nelems
    double nelems = mxGetScalar(prhs[0]);

    // Create a cell array for the output
    plhs[0] = mxCreateCellMatrix(nelems, 1);

    // Loop through the elements
    for (int i = 0; i < (int)nelems; i++) {
        mxArray *xCell = mxGetCell(prhs[1], i);

        mwSize numRows = mxGetM(xCell);
        mwSize numCols = mxGetN(xCell);

        // Create the zero matrix for the current element
        mxArray *uCell = mxCreateDoubleMatrix(numRows, numCols, mxREAL);

        // Assign the zero matrix to the output cell array
        mxSetCell(plhs[0], i, uCell);
    }
}
