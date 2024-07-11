#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the proper number of input and output arguments
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MyToolbox:proj:nrhs", "Four inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:proj:nlhs", "One output required.");
    }

    // Get the scalar nelems
    double nelems = mxGetScalar(prhs[0]);
    double nob = mxGetScalar(prhs[1]);

    // Create a cell array for the output
    plhs[0] = mxCreateCellMatrix(nelems, 1);

    // Loop through the elements
    for (int i = 0; i < (int)nelems; i++) {
        mxArray *xCell = mxGetCell(prhs[2], i);
        mxArray *uCell = mxGetCell(prhs[3], i);

        double *X = mxGetPr(xCell);
        double *U = mxGetPr(uCell);

        mwSize numRows = mxGetM(xCell);
        mwSize numCols = mxGetN(xCell);

        // Create the output matrix for the current element
        mxArray *vCell = mxCreateDoubleMatrix(numRows, numCols, mxREAL);
        double *V = mxGetPr(vCell);

        if (i < (int)nob) {
            // Compute the sum of X.*U
            double sumXU = 0;
            for (int j = 0; j < numRows * numCols; j++) {
                sumXU += X[j] * U[j];
            } 
           // Compute V = U - X.*sumXU
           for (int j = 0; j < numRows * numCols; j++) {
                V[j] = U[j] - X[j] * sumXU;
            }
        } else {
            for (int j = 0; j < numRows * numCols; j++) {
                V[j] = U[j];
            }
        }

        // Assign the computed matrix to the output cell array
        mxSetCell(plhs[0], i, vCell);
    }
}
