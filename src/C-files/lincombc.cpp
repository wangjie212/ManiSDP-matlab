#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the proper number of input and output arguments
    if (nrhs != 3 && nrhs != 5) {
        mexErrMsgIdAndTxt("MyToolbox:lincomb:nrhs", "Three or five inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:lincomb:nlhs", "One output required.");
    }

    // Get the scalar nelems
    double nelems = mxGetScalar(prhs[0]);

    // Create a cell array for the output
    plhs[0] = mxCreateCellMatrix(nelems, 1);

    // Retrieve the scalar multipliers
    double a1 = mxGetScalar(prhs[1]);
    
    if (nrhs == 3) { // Case with three input arguments
        for (int i = 0; i < (int)nelems; i++) {
            mxArray *u1Cell = mxGetCell(prhs[2], i);
            double *U1 = mxGetPr(u1Cell);

            mwSize numRows = mxGetM(u1Cell);
            mwSize numCols = mxGetN(u1Cell);

            // Create the output matrix for the current element
            mxArray *vCell = mxCreateDoubleMatrix(numRows, numCols, mxREAL);
            double *V = mxGetPr(vCell);

            // Compute V = a1 * u1{i}
            for (int j = 0; j < numRows * numCols; j++) {
                V[j] = a1 * U1[j];
            }

            // Assign the computed matrix to the output cell array
            mxSetCell(plhs[0], i, vCell);
        }
    } else if (nrhs == 5) { // Case with five input arguments
        double a2 = mxGetScalar(prhs[3]);

        for (int i = 0; i < (int)nelems; i++) {
            mxArray *u1Cell = mxGetCell(prhs[2], i);
            mxArray *u2Cell = mxGetCell(prhs[4], i);
            
            double *U1 = mxGetPr(u1Cell);
            double *U2 = mxGetPr(u2Cell);

            mwSize numRows = mxGetM(u1Cell);
            mwSize numCols = mxGetN(u1Cell);

            // Create the output matrix for the current element
            mxArray *vCell = mxCreateDoubleMatrix(numRows, numCols, mxREAL);
            double *V = mxGetPr(vCell);

            // Compute V = a1 * u1{i} + a2 * u2{i}
            for (int j = 0; j < numRows * numCols; j++) {
                V[j] = a1 * U1[j] + a2 * U2[j];
            }

            // Assign the computed matrix to the output cell array
            mxSetCell(plhs[0], i, vCell);
        }
    }
}
