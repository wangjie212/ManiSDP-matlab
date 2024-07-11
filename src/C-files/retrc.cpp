#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the proper number of input and output arguments
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:retr:nrhs", "Three inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:retr:nlhs", "One output required.");
    }

    // Get the scalar nelems
    double nelems = mxGetScalar(prhs[0]);

    // Get the cell arrays x and u
    mxArray *x = (mxArray *)prhs[1];
    mxArray *u = (mxArray *)prhs[2];

    // Create a cell array for the output
    plhs[0] = mxCreateCellMatrix(nelems, 1);

    // Loop through the elements
    for (int i = 0; i < (int)nelems; i++) {
        mxArray *xCell = mxGetCell(x, i);
        mxArray *uCell = mxGetCell(u, i);
        mwSize rows = mxGetM(xCell);
        mwSize cols = mxGetN(xCell);

        // Create yCell for the current element
        mxArray *yCell = mxCreateDoubleMatrix(rows, cols, mxREAL);
        double *Y = mxGetPr(yCell);
        double *X = mxGetPr(xCell);
        double *U = mxGetPr(uCell);

        for (int col = 0; col < cols; col++) {
            double sumSquares = 0;
            for (int row = 0; row < rows; row++) {
                int idx = row + col * rows;
                Y[idx] = X[idx] + U[idx];
                sumSquares += Y[idx] * Y[idx];
            }
            double normFactor = sqrt(sumSquares);
            for (int row = 0; row < rows; row++) {
                Y[row + col * rows] /= normFactor;
            }
        }

        // Assign the result to the output cell array
        mxSetCell(plhs[0], i, yCell);
    }
}
