#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

double gaussRand() {
    static double V1, V2, S;
    static int phase = 0;
    double X;
    srand(3232);
    if (phase == 0) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while (S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else {
        X = V2 * sqrt(-2 * log(S) / S);
    }

    phase = 1 - phase;
    return X;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the proper number of input and output arguments
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MyToolbox:rand_custom:nrhs", "Four inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:rand_custom:nlhs", "One output required.");
    }

    // Get the scalar nelems
    double nelems = mxGetScalar(prhs[0]);
    double nob = mxGetScalar(prhs[1]);

    // Retrieve pset and nset
    double *pset = mxGetPr(prhs[2]);
    double *nset = mxGetPr(prhs[3]);

    // Create a cell array for the output
    plhs[0] = mxCreateCellMatrix(nelems, 1);

    // Seed the random number generator
    srand((unsigned int)time(NULL));

    // Loop through the elements
    for (int i = 0; i < (int)nelems; i++) {
        mwSize p = (mwSize)pset[i];
        mwSize n = (mwSize)nset[i];

        // Create the random matrix for the current element
        mxArray *xCell = mxCreateDoubleMatrix(p, n, mxREAL);
        double *X = mxGetPr(xCell);

        for (int j = 0; j < p * n; j++) {
            X[j] = gaussRand(); // Random number from standard Gaussian distribution
        }
        
        if (i < (int)nob) {
           // Normalize the columns by the square root of the sum of squares
           for (int col = 0; col < n; col++) {
               double sumSquares = 0;
               for (int row = 0; row < p; row++) {
                   sumSquares += X[row + col * p] * X[row + col * p];
               }
               double normFactor = sqrt(sumSquares);
               for (int row = 0; row < p; row++) {
                   X[row + col * p] /= normFactor;
               }
           }
        }

        // Assign the random matrix to the output cell array
        mxSetCell(plhs[0], i, xCell);
    }
}
