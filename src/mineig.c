#include "mex.h"
#include "lapacke.h"
#include "lapacke_config.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])

{
    double *A, *W, *Z;
    lapack_int *ifail = NULL, n, M;
               
    /* Assign pointers to the various parameters */
    A = mxGetPr(prhs[0]);
    n = mxGetN(prhs[0]);
    
    /* Create a matrix for the return argument */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);

    Z = mxGetPr(plhs[0]);
    W = mxGetPr(plhs[1]);
   
    /* Do the actual computations in a subroutine */
    // LAPACKE_dsyevx(LAPACK_COL_MAJOR, 'V', 'V', 'U', (lapack_int)n, A, (lapack_int)n, 0.0, 100000.0, (lapack_int)1, (lapack_int)1, 0.0, (lapack_int*)M, W, Z, (lapack_int)n, ifail);
    LAPACKE_dsyevr(LAPACK_COL_MAJOR, 'V', 'I', 'U', n, A, n, 0.000001, 100000.0, 1, 1, -1.0,   &M, W, Z, n, ifail);
    return;
}