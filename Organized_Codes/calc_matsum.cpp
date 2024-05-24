// summing rows of a given matrix

#include <stdio.h>
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double *input, *output;
    int i, j, n_rows, n_cols;
    double sum;

    //if ((nrhs != 1) || (nlhs > 1)) {
    //    printf("Error! Expecting exactly 1 rhs and up to 1 lhs argument!\n");
    //    return;
    //}
        
    input = mxGetPr(prhs[0]);
    n_rows = mxGetM(prhs[0]);
    n_cols = mxGetN(prhs[0]);


    if (nlhs == 1) {
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        output = mxGetPr(plhs[0]);
    }
        
    sum = 0;
    for(i=0; i<n_rows; i++) {
        for(j=0; j<n_cols; j++) 
            sum += input[(i*n_cols)+j];
            //printf("value = %f\\", input[1]);      
    }
    if (nlhs == 1)
        *output = sum;
} /* end mexFunction */