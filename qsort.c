#include "matrix.h"

static const mxDouble *array = NULL;
int compare(const void *pa, const void *pb);
int compare(const void *pa, const void *pb) {
    const mxDouble a = array[(size_t) * (const mxDouble *)pa];
    const mxDouble b = array[(size_t) * (const mxDouble *)pb];
    return a == b ? 0 : (a < b ? 1 : -1);
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nlhs == 1 && nrhs == 1) {
        if (mxGetNumberOfDimensions(prhs[0]) != 2) {
            return;
        }
        const size_t m = mxGetM(prhs[0]);
        const size_t n = mxGetN(prhs[0]);
        const size_t numel = mxGetNumberOfElements(prhs[0]);
        plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
        mxDouble *index = mxGetDoubles(plhs[0]);
        for (size_t i = 0; i < numel; i++) {
            index[i] = (mxDouble)i;
        }
        array = mxGetDoubles(prhs[0]);
        qsort(mxGetDoubles(plhs[0]), mxGetNumberOfElements(plhs[0]), mxGetElementSize(plhs[0]), compare);
        for (size_t i = 0; i < numel; i++) {
            index[i] += 1;
        }
    }
}
