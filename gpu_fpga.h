// extern int workers, gangs, size;
// extern float *A, *B, *D, *E;

void funcFPGA(
    float* X_result,
    const float* VAL,
    const int* COL_IND,
    const int* ROW_PTR,
    const float* B,
    const int N,
    const int K,
    const int VAL_SIZE
);
