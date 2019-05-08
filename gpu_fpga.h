// extern int workers, gangs, size;
// extern float *A, *B, *D, *E;

void funcFPGA(
    float* restrict X_result,
    const float* restrict VAL,
    const int* restrict COL_IND,
    const int* restrict ROW_PTR,
    const float* restrict B,
    const int N,
    const int K,
    const int VAL_SIZE
);
