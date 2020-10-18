// extern int workers, gangs, size;
// extern float *A, *B, *D, *E;
// #include "openacc.h"

enum DEVICE{ GPU,
FPGA };

void funcFPGA(
    float* X_result,
    float* VAL,
    int* COL_IND,
    int* ROW_PTR,
    float* B,
    int N,
    int K,
    int VAL_SIZE
);

// void initFPGA(
//     float* restrict X_result,
//     float* restrict VAL,
//     int* restrict COL_IND,
//     int* restrict ROW_PTR,
//     float* restrict B,
//     int N,
//     int K,
//     int VAL_SIZE
// );

// void shutdownFPGA(
//     float* restrict X_result,
//     float* restrict VAL,
//     int* restrict COL_IND,
//     int* restrict ROW_PTR,
//     float* restrict B,
//     int N,
//     int K,
//     int VAL_SIZE
// );

void matmul(float *a, float *b, float *c, int N, float *e, float *f);
void MatrixMultiplication_openmp(float *a, float *b, float *c, int N);
void h_matrix_vector_malti(float *a, float *b, float *c, int N);
void verify_gpu(float *h_c, float *c_CPU, unsigned long N);

void verify_fpga(
    float* FPGA_calc_result,
    float* VAL,
    int* COL_IND,
    int* ROW_PTR,
    float* B,
    int N,
    int K,
    int VAL_SIZE
	  );
