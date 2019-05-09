#include <iostream>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <cassert>
#include <cstring>
extern "C"{
	#include <bebop/smc/sparse_matrix.h>
	#include <bebop/smc/sparse_matrix_ops.h>
	#include <bebop/smc/csr_matrix.h>
}
#include "gpu_fpga.h"

void matmul(float *a, float *b, float *c, int N)
{
	int i, j, k;

// <<<dim3(numblock, numblock), dim3(numthread, numthread)>>>
#pragma acc data copyout(c[:N*N]) copyin(a[:N*N], b[:N*N])
#pragma acc kernels
	{
#pragma acc loop gang(N/32) vector(32) independent
		for (i = 0; i < N; ++i)
		{
#pragma acc loop gang(N/32) vector(32) independent
			for (j = 0; j < N; ++j)
			{
				float sum = 0.0;
#pragma acc loop reduction(+:sum)
				for (k = 0; k < N; ++k)
					sum += a[i * N + k] * b[k * N + j];
				c[i * N + j] = sum;
			}
		}
	}
}

void matrix_vector_malti(float *a, float *b, float *c, int N)
{
	int i, j;

// <<<dim3(numblock), dim3(numthread)>>>
#pragma acc data copyout(c[:N]) copyin(a[:N*N], b[:N])
#pragma acc kernels
	{
#pragma acc loop independent gang(N/32) vector(32)
		for (i = 0; i < N; ++i)
		{
			float sum = 0.0;
#pragma acc loop reduction(+:sum)
			for (j = 0; j < N; ++j)
				sum += a[i * N + j] * b[j];
			c[i] = sum;
		}
	}
}

void MatrixMultiplication_openmp(float *a, float *b, float *c, int N)
{
	int i, j, k;
	int chunk;
// #ifdef _OPENMP
	// omp_set_num_threads(numstream);
	if (omp_get_thread_num() == 0)
	{
		printf("Number of OpenMP threads %d\n", omp_get_max_threads());
		chunk = N / omp_get_max_threads();
	}
// #endif

#pragma omp parallel shared(a, b, c, chunk) private(i, j, k)
	{
#pragma omp for
		for (i = 0; i < N; ++i)
		{
			for (j = 0; j < N; ++j)
			{
				float sum = 0.0;
				for (k = 0; k < N; ++k)
					sum += a[i * N + k] * b[k * N + j];
				c[i * N + j] = sum;
			}
		}
	}
}

void h_matrix_vector_malti(float *a, float *b, float *c, int N)
{
	int i, j;
	int chunk;
// #ifdef _OPENMP
	if (omp_get_thread_num() == 0)
	{
		printf("Number of OpenMP threads %d\n", omp_get_max_threads());
		chunk = N / omp_get_max_threads();
	}
// #endif

#pragma omp parallel shared(a, b, c, chunk) private(i, j)
	{
#pragma omp for
		for (i = 0; i < N; ++i)
		{
			float sum = 0.0;
			for (j = 0; j < N; ++j)
				sum += a[i * N + j] * b[j];
			c[i] = sum;
		}
	}
}

void verify_gpu(float *h_c, float *c_CPU, unsigned long N)
{
	double cpu_sum = 0.0;
	double gpu_sum = 0.0;
	double rel_err = 0.0;

#pragma omp parallel for reduction(+:cpu_sum, gpu_sum)
	for (unsigned long i = 0; i < N; ++i)
	{
		// printf("(CPU) %f\n", c_CPU[i]);
		// printf("(GPU) %f\n", h_c[i]);
		cpu_sum += (double)c_CPU[i] * c_CPU[i];
		gpu_sum += (double)h_c[i] * h_c[i];
	}

	cpu_sum = sqrt(cpu_sum);
	gpu_sum = sqrt(gpu_sum);
	if (cpu_sum > gpu_sum)
	{
		rel_err = (cpu_sum - gpu_sum) / cpu_sum;
	}
	else
	{
		rel_err = (gpu_sum - cpu_sum) / cpu_sum;
	}

	if (rel_err < 1e-6)
	{
		printf("Verification Successful err = %e\n", rel_err);
	}
	else
	{
		printf("Verification Fail err = %e\n", rel_err);
	}
	printf("ResultGPU = %lf\n", gpu_sum);
	printf("ResultCPU = %lf\n", cpu_sum);
}

void verify_fpga(
    float* FPGA_calc_result,
    float* VAL,
    int* COL_IND,
    int* ROW_PTR,
    float* B,
    int N,
    int K,
    int VAL_SIZE
	  )
{
	// float *x = new float[N], *r = new float[N], *p = new float[N], *y = new float[N], alfa, beta;
	// float *VAL_local = new float[VAL_SIZE];
	// int *COL_IND_local = new int[VAL_SIZE], *ROW_PTR_local = new int[N + 1];
	// float temp_sum, temp_pap, temp_rr1, temp_rr2;
  int error = N;

	float x[N], r[N], p[N], y[N], alfa, beta;
	float VAL_local[VAL_SIZE];
	int COL_IND_local[VAL_SIZE], ROW_PTR_local[N + 1];
	float temp_sum, temp_pap, temp_rr1, temp_rr2, sum = 0, sum_cpu = 0;

  std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

	temp_rr1 = 0.0f;
	for(int i = 0; i < N; ++i){
		ROW_PTR_local[i] = ROW_PTR[i];
		x[i] = 0.0f;
		r[i] = B[i];
		p[i] = B[i];
		temp_rr1 += r[i] * r[i];
	}
	ROW_PTR_local[N] = ROW_PTR[N];

	for(int i = 0; i < VAL_SIZE; ++i){
		COL_IND_local[i] = COL_IND[i];
		VAL_local[i] = VAL[i];
	}

	for(int i = 0; i < K; ++i){
		temp_pap = 0.0f;
		for(int j = 0; j < N; ++j){
			temp_sum = 0.0f;
			for(int l = ROW_PTR_local[j]; l < ROW_PTR_local[j + 1]; ++l){
				temp_sum += p[COL_IND_local[l]] * VAL_local[l];
			}
			y[j] = temp_sum;
			temp_pap += p[j] * temp_sum;
		}

		alfa = temp_rr1 / temp_pap;

		temp_rr2 = 0.0f;
		for(int j = 0; j < N; ++j){
			x[j] += alfa * p[j];
			r[j] -= alfa * y[j];
			temp_rr2 += r[j] * r[j];
		}

		beta = temp_rr2 / temp_rr1;

		for(int j = 0; j < N; ++j){
			p[j] = r[j] + beta * p[j];
		}
		temp_rr1 = temp_rr2;

	}

  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();

// if (fetestexcept(FE_INVALID)) {
//    puts("浮動小数点例外が発生しました");
// }
// if (fetestexcept(FE_DIVBYZERO)) {
//    puts("ゼロ除算が発生しました");
// }
// if (fetestexcept(FE_OVERFLOW)) {
//    puts("オーバーフローが発生しました");
// }
// if (fetestexcept(FE_UNDERFLOW)) {
//    puts("アンダーフローが発生しました");
// }
// if (fetestexcept(FE_INEXACT)) {
//    puts("不正確な結果が発生しました");
// }

	for(int j = 0; j < N; ++j){
    // std::cout << "FPGA" << FPGA_calc_result[j] << ", CPU"<< x[j] << std::endl;
		if(FPGA_calc_result[j] != x[j]) {
      error = j;
      // break;
    }
    sum += FPGA_calc_result[j];
    sum_cpu += x[j];
	}

  if (error == N) {
    std::cout << std::string(30, '-') << std::endl;
    std::cout << "FPGA Verification: PASS" << std::endl;
    std::cout << "ResultFPGA = " << sum << std::endl;
  } else {
    std::cout << "Error! FPGA Verification failed..." << error << std::endl;
    std::cout << "ResultFPGA = " << sum << std::endl;
    std::cout << "ResultCPU  = " << sum_cpu << std::endl;
   }
  std::cout << "CG CPU elapsed time: " << std::fixed << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << " usec" << std::endl;
}

int main(int argc, char *argv[])
{
	struct sparse_matrix_t* A_ = load_sparse_matrix(MATRIX_MARKET, "bcsstk17.mtx");
	assert(A_ != NULL);
	int errcode = sparse_matrix_convert(A_, CSR);
	if (errcode != 0)
	{
		fprintf(stderr, "*** Conversion failed! ***\n");
		// Note: Don't call destroy_sparse_matrix (A_) unless you 
		// can call free on val, ind and ptr.
		free(A_);
		exit(EXIT_FAILURE);
	}

  struct csr_matrix_t* A = (struct csr_matrix_t*) A_->repr;
  assert (A);
  assert (A->nnz == (A->rowptr[A->m] - A->rowptr[0]));

	// check command line arguments
	///////////////////////////////////////////
	if (argc == 1)
	{
		std::cout << "usage: ./host <numdata_h> <valsize> <numtry>"   << std::endl;
		exit(0);
	}
	if (argc != 4)
	{
		std::cerr << "Error! The number of arguments is wrong."       << std::endl;
		exit(1);
	}

  const int  numdata_h = A->n; // std::stoull(std::string(argv[2]));
	int N = numdata_h;
  const int  valsize   = A->nnz;
	int VAL_SIZE = valsize;
	const int numtry = std::stoull(std::string(argv[3]));
	const unsigned long numbyte = numdata_h * numdata_h * sizeof(float); // this sample uses "float"

	printf("numdata_h: %d, valsize: %d, numtry: %d\n", numdata_h, valsize, numtry);

	// host memory settings
	///////////////////////////////////////////

	/***** GPU *****/
	// static const int numthread = 16;
	// const int numblock = (numdata_h % numthread) ? (numdata_h / numthread) + 1 : (numdata_h / numthread);
	float *h_a, *h_b, *h_c, *c_CPU, *h_vec_b, *h_vec_mul, *vec_b_CPU;

	h_a = (float *)malloc(numbyte);
	h_b = (float *)malloc(numbyte);
	h_c = (float *)malloc(numbyte);
	c_CPU = (float *)malloc(numbyte);
	vec_b_CPU = (float *)malloc(numdata_h * sizeof(float));
	h_vec_mul = (float *)malloc(numdata_h * sizeof(float));
	h_vec_b = (float *)malloc(numdata_h * sizeof(float));

	for (int i = 0; i < numdata_h; ++i)
	{
		for (int j = 0; j < numdata_h; ++j)
		{
			h_a[i * numdata_h + j] = (j + 1) / 2 * 0.0001f;
			h_b[i * numdata_h + j] = 0.5f;
			h_c[i * numdata_h + j] = 0.0f;
			c_CPU[i * numdata_h + j] = 0.0f;
		}
		h_vec_b[i] = 0.0f;
		h_vec_mul[i] = 0.01f;
		vec_b_CPU[i] = 0.0f;
	}

	/***** FPGA *****/
	int K = numtry;
	float *FPGA_calc_result; // X_result;
	float *VAL;
	int *COL_IND;
	int *ROW_PTR;
	float *B;

  posix_memalign((void **)&FPGA_calc_result, 64, N * sizeof(float));
  posix_memalign((void **)&VAL, 64, VAL_SIZE * sizeof(float));
  posix_memalign((void **)&COL_IND, 64, VAL_SIZE * sizeof(int));
  posix_memalign((void **)&ROW_PTR, 64, (N+1) * sizeof(int));
  posix_memalign((void **)&B, 64, N * sizeof(float));

  double *VAL_temp;
  posix_memalign((void **)&VAL_temp, 64, VAL_SIZE * sizeof(double));
   
  memcpy(VAL_temp, A->values, VAL_SIZE * sizeof (double));
  memcpy(COL_IND, A->colidx, VAL_SIZE * sizeof (int));
  memcpy(ROW_PTR, A->rowptr, (N+1) * sizeof (int));
  for (int i = 0; i < VAL_SIZE; ++i)
  {
        VAL[i] = (float)VAL_temp[i];
  }

	// device memory settings
	///////////////////////////////////////////

	// main routine
	///////////////////////////////////////////

	/***** GPU *****/
	std::chrono::system_clock::time_point start_gpu = std::chrono::system_clock::now();

	matmul(h_a, h_b, h_c, numdata_h);
	matrix_vector_malti(h_c, h_vec_mul, h_vec_b, numdata_h);

	std::chrono::system_clock::time_point end_gpu = std::chrono::system_clock::now();

	/***** FPGA *****/
	for (int j = 0; j < N; ++j)
	{
		FPGA_calc_result[j] = 0;
		// ROW_PTR[j] = A->rowptr[j];
		B[j] = h_vec_b[j] - VAL[j] * 1; //000000.0; // b - Ax
	}
	ROW_PTR[N] = N;

	std::chrono::system_clock::time_point start_fpga = std::chrono::system_clock::now();

	funcFPGA(FPGA_calc_result, VAL, COL_IND, ROW_PTR, B, N, K, VAL_SIZE);
	
	std::chrono::system_clock::time_point end_fpga = std::chrono::system_clock::now();

	std::cout << "GPU  elapsed time: " << std::fixed << std::chrono::duration_cast<std::chrono::microseconds>(end_gpu-start_gpu).count() << " usec" << std::endl;
	std::cout << std::string(30, '-') << std::endl;

	std::cout << "FPGA elapsed time: " << std::fixed << std::chrono::duration_cast<std::chrono::microseconds>(end_fpga-start_fpga).count() << " usec" << std::endl;
	std::cout << std::string(30, '-') << std::endl;

	// verification
	///////////////////////////////////////////
	// MatrixMultiplication_openmp(h_a, h_b, c_CPU, numdata_h);    // 本番はコメントアウトして良い
	// h_matrix_vector_malti(c_CPU, h_vec_mul, vec_b_CPU, numdata_h);    // 本番はコメントアウトして良い

	// verify_gpu(h_c, c_CPU, numdata_h*numdata_h); // 行列積チェック
	verify_gpu(h_vec_b, vec_b_CPU, numdata_h); // h_vec_b チェック

	verify_fpga(FPGA_calc_result, VAL, COL_IND, ROW_PTR, B, N, K, VAL_SIZE);

	// cleanup
	///////////////////////////////////////////
  // destroy_sparse_matrix(A_);
	// destroy_csr_matrix(A);
	

	return 0;
}
