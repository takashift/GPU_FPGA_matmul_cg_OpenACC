#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <iostream>
// #include <cstdlib>
// #include <cmath>
#include <omp.h>
#include <string.h>
#include <time.h>
// #include <chrono>
// #include <cassert>
// #include <cstring>
#include "gpu_fpga.h"
// #include <openacc.h> 
#include "/home/tsunashima/OpenARC/include/openacc.h" // Omniのopenacc.hとの干渉防止の為にフルパス指定

#define BLOCK_SIZE 10974 //38000
#define V_SIZE 219812 //220000
#define SPLIT_SIZE 7

float time_diff(struct timespec *start, struct timespec *end) {
	return end->tv_sec - start->tv_sec + (end->tv_nsec - start->tv_nsec) * 1E-9;
}


void funcFPGA(
    float* restrict X_result,
    float* restrict VAL,
    int* restrict COL_IND,
    int* restrict ROW_PTR,
    float* restrict B,
    int N,
    int K,
    int VAL_SIZE
    )
{
#pragma accomn target_dev(FPGA)
// ここに波括弧があるとacc_init()がエラーになる

// acc_init(acc_device_altera_emulator);

int N_init1 = BLOCK_SIZE/SPLIT_SIZE;
int N_init2 = 2*BLOCK_SIZE/SPLIT_SIZE;
int N_init3 = 3*BLOCK_SIZE/SPLIT_SIZE;
int N_init4 = 4*BLOCK_SIZE/SPLIT_SIZE;
int N_init5 = 5*BLOCK_SIZE/SPLIT_SIZE;
int N_init6 = 6*BLOCK_SIZE/SPLIT_SIZE;
int N_end   = 1*BLOCK_SIZE/SPLIT_SIZE + BLOCK_SIZE%SPLIT_SIZE;

// #define N_init1 BLOCK_SIZE/SPLIT_SIZE
// #define N_init2 2*BLOCK_SIZE/SPLIT_SIZE
// #define N_init3 3*BLOCK_SIZE/SPLIT_SIZE
// #define N_init4 4*BLOCK_SIZE/SPLIT_SIZE
// #define N_init5 5*BLOCK_SIZE/SPLIT_SIZE
// #define N_init6 6*BLOCK_SIZE/SPLIT_SIZE
// #define N_end   1*BLOCK_SIZE/SPLIT_SIZE + BLOCK_SIZE%SPLIT_SIZE

#pragma acc enter data create(VAL[0:VAL_SIZE], COL_IND[0:VAL_SIZE], ROW_PTR[0:N+1], B[0:N], N_init1, N_init2, N_init3, N_init4, N_init5, N_init5, N_init6, N_end, N, K, VAL_SIZE) create(X_result[0:N])
// #pragma acc enter data create(VAL[0:VAL_SIZE], COL_IND[0:VAL_SIZE], ROW_PTR[0:N+1], B[0:N], N, K, VAL_SIZE) create(X_result[0:N])

struct timespec fpga_copyin_start;
struct timespec fpga_copyin_end;
struct timespec fpga_calc_start;
struct timespec fpga_calc_end;
struct timespec fpga_copyout_start;
struct timespec fpga_copyout_end;
float pipe_temp_pap_from0to1[1];
float pipe_temp_pap_from0to2[1];
float pipe_temp_pap_from0to3[1];
float pipe_temp_pap_from0to4[1];
float pipe_temp_pap_from0to5[1];
float pipe_temp_pap_from0to6[1];
float pipe_temp_pap_from1to0[1];
float pipe_temp_pap_from1to2[1];
float pipe_temp_pap_from1to3[1];
float pipe_temp_pap_from1to4[1];
float pipe_temp_pap_from1to5[1];
float pipe_temp_pap_from1to6[1];
float pipe_temp_pap_from2to0[1];
float pipe_temp_pap_from2to1[1];
float pipe_temp_pap_from2to3[1];
float pipe_temp_pap_from2to4[1];
float pipe_temp_pap_from2to5[1];
float pipe_temp_pap_from2to6[1];
float pipe_temp_pap_from3to0[1];
float pipe_temp_pap_from3to1[1];
float pipe_temp_pap_from3to2[1];
float pipe_temp_pap_from3to4[1];
float pipe_temp_pap_from3to5[1];
float pipe_temp_pap_from3to6[1];
float pipe_temp_pap_from4to0[1];
float pipe_temp_pap_from4to1[1];
float pipe_temp_pap_from4to2[1];
float pipe_temp_pap_from4to3[1];
float pipe_temp_pap_from4to5[1];
float pipe_temp_pap_from4to6[1];
float pipe_temp_pap_from5to0[1];
float pipe_temp_pap_from5to1[1];
float pipe_temp_pap_from5to2[1];
float pipe_temp_pap_from5to3[1];
float pipe_temp_pap_from5to4[1];
float pipe_temp_pap_from5to6[1];
float pipe_temp_pap_from6to0[1];
float pipe_temp_pap_from6to1[1];
float pipe_temp_pap_from6to2[1];
float pipe_temp_pap_from6to3[1];
float pipe_temp_pap_from6to4[1];
float pipe_temp_pap_from6to5[1];

float pipe_temp_rr2_from0to1[1];
float pipe_temp_rr2_from0to2[1];
float pipe_temp_rr2_from0to3[1];
float pipe_temp_rr2_from0to4[1];
float pipe_temp_rr2_from0to5[1];
float pipe_temp_rr2_from0to6[1];
float pipe_temp_rr2_from1to0[1];
float pipe_temp_rr2_from1to2[1];
float pipe_temp_rr2_from1to3[1];
float pipe_temp_rr2_from1to4[1];
float pipe_temp_rr2_from1to5[1];
float pipe_temp_rr2_from1to6[1];
float pipe_temp_rr2_from2to0[1];
float pipe_temp_rr2_from2to1[1];
float pipe_temp_rr2_from2to3[1];
float pipe_temp_rr2_from2to4[1];
float pipe_temp_rr2_from2to5[1];
float pipe_temp_rr2_from2to6[1];
float pipe_temp_rr2_from3to0[1];
float pipe_temp_rr2_from3to1[1];
float pipe_temp_rr2_from3to2[1];
float pipe_temp_rr2_from3to4[1];
float pipe_temp_rr2_from3to5[1];
float pipe_temp_rr2_from3to6[1];
float pipe_temp_rr2_from4to0[1];
float pipe_temp_rr2_from4to1[1];
float pipe_temp_rr2_from4to2[1];
float pipe_temp_rr2_from4to3[1];
float pipe_temp_rr2_from4to5[1];
float pipe_temp_rr2_from4to6[1];
float pipe_temp_rr2_from5to0[1];
float pipe_temp_rr2_from5to1[1];
float pipe_temp_rr2_from5to2[1];
float pipe_temp_rr2_from5to3[1];
float pipe_temp_rr2_from5to4[1];
float pipe_temp_rr2_from5to6[1];
float pipe_temp_rr2_from6to0[1];
float pipe_temp_rr2_from6to1[1];
float pipe_temp_rr2_from6to2[1];
float pipe_temp_rr2_from6to3[1];
float pipe_temp_rr2_from6to4[1];
float pipe_temp_rr2_from6to5[1];

float pipe_p_from0to1[1];
float pipe_p_from0to2[1];
float pipe_p_from0to3[1];
float pipe_p_from0to4[1];
float pipe_p_from0to5[1];
float pipe_p_from0to6[1];
float pipe_p_from1to0[1];
float pipe_p_from1to2[1];
float pipe_p_from1to3[1];
float pipe_p_from1to4[1];
float pipe_p_from1to5[1];
float pipe_p_from1to6[1];
float pipe_p_from2to0[1];
float pipe_p_from2to1[1];
float pipe_p_from2to3[1];
float pipe_p_from2to4[1];
float pipe_p_from2to5[1];
float pipe_p_from2to6[1];
float pipe_p_from3to0[1];
float pipe_p_from3to1[1];
float pipe_p_from3to2[1];
float pipe_p_from3to4[1];
float pipe_p_from3to5[1];
float pipe_p_from3to6[1];
float pipe_p_from4to0[1];
float pipe_p_from4to1[1];
float pipe_p_from4to2[1];
float pipe_p_from4to3[1];
float pipe_p_from4to5[1];
float pipe_p_from4to6[1];
float pipe_p_from5to0[1];
float pipe_p_from5to1[1];
float pipe_p_from5to2[1];
float pipe_p_from5to3[1];
float pipe_p_from5to4[1];
float pipe_p_from5to6[1];
float pipe_p_from6to0[1];
float pipe_p_from6to1[1];
float pipe_p_from6to2[1];
float pipe_p_from6to3[1];
float pipe_p_from6to4[1];
float pipe_p_from6to5[1];

clock_gettime(CLOCK_REALTIME, &fpga_copyin_start);
#pragma acc update device(VAL[0:VAL_SIZE], COL_IND[0:VAL_SIZE], ROW_PTR[0:N+1], B[0:N], N_init1, N_init2, N_init3, N_init4, N_init5, N_init5, N_init6, N_end, N, K, VAL_SIZE)
// #pragma acc update device(VAL[0:VAL_SIZE], COL_IND[0:VAL_SIZE], ROW_PTR[0:N+1], B[0:N], N, K, VAL_SIZE)
clock_gettime(CLOCK_REALTIME, &fpga_copyin_end);
printf("Host to FPGA time: %lf sec\n", time_diff(&fpga_copyin_start, &fpga_copyin_end));


#pragma acc data present(VAL[0:VAL_SIZE], COL_IND[0:VAL_SIZE], ROW_PTR[0:N+1], B[0:N], N_init1, N_init2, N_init3, N_init4, N_init5, N_init5, N_init6, N_end, N, K, VAL_SIZE, X_result[0:N]) pipe(pipe_temp_pap_from0to1, pipe_temp_pap_from0to2, pipe_temp_pap_from0to3, pipe_temp_pap_from0to4, pipe_temp_pap_from0to5, pipe_temp_pap_from0to6, pipe_temp_rr2_from0to1, pipe_temp_rr2_from0to2, pipe_temp_rr2_from0to3, pipe_temp_rr2_from0to4, pipe_temp_rr2_from0to5, pipe_temp_rr2_from0to6, pipe_p_from0to1, pipe_p_from0to2, pipe_p_from0to3, pipe_p_from0to4, pipe_p_from0to5, pipe_p_from0to6, pipe_temp_pap_from1to0, pipe_temp_pap_from1to2, pipe_temp_pap_from1to3, pipe_temp_pap_from1to4, pipe_temp_pap_from1to5, pipe_temp_pap_from1to6, pipe_temp_rr2_from1to0, pipe_temp_rr2_from1to2, pipe_temp_rr2_from1to3, pipe_temp_rr2_from1to4, pipe_temp_rr2_from1to5, pipe_temp_rr2_from1to6, pipe_p_from1to0, pipe_p_from1to2, pipe_p_from1to3, pipe_p_from1to4, pipe_p_from1to5, pipe_p_from1to6, pipe_temp_pap_from2to0, pipe_temp_pap_from2to1, pipe_temp_pap_from2to3, pipe_temp_pap_from2to4, pipe_temp_pap_from2to5, pipe_temp_pap_from2to6, pipe_temp_rr2_from2to0, pipe_temp_rr2_from2to1, pipe_temp_rr2_from2to3, pipe_temp_rr2_from2to4, pipe_temp_rr2_from2to5, pipe_temp_rr2_from2to6, pipe_p_from2to0, pipe_p_from2to1, pipe_p_from2to3, pipe_p_from2to4, pipe_p_from2to5, pipe_p_from2to6, pipe_temp_pap_from3to0, pipe_temp_pap_from3to1, pipe_temp_pap_from3to2, pipe_temp_pap_from3to4, pipe_temp_pap_from3to5, pipe_temp_pap_from3to6, pipe_temp_rr2_from3to0, pipe_temp_rr2_from3to1, pipe_temp_rr2_from3to2, pipe_temp_rr2_from3to4, pipe_temp_rr2_from3to5, pipe_temp_rr2_from3to6, pipe_p_from3to0, pipe_p_from3to1, pipe_p_from3to2, pipe_p_from3to4, pipe_p_from3to5, pipe_p_from3to6, pipe_p_from3to6, pipe_temp_pap_from4to0, pipe_temp_pap_from4to1, pipe_temp_pap_from4to2, pipe_temp_pap_from4to3, pipe_temp_pap_from4to5, pipe_temp_pap_from4to6, pipe_temp_rr2_from4to0, pipe_temp_rr2_from4to1, pipe_temp_rr2_from4to2, pipe_temp_rr2_from4to3, pipe_temp_rr2_from4to5, pipe_temp_rr2_from4to6, pipe_p_from4to0, pipe_p_from4to1, pipe_p_from4to2, pipe_p_from4to3, pipe_p_from4to5, pipe_p_from4to6, pipe_temp_pap_from5to0, pipe_temp_pap_from5to1, pipe_temp_pap_from5to2, pipe_temp_pap_from5to3, pipe_temp_pap_from5to4, pipe_temp_pap_from5to6, pipe_temp_rr2_from5to0, pipe_temp_rr2_from5to1, pipe_temp_rr2_from5to2, pipe_temp_rr2_from5to3, pipe_temp_rr2_from5to4, pipe_temp_rr2_from5to6, pipe_p_from5to0, pipe_p_from5to1, pipe_p_from5to2, pipe_p_from5to3, pipe_p_from5to4, pipe_p_from5to6, pipe_temp_pap_from6to0, pipe_temp_pap_from6to1, pipe_temp_pap_from6to2, pipe_temp_pap_from6to3, pipe_temp_pap_from6to4, pipe_temp_pap_from6to5, pipe_temp_rr2_from6to0, pipe_temp_rr2_from6to1, pipe_temp_rr2_from6to2, pipe_temp_rr2_from6to3, pipe_temp_rr2_from6to4, pipe_temp_rr2_from6to5, pipe_p_from6to0, pipe_p_from6to1, pipe_p_from6to2, pipe_p_from6to3, pipe_p_from6to4, pipe_p_from6to5)
// #pragma acc data present(VAL[0:VAL_SIZE], COL_IND[0:VAL_SIZE], ROW_PTR[0:N+1], B[0:N], N, K, VAL_SIZE, X_result[0:N]) pipe(pipe_temp_pap_from0to1, pipe_temp_pap_from0to2, pipe_temp_pap_from0to3, pipe_temp_pap_from0to4, pipe_temp_pap_from0to5, pipe_temp_pap_from0to6, pipe_temp_rr2_from0to1, pipe_temp_rr2_from0to2, pipe_temp_rr2_from0to3, pipe_temp_rr2_from0to4, pipe_temp_rr2_from0to5, pipe_temp_rr2_from0to6, pipe_p_from0to1, pipe_p_from0to2, pipe_p_from0to3, pipe_p_from0to4, pipe_p_from0to5, pipe_p_from0to6, pipe_temp_pap_from1to0, pipe_temp_pap_from1to2, pipe_temp_pap_from1to3, pipe_temp_pap_from1to4, pipe_temp_pap_from1to5, pipe_temp_pap_from1to6, pipe_temp_rr2_from1to0, pipe_temp_rr2_from1to2, pipe_temp_rr2_from1to3, pipe_temp_rr2_from1to4, pipe_temp_rr2_from1to5, pipe_temp_rr2_from1to6, pipe_p_from1to0, pipe_p_from1to2, pipe_p_from1to3, pipe_p_from1to4, pipe_p_from1to5, pipe_p_from1to6, pipe_temp_pap_from2to0, pipe_temp_pap_from2to1, pipe_temp_pap_from2to3, pipe_temp_pap_from2to4, pipe_temp_pap_from2to5, pipe_temp_pap_from2to6, pipe_temp_rr2_from2to0, pipe_temp_rr2_from2to1, pipe_temp_rr2_from2to3, pipe_temp_rr2_from2to4, pipe_temp_rr2_from2to5, pipe_temp_rr2_from2to6, pipe_p_from2to0, pipe_p_from2to1, pipe_p_from2to3, pipe_p_from2to4, pipe_p_from2to5, pipe_p_from2to6, pipe_temp_pap_from3to0, pipe_temp_pap_from3to1, pipe_temp_pap_from3to2, pipe_temp_pap_from3to4, pipe_temp_pap_from3to5, pipe_temp_pap_from3to6, pipe_temp_rr2_from3to0, pipe_temp_rr2_from3to1, pipe_temp_rr2_from3to2, pipe_temp_rr2_from3to4, pipe_temp_rr2_from3to5, pipe_temp_rr2_from3to6, pipe_p_from3to0, pipe_p_from3to1, pipe_p_from3to2, pipe_p_from3to4, pipe_p_from3to5, pipe_p_from3to6, pipe_p_from3to6, pipe_temp_pap_from4to0, pipe_temp_pap_from4to1, pipe_temp_pap_from4to2, pipe_temp_pap_from4to3, pipe_temp_pap_from4to5, pipe_temp_pap_from4to6, pipe_temp_rr2_from4to0, pipe_temp_rr2_from4to1, pipe_temp_rr2_from4to2, pipe_temp_rr2_from4to3, pipe_temp_rr2_from4to5, pipe_temp_rr2_from4to6, pipe_p_from4to0, pipe_p_from4to1, pipe_p_from4to2, pipe_p_from4to3, pipe_p_from4to5, pipe_p_from4to6, pipe_temp_pap_from5to0, pipe_temp_pap_from5to1, pipe_temp_pap_from5to2, pipe_temp_pap_from5to3, pipe_temp_pap_from5to4, pipe_temp_pap_from5to6, pipe_temp_rr2_from5to0, pipe_temp_rr2_from5to1, pipe_temp_rr2_from5to2, pipe_temp_rr2_from5to3, pipe_temp_rr2_from5to4, pipe_temp_rr2_from5to6, pipe_p_from5to0, pipe_p_from5to1, pipe_p_from5to2, pipe_p_from5to3, pipe_p_from5to4, pipe_p_from5to6, pipe_temp_pap_from6to0, pipe_temp_pap_from6to1, pipe_temp_pap_from6to2, pipe_temp_pap_from6to3, pipe_temp_pap_from6to4, pipe_temp_pap_from6to5, pipe_temp_rr2_from6to0, pipe_temp_rr2_from6to1, pipe_temp_rr2_from6to2, pipe_temp_rr2_from6to3, pipe_temp_rr2_from6to4, pipe_temp_rr2_from6to5, pipe_p_from6to0, pipe_p_from6to1, pipe_p_from6to2, pipe_p_from6to3, pipe_p_from6to4, pipe_p_from6to5)
{
clock_gettime(CLOCK_REALTIME, &fpga_calc_start);



// calc_CG_FPGA_slave(X_result, VAL, COL_IND, ROW_PTR, B, N/SPLIT_SIZE, N, K, VAL_SIZE);
#pragma acc serial async(0) pipein(pipe_temp_pap_from1to0, pipe_temp_rr2_from1to0, pipe_p_from1to0, pipe_temp_pap_from2to0, pipe_temp_rr2_from2to0, pipe_p_from2to0, pipe_temp_pap_from3to0, pipe_temp_rr2_from3to0, pipe_p_from3to0, pipe_temp_pap_from4to0, pipe_temp_rr2_from4to0, pipe_p_from4to0, pipe_temp_pap_from5to0, pipe_temp_rr2_from5to0, pipe_p_from5to0, pipe_temp_pap_from6to0, pipe_temp_rr2_from6to0, pipe_p_from6to0) pipeout(pipe_temp_pap_from0to1, pipe_temp_pap_from0to2, pipe_temp_pap_from0to3, pipe_temp_pap_from0to4, pipe_temp_pap_from0to5, pipe_temp_pap_from0to6, pipe_temp_rr2_from0to1, pipe_temp_rr2_from0to2, pipe_temp_rr2_from0to3, pipe_temp_rr2_from0to4, pipe_temp_rr2_from0to5, pipe_temp_rr2_from0to6, pipe_p_from0to1, pipe_p_from0to2, pipe_p_from0to3, pipe_p_from0to4, pipe_p_from0to5, pipe_p_from0to6)
{
	// V_SIZEの２つはエミュレーションのときはserialディレクティブの外に出す
	float VAL_local[V_SIZE];
	int COL_IND_local[V_SIZE];
	float x[BLOCK_SIZE/SPLIT_SIZE], r[BLOCK_SIZE/SPLIT_SIZE], p[BLOCK_SIZE], y[BLOCK_SIZE/SPLIT_SIZE];
	// pはBLOCK_SIZEのまま、VAL_localとCOL_IND_localも多分行ごとに数が異なるからそのまま
	float alfa, beta;//, x[BLOCK_SIZE], r[BLOCK_SIZE], p[BLOCK_SIZE], y[BLOCK_SIZE];
	int ROW_PTR_local[(BLOCK_SIZE + 1)];
	float temp_sum=0.0f, temp_pap, temp_rr1, temp_rr2;

	temp_rr1 = 0.0f;
#pragma acc loop independent //reduction(+:temp_rr1)
	for(int i = 0; i < N; ++i){
		ROW_PTR_local[i] = ROW_PTR[i];
		p[i] = B[i];
		temp_rr1 += p[i] * p[i];
	}

	for(int i = 0; i < N_init1; ++i) {
		x[i] = 0.0f;
		r[i] = B[i];
	}

	ROW_PTR_local[N_init1] = ROW_PTR[N_init1];

#pragma acc loop independent
	for(int i = 0; i < VAL_SIZE; ++i){
		COL_IND_local[i] = COL_IND[i];
		VAL_local[i] = VAL[i];
	}


	for(int i = 0; i < K; ++i){
		temp_pap = 0.0f;
#pragma acc loop reduction(+:temp_sum)
		for(int j = 0; j < N_init1; ++j){
			temp_sum = 0.0f;
			for(int l = ROW_PTR_local[j]; l < ROW_PTR_local[j + 1]; ++l){
				temp_sum += VAL_local[l] * p[COL_IND_local[l]];
			}
			y[j] = temp_sum;
			temp_pap += p[j] * temp_sum;
		}

// #pragma acc loop //reduction(+:temp_pap)
// 		for(int j = 0; j < N_init1; ++j){
// 			temp_pap += p[j] * y[j];
// 		}

// チャネル必要(temp_pap送信)
pipe_temp_pap_from0to1[0] = temp_pap;
pipe_temp_pap_from0to2[0] = temp_pap;
pipe_temp_pap_from0to3[0] = temp_pap;
pipe_temp_pap_from0to4[0] = temp_pap;
pipe_temp_pap_from0to5[0] = temp_pap;
pipe_temp_pap_from0to6[0] = temp_pap;
// チャネル必要(temp_papを合計)
temp_pap += pipe_temp_pap_from1to0[0];
temp_pap += pipe_temp_pap_from2to0[0];
temp_pap += pipe_temp_pap_from3to0[0];
temp_pap += pipe_temp_pap_from4to0[0];
temp_pap += pipe_temp_pap_from5to0[0];
temp_pap += pipe_temp_pap_from6to0[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		alfa = temp_rr1 / temp_pap;

		temp_rr2 = 0.0f;
#pragma acc loop reduction(+:temp_rr2)
		for(int j = 0; j < N_init1; ++j){
			x[j] += alfa * p[j];
			r[j] -= alfa * y[j];
			temp_rr2 += r[j] * r[j];
		}
// チャネル必要(temp_rr2送信)
pipe_temp_rr2_from0to1[0] = temp_rr2;
pipe_temp_rr2_from0to2[0] = temp_rr2;
pipe_temp_rr2_from0to3[0] = temp_rr2;
pipe_temp_rr2_from0to4[0] = temp_rr2;
pipe_temp_rr2_from0to5[0] = temp_rr2;
pipe_temp_rr2_from0to6[0] = temp_rr2;
// チャネル必要(temp_rr2を合計)
temp_rr2 += pipe_temp_rr2_from1to0[0];
temp_rr2 += pipe_temp_rr2_from2to0[0];
temp_rr2 += pipe_temp_rr2_from3to0[0];
temp_rr2 += pipe_temp_rr2_from4to0[0];
temp_rr2 += pipe_temp_rr2_from5to0[0];
temp_rr2 += pipe_temp_rr2_from6to0[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		beta = temp_rr2 / temp_rr1;

#pragma acc loop independent
		for(int j = 0; j < N_init1; ++j){
			p[j] = r[j] + beta * p[j];
		}

// カーネルごとの送信の順番を揃える必要がある
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	pipe_p_from0to1[0] = p[j];
	pipe_p_from0to2[0] = p[j];
	pipe_p_from0to3[0] = p[j];
	pipe_p_from0to4[0] = p[j];
	pipe_p_from0to5[0] = p[j];
	pipe_p_from0to6[0] = p[j];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init1] = pipe_p_from1to0[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init2] = pipe_p_from2to0[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init3] = pipe_p_from3to0[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init4] = pipe_p_from4to0[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init5] = pipe_p_from5to0[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_end; ++j){
	p[j+N_init6] = pipe_p_from6to0[0];
}

		temp_rr1 = temp_rr2;
	}

#pragma acc loop independent
	for(int j = 0; j < N_init1; ++j){
		X_result[j] = x[j];
	}
}














#pragma acc serial async(1) pipein(pipe_temp_pap_from0to1, pipe_temp_rr2_from0to1, pipe_p_from0to1, pipe_temp_pap_from2to1, pipe_temp_rr2_from2to1, pipe_p_from2to1, pipe_temp_pap_from3to1, pipe_temp_rr2_from3to1, pipe_p_from3to1, pipe_temp_pap_from4to1, pipe_temp_rr2_from4to1, pipe_p_from4to1, pipe_temp_pap_from5to1, pipe_temp_rr2_from5to1, pipe_p_from5to1, pipe_temp_pap_from6to1, pipe_temp_rr2_from6to1, pipe_p_from6to1) pipeout(pipe_temp_pap_from1to0, pipe_temp_pap_from1to2, pipe_temp_pap_from1to3, pipe_temp_pap_from1to4, pipe_temp_pap_from1to5, pipe_temp_pap_from1to6, pipe_temp_rr2_from1to0, pipe_temp_rr2_from1to2, pipe_temp_rr2_from1to3, pipe_temp_rr2_from1to4, pipe_temp_rr2_from1to5, pipe_temp_rr2_from1to6, pipe_p_from1to0, pipe_p_from1to2, pipe_p_from1to3, pipe_p_from1to4, pipe_p_from1to5, pipe_p_from1to6)
{
	float alfa, beta;
	float x[BLOCK_SIZE/SPLIT_SIZE], r[BLOCK_SIZE/SPLIT_SIZE], p[BLOCK_SIZE], y[BLOCK_SIZE/SPLIT_SIZE];
	float VAL_local[V_SIZE];    // エミュレーションのときはコメントアウト
	int COL_IND_local[V_SIZE];  // エミュレーションのときはコメントアウト
	int ROW_PTR_local[(BLOCK_SIZE + 1)];
	float temp_sum=0.0f, temp_pap, temp_rr1, temp_rr2;

	temp_rr1 = 0.0f;
#pragma acc loop independent //reduction(+:temp_rr1)
	for(int i = 0; i < N; ++i){
		ROW_PTR_local[i] = ROW_PTR[i];
		p[i] = B[i];
		temp_rr1 += p[i] * p[i];
	}

	for(int i = 0; i < N_init1; ++i) {
		x[i] = 0.0f;
		r[i] = B[i+N_init1];
	}

	ROW_PTR_local[N] = ROW_PTR[N];

#pragma acc loop independent
	for(int i = 0; i < VAL_SIZE; ++i){
		COL_IND_local[i] = COL_IND[i];
		VAL_local[i] = VAL[i];
	}

	for(int i = 0; i < K; ++i){
		temp_pap = 0.0f;
#pragma acc loop reduction(+:temp_sum)
		for(int j = 0; j < N_init1; ++j){
			temp_sum = 0.0f;
			for(int l = ROW_PTR_local[j+N_init1]; l < ROW_PTR_local[j + N_init1 + 1]; ++l){
				temp_sum += VAL_local[l] * p[COL_IND_local[l]];
			}
			y[j] = temp_sum;
			temp_pap += p[j+N_init1] * temp_sum;
		}

// #pragma acc loop //reduction(+:temp_pap)
// 		for(int j = N_init1; j < N_init2; ++j){
// 			temp_pap += p[j] * y[j];
// 		}

// チャネル必要(temp_pap送信)
pipe_temp_pap_from1to0[0] = temp_pap;
pipe_temp_pap_from1to2[0] = temp_pap;
pipe_temp_pap_from1to3[0] = temp_pap;
pipe_temp_pap_from1to4[0] = temp_pap;
pipe_temp_pap_from1to5[0] = temp_pap;
pipe_temp_pap_from1to6[0] = temp_pap;
// チャネル必要(合計を受信)
temp_pap += pipe_temp_pap_from0to1[0];
temp_pap += pipe_temp_pap_from2to1[0];
temp_pap += pipe_temp_pap_from3to1[0];
temp_pap += pipe_temp_pap_from4to1[0];
temp_pap += pipe_temp_pap_from5to1[0];
temp_pap += pipe_temp_pap_from6to1[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		alfa = temp_rr1 / temp_pap;

		temp_rr2 = 0.0f;
#pragma acc loop reduction(+:temp_rr2)
		for(int j = 0; j < N_init1; ++j){
			x[j] += alfa * p[j+N_init1];
			r[j] -= alfa * y[j];
			temp_rr2 += r[j] * r[j];
		}
// チャネル必要(temp_rr2送信)
pipe_temp_rr2_from1to0[0] = temp_rr2;
pipe_temp_rr2_from1to2[0] = temp_rr2;
pipe_temp_rr2_from1to3[0] = temp_rr2;
pipe_temp_rr2_from1to4[0] = temp_rr2;
pipe_temp_rr2_from1to5[0] = temp_rr2;
pipe_temp_rr2_from1to6[0] = temp_rr2;
// チャネル必要(合計を受信)
temp_rr2 += pipe_temp_rr2_from0to1[0];
temp_rr2 += pipe_temp_rr2_from2to1[0];
temp_rr2 += pipe_temp_rr2_from3to1[0];
temp_rr2 += pipe_temp_rr2_from4to1[0];
temp_rr2 += pipe_temp_rr2_from5to1[0];
temp_rr2 += pipe_temp_rr2_from6to1[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		beta = temp_rr2 / temp_rr1;

// #pragma acc loop independent
		for(int j = 0; j < N_init1; ++j){
			p[j+N_init1] = r[j] + beta * p[j+N_init1];
		}

// カーネルごとの送信の順番を揃える必要がある
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j] = pipe_p_from0to1[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	pipe_p_from1to0[0] = p[j+N_init1];
	pipe_p_from1to2[0] = p[j+N_init1];
	pipe_p_from1to3[0] = p[j+N_init1];
	pipe_p_from1to4[0] = p[j+N_init1];
	pipe_p_from1to5[0] = p[j+N_init1];
	pipe_p_from1to6[0] = p[j+N_init1];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init2] = pipe_p_from2to1[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init3] = pipe_p_from3to1[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init4] = pipe_p_from4to1[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init5] = pipe_p_from5to1[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_end; ++j){
	p[j+N_init6] = pipe_p_from6to1[0];
}

		temp_rr1 = temp_rr2;
	}

// #pragma acc loop independent
	for(int j = 0; j < N_init1; ++j){
		X_result[j+N_init1] = x[j];
	}
}












#pragma acc serial async(2) pipein(pipe_temp_pap_from0to2, pipe_temp_rr2_from0to2, pipe_p_from0to2, pipe_temp_pap_from1to2, pipe_temp_rr2_from1to2, pipe_p_from1to2, pipe_temp_pap_from3to2, pipe_temp_rr2_from3to2, pipe_p_from3to2, pipe_temp_pap_from4to2, pipe_temp_rr2_from4to2, pipe_p_from4to2, pipe_temp_pap_from5to2, pipe_temp_rr2_from5to2, pipe_p_from5to2, pipe_temp_pap_from6to2, pipe_temp_rr2_from6to2, pipe_p_from6to2) pipeout(pipe_temp_pap_from2to0, pipe_temp_pap_from2to1, pipe_temp_pap_from2to3, pipe_temp_pap_from2to4, pipe_temp_pap_from2to5, pipe_temp_pap_from2to6, pipe_temp_rr2_from2to0, pipe_temp_rr2_from2to1, pipe_temp_rr2_from2to3, pipe_temp_rr2_from2to4, pipe_temp_rr2_from2to5, pipe_temp_rr2_from2to6, pipe_p_from2to0, pipe_p_from2to1, pipe_p_from2to3, pipe_p_from2to4, pipe_p_from2to5, pipe_p_from2to6)
{
	float alfa, beta;
	float x[BLOCK_SIZE/SPLIT_SIZE], r[BLOCK_SIZE/SPLIT_SIZE], p[BLOCK_SIZE], y[BLOCK_SIZE/SPLIT_SIZE];
	float VAL_local[V_SIZE];    // エミュレーションのときはコメントアウト
	int COL_IND_local[V_SIZE];  // エミュレーションのときはコメントアウト
	int ROW_PTR_local[(BLOCK_SIZE + 1)];
	float temp_sum=0.0f, temp_pap, temp_rr1, temp_rr2;

	temp_rr1 = 0.0f;
#pragma acc loop independent //reduction(+:temp_rr1)
	for(int i = 0; i < N; ++i){
		ROW_PTR_local[i] = ROW_PTR[i];
		p[i] = B[i];
		temp_rr1 += p[i] * p[i];
	}

	for(int i = 0; i < N_init1; ++i) {
		x[i] = 0.0f;
		r[i] = B[i+N_init2];
	}

	ROW_PTR_local[N] = ROW_PTR[N];

#pragma acc loop independent
	for(int i = 0; i < VAL_SIZE; ++i){
		COL_IND_local[i] = COL_IND[i];
		VAL_local[i] = VAL[i];
	}

	for(int i = 0; i < K; ++i){
		temp_pap = 0.0f;
#pragma acc loop reduction(+:temp_sum)
		for(int j = 0; j < N_init1; ++j){
			temp_sum = 0.0f;
			for(int l = ROW_PTR_local[j+N_init2]; l < ROW_PTR_local[j +N_init2+ 1]; ++l){
				temp_sum += VAL_local[l] * p[COL_IND_local[l]];
			}
			y[j] = temp_sum;
			temp_pap += p[j+N_init2] * temp_sum;
		}

// #pragma acc loop //reduction(+:temp_pap)
// 		for(int j = N_init2; j < N_init3; ++j){
// 			temp_pap += p[j] * y[j];
// 		}

// チャネル必要(temp_pap送信)
pipe_temp_pap_from2to0[0] = temp_pap;
pipe_temp_pap_from2to1[0] = temp_pap;
pipe_temp_pap_from2to3[0] = temp_pap;
pipe_temp_pap_from2to4[0] = temp_pap;
pipe_temp_pap_from2to5[0] = temp_pap;
pipe_temp_pap_from2to6[0] = temp_pap;
// チャネル必要(合計を受信)
temp_pap += pipe_temp_pap_from0to2[0];
temp_pap += pipe_temp_pap_from1to2[0];
temp_pap += pipe_temp_pap_from3to2[0];
temp_pap += pipe_temp_pap_from4to2[0];
temp_pap += pipe_temp_pap_from5to2[0];
temp_pap += pipe_temp_pap_from6to2[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		alfa = temp_rr1 / temp_pap;

		temp_rr2 = 0.0f;
#pragma acc loop reduction(+:temp_rr2)
		for(int j = 0; j < N_init1; ++j){
			x[j] += alfa * p[j+N_init2];
			r[j] -= alfa * y[j];
			temp_rr2 += r[j] * r[j];
		}
// チャネル必要(temp_rr2送信)
pipe_temp_rr2_from2to0[0] = temp_rr2;
pipe_temp_rr2_from2to1[0] = temp_rr2;
pipe_temp_rr2_from2to3[0] = temp_rr2;
pipe_temp_rr2_from2to4[0] = temp_rr2;
pipe_temp_rr2_from2to5[0] = temp_rr2;
pipe_temp_rr2_from2to6[0] = temp_rr2;
// チャネル必要(合計を受信)
temp_rr2 += pipe_temp_rr2_from0to2[0];
temp_rr2 += pipe_temp_rr2_from1to2[0];
temp_rr2 += pipe_temp_rr2_from3to2[0];
temp_rr2 += pipe_temp_rr2_from4to2[0];
temp_rr2 += pipe_temp_rr2_from5to2[0];
temp_rr2 += pipe_temp_rr2_from6to2[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		beta = temp_rr2 / temp_rr1;

// #pragma acc loop independent
		for(int j = 0; j < N_init1; ++j){
			p[j+N_init2] = r[j] + beta * p[j+N_init2];
		}

// カーネルごとの送信の順番を揃える必要がある
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j] = pipe_p_from0to2[0];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init1] = pipe_p_from1to2[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	pipe_p_from2to0[0] = p[j+N_init2];
	pipe_p_from2to1[0] = p[j+N_init2];
	pipe_p_from2to3[0] = p[j+N_init2];
	pipe_p_from2to4[0] = p[j+N_init2];
	pipe_p_from2to5[0] = p[j+N_init2];
	pipe_p_from2to6[0] = p[j+N_init2];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init3] = pipe_p_from3to2[0];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init4] = pipe_p_from4to2[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init5] = pipe_p_from5to2[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_end; ++j){
	p[j+N_init6] = pipe_p_from6to2[0];
}

		temp_rr1 = temp_rr2;
	}

// #pragma acc loop independent
	for(int j = 0; j < N_init1; ++j){
		X_result[j+N_init2] = x[j];
	}
}












#pragma acc serial async(3) pipein(pipe_temp_pap_from0to3, pipe_temp_rr2_from0to3, pipe_p_from0to3, pipe_temp_pap_from1to3, pipe_temp_rr2_from1to3, pipe_p_from1to3, pipe_temp_pap_from2to3, pipe_temp_rr2_from2to3, pipe_p_from2to3, pipe_temp_pap_from4to3, pipe_temp_rr2_from4to3, pipe_p_from4to3, pipe_temp_pap_from5to3, pipe_temp_rr2_from5to3, pipe_p_from5to3, pipe_temp_pap_from6to3, pipe_temp_rr2_from6to3, pipe_p_from6to3) pipeout(pipe_temp_pap_from3to0, pipe_temp_pap_from3to1, pipe_temp_pap_from3to2, pipe_temp_pap_from3to4, pipe_temp_pap_from3to5, pipe_temp_pap_from3to6, pipe_temp_rr2_from3to0, pipe_temp_rr2_from3to1, pipe_temp_rr2_from3to2, pipe_temp_rr2_from3to4, pipe_temp_rr2_from3to5, pipe_temp_rr2_from3to6, pipe_p_from3to0, pipe_p_from3to1, pipe_p_from3to2, pipe_p_from3to4, pipe_p_from3to5, pipe_p_from3to6)
{
	float alfa, beta;
	float x[BLOCK_SIZE/SPLIT_SIZE], r[BLOCK_SIZE/SPLIT_SIZE], p[BLOCK_SIZE], y[BLOCK_SIZE/SPLIT_SIZE];
	float VAL_local[V_SIZE];    // エミュレーションのときはコメントアウト
	int COL_IND_local[V_SIZE];  // エミュレーションのときはコメントアウト
	int ROW_PTR_local[(BLOCK_SIZE + 1)];
	float temp_sum=0.0f, temp_pap, temp_rr1, temp_rr2;

	temp_rr1 = 0.0f;
#pragma acc loop independent //reduction(+:temp_rr1)
	for(int i = 0; i < N; ++i){
		ROW_PTR_local[i] = ROW_PTR[i];
		p[i] = B[i];
		temp_rr1 += p[i] * p[i];
	}

	for(int i = 0; i < N_init1; ++i) {
		x[i] = 0.0f;
		r[i] = B[i+N_init3];
	}

	ROW_PTR_local[N] = ROW_PTR[N];

#pragma acc loop independent
	for(int i = 0; i < VAL_SIZE; ++i){
		COL_IND_local[i] = COL_IND[i];
		VAL_local[i] = VAL[i];
	}

	for(int i = 0; i < K; ++i){
		temp_pap = 0.0f;
#pragma acc loop reduction(+:temp_sum)
		for(int j = 0; j < N_init1; ++j){
			temp_sum = 0.0f;
			for(int l = ROW_PTR_local[j+N_init3]; l < ROW_PTR_local[j +N_init3+ 1]; ++l){
				temp_sum += VAL_local[l] * p[COL_IND_local[l]];
			}
			y[j] = temp_sum;
			temp_pap += p[j+N_init3] * temp_sum;
		}

// #pragma acc loop //reduction(+:temp_pap)
// 		for(int j = N_init3; j < N_init4; ++j){
// 			temp_pap += p[j] * y[j];
// 		}

// チャネル必要(temp_pap送信)
pipe_temp_pap_from3to0[0] = temp_pap;
pipe_temp_pap_from3to1[0] = temp_pap;
pipe_temp_pap_from3to2[0] = temp_pap;
pipe_temp_pap_from3to4[0] = temp_pap;
pipe_temp_pap_from3to5[0] = temp_pap;
pipe_temp_pap_from3to6[0] = temp_pap;
// チャネル必要(合計を受信)
temp_pap += pipe_temp_pap_from0to3[0];
temp_pap += pipe_temp_pap_from1to3[0];
temp_pap += pipe_temp_pap_from2to3[0];
temp_pap += pipe_temp_pap_from4to3[0];
temp_pap += pipe_temp_pap_from5to3[0];
temp_pap += pipe_temp_pap_from6to3[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		alfa = temp_rr1 / temp_pap;

		temp_rr2 = 0.0f;
#pragma acc loop reduction(+:temp_rr2)
		for(int j = 0; j < N_init1; ++j){
			x[j] += alfa * p[j+N_init3];
			r[j] -= alfa * y[j];
			temp_rr2 += r[j] * r[j];
		}
// チャネル必要(temp_rr2送信)
pipe_temp_rr2_from3to0[0] = temp_rr2;
pipe_temp_rr2_from3to1[0] = temp_rr2;
pipe_temp_rr2_from3to2[0] = temp_rr2;
pipe_temp_rr2_from3to4[0] = temp_rr2;
pipe_temp_rr2_from3to5[0] = temp_rr2;
pipe_temp_rr2_from3to6[0] = temp_rr2;
// チャネル必要(合計を受信)
temp_rr2 += pipe_temp_rr2_from0to3[0];
temp_rr2 += pipe_temp_rr2_from1to3[0];
temp_rr2 += pipe_temp_rr2_from2to3[0];
temp_rr2 += pipe_temp_rr2_from4to3[0];
temp_rr2 += pipe_temp_rr2_from5to3[0];
temp_rr2 += pipe_temp_rr2_from6to3[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		beta = temp_rr2 / temp_rr1;

// #pragma acc loop independent
		for(int j = 0; j < N_init1; ++j){
			p[j+N_init3] = r[j] + beta * p[j+N_init3];
		}

// カーネルごとの送信の順番を揃える必要がある
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j] = pipe_p_from0to3[0];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init1] = pipe_p_from1to3[0];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init2] = pipe_p_from2to3[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	pipe_p_from3to0[0] = p[j+N_init3];
	pipe_p_from3to1[0] = p[j+N_init3];
	pipe_p_from3to2[0] = p[j+N_init3];
	pipe_p_from3to4[0] = p[j+N_init3];
	pipe_p_from3to5[0] = p[j+N_init3];
	pipe_p_from3to6[0] = p[j+N_init3];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init4] = pipe_p_from4to3[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init5] = pipe_p_from5to3[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_end; ++j){
	p[j+N_init6] = pipe_p_from6to3[0];
}

		temp_rr1 = temp_rr2;
	}

// #pragma acc loop independent
	for(int j = 0; j < N_init1; ++j){
		X_result[j+N_init3] = x[j];
	}
}












#pragma acc serial async(4) pipein(pipe_temp_pap_from0to4, pipe_temp_rr2_from0to4, pipe_p_from0to4, pipe_temp_pap_from1to4, pipe_temp_rr2_from1to4, pipe_p_from1to4, pipe_temp_pap_from2to4, pipe_temp_rr2_from2to4, pipe_p_from2to4, pipe_temp_pap_from3to4, pipe_temp_rr2_from3to4, pipe_p_from3to4, pipe_temp_pap_from5to4, pipe_temp_rr2_from5to4, pipe_p_from5to4, pipe_temp_pap_from6to4, pipe_temp_rr2_from6to4, pipe_p_from6to4) pipeout(pipe_temp_pap_from4to0, pipe_temp_pap_from4to1, pipe_temp_pap_from4to2, pipe_temp_pap_from4to3, pipe_temp_pap_from4to5, pipe_temp_pap_from4to6, pipe_temp_rr2_from4to0, pipe_temp_rr2_from4to1, pipe_temp_rr2_from4to2, pipe_temp_rr2_from4to3, pipe_temp_rr2_from4to5, pipe_temp_rr2_from4to6, pipe_p_from4to0, pipe_p_from4to1, pipe_p_from4to2, pipe_p_from4to3, pipe_p_from4to5, pipe_p_from4to6)
{
	float alfa, beta;
	float x[BLOCK_SIZE/SPLIT_SIZE], r[BLOCK_SIZE/SPLIT_SIZE], p[BLOCK_SIZE], y[BLOCK_SIZE/SPLIT_SIZE];
	float VAL_local[V_SIZE];    // エミュレーションのときはコメントアウト
	int COL_IND_local[V_SIZE];  // エミュレーションのときはコメントアウト
	int ROW_PTR_local[(BLOCK_SIZE + 1)];
	float temp_sum=0.0f, temp_pap, temp_rr1, temp_rr2;

	temp_rr1 = 0.0f;
#pragma acc loop independent //reduction(+:temp_rr1)
	for(int i = 0; i < N; ++i){
		ROW_PTR_local[i] = ROW_PTR[i];
		p[i] = B[i];
		temp_rr1 += p[i] * p[i];
	}

	for(int i = 0; i < N_init1; ++i) {
		x[i] = 0.0f;
		r[i] = B[i+N_init4];
	}

	ROW_PTR_local[N] = ROW_PTR[N];

#pragma acc loop independent
	for(int i = 0; i < VAL_SIZE; ++i){
		COL_IND_local[i] = COL_IND[i];
		VAL_local[i] = VAL[i];
	}

	for(int i = 0; i < K; ++i){
		temp_pap = 0.0f;
#pragma acc loop reduction(+:temp_sum)
		for(int j = 0; j < N_init1; ++j){
			temp_sum = 0.0f;
			for(int l = ROW_PTR_local[j+N_init4]; l < ROW_PTR_local[j +N_init4+ 1]; ++l){
				temp_sum += VAL_local[l] * p[COL_IND_local[l]];
			}
			y[j] = temp_sum;
			temp_pap += p[j+N_init4] * temp_sum;
		}

// #pragma acc loop //reduction(+:temp_pap)
// 		for(int j = N_init4; j < N_init5; ++j){
// 			temp_pap += p[j] * y[j];
// 		}

// チャネル必要(temp_pap送信)
pipe_temp_pap_from4to0[0] = temp_pap;
pipe_temp_pap_from4to1[0] = temp_pap;
pipe_temp_pap_from4to2[0] = temp_pap;
pipe_temp_pap_from4to3[0] = temp_pap;
pipe_temp_pap_from4to5[0] = temp_pap;
pipe_temp_pap_from4to6[0] = temp_pap;
// チャネル必要(合計を受信)
temp_pap += pipe_temp_pap_from0to4[0];
temp_pap += pipe_temp_pap_from1to4[0];
temp_pap += pipe_temp_pap_from2to4[0];
temp_pap += pipe_temp_pap_from3to4[0];
temp_pap += pipe_temp_pap_from5to4[0];
temp_pap += pipe_temp_pap_from6to4[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		alfa = temp_rr1 / temp_pap;

		temp_rr2 = 0.0f;
#pragma acc loop reduction(+:temp_rr2)
		for(int j = 0; j < N_init1; ++j){
			x[j] += alfa * p[j+N_init4];
			r[j] -= alfa * y[j];
			temp_rr2 += r[j] * r[j];
		}
// チャネル必要(temp_rr2送信)
pipe_temp_rr2_from4to0[0] = temp_rr2;
pipe_temp_rr2_from4to1[0] = temp_rr2;
pipe_temp_rr2_from4to2[0] = temp_rr2;
pipe_temp_rr2_from4to3[0] = temp_rr2;
pipe_temp_rr2_from4to5[0] = temp_rr2;
pipe_temp_rr2_from4to6[0] = temp_rr2;
// チャネル必要(合計を受信)
temp_rr2 += pipe_temp_rr2_from0to4[0];
temp_rr2 += pipe_temp_rr2_from1to4[0];
temp_rr2 += pipe_temp_rr2_from2to4[0];
temp_rr2 += pipe_temp_rr2_from3to4[0];
temp_rr2 += pipe_temp_rr2_from5to4[0];
temp_rr2 += pipe_temp_rr2_from6to4[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		beta = temp_rr2 / temp_rr1;

// #pragma acc loop independent
		for(int j = 0; j < N_init1; ++j){
			p[j+N_init4] = r[j] + beta * p[j+N_init4];
		}

// カーネルごとの送信の順番を揃える必要がある
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j] = pipe_p_from0to4[0];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init1] = pipe_p_from1to4[0];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init2] = pipe_p_from2to4[0];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init3] = pipe_p_from3to4[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	pipe_p_from4to0[0] = p[j+N_init4];
	pipe_p_from4to1[0] = p[j+N_init4];
	pipe_p_from4to2[0] = p[j+N_init4];
	pipe_p_from4to3[0] = p[j+N_init4];
	pipe_p_from4to5[0] = p[j+N_init4];
	pipe_p_from4to6[0] = p[j+N_init4];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init5] = pipe_p_from5to4[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_end; ++j){
	p[j+N_init6] = pipe_p_from6to4[0];
}

		temp_rr1 = temp_rr2;
	}

// #pragma acc loop independent
	for(int j = 0; j < N_init1; ++j){
		X_result[j+N_init4] = x[j];
	}
}












#pragma acc serial async(5) pipein(pipe_temp_pap_from0to5, pipe_temp_rr2_from0to5, pipe_p_from0to5, pipe_temp_pap_from1to5, pipe_temp_rr2_from1to5, pipe_p_from1to5, pipe_temp_pap_from2to5, pipe_temp_rr2_from2to5, pipe_p_from2to5, pipe_temp_pap_from3to5, pipe_temp_rr2_from3to5, pipe_p_from3to5, pipe_temp_pap_from4to5, pipe_temp_rr2_from4to5, pipe_p_from4to5, pipe_temp_pap_from6to5, pipe_temp_rr2_from6to5, pipe_p_from6to5) pipeout(pipe_temp_pap_from5to0, pipe_temp_pap_from5to1, pipe_temp_pap_from5to2, pipe_temp_pap_from5to3, pipe_temp_pap_from5to4, pipe_temp_pap_from5to6, pipe_temp_rr2_from5to0, pipe_temp_rr2_from5to1, pipe_temp_rr2_from5to2, pipe_temp_rr2_from5to3, pipe_temp_rr2_from5to4, pipe_temp_rr2_from5to6, pipe_p_from5to0, pipe_p_from5to1, pipe_p_from5to2, pipe_p_from5to3, pipe_p_from5to4, pipe_p_from5to6)
{
	float alfa, beta;
	float x[BLOCK_SIZE/SPLIT_SIZE], r[BLOCK_SIZE/SPLIT_SIZE], p[BLOCK_SIZE], y[BLOCK_SIZE/SPLIT_SIZE];
	float VAL_local[V_SIZE];    // エミュレーションのときはコメントアウト
	int COL_IND_local[V_SIZE];  // エミュレーションのときはコメントアウト
	int ROW_PTR_local[(BLOCK_SIZE + 1)];
	float temp_sum=0.0f, temp_pap, temp_rr1, temp_rr2;

	temp_rr1 = 0.0f;
#pragma acc loop independent //reduction(+:temp_rr1)
	for(int i = 0; i < N; ++i){
		ROW_PTR_local[i] = ROW_PTR[i];
		p[i] = B[i];
		temp_rr1 += p[i] * p[i];
	}

	for(int i = 0; i < N_init1; ++i) {
		x[i] = 0.0f;
		r[i] = B[i+N_init5];
	}

	ROW_PTR_local[N] = ROW_PTR[N];

#pragma acc loop independent
	for(int i = 0; i < VAL_SIZE; ++i){
		COL_IND_local[i] = COL_IND[i];
		VAL_local[i] = VAL[i];
	}

	for(int i = 0; i < K; ++i){
		temp_pap = 0.0f;
#pragma acc loop reduction(+:temp_sum)
		for(int j = 0; j < N_init1; ++j){
			temp_sum = 0.0f;
			for(int l = ROW_PTR_local[j+N_init5]; l < ROW_PTR_local[j +N_init5+ 1]; ++l){
				temp_sum += VAL_local[l] * p[COL_IND_local[l]];
			}
			y[j] = temp_sum;
			temp_pap += p[j+N_init5] * temp_sum;
		}

// #pragma acc loop //reduction(+:temp_pap)
// 		for(int j = N_init5; j < N_init5; ++j){
// 			temp_pap += p[j] * y[j];
// 		}

// チャネル必要(temp_pap送信)
pipe_temp_pap_from5to0[0] = temp_pap;
pipe_temp_pap_from5to1[0] = temp_pap;
pipe_temp_pap_from5to2[0] = temp_pap;
pipe_temp_pap_from5to3[0] = temp_pap;
pipe_temp_pap_from5to4[0] = temp_pap;
pipe_temp_pap_from5to6[0] = temp_pap;
// チャネル必要(合計を受信)
temp_pap += pipe_temp_pap_from0to5[0];
temp_pap += pipe_temp_pap_from1to5[0];
temp_pap += pipe_temp_pap_from2to5[0];
temp_pap += pipe_temp_pap_from3to5[0];
temp_pap += pipe_temp_pap_from4to5[0];
temp_pap += pipe_temp_pap_from6to5[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		alfa = temp_rr1 / temp_pap;

		temp_rr2 = 0.0f;
#pragma acc loop reduction(+:temp_rr2)
		for(int j = 0; j < N_init1; ++j){
			x[j] += alfa * p[j+N_init5];
			r[j] -= alfa * y[j];
			temp_rr2 += r[j] * r[j];
		}
// チャネル必要(temp_rr2送信)
pipe_temp_rr2_from5to0[0] = temp_rr2;
pipe_temp_rr2_from5to1[0] = temp_rr2;
pipe_temp_rr2_from5to2[0] = temp_rr2;
pipe_temp_rr2_from5to3[0] = temp_rr2;
pipe_temp_rr2_from5to4[0] = temp_rr2;
pipe_temp_rr2_from5to6[0] = temp_rr2;
// チャネル必要(合計を受信)
temp_rr2 += pipe_temp_rr2_from0to5[0];
temp_rr2 += pipe_temp_rr2_from1to5[0];
temp_rr2 += pipe_temp_rr2_from2to5[0];
temp_rr2 += pipe_temp_rr2_from3to5[0];
temp_rr2 += pipe_temp_rr2_from4to5[0];
temp_rr2 += pipe_temp_rr2_from6to5[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		beta = temp_rr2 / temp_rr1;

// #pragma acc loop independent
		for(int j = 0; j < N_init1; ++j){
			p[j+N_init5] = r[j] + beta * p[j+N_init5];
		}

// カーネルごとの送信の順番を揃える必要がある
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j] = pipe_p_from0to5[0];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init1] = pipe_p_from1to5[0];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init2] = pipe_p_from2to5[0];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init3] = pipe_p_from3to5[0];
}
for(int j = 0; j < N_init1; ++j){
	p[j+N_init4] = pipe_p_from4to5[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	pipe_p_from5to0[0] = p[j+N_init5];
	pipe_p_from5to1[0] = p[j+N_init5];
	pipe_p_from5to2[0] = p[j+N_init5];
	pipe_p_from5to3[0] = p[j+N_init5];
	pipe_p_from5to4[0] = p[j+N_init5];
	pipe_p_from5to6[0] = p[j+N_init5];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_end; ++j){
	p[j+N_init6] = pipe_p_from6to5[0];
}

		temp_rr1 = temp_rr2;
	}

// #pragma acc loop independent
	for(int j = 0; j < N_init1; ++j){
		X_result[j+N_init5] = x[j];
	}
}












#pragma acc serial async(6) pipein(pipe_temp_pap_from0to6, pipe_temp_rr2_from0to6, pipe_p_from0to6, pipe_temp_pap_from1to6, pipe_temp_rr2_from1to6, pipe_p_from1to6, pipe_temp_pap_from2to6, pipe_temp_rr2_from2to6, pipe_p_from2to6, pipe_temp_pap_from3to6, pipe_temp_rr2_from3to6, pipe_p_from3to6, pipe_temp_pap_from4to6, pipe_temp_rr2_from4to6, pipe_p_from4to6, pipe_temp_pap_from5to6, pipe_temp_rr2_from5to6, pipe_p_from5to6) pipeout(pipe_temp_pap_from6to0, pipe_temp_pap_from6to1, pipe_temp_pap_from6to2, pipe_temp_pap_from6to3, pipe_temp_pap_from6to4, pipe_temp_pap_from6to5, pipe_temp_rr2_from6to0, pipe_temp_rr2_from6to1, pipe_temp_rr2_from6to2, pipe_temp_rr2_from6to3, pipe_temp_rr2_from6to4, pipe_temp_rr2_from6to5, pipe_p_from6to0, pipe_p_from6to1, pipe_p_from6to2, pipe_p_from6to3, pipe_p_from6to4, pipe_p_from6to5)
{
	float alfa, beta;
	float x[BLOCK_SIZE/SPLIT_SIZE+BLOCK_SIZE%SPLIT_SIZE], r[BLOCK_SIZE/SPLIT_SIZE+BLOCK_SIZE%SPLIT_SIZE], p[BLOCK_SIZE], y[BLOCK_SIZE/SPLIT_SIZE+BLOCK_SIZE%SPLIT_SIZE];
	float VAL_local[V_SIZE];    // エミュレーションのときはコメントアウト
	int COL_IND_local[V_SIZE];  // エミュレーションのときはコメントアウト
	int ROW_PTR_local[(BLOCK_SIZE + 1)];
	float temp_sum=0.0f, temp_pap, temp_rr1, temp_rr2;

	temp_rr1 = 0.0f;
#pragma acc loop independent //reduction(+:temp_rr1)
	for(int i = 0; i < N; ++i){
		ROW_PTR_local[i] = ROW_PTR[i];
		p[i] = B[i];
		temp_rr1 += p[i] * p[i];
	}

	for(int i = 0; i < N_end; ++i) {
		x[i] = 0.0f;
		r[i] = B[i+N_init6];
	}

	ROW_PTR_local[N] = ROW_PTR[N];

#pragma acc loop independent
	for(int i = 0; i < VAL_SIZE; ++i){
		COL_IND_local[i] = COL_IND[i];
		VAL_local[i] = VAL[i];
	}

	for(int i = 0; i < K; ++i){
		temp_pap = 0.0f;
#pragma acc loop reduction(+:temp_sum)
		for(int j = 0; j < N_end; ++j){
			temp_sum = 0.0f;
			for(int l = ROW_PTR_local[j+N_init6]; l < ROW_PTR_local[j +N_init6+ 1]; ++l){
				temp_sum += VAL_local[l] * p[COL_IND_local[l]];
			}
			y[j] = temp_sum;
			temp_pap += p[j+N_init6] * temp_sum;
		}

// #pragma acc loop //reduction(+:temp_pap)
// 		for(int j = N_init6; j < N; ++j){
// 			temp_pap += p[j] * y[j];
// 		}

// チャネル必要(temp_pap送信)
pipe_temp_pap_from6to0[0] = temp_pap;
pipe_temp_pap_from6to1[0] = temp_pap;
pipe_temp_pap_from6to2[0] = temp_pap;
pipe_temp_pap_from6to3[0] = temp_pap;
pipe_temp_pap_from6to4[0] = temp_pap;
pipe_temp_pap_from6to5[0] = temp_pap;
// チャネル必要(合計を受信)
temp_pap += pipe_temp_pap_from0to6[0];
temp_pap += pipe_temp_pap_from1to6[0];
temp_pap += pipe_temp_pap_from2to6[0];
temp_pap += pipe_temp_pap_from3to6[0];
temp_pap += pipe_temp_pap_from4to6[0];
temp_pap += pipe_temp_pap_from5to6[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		alfa = temp_rr1 / temp_pap;

		temp_rr2 = 0.0f;
#pragma acc loop reduction(+:temp_rr2)
		for(int j = 0; j < N_end; ++j){
			x[j] += alfa * p[j+N_init6];
			r[j] -= alfa * y[j];
			temp_rr2 += r[j] * r[j];
		}
// チャネル必要(temp_rr2送信)
pipe_temp_rr2_from6to0[0] = temp_rr2;
pipe_temp_rr2_from6to1[0] = temp_rr2;
pipe_temp_rr2_from6to2[0] = temp_rr2;
pipe_temp_rr2_from6to3[0] = temp_rr2;
pipe_temp_rr2_from6to4[0] = temp_rr2;
pipe_temp_rr2_from6to5[0] = temp_rr2;
// チャネル必要(合計を受信)
temp_rr2 += pipe_temp_rr2_from0to6[0];
temp_rr2 += pipe_temp_rr2_from1to6[0];
temp_rr2 += pipe_temp_rr2_from2to6[0];
temp_rr2 += pipe_temp_rr2_from3to6[0];
temp_rr2 += pipe_temp_rr2_from4to6[0];
temp_rr2 += pipe_temp_rr2_from5to6[0];

// 片方だけで計算しようかと思ったけど、パイプのほうが重そうなので全カーネルで計算
		beta = temp_rr2 / temp_rr1;

// #pragma acc loop independent
		for(int j = 0; j < N_end; ++j){
			p[j+N_init6] = r[j] + beta * p[j+N_init6];
		}

// カーネルごとの送信の順番を揃える必要がある
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j] = pipe_p_from0to6[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init1] = pipe_p_from1to6[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init2] = pipe_p_from2to6[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init3] = pipe_p_from3to6[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init4] = pipe_p_from4to6[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_init1; ++j){
	p[j+N_init5] = pipe_p_from5to6[0];
}
// チャネル必要（各カーネルで担当外の要素の更新が必要）
for(int j = 0; j < N_end; ++j){
	pipe_p_from6to0[0] = p[j+N_init6];
	pipe_p_from6to1[0] = p[j+N_init6];
	pipe_p_from6to2[0] = p[j+N_init6];
	pipe_p_from6to3[0] = p[j+N_init6];
	pipe_p_from6to4[0] = p[j+N_init6];
	pipe_p_from6to5[0] = p[j+N_init6];
}

		temp_rr1 = temp_rr2;
	}

// #pragma acc loop independent
	for(int j = 0; j < N_end; ++j){
		X_result[j+N_init6] = x[j];
	}
}












#pragma acc wait

clock_gettime(CLOCK_REALTIME, &fpga_calc_end);
printf("FPGA caluctation time: %lf sec\n", time_diff(&fpga_calc_start, &fpga_calc_end));

clock_gettime(CLOCK_REALTIME, &fpga_copyout_start);
#pragma acc update host(X_result[0:N])
clock_gettime(CLOCK_REALTIME, &fpga_copyout_end);
printf("FPGA to Host time: %lf sec\n", time_diff(&fpga_copyout_start, &fpga_copyout_end));
}

#pragma acc exit data delete(VAL[0:VAL_SIZE], COL_IND[0:VAL_SIZE], ROW_PTR[0:N+1], B[0:N], N, K, VAL_SIZE) delete(X_result[0:N])
	// acc_shutdown(acc_device_altera);
}



void matmul(float *a, float *b, float *c, int N, float *e, float *f)
{
	struct timespec gpu_copyin_start;
	struct timespec gpu_copyin_end;
	struct timespec gpu_calc_start;
	struct timespec gpu_calc_end;
	struct timespec gpu_copyout_start;
	struct timespec gpu_copyout_end;

#pragma accomn target_dev(GPU)
{
	int i, j, k;
	float sum;
#pragma acc enter data create(N, a[:N*N], b[:N*N], c[:N*N], e[:N], f[:N])
	// std::chrono::system_clock::time_point start_gpu;

#pragma acc data present(N, c, a, b, f, e)
{
// start_gpu = std::chrono::system_clock::now();

clock_gettime(CLOCK_REALTIME, &gpu_copyin_start);
#pragma acc update device(N, a[:N*N], b[:N*N], e[:N])
clock_gettime(CLOCK_REALTIME, &gpu_copyin_end);
printf("Host to GPU time: %lf sec\n", time_diff(&gpu_copyin_start, &gpu_copyin_end));

clock_gettime(CLOCK_REALTIME, &gpu_calc_start);
#pragma acc kernels
	{
#pragma acc loop gang vector(32) independent
		for (i = 0; i < N; ++i)
		{
#pragma acc loop gang vector(32) independent
			for (j = 0; j < N; ++j)
			{
				sum = 0.0;
#pragma acc loop independent reduction(+:sum)
				for (k = 0; k < N; ++k)
					sum += a[i * N + k] * b[k * N + j];
				c[i * N + j] = sum;
			}
		}

//matrix_vector_malti

#pragma acc loop gang vector independent
		for (i = 0; i < N; ++i)
		{
			sum = 0.0;
#pragma acc loop independent reduction(+:sum)
			for (j = 0; j < N; ++j)
				sum += c[i * N + j] * e[j];
			f[i] = sum;
		}
	}
clock_gettime(CLOCK_REALTIME, &gpu_calc_end);
printf("GPU caluctation time: %lf sec\n", time_diff(&gpu_calc_start, &gpu_calc_end));

clock_gettime(CLOCK_REALTIME, &gpu_copyout_start);
#pragma acc update host(f[:N])
clock_gettime(CLOCK_REALTIME, &gpu_copyout_end);
printf("GPU to Host time: %lf sec\n", time_diff(&gpu_copyout_start, &gpu_copyout_end));

}
}

// return start_gpu;
}

// void matrix_vector_malti(float *a, float *b, float *c, int N)
// {
// 	int i, j;

// // <<<dim3(numblock), dim3(numthread)>>>
// #pragma acc update device(a[:N*N], b[:N])

// #pragma acc kernels present(c, a, b)
// 	{
// #pragma acc loop independent gang(N/32) vector(32)
// 		for (i = 0; i < N; ++i)
// 		{
// 			float sum = 0.0;
// #pragma acc loop reduction(+:sum)
// 			for (j = 0; j < N; ++j)
// 				sum += a[i * N + j] * b[j];
// 			c[i] = sum;
// 		}
// 	}

// #pragma acc update self(c[:N])
// }

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
	double diff_sum = 0.0;
	double diff_temp;
	double diff_pow;
	double diff[N];
	double err_max = 0.0;
	double rel_err = 0.0;
	unsigned long i;
	int error = 0;

#pragma omp parallel for reduction(+:cpu_sum, diff_sum) private(diff_temp, diff_pow)
	for (i = 0; i < N; ++i)
	{
		// printf("(CPU) %f\n", c_CPU[i]);
		// printf("(GPU) %f\n", h_c[i]);
		// if(h_c[i] == c_CPU[i])
		// 	diff[i] = 0.0;
		// else {
			cpu_sum  += c_CPU[i] * c_CPU[i];
			diff_temp = h_c[i] - c_CPU[i];
			diff_pow = diff_temp * diff_temp;
			diff[i]   = sqrt(diff_pow);
			diff_sum += diff_pow;
		// }
	}

	for (i = 0; i < N; ++i) {
		if(diff[i] != 0.0)
		    error++;
		if(diff[i]/sqrt(c_CPU[i]*c_CPU[i]) > err_max) {
			err_max = diff[i]/sqrt(c_CPU[i]*c_CPU[i]);
		}
	}
	printf("GPU err_max: %e\n"
		"num error: %d\n"
	, err_max, error);
	
	// gpu_sum = sqrt(gpu_sum);
	// diff_sum = sqrt(diff_sum);

	rel_err = sqrt(diff_sum) / sqrt(cpu_sum);

	if (rel_err < 1e-6)
	{
		printf("GPU Verification: PASS\n"
		"Verification Successful err = %e\n", rel_err);
	}
	else
	{
		printf("Error! GPU Verification failed...\n"
		"Verification Fail err = %e\n", rel_err);
	}
	// printf("ResultGPU = %lf\n", gpu_sum);
	// printf("ResultCPU = %lf\n", cpu_sum);
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
  int error = 0;

	float x[N], r[N], p[N], y[N], alfa, beta;
	float VAL_local[VAL_SIZE];
	int COL_IND_local[VAL_SIZE], ROW_PTR_local[N + 1];
	float temp_sum, temp_pap, temp_rr1, temp_rr2, cpu_sum = 0; //sum = 0

//   std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
struct timespec host_CG_start;
struct timespec host_CG_end;
clock_gettime(CLOCK_REALTIME, &host_CG_start);

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
				temp_sum += VAL_local[l] * p[COL_IND_local[l]];
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
	
	for(int j = 0; j < 10974; ++j){
		B[j] = x[j];  // 無駄な処理だけど、FPGAと実行時間をフェアにするために実行
	}
clock_gettime(CLOCK_REALTIME, &host_CG_end);
printf("Host CG calc time: %lf sec\n", time_diff(&host_CG_start, &host_CG_end));

	for(int j = 0; j < 10974; ++j){
		r[j] = B[j];  // コンパイラがloadが無いとメモリコピーを無視する可能性があるので対策
	}

//   std::chrono::system_clock::time_point end = std::chrono::system_clock::now();

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

	float diff_sum = 0.0f;
	float diff_temp;
	float diff_pow;
	float diff[N];
	float err_max = 0.0f;
	float rel_err = 0.0f;
	float fpga_sum = 0.0f;

	for(int i = 0; i < N; ++i) {
    // std::cout << "FPGA" << FPGA_calc_result[j] << ", CPU"<< x[j] << std::endl;
		// if(FPGA_calc_result[i] == x[i])
			// diff[i] = 0.0;
		// else {
			cpu_sum  += x[i] * x[i];
			diff_temp = FPGA_calc_result[i] - x[i];
			diff_pow = diff_temp * diff_temp;
			diff[i]   = sqrtf(diff_pow);
			diff_sum += diff_pow;
			fpga_sum += FPGA_calc_result[i] * FPGA_calc_result[i];
		// }
    }
    // sum += FPGA_calc_result[j];
    // cpu_sum += x[j];

printf("fpga_sum: %lf\n", fpga_sum);
printf("sqrtf(fpga_sum): %lf\n", sqrtf(fpga_sum));
printf("cpusum: %lf\n", cpu_sum);
printf("sqrtf(cpusum): %lf\n", sqrtf(cpu_sum));
printf("diff_sum: %lf\n", diff_sum);
printf("sqrtf(diff_sum): %lf\n", sqrtf(diff_sum));

	for (int i = 0; i < N; ++i) {
		if(diff[i] != 0.0)
		    error++;
		if(diff[i]/sqrtf(x[i]*x[i]) > err_max) {
			err_max = diff[i]/sqrtf(x[i]*x[i]);
		}
	}
	printf("-----------------------------------------------------------------\n"
		"FPGA err_max: %e\n"
		"num error: %d\n"
	, err_max, error);
	
	rel_err = sqrtf(diff_sum) / sqrtf(cpu_sum);

	if (rel_err < 1e-6)
	{
		printf("FPGA Verification: PASS\n"
		"Verification Successful err = %e\n"	
		, rel_err);
		// "ResultFPGA = %f\n"
		// , sum);	}
	}
	else
	{
		printf("Error! FPGA Verification failed...\n"
			"Verification Fail err = %e\n", rel_err);
		// "ResultFPGA = %f\n"
		// "ResultCPU  = %f\n"
		// "num error: %d\n"
		// , sum, cpu_sum, error);
	}

//   std::cout << "CG CPU elapsed time: " << std::fixed << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() << " usec" << std::endl;
}


int main(int argc, char *argv[])
{
	struct sparse_matrix_t* A_ = load_sparse_matrix(MATRIX_MARKET, "bcsstk17.mtx");
	// assert(A_ != NULL);
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
//   assert (A);
//   assert (A->nnz == (A->rowptr[A->m] - A->rowptr[0]));

	// check command line arguments
	///////////////////////////////////////////
	// if (argc == 1)
	// {
	// 	// std::cout << "usage: ./host <numdata_h> <valsize> <numtry>"   << std::endl;
	// 	exit(0);
	// }
	// if (argc != 4)
	// {
	// 	// std::cerr << "Error! The number of arguments is wrong."       << std::endl;
	// 	exit(1);
	// }

  const int  numdata_h = A->n; // std::stoull(std::string(argv[2]));
	int N = numdata_h;
  const int  valsize   = A->nnz;
	int VAL_SIZE = valsize;
	const int numtry = 1000;
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

	// 初期化？
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

	matmul(h_a, h_b, h_c, numdata_h, h_vec_mul, h_vec_b);
	// matrix_vector_malti(h_c, h_vec_mul, h_vec_b, numdata_h);

	// std::chrono::system_clock::time_point end_gpu = std::chrono::system_clock::now();

// #pragma acc exit data delete(h_a, h_b, h_c)
// #pragma acc exit data delete(h_vec_mul, h_vec_b)

	/***** FPGA *****/
	for (int j = 0; j < N; ++j)
	{
		FPGA_calc_result[j] = 0;
		// ROW_PTR[j] = A->rowptr[j];
		B[j] = h_vec_b[j] - VAL[j] * 1; //000000.0; // b - Ax
	}
	// ROW_PTR[N] = N;

		// initFPGA(FPGA_calc_result, VAL, COL_IND, ROW_PTR, B, N, K, VAL_SIZE);

		// std::chrono::system_clock::time_point start_fpga = std::chrono::system_clock::now();

		// sendDataToFPGA(FPGA_calc_result, VAL, COL_IND, ROW_PTR, B, N, K, VAL_SIZE);

	funcFPGA(FPGA_calc_result, VAL, COL_IND, ROW_PTR, B, N, K, VAL_SIZE);

		// recvDataFromFPGA(FPGA_calc_result, VAL, COL_IND, ROW_PTR, B, N, K, VAL_SIZE);
		
		// std::chrono::system_clock::time_point end_fpga = std::chrono::system_clock::now();

		// shutdownFPGA(FPGA_calc_result, VAL, COL_IND, ROW_PTR, B, N, K, VAL_SIZE);

		// std::cout << "GPU  elapsed time: " << std::fixed << std::chrono::duration_cast<std::chrono::microseconds>(end_gpu-start_gpu).count() << " usec" << std::endl;
		// std::cout << std::string(30, '-') << std::endl;

		// std::cout << "FPGA elapsed time: " << std::fixed << std::chrono::duration_cast<std::chrono::microseconds>(end_fpga-start_fpga).count() << " usec" << std::endl;
		// std::cout << std::string(30, '-') << std::endl;

	// verification
	///////////////////////////////////////////
	// MatrixMultiplication_openmp(h_a, h_b, c_CPU, numdata_h);    // 本番はコメントアウトして良い
	// h_matrix_vector_malti(c_CPU, h_vec_mul, vec_b_CPU, numdata_h);    // 本番はコメントアウトして良い

	// verify_gpu(h_c, c_CPU, numdata_h*numdata_h); // 行列積チェック（予備のコード）
	// verify_gpu(h_vec_b, vec_b_CPU, numdata_h); // h_vec_b チェック

	// FPGAの結果検証ではCPUでの計算にもGPUの結果を用いている（GPUの結果検証とは切り離している）
	verify_fpga(FPGA_calc_result, VAL, COL_IND, ROW_PTR, B, N, K, VAL_SIZE);

	// cleanup
	///////////////////////////////////////////
  // destroy_sparse_matrix(A_);
	// destroy_csr_matrix(A);
	

	return 0;
}
