# GPU_FPGA_matmul_cg_OpenACC
The matmul-CG program for MHOAT: Multi-Hybrid OpenACC Translator.  
GPU and FPGA cooperative computing in a single computing process with MHOAT and backend compier: PGI(nvc) and OpenARC.

- Calculation on GPU:  matrix-matrix multiplication and matrix-vector multiprication
- Calculation on FPGA: Conjugate Gradient solver with GPU’s result vector

CUDA+OpenCL ver.: https://github.com/takashift/cuda_opencl_matmul_cg