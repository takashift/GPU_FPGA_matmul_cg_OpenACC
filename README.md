# GPU_FPGA_matmul_cg_OpenACC
PGI Compiler と OpenARC を組み合わせて、OpenACC のプログラム同士をリンクコンパイルすることで一つのプロセスから GPU と FPGA を両方使うプログラム。  
GPUで行列積とベクトル化、FPGAでCG法の計算を行う。
HPC169で発表。