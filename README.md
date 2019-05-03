# gpu_fpga_hybrid_matmul_cg
PGI Compiler と OpenARC を組み合わせて、OpenACC のプログラム同士をリンクコンパイルすることで GPU と FPGA を両方使うプログラム。
GPUで行列積とベクトル化、FPGAでCG法の計算を行う。