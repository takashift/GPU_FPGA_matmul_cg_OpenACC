# GPU_FPGA_matmul_cg_OpenACC
PGI Compiler と OpenARC を組み合わせて、OpenACC のプログラム同士をリンクコンパイルすることで一つのプロセスから GPU と FPGA を両方使うプログラム。

GPUで行列積とベクトル化、FPGAでCG法の計算を行う。

OmniCompiler に実装予定の OpenACC to OpenACC の ソース to ソース コンパイラの出力を想定して記述。

HPC169で発表。

# 実行手順
OpenARC V0.14時点

ジョブスクリプト内のパスは適宜変更。
- PGIについて補足
    - FPGA側はCであるがOpenARCで変換後にはC++となるため、GPU側の方はマングリング問題を回避するためにpgc++でコンパイルする。
    - -ta=tesla:cc60はTesla P100向け

Makefileの以下については各自変更する。
- BENCHMARK
- CSRCS
- AOCL_BOARD
- CLIBS2（ライブラリパスのオプション設定）
