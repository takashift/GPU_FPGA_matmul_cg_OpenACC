# GPU_FPGA_matmul_cg_OpenACC
PGI Compiler と OpenARC を組み合わせて、OpenACC のプログラム同士をリンクコンパイルすることで一つのプロセスから GPU と FPGA を両方使うプログラム。

GPUで行列積とベクトル化、FPGAでCG法の計算を行う。

HPC169で発表。

# 実行手順
OpenARC V0.14時点

ジョブスクリプト内のパスは適宜変更。
- PGIについて補足
    - FPGA側はCであるがOpenARCで変換後にはC++となるため、GPU側の方はマングリング問題を回避するためにC++で記述(pgc++でコンパイル)する。
    - -ta=tesla:cc60はTesla P100向け

Makefileの以下については各自変更する。
- BENCHMARK
- CSRCS
- AOCL_BOARD
- CLIBS2（ライブラリパスのオプション設定）

また、OpenARC-master/にはmake.templateのCOMPILE_HOSTを以下のように変更したファイルを用意する。

### make.template.obj
```
COMPILE_HOST: $($(CETUS_OUTPUT)/$(OCXXSRCS)) $(CXXSRCS) makedirectories
        $(COPY_CXXSRCS)
        $(PRECMD)
        cd $(CETUS_OUTPUT); $(CXX) $(DEFSET_ACC) $(CFLAGS2) -I ../ -o $(BENCHMARK_OPENACC).o -c -g $(OCXXSRCS) $(CXXSRCS) $(CLIBS2_BASE) $(CLIBS2); cp $(KERNEL_FILE) ../$(TARGET); if [ -f "$(OPENARCLIB)/binBuilder_$(OPENARCLIB_SUFFIX)" ]; then cp $(OPENARCLIB)/binBuilder_$(OPENARCLIB_SUFFIX) ../$(TARGET); fi; cp $(OPENARCLIB)/Timer ../$(TARGET); cd ../
        cd $(CETUS_OUTPUT); grep resilience $(KERNEL_FILE) > /dev/null && cp $(OPENARCLIB)/$(RESILIENCE_FILE) $(OPENARCLIB)/resilience.h ../$(TARGET); cd ../
```

### make.template.pgi
```
COMPILE_HOST: $($(CETUS_OUTPUT)/$(OCXXSRCS)) $(CXXSRCS) makedirectories
        $(COPY_CXXSRCS)
        $(PRECMD)
        cd $(CETUS_OUTPUT); pgc++ -acc -ta=tesla:cc60 -g -std=c++11 -mp -Minfo=accel $(DEFSET_ACC) $(CFLAGS2) -I ../ -o ../$(BENCHMARK_OPENACC) matmul_gpu.o gpu_fpga_matmul_cg_ACC.o $(CLIBS2_BASE) $(CLIBS2); cp $(KERNEL_FILE) ../$(TARGET); if [ -f "$(OPENARCLIB)/binBuilder_$(OPENARCLIB_SUFFIX)" ]; then cp $(OPENARCLIB)/binBuilder_$(OPENARCLIB_SUFFIX) ../$(TARGET); fi; cp $(OPENARCLIB)/Timer ../$(TARGET); cd ../
        cd $(CETUS_OUTPUT); grep resilience $(KERNEL_FILE) > /dev/null && cp $(OPENARCLIB)/$(RESILIENCE_FILE) $(OPENARCLIB)/resilience.h ../$(TARGET); cd ../
```

## コンパイル方法
```
./compile_kernel.sh
sbatch -p syn3 ./aoc.sh
./compile_host.sh
```

## 実行方法
```
srun -l run.sh
```