#!/bin/sh

module load quartus/17.1.2.304 aocl/a10pl4_4
module load cuda/9.2.148
module load pgi/18.10
export OPENARC_ARCH=3
export ACC_DEVICE_TYPE=acc_device_not_host
export ACC_DEVICE_NUM=0
export OPENARCRT_UNIFIEDMEM=0
export OPENARC_FPGA=ARRIA_X
#export ACC_DEVICE_TYPE=RADEON  # 必要ないかも知れない？
#export openarc=/home/tsunashima/openarc-fpga-master
#./O2GBuild.script
#make
#cd cetus_output
#aoc openarc_kernel.aoco
#mv openarc_kernel.aocx ../bin/
#cd ../bin
#aocl diagnose
#./each_block2_gpu_ACC

bin/gpu_fpga_matmul_cg_ACC 1000 1000 1000

