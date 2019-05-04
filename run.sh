#!/bin/sh

module load quartus/17.1.2.304 aocl/a10pl4_4
module load cuda/9.2.148
module load pgi/18.10
export OPENARC_ARCH=3

#export openarc=/home/tsunashima/openarc-fpga-master
#./O2GBuild.script
#make
#cd cetus_output
#aoc openarc_kernel.aoco
#mv openarc_kernel.aocx ../bin/
#cd ../bin
#aocl diagnose
#./each_block2_gpu_ACC

cd bin
./gpu_fpga_matmul_cg_ACC 1000 1000 1000

