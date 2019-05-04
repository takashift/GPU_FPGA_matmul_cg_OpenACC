#!/bin/sh

module load quartus/17.1.2.304 aocl/a10pl4_4
#pgi
module load cuda/9.2.148
module load pgi/18.10

export OPENARC_ARCH=3
export ACC_DEVICE_TYPE=acc_device_not_host
export ACC_DEVICE_NUM=0
export OPENARCRT_UNIFIEDMEM=0
export OPENARC_FPGA=ARRIA_X
#export ACC_DEVICE_TYPE=RADEON  # 必要ないかも知れない？
export openarc=/home/tsunashima/openarc-fpga-master
#./O2GBuild.script
cp $openarc/make.template.obj $openarc/make.template
make COMPILE_HOST
cd cetus_output/
pgc++ -acc -c -ta=tesla:cc60 -O3 -std=c++11 -mp -Minfo=accel ../matmul_gpu.cpp
cd ..
cp $openarc/make.template.pgi $openarc/make.template
make COMPILE_HOST
cp $openarc/make.template.old $openarc/make.template

#aoc openarc_kernel.aoco
#cp openarc_kernel.aocx ../bin/
#cd ../bin
#./each_block1_ACC
