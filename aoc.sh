#!/bin/sh

module load quartus/17.1.2.304 aocl/a10pl4_4

cd cetus_output
aoc openarc_kernel.aoco
cp openarc_kernel.aocx ../bin/
