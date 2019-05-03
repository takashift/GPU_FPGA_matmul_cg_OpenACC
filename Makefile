include ../../../../make.header

########################
# Set the program name #
########################
BENCHMARK = gpu_fpga_matmul_cg

########################################
# Set the input C source files (CSRCS) #
########################################
CSRCS = cg_fpga.c

#########################################
# Set macros used for the input program #
#########################################
#_N_ ?= 512
#DEFSET_CPU = -D_N_=$(_N_)
#DEFSET_ACC = -D_N_=$(_N_)

#########################################################
# Makefile options that the user can overwrite          #
# OMP: set to 1 to use OpenMP (default: 0)              # 
# MODE: set to profile to use a built-in profiling tool #
#       (default: normal)                               #
#       If this is set to profile, the runtime system   #
#       will print profiling results according to the   #
#       verbosity level set by OPENARCRT_VERBOSITY      #
#       environment variable.                           # 
#########################################################
OMP ?= 0
MODE ?= normal
AOCL_BOARD ?= a10pl4_dd4gb_gx115_m512
AOCL_FLAGS ?= -v -c -report

#########################################################
# Use the following macros to give program-specific     #
# compiler flags and libraries                          #
# - CFLAGS1 and CLIBS1 to compile the input C program   #
# - CFLAGS2 and CLIBS2 to compile the OpenARC-generated #
#   output C++ program                                  # 
#########################################################
#CFLAGS1 =  
#CFLAGS2 =  
#CLIBS1 = 
#CLIBS2 = 

################################################
# TARGET is where the output binary is stored. #
################################################
#TARGET ?= ./bin

include ../../../../make.template
