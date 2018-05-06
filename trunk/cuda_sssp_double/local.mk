###############################################################################
#
# File      : common.mk
# Author    : Alex Stivala (astivala)
# Created   : July 2008
#
# $Id: Makefile 1690 2008-07-16 04:35:47Z astivala $
#
# Definitions of compilers and compile options for all Makefiles.
# Used for compiling with MPI and POSIX threads, and CUDA.
# Use GNU make.
#
# set MODE=DEVICEEMU to build with device emulation and debugging on. 
# default is to build with optimizations on and no debug or profile, for GPU
# set MODE=GPUDEBUG to build with debugging on GPU with cuda-gdb enabled.
# Otherwise,
# set MODE=DEBUG to build with debugging and verbose printing on host,
# default is to build with optimizations on and no debug or profile, for GPU
# execution
#
###############################################################################

CC           = gcc
LD           = gcc

# must ensure that MPICC uses the same compiler CC 
# e.g. for gcc must "module load openmpi-gcc" to get OpenMPI using gcc
# on tango.vpac.org
MPICC        = mpicc
MPILD        = mpicc

# Requires the CUDA SDK 
CPLUSPLUS    = nvcc
NVCC         = nvcc
CUDA_LD	     = nvcc
CUDA_INCLUDE = /usr/local/cuda/4.2/cuda/include
CUDA_SDK_ROOT = /usr/local/cuda/4.2/cuda/sdk
CUDA_SDK_INCLUDE = $(HOME)/cuda/sdk/C/common/inc
CUDA_SDK_UTIL_LIB = $(HOME)/cuda/sdk/C/lib  # libcutil must be built by user


# requires Boost c++ library
BOOST_INCLUDE = $(HOME)/boost_1_51_0

# NVIDIA CUDA compute capability - we require double precision
# so need compute capability 1.3 at least
GPU_ARCH = sm_13

CUDA_CPPFLAGS   = -I$(CUDA_SDK_INCLUDE) -DCUDA -I$(BOOST_INCLUDE)
HOST_CPPFLAGS   = 
CDEBUG     = -g #-DDEBUG
#OPTFLAGS = -O0 # cannot use -O 1 or 2 or 3 on atreyu gcc 4.5.1 causes bellmanford (only) to be 3 orders of magnitude slower !?!?!?
OPTFLAGS = -O3
NVCCOPTFLAGS = --ptxas-options=-v --use_fast_math -O3  
CFLAGS     = $(OPTFLAGS)

CPPFLAGS = -DUNIX -I$(BOOST_INCLUDE)
LDFLAGS =   #-pg

CDEBUG = -g -pg -DDEBUG
CFLAGS     += -Wall  
#              the following warnings are not implied by -Wall
CFLAGS     += -Wextra -Wfloat-equal  \
              -Wundef -Wshadow \
              -Wpointer-arith -Wcast-qual -Wcast-align\
              -Wwrite-strings \
              -Wmissing-declarations #-Wunreachable-code

ifeq ($(MODE), DEVICEEMU)
  NVCFLAGS   = $(CDEBUG) -deviceemu
else
  ifeq ($(MODE),GPUDEBUG)
    CFLAGS += -g -G
    NVCFLAGS = -g -G
  else
    NVCFLAGS   = $(NVCCOPTFLAGS)
  endif
  ifeq ($(MODE),DEBUG)
    CFLAGS = $(CDEBUG)
    NVCFLAGS = $(CDEBUG) -G
  endif
endif
#NVCFLAGS += -pg


PTHREAD_CFLAGS = -pthread
PTHREAD_LDFLAGS = -pthread

CUDA_LDLIBPATH  = -L$(CUDA_SDK_UTIL_LIB)
CUDA_LDLIBS     = -lcutil_x86_64
CUDA_LDFLAGS = -arch $(GPU_ARCH)

LDLIBPATH  = 
LDLIBS     = -lm



# Program to build TAGS file for EMACS
MAKETAGS   = etags

MAKEDEPEND = gcc -MM -x c++ $(CPPFLAGS)

