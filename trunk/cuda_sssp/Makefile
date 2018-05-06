###############################################################################
#
# File      : Makefile
# Author    : Alex Stivala (astivala)
# Created   : February 2009
#
# $Id: Makefile 227 2011-04-18 04:52:01Z astivala $
#
# Makefile for building CUDA programs
#
#
# The default target "all" will compile and link the program.
#
# The "clean" target will remove object and executable files.
#
# set MODE=DEVICEEMU to build with device emulation and debugging on. 
# default is to build with optimizations on and no debug or profile, for GPU
# set MODE=GPUDEBUG to build with debugging on GPU with cuda-gdb enabled.
# Otherwise,
# set MODE=DEBUG to build with debugging and verbose printing on host,
# default is to build with optimizations on and no debug or profile, for GPU
# execution.
#
#
###############################################################################


# Requires the CUDA SDK 
CUDA_INCLUDE = /usr/local/cuda/4.2/cuda/include
CUDA_SDK_ROOT = /usr/local/cuda/4.2/cuda/sdk
CUDA_SDK_INCLUDE = $(CUDA_SDK_ROOT)/C/common/inc
CUDA_SDK_UTIL_LIB = $(HOME)/cuda/sdk/C/lib  # libcutil must be built by user

CUDA_KERNELS  = all_kernels.cu
CUDA_HOSTCODE = sssp.cu harish_host.cu okuyama_host.cu pape_cuda_host.cu
CPPSRCS         = sssp_gold.cpp 
CSRCS          = sssp_pape.c
HDRS          = harish_kernels.h okuyama_kernels.h sssp_gold.h sssp.h \
               harish_host.h harish_kernels.h okuyama_kernels.h okuyama_host.h \
               atomic_devfunc.h sssp_pape.h pape_kernels.h pape_cuda_host.h
CUDA_OBJS     = $(CUDA_KERNELS:.cu=.o)
CUDA_HOST_OBJS= $(CUDA_HOSTCODE:.cu=.o)
CPPOBJS         = $(CPPSRCS:.cpp=.o)
COBJS          = $(CSRCS:.c=.o)
OBJS          = $(CUDA_OBJS) $(CPPOBJS) $(CUDA_HOST_OBJS) $(COBJS)

SRCS	= $(CPPSRCS) $(CUDA_HOSTCODE)

CC           = gcc
CPLUSPLUS    = nvcc
NVCC         = nvcc
#CPPFLAGS   = -D_XOPEN_SOURCE=1 -D_XOPEN_SOURCE_EXTENDED=1


CUDA_CPPFLAGS   = -I$(CUDA_SDK_INCLUDE) -DCUDA -I$(HOME)/boost_1_51_0
HOST_CPPFLAGS   = 
CDEBUG     = -g #-DDEBUG
OPTFLAGS   = -O3
NVCCOPTFLAGS = --ptxas-options=-v --use_fast_math -O3  
CFLAGS     = $(OPTFLAGS)
#CFLAGS     = $(OPTFLAGS) -Wall -pg
##              the following warnings are not implied by -Wall
#CFLAGS     += -Wextra -Wfloat-equal  \
#              -Wundef -Wshadow \
#              -Wpointer-arith -Wcast-qual -Wcast-align\
#              -Wwrite-strings \
#              -Wmissing-declarations -Wunreachable-code
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


# NVIDIA CUDA compute capability - we do not require double precision
# but now we want to compile with cpabilty 2.0 for Fermi architecture
#GPU_ARCH = sm_20
# we need to use atomicCAS on 64 bits requiring compute capability 
# at least 1.2
# but have found that if we have 1.2 or 1.3, the we require -G0 or -G1 
# or the harish/okuyama code (atomic functions) doesn't work,
# using arch 2.0 solves this problem (don't know why though)
GPU_ARCH = sm_20
#GPU_ARCH = sm_12
#GPU_ARCH = sm_13


LD         = nvcc
LDFLAGS    = -arch $(GPU_ARCH)
ifeq ($(MODE),DEVICEEMU)
  LDFLAGS += -deviceemu
else
  ifeq ($(MODE),GPUDEBUG)
    LDFLAGS += -g -G
  endif
  ifeq ($(MODE),DEBUG)
    LDFLAGS += -g -G
   endif
endif

#LDFLAGS += -pg

LDLIBPATH  = -L$(CUDA_SDK_UTIL_LIB)
LDLIBS     = -lcutil_x86_64

MAKEDEPEND = gcc -MM $(CPPFLAGS)
DEPENDFILE = .depend

# Program to build TAGS file for EMACS
MAKETAGS   = etags

all: sssp

sssp: $(OBJS)
	$(LD) -o sssp $(OBJS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)


depend: $(SRCS) $(HDRS)
	$(MAKEDEPEND) $(SRCS) $(HDRS) > $(DEPENDFILE)

# make the TAGS file for EMACS
tags:
	$(MAKETAGS) $(SRCS) $(HDRS)


clean:
	$(RM) $(OBJS)
	$(RM) sssp


all_kernels.o: all_kernels.cu atomic_devfunc.h atomic_devfunc_types.h  harish_kernels.cu okuyama_kernels.cu pape_kernels.cu

%.o : %.cu
	$(NVCC) -c $(NVCFLAGS) $(CUDA_CPPFLAGS) -arch $(GPU_ARCH) $< -o $@

%.o : %.cpp
	$(CPLUSPLUS) -c $(CFLAGS) $(CUDA_CPPFLAGS) -o $@ $<

%.o : %.c
	$(CC) -c $(CFLAGS) $(CUDA_CPPFLAGS) -o $@ $<



include $(DEPENDFILE)


