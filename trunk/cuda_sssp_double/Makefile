###############################################################################
#
# File      : Makefile
# Author    : Alex Stivala (astivala)
# Created   : February 2009
#
# $Id: Makefile 844 2011-11-10 02:58:03Z astivala $
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

include common.mk
-include local.mk

# location of modified Sparse BLAS reference implementation
SPBLAS_DIR = spblas_0_8

#location of modified CUSP library
CUSP_DIR = cusp

CPPFLAGS += -I$(SPBLAS_DIR)
CUDA_LDLIBPATH += -L$(SPBLAS_DIR)
CUDA_LDLIBS += -lsparseblas -lm
CUDA_CPPFLAGS += -I$(CUSP_DIR)

CUDA_KERNELS  = all_kernels.cu
CUDA_HOSTCODE = sssp.cu pape_cuda_host.cu  bellmanford_cusp.cu
CPPSRCS         = sssp_gold.cpp 
CSRCS          = sssp_pape.c bellmanford_spblas.c
HDRS          = sssp_gold.h sssp.h 
CUDA_OBJS     = $(CUDA_KERNELS:.cu=.o)
CUDA_HOST_OBJS= $(CUDA_HOSTCODE:.cu=.o) bellmanford_cusp_gpu.o
CPPOBJS         = $(CPPSRCS:.cpp=.o)
COBJS          = $(CSRCS:.c=.o)
OBJS          = $(CUDA_OBJS) $(CPPOBJS) $(CUDA_HOST_OBJS) $(COBJS)

SRCS	= $(CPPSRCS) $(CUDA_HOSTCODE)


DEPENDFILE = .depend


all: sssp

sssp: $(OBJS)
	$(CUDA_LD) -o sssp $(OBJS) $(CUDA_LDFLAGS) $(CUDA_LDLIBPATH) $(CUDA_LDLIBS)


depend: $(SRCS) $(HDRS)
	$(MAKEDEPEND) $(SRCS) $(HDRS) > $(DEPENDFILE)

# make the TAGS file for EMACS
tags:
	$(MAKETAGS) $(SRCS) $(HDRS)


clean:
	$(RM) $(OBJS)
	$(RM) sssp


all_kernels.o: all_kernels.cu pape_kernels.cu

# tricky: bellmanford_cusp.cu is compiled to two difference objcts,
# bellmanford_cusp.o with CUSP_DEVICE undefined to make host (no device) versions
# and bellmanford_cusp_gpu.o to make version using GPU device
bellmanford_cusp_gpu.o: bellmanford_cusp.cu
	$(NVCC) -c $(NVCFLAGS) $(CUDA_CPPFLAGS) -DCUSP_DEVICE -arch $(GPU_ARCH) $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

%.o : %.cu
	$(NVCC) -c $(NVCFLAGS) $(CUDA_CPPFLAGS) -arch $(GPU_ARCH) $< -o $@

%.o : %.cpp
	$(CPLUSPLUS) -c $(NVCFLAGS) $(CUDA_CPPFLAGS) -o $@ $<




include $(DEPENDFILE)


