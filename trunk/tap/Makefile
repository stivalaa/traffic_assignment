##############################################################################
#
# File      : Makefile
# Author    : Alex Stivala (astivala)
# Created   : February 2009
#
# $Id: Makefile 776 2011-10-03 04:27:19Z astivala $
#
# The default target "all" will compile and link the program.
#
# The "clean" target will remove object and executable files.
#
# set MODE=DEBUG to build with debugging and verbose printing on host,
# default is to build with optimizations on and no debug or profile
#
#
#
##############################################################################

include common.mk
-include local.mk


DEPENDFILE = .depend

COMMON_C_OBJS =  parsetapfiles.o tap_functions.o utils.o
COMMON_C_SRCS = $(COMMON_C_OBJS:.o=.c)
PAPE_C_OBJS =    sssp_pape.o
PAPE_C_SRCS = $(PAPE_C_OBJS:.o=.c)
GREEDY_C_OBJS =    tap_greedy.o
GREEDY_C_SRCS = $(GREEDY_C_OBJS:.o=.c)
FW_C_OBJS =    tap_frankwolfe.o
FW_C_SRCS = $(FW_C_OBJS:.o=.c)
FW_PTHREAD_C_OBJS = tap_frankwolfe_pthread.o
FW_PTHREAD_C_SRCS = $(FW_PTHREAD_C_OBJS:.o=.c)
FW_PTHREAD_MAIN_C_OBJS =  tap_main.o
FW_PTHREAD_MAIN_C_SRCS = $(FW_PTHREAD_MAIN_C_OBJS:.o=.c)
FW_MPI_C_OBJS = tap_mpi_main.o
FW_MPI_C_SRCS = $(FW_MPI_C_OBJS:.o=.c)

FW_CUDA_MAIN_C_OBJS = tap_cuda_main.o
FW_CUDA_MAIN_C_SRCS = $(FW_CUDA_MAIN_C_OBJS:.o=.cu)
FW_CUDA_C_OBJS = tap_frankwolfe_cuda.o
FW_CUDA_C_SRCS = $(FW_CUDA_C_OBJS:.o=.cu)

CUDA_KERNELS  = pape_kernels.cu volume_update_kernels.cu
CUDA_HOSTCODE = pape_cuda_host.cu volume_update_host.cu

CUDA_OBJS     = $(CUDA_KERNELS:.cu=.o)
CUDA_HOST_OBJS= $(CUDA_HOSTCODE:.cu=.o)

OBJS = $(COMMON_C_OBJS) $(PAPE_C_OBJS) $(FW_C_OBJS) $(FW_PTHREAD_C_OBJS) $(GREEDY_C_OBJS) $(FW_PTHREAD_MAIN_C_OBJS) $(FW_MPI_C_OBJS) $(CUDA_OBJS) $(CUDA_HOST_OBJS) $(FW_CUDA_MAIN_C_OBJS) $(FW_CUDA_C_OBJS)
SRCS = $(COMMON_C_SRCS) $(PAPE_C_SRCS) $(FW_C_SRCS) $(FW_PTHREAD_C_SRCS) $(GREEDY_C_SRCS) $(FW_PTHREAD_MAIN_C_SRCS) $(FW_MPI_C_SRCS) $(CUDA_KERNELS) $(CUDA_HOSTCODE) $(FW_CUDA_MAIN_C_SRCS) $(FW_CUDA_C_SRCS)

all: tap_greedy tap_frankwolfe tap_frankwolfe_pthread tap_frankwolfe_mpi tap_frankwolfe_cuda

tap_greedy: $(COMMON_C_OBJS)  $(PAPE_C_OBJS) $(GREEDY_C_OBJS)
	$(LD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

tap_frankwolfe: $(COMMON_C_OBJS)  $(PAPE_C_OBJS) $(FW_C_OBJS)
	$(LD) $(LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS)

tap_frankwolfe_pthread: $(COMMON_C_OBJS)  $(PAPE_C_OBJS)  $(FW_PTHREAD_C_OBJS) $(FW_PTHREAD_MAIN_C_OBJS)
	$(LD) $(LDFLAGS) $(PTHREAD_LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS) 

tap_frankwolfe_mpi: $(COMMON_C_OBJS)  $(PAPE_C_OBJS)  $(FW_PTHREAD_C_OBJS) $(FW_MPI_C_OBJS)
	$(MPILD) $(LDFLAGS) $(PTHREAD_LDFLAGS) $(LDLIBPATH)  -o $@ $^  $(LDLIBS) 

tap_frankwolfe_cuda: $(FW_CUDA_MAIN_C_OBJS) $(FW_CUDA_C_OBJS) $(CUDA_HOST_OBJS) $(CUDA_OBJS)  $(COMMON_C_OBJS)  $(PAPE_C_OBJS)
	$(CUDA_LD) $(CUDA_LDFLAGS)  $(CUDA_LDLIBPATH) -o $@ $^ $(CUDA_LDLIBS)

tags:
	$(MAKETAGS) $(SRCS)

depend: $(SRCS)
	$(MAKEDEPEND) $(SRCS) > $(DEPENDFILE)

clean:
	rm -f $(OBJS)
	rm -f tap_greedy tap_frankwolfe tap_frankwolfe_pthread


tap_frankwolfe_pthread.o: tap_frankwolfe_pthread.c
	$(CC) $(PTHREAD_CFLAGS) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

tap_mpi_main.o: tap_mpi_main.c
	$(MPICC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

%.o : %.cu
	$(NVCC) -c $(NVCFLAGS) $(CUDA_CPPFLAGS) -arch $(GPU_ARCH) $< -o $@



###############################################################################
#
# scripts for generating results stats
# 
###############################################################################


R       = R --vanilla --slave

stats: error_stats.txt

error_stats.txt: compare_individual_with_subsets.rtab	\
                 compare_pairwise_oneonly.rtab \
	         compare_pairwise_with_subsets.rtab \
	         get_relative_error_stats.r
	$(R) -f get_relative_error_stats.r > $@

compare_individual_with_subsets.rtab: ../scripts/compare_pairwise_with_subsets.py  ChicagoRegional_mods.txt results/ChicagoRegional_fw_out.txt results/ChicagoRegional_mods_goliath.err /dev/null results/ChicagoRegional_mods_subsets.err
	../scripts/compare_pairwise_with_subsets.py -p ChicagoRegional_mods.txt results/ChicagoRegional_fw_out.txt results/ChicagoRegional_mods_goliath.err /dev/null results/ChicagoRegional_mods_subsets.err | sed 's/froz/"froz/g' | sed 's/[)]/\)"/g'  > $@

compare_pairwise_with_subsets.rtab: ../scripts/compare_pairwise_with_subsets.py  ChicagoRegional_mods.txt results/ChicagoRegional_fw_out.txt results/ChicagoRegional_mods_goliath.err results/ChicagoRegional_pairwise.err results/ChicagoRegional_mods_subsets.err make_subsets_rtab.sh
	./make_subsets_rtab.sh > $@

compare_pairwise_oneonly.rtab: ../scripts/compare_pairwise_with_subsets.py   ChicagoRegional_mods.txt results/ChicagoRegional_fw_out.txt results/ChicagoRegional_mods_goliath.err ChicagoRegional_only1pairwise.err results/ChicagoRegional_mods_subsets.err
	../scripts/compare_pairwise_with_subsets.py -i  ChicagoRegional_only1pairwise.err ChicagoRegional_mods.txt results/ChicagoRegional_fw_out.txt results/ChicagoRegional_mods_goliath.err results/ChicagoRegional_pairwise.err results/ChicagoRegional_mods_subsets.err | sed 's/froz/"froz/g' | sed 's/[)]/\)"/g'  > $@

compare_pairwise_individual.rtab: ../scripts/compare_pairwise_with_subsets.py  ChicagoRegional_mods.txt results/ChicagoRegional_fw_out.txt results/ChicagoRegional_mods_goliath.err results/ChicagoRegional_mods_subsets.err
	../scripts/compare_pairwise_with_subsets.py  ChicagoRegional_mods.txt results/ChicagoRegional_fw_out.txt results/ChicagoRegional_mods_goliath.err /dev/null results/ChicagoRegional_mods_subsets.err | sed 's/froz/"froz/g' | sed 's/[)]/\)"/g'  > $@


ChicagoRegional_xy.csv: make_xy_csv_file_from_rtabs.sh compare_pairwise_linkflow_summary.rtab
	./make_xy_csv_file_from_rtabs.sh > $@


ChicagoRegional_centroid_distance.rtab: compare_pairwise_linkflow_summary.rtab
	grep -v '^#' compare_pairwise_linkflow_summary.rtab |awk 'NR > 1' | cut -f1,2,7|sort -k3,3nr  > $@

###############################################################################

include $(DEPENDFILE)


