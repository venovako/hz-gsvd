# Assumes HDF5, compiled with GNU compilers.
RM=rm -rf
AR=ar

#gfortran
FC=h5fc
# assumes that h5fc sets at least: -fdefault-integer-8
FORFLAGS=-DUSE_GNU -fPIC -fexceptions -fasynchronous-unwind-tables -fno-omit-frame-pointer -fvect-cost-model=unlimited -ffree-line-length-none -fopenmp -fstack-arrays
ARFLAGS=rsv

ifndef MARCH
ifeq ($(shell uname -m),ppc64le)
MARCH=mcpu=native -mpower8-fusion -mtraceback=full
else # !ppc64le
MARCH=march=native
endif # ?ppc64le
endif # !MARCH
CPUFLAGS=-$(MARCH)
FPUFLAGS=-ffp-contract=fast

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) $(CPUFLAGS) -fno-math-errno
DBGFLAGS=-DNDEBUG
else # DEBUG
OPTFLAGS=-Og $(CPUFLAGS)
DBGFLAGS=-g -fcheck=all -finit-local-zero -finit-real=snan -finit-derived
FPUFLAGS += #-ffpe-trap=invalid,zero,overflow
endif # ?NDEBUG

MKL_D=-DMKL_ILP64
ifdef MKLROOT
ifeq ($(USE_MKL),sequential)
MKL_D += -DUSE_MKL -DUSE_LAPACK_SEQ -DMKL_DIRECT_CALL -DMKL_DIRECT_CALL_SEQ
MKL_l=-Wl,--start-group ${MKLROOT}/lib/libmkl_gf_ilp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
else # USE_MKL=parallel
MKL_D += -DUSE_MKL -DUSE_LAPACK_PAR -DMKL_DIRECT_CALL
MKL_l=-Wl,--start-group ${MKLROOT}/lib/libmkl_gf_ilp64.a ${MKLROOT}/lib/libmkl_gnu_thread.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
endif # ?USE_MKL
MKL_I=-I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
else # !MKLROOT
MKL_D += -DUSE_LAPACK_SEQ
MKL_l=-L$(HOME)/lapack-ilp64 -ltmglib -llapack -lrefblas
endif # ?MKLROOT

LIBFLAGS=-I../HZ $(MKL_D) $(MKL_I)
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz $(MKL_l) $(shell if [ -L /usr/lib64/libmemkind.so ]; then echo '-lmemkind'; fi)
