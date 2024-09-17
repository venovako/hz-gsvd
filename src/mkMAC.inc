# Assumes installations of gfortran and HDF5, and an FMA-enabled CPU.
RM=rm -rf
AR=ar

#gfortran
FC=h5fc
# assumes that h5fc sets at least: -fdefault-integer-8
FORFLAGS=-DUSE_GNU -fPIC -fexceptions -fasynchronous-unwind-tables -fno-omit-frame-pointer -fvect-cost-model=unlimited -ffree-line-length-none -fopenmp -fstack-arrays
ARFLAGS=rsv

ifndef CPU
CPU=native
endif # !CPU
CPUFLAGS=-march=$(CPU) -Wa,-q
FPUFLAGS=-ffp-contract=fast

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) $(CPUFLAGS)
DBGFLAGS=-DNDEBUG
else # DEBUG
OPTFLAGS=-O0 $(CPUFLAGS)
DBGFLAGS=-g -fcheck=all -finit-local-zero -finit-real=snan -finit-derived
FPUFLAGS += -ffpe-trap=invalid,zero,overflow
endif # ?NDEBUG

ifeq ($(USE_MKL),sequential)
MKL_D=-DUSE_MKL -DMKL_ILP64 -DUSE_LAPACK_SEQ -DMKL_DIRECT_CALL -DMKL_DIRECT_CALL_SEQ
MKL_l=-lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
else # USE_MKL=parallel
MKL_D=-DUSE_MKL -DMKL_ILP64 -DUSE_LAPACK_PAR -DMKL_DIRECT_CALL
MKL_l=-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lgomp -lpthread -lm -ldl
endif # ?USE_MKL
MKL_I=-I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
MKL_L=-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -L${MKLROOT}/../compiler/lib -Wl,-rpath,${MKLROOT}/../compiler/lib

LIBFLAGS=-I../HZ $(MKL_D) $(MKL_I)
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz $(MKL_L) $(MKL_l)
