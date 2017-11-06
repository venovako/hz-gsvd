# Assumes HDF5, compiled with GNU compilers.
RM=rm -rf
AR=ar

FC=h5fc #gfortran

FORFLAGS=-DUSE_GNU -ffree-line-length-none -fopenmp -fexceptions -frecursive
ARFLAGS=rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -march=native
DBGFLAGS=-DNDEBUG
FPUFLAGS= #-DHAVE_FMA
else # DEBUG
OPTFLAGS=-O0 -march=native
DBGFLAGS=-g -fcheck=all -finit-local-zero -finit-real=snan #-finit-derived
FPUFLAGS=-ffpe-trap=invalid,zero,overflow #-DHAVE_FMA
endif # ?NDEBUG

ifeq ($(USE_MKL),sequential)
MKL_D=-DUSE_MKL -DUSE_LAPACK_SEQ -DMKL_DIRECT_CALL -DMKL_DIRECT_CALL_SEQ
MKL_l=-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
else # USE_MKL=parallel
MKL_D=-DUSE_MKL -DUSE_LAPACK_PAR -DMKL_DIRECT_CALL
MKL_l=-lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
endif # ?USE_MKL
MKL_I=-I${MKLROOT}/include
MKL_L=-L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -L${MKLROOT}/../compiler/lib/intel64 -Wl,-rpath=${MKLROOT}/../compiler/lib/intel64 -Wl,--no-as-needed

LIBFLAGS=-I../HZ $(MKL_D) $(MKL_I)
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz $(MKL_L) $(MKL_l)
