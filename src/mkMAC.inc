# Assumes Homebrew installations of gfortran and HDF5, and an FMA-enabled CPU.
RM=rm -rf
AR=ar

#gfortran
ifdef USE_MPI
FC=h5pfc
else # no MPI
FC=h5fc
endif # ?USE_MPI

FORFLAGS=-DUSE_GNU -ffree-line-length-none -fopenmp -fexceptions -frecursive
ARFLAGS=rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -march=native -Wa,-q
DBGFLAGS=-DNDEBUG
FPUFLAGS=-DHAVE_FMA
else # DEBUG
OPTFLAGS=-O0 -march=native -Wa,-q
DBGFLAGS=-g -fcheck=all -finit-local-zero -finit-real=snan -finit-derived
FPUFLAGS=-DHAVE_FMA -ffpe-trap=invalid,zero,overflow
endif # ?NDEBUG

ifeq ($(USE_MKL),sequential)
MKL_D=-DUSE_MKL -DUSE_LAPACK_SEQ -DMKL_DIRECT_CALL -DMKL_DIRECT_CALL_SEQ
MKL_l=-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl
else # USE_MKL=parallel
MKL_D=-DUSE_MKL -DUSE_LAPACK_PAR -DMKL_DIRECT_CALL
MKL_l=-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl
endif # ?USE_MKL
MKL_I=-I${MKLROOT}/include
MKL_L=-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -L${MKLROOT}/../compiler/lib -Wl,-rpath,${MKLROOT}/../compiler/lib

LIBFLAGS=-I../HZ $(MKL_D) $(MKL_I)
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz $(MKL_L) $(MKL_l)
