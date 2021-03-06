# Assumes HDF5, compiled with GNU compilers.
RM=rm -rf
AR=ar

#gfortran
FC=h5fc

FORFLAGS=-DUSE_GNU -fdefault-integer-8 -ffree-line-length-none -fopenmp -fexceptions -frecursive -fstack-arrays
ARFLAGS=rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -march=native -fgcse-las -fgcse-sm -fipa-pta -ftree-loop-distribution -ftree-loop-im -ftree-loop-ivcanon -fivopts -fvect-cost-model=unlimited -fvariable-expansion-in-unroller
DBGFLAGS=-DNDEBUG -fopt-info-optimized-vec -pedantic
FPUFLAGS=-ffp-contract=fast
else # DEBUG
OPTFLAGS=-O0 -march=native
DBGFLAGS=-g -fcheck=all -finit-local-zero -finit-real=snan -finit-derived -pedantic
FPUFLAGS=-ffp-contract=fast -ffpe-trap=invalid,zero,overflow
endif # ?NDEBUG

ifeq ($(USE_MKL),sequential)
MKL_D=-DUSE_MKL -DMKL_ILP64 -DUSE_LAPACK_SEQ -DMKL_DIRECT_CALL -DMKL_DIRECT_CALL_SEQ
MKL_l=-lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
else # USE_MKL=parallel
MKL_D=-DUSE_MKL -DMKL_ILP64 -DUSE_LAPACK_PAR -DMKL_DIRECT_CALL
MKL_l=-lmkl_gf_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
endif # ?USE_MKL
MKL_I=-I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
MKL_L=-L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -L${MKLROOT}/../compiler/lib/intel64 -Wl,-rpath=${MKLROOT}/../compiler/lib/intel64 -Wl,--no-as-needed

LIBFLAGS=-I../HZ $(MKL_D) $(MKL_I)
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz $(MKL_L) $(MKL_l)
