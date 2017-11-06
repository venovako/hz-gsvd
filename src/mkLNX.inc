# Assumes HDF5, compiled with Intel compilers, and an FMA-enabled CPU.
RM=rm -rf
AR=xiar

#ifort
ifdef USE_MPI
FC=h5pfc
else # no MPI
FC=h5fc
endif # ?USE_MPI

FORFLAGS=-DUSE_INTEL -qopenmp -fexceptions -reentrancy threaded -standard-semantics
ARFLAGS=-qnoipo -lib rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -xHost
DBGFLAGS=-DNDEBUG -fno-omit-frame-pointer -qopt-report=5 -traceback
FPUFLAGS=-DHAVE_FMA -fp-model source -fma -no-ftz -no-complex-limited-range -prec-div -prec-sqrt
else # DEBUG
OPTFLAGS=-O0 -xHost
DBGFLAGS=-fno-omit-frame-pointer -g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -debug-parameters all -check all -warn all -traceback
FPUFLAGS=-DHAVE_FMA -fp-model strict -assume ieee_fpe_flags -fma -fp-stack-check -no-ftz -no-complex-limited-range -prec-div -prec-sqrt
endif # ?NDEBUG

ifeq ($(USE_MKL),sequential)
MKL_D=-DUSE_MKL -DUSE_LAPACK_SEQ -DMKL_DIRECT_CALL -DMKL_DIRECT_CALL_SEQ
MKL=sequential
else # USE_MKL=parallel
MKL_D=-DUSE_MKL -DUSE_LAPACK_PAR -DMKL_DIRECT_CALL
MKL=parallel
endif # ?USE_MKL

LIBFLAGS=-I../HZ $(MKL_D) -mkl=$(MKL) 
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz
