# Assumes HDF5, compiled with Intel compilers.
RM=rm -rf
AR=xiar

FC=h5fc #ifort

FORFLAGS=-qopenmp -fexceptions -reentrancy threaded -standard-semantics
ARFLAGS=-qnoipo -lib rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -xHost
DBGFLAGS=-DNDEBUG -fno-omit-frame-pointer -qopt-report=5 -traceback
FPUFLAGS=-fp-model source -fma -no-ftz -no-complex-limited-range -prec-div -prec-sqrt
else # DEBUG
OPTFLAGS=-O0 -xHost
DBGFLAGS=-fno-omit-frame-pointer -g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -debug-parameters all -check all -warn all -traceback
FPUFLAGS=-fp-model strict -assume ieee_fpe_flags -fma -fp-stack-check -no-ftz -no-complex-limited-range -prec-div -prec-sqrt
endif # ?NDEBUG

LIBFLAGS=-I../HZ -DUSE_MKL -mkl #-DMKL_DIRECT_CALL
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz
