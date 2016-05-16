HDF5=$(HOME)/hdf5_intel-nompi
H5FC=$(HDF5)/bin/h5fc #h5pfc

RM=rm -rf

FC=$(H5FC)
AR=xiar

FORFLAGS=-DUSE_KNC -mmic -qopenmp -fexceptions -reentrancy threaded -standard-semantics
ARFLAGS=-qnoipo -lib rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG)
DBGFLAGS=-DNDEBUG -fno-omit-frame-pointer -qopt-report=5 -traceback
FPUFLAGS=-fp-model source -fma -no-ftz -no-complex-limited-range -prec-div -prec-sqrt
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-fno-omit-frame-pointer -g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -debug-parameters all -check all -warn all -traceback
FPUFLAGS=-fp-model strict -assume ieee_fpe_flags -fma -fp-stack-check -no-ftz -no-complex-limited-range -prec-div -prec-sqrt
endif # ?NDEBUG

LIBFLAGS=-DUSE_MKL -I../HZ
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz
