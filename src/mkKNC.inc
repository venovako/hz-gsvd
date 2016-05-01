HDF5=$(HOME)/hdf5_intel-nompi
H5FC=$(HDF5)/bin/h5fc #h5pfc
H5CC=$(HDF5)/bin/h5cc #h5pcc

RM=rm -rf

FC=$(H5FC)
AR=xiar

FORFLAGS=-DHAVE_KNC -qopenmp -fexceptions -fno-omit-frame-pointer -reentrancy threaded -standard-semantics -traceback
C99FLAGS=-DHAVE_KNC -qopenmp -fexceptions -fno-omit-frame-pointer -std=gnu99 -traceback
ARFLAGS=-qnoipo -lib rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG)
DBGFLAGS=-DNDEBUG -qopt-report=5
FPUFLAGS=-fp-model source -fma -no-ftz -no-complex-limited-range -prec-div -prec-sqrt
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -debug-parameters all -check all -warn all
FPUFLAGS=-fp-model strict -assume ieee_fpe_flags -fma -fp-stack-check -no-ftz -no-complex-limited-range -prec-div -prec-sqrt
endif # NDEBUG

LIBFLAGS=-DUSE_MKL -I../HZ
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz
