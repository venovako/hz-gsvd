HDF5=$(HOME)/hdf5_gnu-nompi
H5FC=$(HDF5)/bin/h5fc #h5pfc

RM=rm -rf

FC=$(H5FC)
AR=ar

FORFLAGS=-DUSE_GNU -fopenmp -fexceptions -frecursive
ARFLAGS=rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -march=native -Wa,-q
DBGFLAGS=-DNDEBUG
FPUFLAGS=
else # DEBUG
OPTFLAGS=-O0 -march=native -Wa,-q
DBGFLAGS=-g -fcheck=all -finit-local-zero -finit-real=snan
FPUFLAGS=-ffpe-trap=invalid,zero,overflow
endif # ?NDEBUG

LIBFLAGS=-I../HZ
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz -L/usr/local/lib -ltmglib -llapack -lrefblas
