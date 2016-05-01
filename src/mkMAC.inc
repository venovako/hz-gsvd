HDF5=$(HOME)/hdf5_gnu-nompi
H5FC=$(HDF5)/bin/h5fc #h5pfc
H5CC=$(HDF5)/bin/h5cc #h5pcc

RM=rm -rf

FC=$(H5FC)
AR=ar

FORFLAGS=-DUSE_GNU -fopenmp -fexceptions -frecursive -march=native -Wa,-q
C99FLAGS=-fexceptions -march=native -ffp-contract=on -integrated-as
ARFLAGS=rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG)
DBGFLAGS=-DNDEBUG
FPUFLAGS=
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-g -fcheck=all -finit-local-zero -finit-real=snan
FPUFLAGS=-ffpe-trap=invalid,zero,overflow
endif # NDEBUG

LIBFLAGS=-I../HZ
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz -L/usr/local/lib -ltmglib -llapack -lrefblas
