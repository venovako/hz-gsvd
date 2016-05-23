RM=rm -rf
FC=gfortran
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

HDF5=$(HOME)/hdf5_gnu-nompi
LIBFLAGS=-I../HZ -I$(HDF5)/include
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz -L$(HDF5)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -L/usr/local/lib -ltmglib -llapack -lrefblas
