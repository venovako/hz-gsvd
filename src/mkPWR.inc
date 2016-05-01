HDF5=$(HOME)/hdf5_ibm-nompi
H5FC=$(HDF5)/bin/h5fc #h5pfc
H5CC=$(HDF5)/bin/h5cc #h5pcc

RM=rm -rf

FC=$(H5FC)
AR=ar

FORFLAGS=-qfloat=nans:subnormals -qmaxmem=-1 -qnosave -qsclk=micro -qsmp=omp -qtune=pwr8:smt8 # -WF,-qfpp
ARFLAGS=rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG)
DBGFLAGS=-WF,-DNDEBUG
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-g
endif # NDEBUG

LIBFLAGS=-I../HZ
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS)
#LAPACK=-L$(HOME)/lapack -ltmglib -llapack -lrefblas
LAPACK=-lessl -L$(HOME)/lapack -ltmglib -lessl -llapack -lessl -lrefblas -lessl
LDFLAGS=-L../HZ -lhz $(LAPACK)
