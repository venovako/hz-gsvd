HDF5=$(HOME)/hdf5_ibm-nompi
H5FC=$(HDF5)/bin/h5fc #h5pfc

RM=rm -rf

FC=$(H5FC)
AR=ar

FORFLAGS=-qfloat=nans:subnormals -qnosave -qsclk=micro -qsmp=omp # -WF,-qfpp
ARFLAGS=rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -qmaxmem=-1 -qtune=pwr8:smt8
DBGFLAGS=-WF,-DNDEBUG
else # DEBUG
OPTFLAGS=-O0 -qmaxmem=-1 -qtune=pwr8:smt8
DBGFLAGS=-g
endif # ?NDEBUG

LIBFLAGS=-I../HZ
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS)
#LAPACK=-L$(HOME)/lapack -ltmglib -llapack -lrefblas
LAPACK=-lessl -L$(HOME)/lapack -ltmglib -lessl -llapack -lessl -lrefblas -lessl
LDFLAGS=-L../HZ -lhz $(LAPACK)
