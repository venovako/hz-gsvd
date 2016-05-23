RM=rm -rf
FC=xlf2008_r
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

HDF5=$(HOME)/hdf5_ibm-nompi
LIBFLAGS=-I../HZ -I$(HDF5)/include
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS)
#LAPACK=-L$(HOME)/lapack -ltmglib -llapack -lrefblas
LAPACK=-lessl -L$(HOME)/lapack -ltmglib -lessl -llapack -lessl -lrefblas -lessl
LDFLAGS=-L../HZ -lhz -L$(HDF5)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 $(LAPACK)
