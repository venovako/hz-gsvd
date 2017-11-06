RM=rm -rf
AR=ar

FC=h5fc #xlf2008_r

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
LAPACK=-lessl -L$(HOME)/lapack -ltmglib -lessl -llapack -lessl -lrefblas -lessl
LDFLAGS=-L../HZ -lhz $(LAPACK)
