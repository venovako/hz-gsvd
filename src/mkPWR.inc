RM=rm -rf

FC=h5fc
AR=ar

FORFLAGS=-qsmp=omp # -WF,-qfpp
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
LDFLAGS=-L../HZ -lhz -L$(HOME)/lapack -ltmglib -llapack -lrefblas
