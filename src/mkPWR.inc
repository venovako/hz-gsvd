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
#LAPACK=-L$(HOME)/lapack -ltmglib -llapack -lrefblas
LAPACK=-lessl -L$(HOME)/lapack -ltmglib -lessl -llapack -lessl -lrefblas -lessl
LDFLAGS=-L../HZ -lhz $(LAPACK)
