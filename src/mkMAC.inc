RM=rm -rf

CC=h5pcc
FC=h5pfc
AR=ar

C99FLAGS=-fexceptions
FORFLAGS=-fopenmp
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
CCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C99FLAGS) $(FPUFLAGS)
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz -L/usr/local/lib -ltmglib -llapack -lrefblas
