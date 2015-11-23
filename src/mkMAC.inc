RM=rm -rf

FC=h5pfc
AR=ar

FORFLAGS=-fopenmp
ARFLAGS=rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG)
DBGFLAGS=-DNDEBUG
FPUFLAGS=-DHAVE_FMA
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-g -fcheck=all -finit-local-zero -finit-real=snan
FPUFLAGS=-DHAVE_FMA -ffpe-trap=invalid,zero,overflow
endif # NDEBUG

LIBFLAGS=-DUSE_OPENBLAS -I../HZ
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz -L/usr/local/lib -ltmglib -llapack -lopenblas
