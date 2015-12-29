RM=rm -rf

CC=h5pcc
FC=h5pfc
AR=xiar

C99FLAGS=-DHAVE_PHI -qopenmp -std=gnu99 -fexceptions -traceback
FORFLAGS=-DHAVE_PHI -qopenmp -standard-semantics -fexceptions -traceback
ARFLAGS=-qnoipo -lib rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -fno-omit-frame-pointer
DBGFLAGS=-DNDEBUG -qopt-report=5
FPUFLAGS=-DUSE_DIV -fp-model source -no-ftz -fma -prec-div -prec-sqrt
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -debug-parameters all -check all -warn all
FPUFLAGS=-DUSE_DIV -fp-model strict -assume ieee_fpe_flags -fp-stack-check -no-ftz -fma -prec-div -prec-sqrt
endif # NDEBUG

LIBFLAGS=-DUSE_MKL -I../HZ
CCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C99FLAGS) $(FPUFLAGS)
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz
