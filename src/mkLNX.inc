RM=rm -rf

FC=h5pfc
AR=xiar

FORFLAGS=-qopenmp -standard-semantics -fexceptions -traceback # -xHost
ARFLAGS=-qnoipo -lib rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG)
DBGFLAGS=-DNDEBUG -qopt-report=5
FPUFLAGS=-fp-model source -DHAVE_FMA -no-ftz -fma -prec-div -prec-sqrt
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -debug-parameters all -check all -warn all
FPUFLAGS=-fp-model strict -assume ieee_fpe_flags -fp-stack-check -DHAVE_FMA -no-ftz -fma -prec-div -prec-sqrt
endif # NDEBUG

LIBFLAGS=-DUSE_MKL -I../HZ
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz
