RM=rm -rf

FC=h5pfc
AR=xiar

FORFLAGS=-openmp -standard-semantics -traceback
ARFLAGS=-qnoipo -lib rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) # -xHost
DBGFLAGS=-DNDEBUG -qopt-report=5
FPUFLAGS=-DHAVE_FMA -fp-model source -no-ftz -fma -prec-div -prec-sqrt
else # DEBUG
OPTFLAGS=-O0 # -xHost
DBGFLAGS=-g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -debug-parameters all -check all -warn all
FPUFLAGS=-DHAVE_FMA -fp-model strict -assume ieee_fpe_flags -no-ftz -fma -prec-div -prec-sqrt -fp-stack-check
endif # NDEBUG

LIBFLAGS=-DUSE_MKL -I../HZ
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz
