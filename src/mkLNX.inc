RM=rm -rf

CC=h5pcc
FC=h5pfc

AR=xiar
ARFLAGS=-qnoipo -lib rsv

CFLAGS=-std=gnu99 -DNDEBUG -O3 -xHost -fp-model source -no-ftz -fma -prec-div -prec-sqrt -fPIC
FORFLAGS=-openmp -xHost -fPIC -standard-semantics -traceback -warn all

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG)
DBGFLAGS=-DNDEBUG -qopt-report=5
FPUFLAGS=-fp-model source -no-ftz -fma -prec-div -prec-sqrt
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-check_mpi -g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -debug-parameters all -check all
FPUFLAGS=-fp-model strict -assume ieee_fpe_flags -no-ftz -fma -prec-div -prec-sqrt -fp-stack-check
endif # NDEBUG

LIBFLAGS=-DUSE_MKL -I../HZ
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz
