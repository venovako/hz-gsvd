# Assumes HDF5, compiled with GNU compilers, and an FMA-enabled CPU.
RM=rm -rf
AR=ar

FC=h5fc #gfortran

FORFLAGS=-DUSE_GNU -fopenmp -fexceptions -frecursive
ARFLAGS=rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -march=native
DBGFLAGS=-DNDEBUG
FPUFLAGS=-DHAVE_FMA
else # DEBUG
OPTFLAGS=-O0 -march=native
DBGFLAGS=-g -fcheck=all -finit-local-zero -finit-real=snan #-finit-derived
FPUFLAGS=-ffpe-trap=invalid,zero,overflow
endif # ?NDEBUG

LIBFLAGS=-I../HZ -DUSE_MKL -I${MKLROOT}/include
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -L${MKLROOT}/../compiler/lib/intel64 -Wl,-rpath=${MKLROOT}/../compiler/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
