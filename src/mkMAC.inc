# Assumes Homebrew installations of gfortran and HDF5, and an FMA-enabled CPU.
RM=rm -rf
AR=ar

FC=h5fc #gfortran

FORFLAGS=-DUSE_GNU -fopenmp -fexceptions -frecursive
ARFLAGS=rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -march=native -Wa,-q
DBGFLAGS=-DNDEBUG
FPUFLAGS=-DHAVE_FMA
else # DEBUG
OPTFLAGS=-O0 -march=native -Wa,-q
DBGFLAGS=-g -fcheck=all -finit-local-zero -finit-real=snan -finit-derived
FPUFLAGS=-DHAVE_FMA -ffpe-trap=invalid,zero,overflow
endif # ?NDEBUG

LIBFLAGS=-I../HZ -DUSE_MKL -I${MKLROOT}/include
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -L${MKLROOT}/../compiler/lib -Wl,-rpath,${MKLROOT}/../compiler/lib -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl
