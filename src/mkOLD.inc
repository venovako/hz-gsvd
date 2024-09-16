# Assumes HDF5, compiled with Intel compilers, and an FMA-enabled CPU.
RM=rm -rf
AR=xiar

#ifort
FC=h5fc

FORFLAGS=-DUSE_INTEL -i8 -qopenmp -reentrancy threaded -standard-semantics
ARFLAGS=-qnoipo -lib rsv

ifndef CPU
CPU=Host
endif # !CPU
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -x$(CPU) -qopt-zmm-usage=high
ifdef PROFILE
DBGFLAGS=-DNDEBUG -DPROFILE -fno-omit-frame-pointer -g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -debug-parameters all -qopt-report=5 -traceback
else # no profiling
DBGFLAGS=-DNDEBUG -fno-omit-frame-pointer -qopt-report=5 -traceback
endif # ?PROFILE
FPUFLAGS=-fp-model precise -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt
else # DEBUG
OPTFLAGS=-O0 -x$(CPU) -qopt-zmm-usage=high
DBGFLAGS=-fno-omit-frame-pointer -g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -debug-parameters all -check all -warn all -traceback
FPUFLAGS=-fp-model strict -assume ieee_fpe_flags -fma -fp-stack-check -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt
endif # ?NDEBUG

ifeq ($(USE_MKL),sequential)
MKL_D=-DUSE_MKL -DMKL_ILP64 -DUSE_LAPACK_SEQ -DMKL_DIRECT_CALL -DMKL_DIRECT_CALL_SEQ
MKL_l=-lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
else # USE_MKL=parallel
MKL_D=-DUSE_MKL -DMKL_ILP64 -DUSE_LAPACK_PAR -DMKL_DIRECT_CALL
MKL_l=-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
endif # ?USE_MKL
MKL_I=-I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
MKL_L=-L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -L${MKLROOT}/../compiler/lib/intel64 -Wl,-rpath=${MKLROOT}/../compiler/lib/intel64 -Wl,--no-as-needed

LIBFLAGS=-I../HZ $(MKL_D) $(MKL_I)
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz $(MKL_L) $(MKL_l)
