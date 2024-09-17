# Assumes HDF5, compiled with Intel compilers, and an FMA-enabled CPU.
RM=rm -rf
AR=xiar

#ifort
FC=h5fc
ifndef CPU
CPU=Host
# COMMON-AVX512 for KNLs
endif # !CPU
ARFLAGS=-qnoipo -lib rsv
# assumes that h5fc sets at least: -i8 -standard-semantics
CPUFLAGS=-DUSE_INTEL -x$(CPU) -fPIC -fexceptions -fasynchronous-unwind-tables -fno-omit-frame-pointer -qopt-multi-version-aggressive -qopt-zmm-usage=high -qopenmp -traceback -reentrancy threaded
FPUFLAGS=-fp-model precise -fma -fprotect-parens -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -qsimd-honor-fp-model -qsimd-serialize-fp-reduction
GDBFLAGS=-g -debug extended -debug inline-debug-info -debug pubnames -debug parallel -debug-parameters all

ifndef CPU
CPU=Host
endif # !CPU
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) $(CPUFLAGS) -fno-math-errno -inline-level=2
DBGFLAGS=-DNDEBUG -qopt-report=5
ifdef PROFILE
DBGFLAGS += -DPROFILE $(GDBFLAGS) 
endif # PROFILE
FPUFLAGS =-fp-model precise -fma -fprotect-parens -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -qsimd-honor-fp-model -qsimd-serialize-fp-reduction
else # DEBUG
OPTFLAGS=-O0 $(CPUFLAGS)
DBGFLAGS=$(GDBFLAGS) -check all -warn all
FPUFLAGS += -assume ieee_fpe_flags -fp-stack-check
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
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz $(MKL_L) $(MKL_l)
