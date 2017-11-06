# Assumes HDF5, compiled with IBM XL compilers, and a POWER8 CPU.
RM=rm -rf
AR=ar

FC=h5fc #xlf2008_r

FORFLAGS=-qnosave -qsclk=micro -qsmp=omp #-WF,-qfpp
ARFLAGS=rsv

ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -qmaxmem=-1 -qtune=pwr8:smt8
DBGFLAGS=-WF,-DNDEBUG
FPUFLAGS=-WF,-DHAVE_FMA -qfloat=nans:subnormals
else # DEBUG
OPTFLAGS=-O0 -qmaxmem=-1 -qtune=pwr8:smt8
DBGFLAGS=-g
FPUFLAGS=-WF,-DHAVE_FMA -qfloat=nans:subnormals
endif # ?NDEBUG

ifeq ($(USE_ESSL),smp)
ESSL_D=-WF,-DUSE_ESSL -WF,-DUSE_LAPACK_PAR
ESSL_l=-lesslsmp -lessl
else # USE_ESSL=seq
ESSL_D=-WF,-DUSE_ESSL -WF,-DUSE_LAPACK_SEQ
ESSL_l=-lessl
endif # ?USE_ESSL
LAPACK=$(ESSL_l) -L$(HOME)/lapack -ltmglib $(ESSL_l) -llapack $(ESSL_l) -lrefblas $(ESSL_l)

LIBFLAGS=-I../HZ $(ESSL_D)
FCFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
LDFLAGS=-L../HZ -lhz $(LAPACK)
