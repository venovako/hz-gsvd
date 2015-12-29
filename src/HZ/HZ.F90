MODULE HZ

  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: D_ZERO  =  0.0D0
  DOUBLE PRECISION, PARAMETER :: D_MZERO = -0.0D0
  DOUBLE PRECISION, PARAMETER :: D_HALF  =  0.5D0
  DOUBLE PRECISION, PARAMETER :: D_ONE   =  1.0D0
  DOUBLE PRECISION, PARAMETER :: D_MONE  = -1.0D0
  DOUBLE PRECISION, PARAMETER :: D_TWO   =  2.0D0
  DOUBLE PRECISION, PARAMETER :: D_MTWO  = -2.0D0
  DOUBLE PRECISION, PARAMETER :: D_FOUR  =  4.0D0

CONTAINS

  PURE DOUBLE PRECISION FUNCTION MY_DFMA(A, B, C)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: A, B, C
    !DIR$ FMA
    MY_DFMA = A * B + C
  END FUNCTION MY_DFMA

! Threading support.
#include "INIT_THRS.F90"
#include "BLAS_SET_NUM_THREADS.F90"
#include "GET_NTHR.F90"

! I/O support.
#include "GET_IOUNIT.F90"
#include "DUMPMTX.F90"

! Alignment support.
#include "GCD.F90"
#include "LCM.F90"
#include "LDALIGN.F90"

! Timing support.
#include "TIMER_START.F90"
#include "TIMER_STOP.F90"
#include "TIMER_PRINT.F90"
#include "TIMER2DBLE.F90"

! Shuffler.
#include "INITSH.F90"
#include "BACKSH.F90"

! Partitioner.
#include "PARTBL.F90"

! Stepper.
#include "STRINI.F90"
#include "MMSTEP.F90"

! Vectorized scaling.
#include "DARR_DIV_SCAL.F90"
#include "DARR_DIV_SCPY.F90"
#include "DARR_MUL_SCAL.F90"
#include "DARR_MUL_SCPY.F90"

! Hari-Zimmermann, Level 0.
#include "ZIMMER0.F90"

! Hari-Zimmermann, Level 1.
#include "ZIMMER1.F90"

! Hari-Zimmermann, Level 2.
#include "ZIMMER2.F90"

END MODULE HZ
