! Intended for 256-bit vectorization.
! Change to DO = 1, N, 8 and DO J = 0, 7 for 512-bit.
PURE SUBROUTINE DARR_DIV_SCAL(N, A, S)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(INOUT) :: A(*)
  DOUBLE PRECISION, INTENT(IN) :: S

  INTEGER :: I, J

  !DIR$ ASSUME_ALIGNED A:64

  DO I = 1, N, 4
     !DIR$ VECTOR ALWAYS, ALIGNED
     DO J = 0, 3
        A(I + J) = A(I + J) / S
     END DO
  END DO

END SUBROUTINE DARR_DIV_SCAL
