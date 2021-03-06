SUBROUTINE TIMER_PRINT(CLK)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: CLK(3)

  INTEGER :: C, Q, R
  INTEGER :: U

  U = GET_IOUNIT('O')
  IF (U .LT. 0) RETURN

  C = CLK(2) - CLK(1)
  Q = C / CLK(3)
  R = MOD(C, CLK(3))

  SELECT CASE (CLK(3))
  CASE (1000)       ! ms
     WRITE (U,1) Q, R
  CASE (1000000)    ! us
     WRITE (U,2) Q, R
  CASE (1000000000) ! ns
     WRITE (U,3) Q, R
  CASE DEFAULT      ! other scale
     WRITE (U,4) Q, R, CLK(3)
  END SELECT

1 FORMAT(I5,'.',I3.3,' s')
2 FORMAT(I5,'.',I6.6,' s')
3 FORMAT(I5,'.',I9.9,' s')
4 FORMAT(I5,'+',I12,'/',I12,' s')

END SUBROUTINE TIMER_PRINT
