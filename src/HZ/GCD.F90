PURE INTEGER FUNCTION GCD(A, B)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: A, B

  INTEGER :: AA, BB, CC

  AA = ABS(A)
  BB = ABS(B)

  IF (AA .LT. BB) THEN
     CC = AA
     AA = BB
     BB = CC
  END IF

  DO WHILE (BB .GT. 0)
     AA = AA - BB
     IF (AA .LT. BB) THEN
        CC = AA
        AA = BB
        BB = CC
     END IF
  END DO

  GCD = AA

END FUNCTION GCD
