SUBROUTINE MY_DZIMMER0(FAST, M, N, NP, F, LDF, G, LDG, V, LDV, MAXCYC, TOL, H, K, SIGMA, NSWEEP, NROT, INFO)

  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: D_SQRT_EPS_2 = 1.05367121277235087D-008
  DOUBLE PRECISION, PARAMETER :: D_TWO_POW_27 = 1.34217728000000000D+008

  LOGICAL, INTENT(IN) :: FAST
  INTEGER, INTENT(IN) :: M, N, NP, LDF, LDG, LDV, MAXCYC
  DOUBLE PRECISION, INTENT(IN) :: TOL
  DOUBLE PRECISION, INTENT(INOUT) :: F(LDF,*), G(LDG,*)
  DOUBLE PRECISION, INTENT(OUT) :: V(LDV,*), H(*), K(*), SIGMA(*)
  INTEGER, INTENT(OUT) :: NSWEEP, INFO
  INTEGER(8), INTENT(OUT) :: NROT(2)

  INTEGER :: ITER, I, P, PP, Q, QQ
  LOGICAL :: INTRAN
  INTEGER(8) :: NROTIN(2)

  DOUBLE PRECISION :: APP, AQQ, APQ, BPQ
  DOUBLE PRECISION :: BPQP, BPQM, XI, ETA

  DOUBLE PRECISION :: CSFP(8)
  !DIR$ ATTRIBUTES ALIGN: 64:: CSFP
  DOUBLE PRECISION :: COSF, SINF, COSP, SINP
  EQUIVALENCE (CSFP(1), COSF), (CSFP(2), SINF), (CSFP(3), COSP), (CSFP(4), SINP)
  DOUBLE PRECISION :: COT2T, TANT, COST, SINT
  EQUIVALENCE (CSFP(5), COT2T), (CSFP(6), TANT), (CSFP(7), COST), (CSFP(8), SINT)

  DOUBLE PRECISION :: FASTR(8), D, FCT, MYTOL
  !DIR$ ATTRIBUTES ALIGN: 64:: FASTR
  EQUIVALENCE (FASTR(6), D), (FASTR(7), FCT), (FASTR(8), MYTOL)

  DOUBLE PRECISION, EXTERNAL :: DNRM2, DDOT
  EXTERNAL :: DLASET, DROTM

  !DIR$ ASSUME_ALIGNED F:64,G:64,V:64, H:64,K:64,SIGMA:64
  !DIR$ ASSUME (MOD(LDF, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDG, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDV, 8) .EQ. 0)

  ! Assume that argument check passed.
  INFO = 0
  NSWEEP = 0
  NROT = 0_8

  ! V = I_N
  CALL DLASET('A', M, N, D_ZERO, D_ONE, V, LDV)

  ! Subnormal column norms should never happen (unless the column is 0-vector), due to sqrt.
  DO Q = 1, N
     D = DNRM2(M, G(1, Q), 1)
     IF (.NOT. (D .EQ. D)) THEN ! QNaN
        INFO = 3 * Q
        RETURN
     ELSE IF (D .GT. HUGE(D)) THEN ! +Inf
        INFO = 3 * Q + 1
        RETURN
     ELSE IF (D .LT. TINY(D)) THEN ! Subnormal
        INFO = 3 * Q + 2
        RETURN
     END IF

     IF (D .NE. D_ONE) THEN
        CALL DARR_DIV_SCAL(M, G(1, Q), D)
        CALL DARR_DIV_SCAL(M, F(1, Q), D)
        D = D_ONE / D
     END IF
     V(Q, Q) = D

     ! Should we rescale A_q such that ||A_q' * D|| is finite & normalized???
  END DO

  FCT = SCALE(EPSILON(FCT), -1) * SQRT(DBLE(M))
  IF (TOL .EQ. D_MONE) THEN ! Compute own TOL.
     MYTOL = FCT
  ELSE IF (TOL .EQ. D_ZERO) THEN ! May be +0 or -0.
     MYTOL = D_ZERO
  ELSE
     MYTOL = MIN(TOL, FCT)
  END IF

  IF (NP .EQ. 0) THEN
     PP = N - 1
  ELSE
     PP = NP
  END IF

  DO ITER = 1, MAXCYC

     ! Sweep initializations.

     NROTIN = 0_8

     ! Sweep = row by row (top to bottom), row = left to right.

     DO P = 1, PP

        IF (NP .EQ. 0) THEN
           QQ = P + 1
        ELSE
           QQ = NP + 1
        END IF

        DO Q = QQ, N

           ! Pivot position [p, q].

           ! No checking for overflow/underflow of norm and dot product computations is performed.
           ! In order for this to work and to avoid DLASSQ-based DNRM2 (probably not vectorized),
           ! it should always be maintained that ||F_q|| <= SQRT(HUGE), by rescaling the columns
           ! of F when/if necessary.  New squares of the column norms (after transformation) are
           ! available before it, because they are already used for the column sorting/swapping.
           ! However, they might be too imprecise to compute the scaling factors directly from them.
           ! They might only suggest whether any scaling should take place.

           APP = DDOT(M, F(1, P), 1, F(1, P), 1)
           AQQ = DDOT(M, F(1, Q), 1, F(1, Q), 1)
           APQ = DDOT(M, F(1, P), 1, F(1, Q), 1)
           BPQ = DDOT(M, G(1, P), 1, G(1, Q), 1)

           IF (ABS(BPQ) .GE. MYTOL) THEN
              INTRAN = .TRUE.
           ELSE
              INTRAN = (ABS(APQ) .GE. (SQRT(APP) * SQRT(AQQ) * MYTOL))
           END IF

           IF (.NOT. INTRAN) CYCLE

           BPQP = SQRT(D_ONE + BPQ)
           BPQM = SQRT(D_ONE - BPQ)
           XI = BPQ / (BPQP + BPQM)
           ETA = BPQ / ((D_ONE + BPQP) * (D_ONE + BPQM))

           IF (ABS(BPQ) .LT. D_SQRT_EPS_2) THEN
              ! FCT = D_ONE
              COT2T = (AQQ - APP) / MY_DFMA(BPQ, -(APP + AQQ), SCALE(APQ, 1))
              D = ABS(COT2T)
              IF (D .GE. D_TWO_POW_27) THEN
                 TANT = SCALE(D_ONE / COT2T, -1)
              ELSE IF (D .LT. D_SQRT_EPS_2) THEN
                 TANT = SIGN(D_ONE / (D + D_ONE), COT2T)
              ELSE
                 TANT = SIGN(D_ONE / (D + SQRT(MY_DFMA(COT2T, COT2T, D_ONE))), COT2T)
              END IF
              IF (ABS(TANT) .LT. D_SQRT_EPS_2) THEN
                 ! COST = D_ONE
                 ! SINT = TANT
                 COSF = MY_DFMA((TANT - ETA), XI, D_ONE)              ! COSF = D_ONE + XI * (TANT  - ETA)
                 SINF = MY_DFMA(MY_DFMA(TANT, ETA, D_ONE), -XI, TANT) ! SINF = TANT  - XI * (D_ONE + ETA * TANT)
                 COSP = MY_DFMA((TANT + ETA), -XI, D_ONE)             ! COSP = D_ONE - XI * (TANT  + ETA)
                 SINP = MY_DFMA(MY_DFMA(TANT, -ETA, D_ONE), XI, TANT) ! SINP = TANT  + XI * (D_ONE - ETA * TANT)
              ELSE
                 COST = MY_RSQRT(MY_DFMA(TANT, TANT, D_ONE))
                 SINT = COST * TANT
                 COSF = MY_DFMA(MY_DFMA(COST, -ETA, SINT), XI, COST)  ! COSF = COST + XI * (SINT - ETA * COST)
                 SINF = MY_DFMA(MY_DFMA(SINT, ETA, COST), -XI, SINT)  ! SINF = SINT - XI * (COST + ETA * SINT)
                 COSP = MY_DFMA(MY_DFMA(COST, ETA, SINT), -XI, COST)  ! COSP = COST - XI * (SINT + ETA * COST)
                 SINP = MY_DFMA(MY_DFMA(SINT, -ETA, COST), XI, SINT)  ! SINP = SINT + XI * (COST - ETA * SINT)
              END IF
           ELSE
              FCT = SQRT(MY_DFMA(BPQ, -BPQ, D_ONE))
              ! COT2T = ((AQQ - APP) * FCT) / (SCALE(APQ, 1) - (APP + AQQ) * BPQ)
              COT2T = ((AQQ - APP) * FCT) / MY_DFMA(BPQ, -(APP + AQQ), SCALE(APQ, 1))
              D = ABS(COT2T)
              ! TANT = SIGN(D_ONE / (ABS(COT2T) + SQRT(FMA(COT2T,COT2T,D_ONE)), COT2T)
              IF (D .GE. D_TWO_POW_27) THEN
                 TANT = SCALE(D_ONE / COT2T, -1)
              ELSE IF (D .LT. D_SQRT_EPS_2) THEN
                 TANT = SIGN(D_ONE / (D + D_ONE), COT2T)
              ELSE
                 TANT = SIGN(D_ONE / (D + SQRT(MY_DFMA(COT2T, COT2T, D_ONE))), COT2T)
              END IF
              IF (ABS(TANT) .LT. D_SQRT_EPS_2) THEN
                 ! COST = D_ONE
                 ! SINT = TANT
                 COSF = MY_DFMA((TANT - ETA), XI, D_ONE)              ! COSF = (D_ONE + XI * (TANT  - ETA))
                 SINF = MY_DFMA(MY_DFMA(TANT, ETA, D_ONE), -XI, TANT) ! SINF = (TANT  - XI * (D_ONE + ETA * TANT))
                 COSP = MY_DFMA((TANT + ETA), -XI, D_ONE)             ! COSP = (D_ONE - XI * (TANT  + ETA))
                 SINP = MY_DFMA(MY_DFMA(TANT, -ETA, D_ONE), XI, TANT) ! SINP = (TANT  + XI * (D_ONE - ETA * TANT))
              ELSE
                 COST = MY_RSQRT(MY_DFMA(TANT, TANT, D_ONE))
                 SINT = COST * TANT
                 COSF = MY_DFMA(MY_DFMA(COST, -ETA, SINT), XI, COST)  ! COSF = (COST + XI * (SINT - ETA * COST))
                 SINF = MY_DFMA(MY_DFMA(SINT, ETA, COST), -XI, SINT)  ! SINF = (SINT - XI * (COST + ETA * SINT))
                 COSP = MY_DFMA(MY_DFMA(COST, ETA, SINT), -XI, COST)  ! COSP = (COST - XI * (SINT + ETA * COST))
                 SINP = MY_DFMA(MY_DFMA(SINT, -ETA, COST), XI, SINT)  ! SINP = (SINT + XI * (COST - ETA * SINT))
              END IF
              !!! DO I = 1, 4 is really needed, but for 512-bit vectorisation I goes up to 8 
              !DIR$ VECTOR ALWAYS, ALIGNED
              DO I = 1, 8
                 CSFP(I) = CSFP(I) / FCT
              END DO
           END IF
           ! Compute the new ~app,~aqq for sorting.
           APP = COSF*COSF*APP - SCALE(COSF*SINP*APQ, 1) + SINP*SINP*AQQ
           AQQ = SINF*SINF*APP + SCALE(SINF*COSP*APQ, 1) + COSP*COSP*AQQ
           ! Perform the rotation.
           IF (APP .GE. AQQ) THEN
              IF ((COSF .NE. D_ONE) .OR. (COSP .NE. D_ONE)) THEN
                 IF ((SINF .EQ. D_MONE) .AND. (SINP .EQ. D_MONE)) THEN
                    FASTR(1) =  D_ONE
                    FASTR(2) =  COSF
                    FASTR(5) =  COSP
                 ELSE
                    FASTR(1) =  D_MONE
                    FASTR(2) =  COSF
                    FASTR(3) =  SINF
                    FASTR(4) = -SINP
                    FASTR(5) =  COSP
                 END IF
              ELSE IF ((SINF .NE. D_ZERO) .OR. (SINP .NE. D_ZERO)) THEN
                 FASTR(1) =  D_ZERO
                 FASTR(3) =  SINF
                 FASTR(4) = -SINP
              ELSE
                 INTRAN = .FALSE.
              END IF
           ELSE ! Swap the columns.
              IF ((SINF .NE. D_ONE) .OR. (SINP .NE. D_MONE)) THEN
                 IF ((COSF .EQ. D_MONE) .AND. (COSP .EQ. D_ONE)) THEN
                    FASTR(1) =  D_ONE
                    FASTR(2) =  SINF
                    FASTR(5) = -SINP
                 ELSE
                    FASTR(1) =  D_MONE
                    FASTR(2) =  SINF
                    FASTR(3) =  COSF
                    FASTR(4) =  COSP
                    FASTR(5) = -SINP
                 END IF
              ELSE IF ((COSF .NE. D_ZERO) .OR. (COSP .NE. D_ZERO)) THEN
                 FASTR(1) =  D_ZERO
                 FASTR(3) =  COSF
                 FASTR(4) =  COSP
              ELSE
                 INTRAN = .FALSE.
              END IF
           END IF

           IF (INTRAN) THEN
              ! Transform the columns of F and G.
              CALL DROTM(M, F(1, P), 1, F(1, Q), 1, FASTR)
              CALL DROTM(M, G(1, P), 1, G(1, Q), 1, FASTR)
              ! Transform the columns of V.
              CALL DROTM(M, V(1, P), 1, V(1, Q), 1, FASTR)

              NROTIN(1) = NROTIN(1) + 1_8
              IF ((ABS(COSF) .NE. D_ONE) .OR. (ABS(COSP) .NE. D_ONE)) NROTIN(2) = NROTIN(2) + 1_8
           END IF

           ! End of loop for pivot position [p, q].

        END DO
     END DO

     ! Convergence Tests.
     IF (NROTIN(1) .GT. 0) THEN
        NROT(1) = NROT(1) + NROTIN(1)
        NROT(2) = NROT(2) + NROTIN(2)
     ELSE
        ! Slow convergence test - no rotations per sweep.
        NSWEEP = ITER
        EXIT
     END IF

     ! End of iterations (sweeps) loop.

  END DO

  IF (NSWEEP .EQ. 0) THEN
     ! No convergence.
     INFO = 1
     NSWEEP = MAXCYC
  END IF

  IF (FAST) THEN
     ! Normalize V.
     DO Q = 1, N
        H(Q) = DNRM2(M, F(1, Q), 1)
        K(Q) = DNRM2(M, G(1, Q), 1)
        FCT = HYPOT(H(Q), K(Q))
        IF (FCT .NE. D_ONE) CALL DARR_DIV_SCAL(M, V(1, Q), FCT)
     END DO
  ELSE
     DO Q = 1, N
        H(Q) = DNRM2(M, F(1, Q), 1)
        IF (H(Q) .NE. D_ONE) CALL DARR_DIV_SCAL(M, F(1, Q), H(Q))
        K(Q) = DNRM2(M, G(1, Q), 1)
        ! Ideally, K(Q) should be equal to 1.
        IF (K(Q) .NE. D_ONE) THEN
           CALL DARR_DIV_SCAL(M, G(1, Q), K(Q))
           SIGMA(Q) = H(Q) / K(Q)
        ELSE
           SIGMA(Q) = H(Q)
        END IF
        FCT = HYPOT(H(Q), K(Q))
        IF (FCT .NE. D_ONE) THEN
           H(Q) = H(Q) / FCT
           K(Q) = K(Q) / FCT
           CALL DARR_DIV_SCAL(M, V(1, Q), FCT)
        END IF
     END DO
  END IF

END SUBROUTINE MY_DZIMMER0

SUBROUTINE DZIMMER0(M, N, F, LDF, G, LDG, V, LDV, MAXCYC, TOL, H, K, SIGMA, NSWEEP, NROT, INFO)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: M, N, LDF, LDG, LDV, MAXCYC
  DOUBLE PRECISION, INTENT(IN) :: TOL
  DOUBLE PRECISION, INTENT(INOUT) :: F(LDF,*), G(LDG,*)
  DOUBLE PRECISION, INTENT(OUT) :: V(LDV,*), H(*), K(*), SIGMA(*)
  INTEGER, INTENT(OUT) :: NSWEEP, INFO
  INTEGER(8), INTENT(OUT) :: NROT(2)

  INTEGER :: MCYCLE, I
  DOUBLE PRECISION :: MYTOL

  !DIR$ ASSUME_ALIGNED F:64,G:64,V:64, H:64,K:64,SIGMA:64
  !DIR$ ASSUME (MOD(LDF, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDG, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDV, 8) .EQ. 0)

  IF (MAXCYC .EQ. -1) THEN
     MCYCLE = HUGE(MAXCYC)
  ELSE
     MCYCLE = MAXCYC
  END IF

  IF (.NOT. (TOL .EQ. TOL)) THEN ! QNaN
     MYTOL = D_MONE
  ELSE IF (TOL .EQ. D_MONE) THEN ! Compute own TOL.
     MYTOL = HUGE(TOL)
  ELSE IF (TOL .EQ. D_ZERO) THEN ! May be +0 or -0.
     MYTOL = ABS(TOL)
  ELSE
     MYTOL = TOL
  END IF

  IF (M .LT. 0) THEN
     INFO = -1
  ELSE IF (N .LT. 0) THEN
     INFO = -2
  ELSE IF (N .GT. M) THEN
     INFO = -2
  ELSE IF (LDF .LT. M) THEN
     INFO = -4
  ELSE IF (LDG .LT. M) THEN
     INFO = -6
  ELSE IF (LDV .LT. M) THEN
     INFO = -8
  ELSE IF (MCYCLE .LT. 0) THEN
     INFO = -9
  ELSE IF (MYTOL .LT. D_ZERO) THEN
     INFO = -10
  ELSE
     INFO = 0
  END IF

  IF (INFO .NE. 0) RETURN
  IF (N .EQ. 0) RETURN

  CALL MY_DZIMMER0(.FALSE., M, N, 0, F, LDF, G, LDG, V, LDV, MCYCLE, MYTOL, H, K, SIGMA, NSWEEP, NROT, INFO)

END SUBROUTINE DZIMMER0
