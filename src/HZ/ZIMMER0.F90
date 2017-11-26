SUBROUTINE MY_SZIMMER0(FAST, M, N, NP, F, LDF, G, LDG, V, LDV, MAXCYC, TOL, H, K, SIGMA, NSWEEP, NROT, INFO)

  IMPLICIT NONE

  REAL, PARAMETER :: S_SQRT_EPS_2 = 2.44140625E-04
  REAL, PARAMETER :: S_TWO_POW125 = 5.79261914E+03

  LOGICAL, INTENT(IN) :: FAST
  INTEGER, INTENT(IN) :: M, N, NP, LDF, LDG, LDV, MAXCYC
  REAL, INTENT(IN) :: TOL
  REAL, INTENT(INOUT) :: F(LDF,N), G(LDG,N)
  REAL, INTENT(OUT) :: V(LDV,N), H(*), K(*), SIGMA(*)
  INTEGER, INTENT(OUT) :: NSWEEP, NROT(2), INFO

  INTEGER :: ITER, I, P, PP, Q, QQ
  LOGICAL :: INTRAN
  INTEGER :: NROTIN(2)

#ifdef USE_KNC
  REAL :: CSFP(16)
#else
  REAL :: CSFP(4)
#endif
  !DIR$ ATTRIBUTES ALIGN: 64:: CSFP
  REAL :: COSF, SINF, COSP, SINP
  EQUIVALENCE (CSFP(1), COSF), (CSFP(2), SINF), (CSFP(3), COSP), (CSFP(4), SINP)

  REAL :: COT2T, TANT, COST, SINT
  REAL :: APP, AQQ, APQ, BPQ, BPQP, BPQM
  REAL :: FASTR(5), D, FCT, MYTOL, XI, ETA

  REAL, EXTERNAL :: SDOT, SNRM2
  EXTERNAL :: SLASET, SROTM, SSCAL

  !DIR$ ASSUME_ALIGNED F:64,G:64,V:64, H:64,K:64,SIGMA:64
  !DIR$ ASSUME (MOD(LDF, 16) .EQ. 0)
  !DIR$ ASSUME (MOD(LDG, 16) .EQ. 0)
  !DIR$ ASSUME (MOD(LDV, 16) .EQ. 0)

  ! Assume that argument check passed.
  INFO = 0
  NSWEEP = 0
  NROT = 0
  CSFP = S_ZERO

  ! V = I_N
  CALL SLASET('A', M, N, S_ZERO, S_ONE, V, LDV)

  ! Subnormal column norms should never happen (unless the column is 0-vector), due to sqrt.
  DO Q = 1, N
     FCT = SNRM2(M, G(1, Q), 1)
     IF (.NOT. (FCT .EQ. FCT)) THEN ! QNaN
        INFO = 3 * Q
        RETURN
     ELSE IF (FCT .GT. HUGE(FCT)) THEN ! +Inf
        INFO = 3 * Q + 1
        RETURN
     ELSE IF (FCT .LT. TINY(FCT)) THEN ! Subnormal
        INFO = 3 * Q + 2
        RETURN
     ELSE IF (FCT .NE. S_ONE) THEN
        D = S_ONE / FCT
        CALL SSCAL(M, D, G(1, Q), 1)
        CALL SSCAL(M, D, F(1, Q), 1)
     ELSE
        D = S_ONE
     END IF
     V(Q, Q) = D
     ! Should we rescale A_q such that ||A_q' * D|| is finite & normalized???
  END DO

  D = EPSILON(D)
  FCT = SCALE(D, -1) * SQRT(REAL(M))
  IF (TOL .EQ. S_MONE) THEN ! Compute own TOL.
     MYTOL = FCT
  ELSE IF (TOL .EQ. S_ZERO) THEN ! May be +0 or -0.
     MYTOL = S_ZERO
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

     NROTIN = 0

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
           ! In order for this to work and to avoid SLASSQ-based SNRM2 (probably not vectorized),
           ! it should always be maintained that ||F_q|| <= SQRT(HUGE), by rescaling the columns
           ! of F when/if necessary.  New squares of the column norms (after transformation) are
           ! available before it, because they are already used for the column sorting/swapping.
           ! However, they might be too imprecise to compute the scaling factors directly from them.
           ! They might only suggest whether any scaling should take place.

           APP = SDOT(M, F(1, P), 1, F(1, P), 1)
           AQQ = SDOT(M, F(1, Q), 1, F(1, Q), 1)
           APQ = SDOT(M, F(1, P), 1, F(1, Q), 1)
           BPQ = SDOT(M, G(1, P), 1, G(1, Q), 1)

           IF (ABS(BPQ) .GE. MYTOL) THEN
              INTRAN = .TRUE.
           ELSE
              INTRAN = (ABS(APQ) .GE. (SQRT(APP) * SQRT(AQQ) * MYTOL))
           END IF

           IF (.NOT. INTRAN) CYCLE

           BPQP = SQRT(S_ONE + BPQ)
           BPQM = SQRT(S_ONE - BPQ)
           XI = BPQ / (BPQP + BPQM)
           ETA = BPQ / ((S_ONE + BPQP) * (S_ONE + BPQM))

           IF (ABS(BPQ) .LT. S_SQRT_EPS_2) THEN
              ! FCT = S_ONE
              COT2T = (AQQ - APP) / SFMA(BPQ, -(APP + AQQ), SCALE(APQ, 1))
              D = ABS(COT2T)
              IF (D .GE. S_TWO_POW125) THEN
                 TANT = SCALE(S_ONE / COT2T, -1)
              ELSE IF (D .LT. S_SQRT_EPS_2) THEN
                 TANT = SIGN(S_ONE / (D + S_ONE), COT2T)
              ELSE
                 TANT = SIGN(S_ONE / (D + SQRT(SFMA(COT2T, COT2T, S_ONE))), COT2T)
              END IF
              IF (ABS(TANT) .LT. S_SQRT_EPS_2) THEN
                 ! COST = S_ONE
                 ! SINT = TANT
                 COSF = SFMA((TANT - ETA), XI, S_ONE)           ! COSF = S_ONE + XI * (TANT  - ETA)
                 SINF = SFMA(SFMA(TANT, ETA, S_ONE), -XI, TANT) ! SINF = TANT  - XI * (S_ONE + ETA * TANT)
                 COSP = SFMA((TANT + ETA), -XI, S_ONE)          ! COSP = S_ONE - XI * (TANT  + ETA)
                 SINP = SFMA(SFMA(TANT, -ETA, S_ONE), XI, TANT) ! SINP = TANT  + XI * (S_ONE - ETA * TANT)
              ELSE
                 COST = SRSQRT(SFMA(TANT, TANT, S_ONE))
                 SINT = COST * TANT
                 COSF = SFMA(SFMA(COST, -ETA, SINT), XI, COST)  ! COSF = COST + XI * (SINT - ETA * COST)
                 SINF = SFMA(SFMA(SINT, ETA, COST), -XI, SINT)  ! SINF = SINT - XI * (COST + ETA * SINT)
                 COSP = SFMA(SFMA(COST, ETA, SINT), -XI, COST)  ! COSP = COST - XI * (SINT + ETA * COST)
                 SINP = SFMA(SFMA(SINT, -ETA, COST), XI, SINT)  ! SINP = SINT + XI * (COST - ETA * SINT)
              END IF
           ELSE
              FCT = SQRT(SFMA(BPQ, -BPQ, S_ONE))
              ! COT2T = ((AQQ - APP) * FCT) / (SCALE(APQ, 1) - (APP + AQQ) * BPQ)
              COT2T = ((AQQ - APP) * FCT) / SFMA(BPQ, -(APP + AQQ), SCALE(APQ, 1))
              D = ABS(COT2T)
              ! TANT = SIGN(S_ONE / (ABS(COT2T) + SQRT(FMA(COT2T,COT2T,S_ONE)), COT2T)
              IF (D .GE. S_TWO_POW125) THEN
                 TANT = SCALE(S_ONE / COT2T, -1)
              ELSE IF (D .LT. S_SQRT_EPS_2) THEN
                 TANT = SIGN(S_ONE / (D + S_ONE), COT2T)
              ELSE
                 TANT = SIGN(S_ONE / (D + SQRT(SFMA(COT2T, COT2T, S_ONE))), COT2T)
              END IF
              IF (ABS(TANT) .LT. S_SQRT_EPS_2) THEN
                 ! COST = S_ONE
                 ! SINT = TANT
                 COSF = SFMA((TANT - ETA), XI, S_ONE)           ! COSF = (S_ONE + XI * (TANT  - ETA))
                 SINF = SFMA(SFMA(TANT, ETA, S_ONE), -XI, TANT) ! SINF = (TANT  - XI * (S_ONE + ETA * TANT))
                 COSP = SFMA((TANT + ETA), -XI, S_ONE)          ! COSP = (S_ONE - XI * (TANT  + ETA))
                 SINP = SFMA(SFMA(TANT, -ETA, S_ONE), XI, TANT) ! SINP = (TANT  + XI * (S_ONE - ETA * TANT))
              ELSE
                 COST = SRSQRT(SFMA(TANT, TANT, S_ONE))
                 SINT = COST * TANT
                 COSF = SFMA(SFMA(COST, -ETA, SINT), XI, COST)  ! COSF = (COST + XI * (SINT - ETA * COST))
                 SINF = SFMA(SFMA(SINT, ETA, COST), -XI, SINT)  ! SINF = (SINT - XI * (COST + ETA * SINT))
                 COSP = SFMA(SFMA(COST, ETA, SINT), -XI, COST)  ! COSP = (COST - XI * (SINT + ETA * COST))
                 SINP = SFMA(SFMA(SINT, -ETA, COST), XI, SINT)  ! SINP = (SINT + XI * (COST - ETA * SINT))
              END IF
#ifdef USE_KNC
              D = S_ONE / FCT
              !DIR$ VECTOR ALWAYS, ALIGNED
              DO I = 1, 8
                 CSFP(I) = CSFP(I) * D
              END DO
#else
              !DIR$ VECTOR ALWAYS, ALIGNED
              DO I = 1, 4
                 CSFP(I) = CSFP(I) / FCT
              END DO
#endif
           END IF
           ! Compute the new ~app,~aqq for sorting.
           APP = COSF*COSF*APP - SCALE(COSF*SINP*APQ, 1) + SINP*SINP*AQQ
           AQQ = SINF*SINF*APP + SCALE(SINF*COSP*APQ, 1) + COSP*COSP*AQQ
           ! Perform the rotation.
           IF (APP .GE. AQQ) THEN
              IF ((COSF .NE. S_ONE) .OR. (COSP .NE. S_ONE)) THEN
                 IF ((SINF .EQ. S_MONE) .AND. (SINP .EQ. S_MONE)) THEN
                    FASTR(1) =  S_ONE
                    FASTR(2) =  COSF
                    FASTR(5) =  COSP
                 ELSE
                    FASTR(1) =  S_MONE
                    FASTR(2) =  COSF
                    FASTR(3) =  SINF
                    FASTR(4) = -SINP
                    FASTR(5) =  COSP
                 END IF
              ELSE IF ((SINF .NE. S_ZERO) .OR. (SINP .NE. S_ZERO)) THEN
                 FASTR(1) =  S_ZERO
                 FASTR(3) =  SINF
                 FASTR(4) = -SINP
              ELSE
                 INTRAN = .FALSE.
              END IF
           ELSE ! Swap the columns.
              IF ((SINF .NE. S_ONE) .OR. (SINP .NE. S_MONE)) THEN
                 IF ((COSF .EQ. S_MONE) .AND. (COSP .EQ. S_ONE)) THEN
                    FASTR(1) =  S_ONE
                    FASTR(2) =  SINF
                    FASTR(5) = -SINP
                 ELSE
                    FASTR(1) =  S_MONE
                    FASTR(2) =  SINF
                    FASTR(3) =  COSF
                    FASTR(4) =  COSP
                    FASTR(5) = -SINP
                 END IF
              ELSE IF ((COSF .NE. S_ZERO) .OR. (COSP .NE. S_ZERO)) THEN
                 FASTR(1) =  S_ZERO
                 FASTR(3) =  COSF
                 FASTR(4) =  COSP
              ELSE
                 INTRAN = .FALSE.
              END IF
           END IF

           IF (INTRAN) THEN
              ! Transform the columns of F and G.
              CALL SROTM(M, F(1, P), 1, F(1, Q), 1, FASTR)
              CALL SROTM(M, G(1, P), 1, G(1, Q), 1, FASTR)
              ! Transform the columns of V.
              CALL SROTM(M, V(1, P), 1, V(1, Q), 1, FASTR)

              NROTIN(1) = NROTIN(1) + 1
              IF ((ABS(COSF) .NE. S_ONE) .OR. (ABS(COSP) .NE. S_ONE)) NROTIN(2) = NROTIN(2) + 1
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
        FCT = HYPOT(SNRM2(M, F(1, Q), 1), SNRM2(M, G(1, Q), 1))
        IF (FCT .NE. S_ONE) THEN
           D = S_ONE / FCT
           CALL SSCAL(M, D, V(1, Q), 1)
        END IF
     END DO
  ELSE
     DO Q = 1, N
        H(Q) = SNRM2(M, F(1, Q), 1)
        IF (H(Q) .NE. S_ONE) THEN
           D = S_ONE / H(Q)
           CALL SSCAL(M, D, F(1, Q), 1)
        END IF
        K(Q) = SNRM2(M, G(1, Q), 1)
        ! Ideally, K(Q) should be equal to 1.
        IF (K(Q) .NE. S_ONE) THEN
           D = S_ONE / K(Q)
           SIGMA(Q) = H(Q) * D
           CALL SSCAL(M, D, G(1, Q), 1)
        ELSE
           SIGMA(Q) = H(Q)
        END IF
        FCT = HYPOT(H(Q), K(Q))
        IF (FCT .NE. S_ONE) THEN
           D = S_ONE / FCT
           H(Q) = H(Q) * D
           K(Q) = K(Q) * D
           CALL SSCAL(M, D, V(1, Q), 1)
        END IF
     END DO
  END IF

END SUBROUTINE MY_SZIMMER0

SUBROUTINE SZIMMER0(M, N, F, LDF, G, LDG, V, LDV, MAXCYC, TOL, H, K, SIGMA, NSWEEP, NROT, INFO)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: M, N, LDF, LDG, LDV, MAXCYC
  REAL, INTENT(IN) :: TOL
  REAL, INTENT(INOUT) :: F(LDF,N), G(LDG,N)
  REAL, INTENT(OUT) :: V(LDV,N), H(*), K(*), SIGMA(*)
  INTEGER, INTENT(OUT) :: NSWEEP, NROT(2), INFO

  INTEGER :: MCYCLE, I
  REAL :: MYTOL

  !DIR$ ASSUME_ALIGNED F:64,G:64,V:64, H:64,K:64,SIGMA:64
  !DIR$ ASSUME (MOD(LDF, 16) .EQ. 0)
  !DIR$ ASSUME (MOD(LDG, 16) .EQ. 0)
  !DIR$ ASSUME (MOD(LDV, 16) .EQ. 0)

  IF (MAXCYC .EQ. -1) THEN
     MCYCLE = HUGE(MAXCYC)
  ELSE
     MCYCLE = MAXCYC
  END IF

  IF (.NOT. (TOL .EQ. TOL)) THEN ! QNaN
     MYTOL = S_MONE
  ELSE IF (TOL .EQ. S_MONE) THEN ! Compute own TOL.
     MYTOL = HUGE(TOL)
  ELSE IF (TOL .EQ. S_ZERO) THEN ! May be +0 or -0.
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
  ELSE IF (MYTOL .LT. S_ZERO) THEN
     INFO = -10
  ELSE
     INFO = 0
  END IF

  IF (INFO .NE. 0) RETURN
  IF (N .EQ. 0) RETURN

  CALL MY_SZIMMER0(.FALSE., M, N, 0, F, LDF, G, LDG, V, LDV, MCYCLE, MYTOL, H, K, SIGMA, NSWEEP, NROT, INFO)

END SUBROUTINE SZIMMER0

SUBROUTINE MY_DZIMMER0(FAST, M, N, NP, F, LDF, G, LDG, V, LDV, MAXCYC, TOL, H, K, SIGMA, NSWEEP, NROT, INFO)

  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: D_SQRT_EPS_2 = 1.05367121277235087D-008
  DOUBLE PRECISION, PARAMETER :: D_TWO_POW_27 = 1.34217728000000000D+008

  LOGICAL, INTENT(IN) :: FAST
  INTEGER, INTENT(IN) :: M, N, NP, LDF, LDG, LDV, MAXCYC
  DOUBLE PRECISION, INTENT(IN) :: TOL
  DOUBLE PRECISION, INTENT(INOUT) :: F(LDF,N), G(LDG,N)
  DOUBLE PRECISION, INTENT(OUT) :: V(LDV,N), H(*), K(*), SIGMA(*)
  INTEGER, INTENT(OUT) :: NSWEEP, NROT(2), INFO

  INTEGER :: ITER, I, P, PP, Q, QQ
  LOGICAL :: INTRAN
  INTEGER :: NROTIN(2)

#ifdef USE_KNC
  DOUBLE PRECISION :: CSFP(8)
#else
  DOUBLE PRECISION :: CSFP(4)
#endif
  !DIR$ ATTRIBUTES ALIGN: 64:: CSFP
  DOUBLE PRECISION :: COSF, SINF, COSP, SINP
  EQUIVALENCE (CSFP(1), COSF), (CSFP(2), SINF), (CSFP(3), COSP), (CSFP(4), SINP)

  DOUBLE PRECISION :: COT2T, TANT, COST, SINT
  DOUBLE PRECISION :: APP, AQQ, APQ, BPQ, BPQP, BPQM
  DOUBLE PRECISION :: FASTR(5), D, FCT, MYTOL, XI, ETA

  DOUBLE PRECISION, EXTERNAL :: DDOT, DNRM2
  EXTERNAL :: DLASET, DROTM, DSCAL

  !DIR$ ASSUME_ALIGNED F:64,G:64,V:64, H:64,K:64,SIGMA:64
  !DIR$ ASSUME (MOD(LDF, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDG, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDV, 8) .EQ. 0)

  ! Assume that argument check passed.
  INFO = 0
  NSWEEP = 0
  NROT = 0
  CSFP = D_ZERO

  ! V = I_N
  CALL DLASET('A', M, N, D_ZERO, D_ONE, V, LDV)

  ! Subnormal column norms should never happen (unless the column is 0-vector), due to sqrt.
  DO Q = 1, N
     FCT = DNRM2(M, G(1, Q), 1)
     IF (.NOT. (FCT .EQ. FCT)) THEN ! QNaN
        INFO = 3 * Q
        RETURN
     ELSE IF (FCT .GT. HUGE(FCT)) THEN ! +Inf
        INFO = 3 * Q + 1
        RETURN
     ELSE IF (FCT .LT. TINY(FCT)) THEN ! Subnormal
        INFO = 3 * Q + 2
        RETURN
     ELSE IF (FCT .NE. D_ONE) THEN
        D = D_ONE / FCT
        CALL DSCAL(M, D, G(1, Q), 1)
        CALL DSCAL(M, D, F(1, Q), 1)
     ELSE
        D = D_ONE
     END IF
     V(Q, Q) = D
     ! Should we rescale A_q such that ||A_q' * D|| is finite & normalized???
  END DO

  D = EPSILON(D)
  FCT = SCALE(D, -1) * SQRT(DBLE(M))
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

     NROTIN = 0

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
              COT2T = (AQQ - APP) / DFMA(BPQ, -(APP + AQQ), SCALE(APQ, 1))
              D = ABS(COT2T)
              IF (D .GE. D_TWO_POW_27) THEN
                 TANT = SCALE(D_ONE / COT2T, -1)
              ELSE IF (D .LT. D_SQRT_EPS_2) THEN
                 TANT = SIGN(D_ONE / (D + D_ONE), COT2T)
              ELSE
                 TANT = SIGN(D_ONE / (D + SQRT(DFMA(COT2T, COT2T, D_ONE))), COT2T)
              END IF
              IF (ABS(TANT) .LT. D_SQRT_EPS_2) THEN
                 ! COST = D_ONE
                 ! SINT = TANT
                 COSF = DFMA((TANT - ETA), XI, D_ONE)           ! COSF = D_ONE + XI * (TANT  - ETA)
                 SINF = DFMA(DFMA(TANT, ETA, D_ONE), -XI, TANT) ! SINF = TANT  - XI * (D_ONE + ETA * TANT)
                 COSP = DFMA((TANT + ETA), -XI, D_ONE)          ! COSP = D_ONE - XI * (TANT  + ETA)
                 SINP = DFMA(DFMA(TANT, -ETA, D_ONE), XI, TANT) ! SINP = TANT  + XI * (D_ONE - ETA * TANT)
              ELSE
                 COST = DRSQRT(DFMA(TANT, TANT, D_ONE))
                 SINT = COST * TANT
                 COSF = DFMA(DFMA(COST, -ETA, SINT), XI, COST)  ! COSF = COST + XI * (SINT - ETA * COST)
                 SINF = DFMA(DFMA(SINT, ETA, COST), -XI, SINT)  ! SINF = SINT - XI * (COST + ETA * SINT)
                 COSP = DFMA(DFMA(COST, ETA, SINT), -XI, COST)  ! COSP = COST - XI * (SINT + ETA * COST)
                 SINP = DFMA(DFMA(SINT, -ETA, COST), XI, SINT)  ! SINP = SINT + XI * (COST - ETA * SINT)
              END IF
           ELSE
              FCT = SQRT(DFMA(BPQ, -BPQ, D_ONE))
              ! COT2T = ((AQQ - APP) * FCT) / (SCALE(APQ, 1) - (APP + AQQ) * BPQ)
              COT2T = ((AQQ - APP) * FCT) / DFMA(BPQ, -(APP + AQQ), SCALE(APQ, 1))
              D = ABS(COT2T)
              ! TANT = SIGN(D_ONE / (ABS(COT2T) + SQRT(FMA(COT2T,COT2T,D_ONE)), COT2T)
              IF (D .GE. D_TWO_POW_27) THEN
                 TANT = SCALE(D_ONE / COT2T, -1)
              ELSE IF (D .LT. D_SQRT_EPS_2) THEN
                 TANT = SIGN(D_ONE / (D + D_ONE), COT2T)
              ELSE
                 TANT = SIGN(D_ONE / (D + SQRT(DFMA(COT2T, COT2T, D_ONE))), COT2T)
              END IF
              IF (ABS(TANT) .LT. D_SQRT_EPS_2) THEN
                 ! COST = D_ONE
                 ! SINT = TANT
                 COSF = DFMA((TANT - ETA), XI, D_ONE)           ! COSF = (D_ONE + XI * (TANT  - ETA))
                 SINF = DFMA(DFMA(TANT, ETA, D_ONE), -XI, TANT) ! SINF = (TANT  - XI * (D_ONE + ETA * TANT))
                 COSP = DFMA((TANT + ETA), -XI, D_ONE)          ! COSP = (D_ONE - XI * (TANT  + ETA))
                 SINP = DFMA(DFMA(TANT, -ETA, D_ONE), XI, TANT) ! SINP = (TANT  + XI * (D_ONE - ETA * TANT))
              ELSE
                 COST = DRSQRT(DFMA(TANT, TANT, D_ONE))
                 SINT = COST * TANT
                 COSF = DFMA(DFMA(COST, -ETA, SINT), XI, COST)  ! COSF = (COST + XI * (SINT - ETA * COST))
                 SINF = DFMA(DFMA(SINT, ETA, COST), -XI, SINT)  ! SINF = (SINT - XI * (COST + ETA * SINT))
                 COSP = DFMA(DFMA(COST, ETA, SINT), -XI, COST)  ! COSP = (COST - XI * (SINT + ETA * COST))
                 SINP = DFMA(DFMA(SINT, -ETA, COST), XI, SINT)  ! SINP = (SINT + XI * (COST - ETA * SINT))
              END IF
#ifdef USE_KNC
              D = D_ONE / FCT
              !DIR$ VECTOR ALWAYS, ALIGNED
              DO I = 1, 8
                 CSFP(I) = CSFP(I) * D
              END DO
#else
              !DIR$ VECTOR ALWAYS, ALIGNED
              DO I = 1, 4
                 CSFP(I) = CSFP(I) / FCT
              END DO
#endif
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

              NROTIN(1) = NROTIN(1) + 1
              IF ((ABS(COSF) .NE. D_ONE) .OR. (ABS(COSP) .NE. D_ONE)) NROTIN(2) = NROTIN(2) + 1
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
        FCT = HYPOT(DNRM2(M, F(1, Q), 1), DNRM2(M, G(1, Q), 1))
        IF (FCT .NE. D_ONE) THEN
           D = D_ONE / FCT
           CALL DSCAL(M, D, V(1, Q), 1)
        END IF
     END DO
  ELSE
     DO Q = 1, N
        H(Q) = DNRM2(M, F(1, Q), 1)
        IF (H(Q) .NE. D_ONE) THEN
           D = D_ONE / H(Q)
           CALL DSCAL(M, D, F(1, Q), 1)
        END IF
        K(Q) = DNRM2(M, G(1, Q), 1)
        ! Ideally, K(Q) should be equal to 1.
        IF (K(Q) .NE. D_ONE) THEN
           D = D_ONE / K(Q)
           SIGMA(Q) = H(Q) * D
           CALL DSCAL(M, D, G(1, Q), 1)
        ELSE
           SIGMA(Q) = H(Q)
        END IF
        FCT = HYPOT(H(Q), K(Q))
        IF (FCT .NE. D_ONE) THEN
           D = D_ONE / FCT
           H(Q) = H(Q) * D
           K(Q) = K(Q) * D
           CALL DSCAL(M, D, V(1, Q), 1)
        END IF
     END DO
  END IF

END SUBROUTINE MY_DZIMMER0

SUBROUTINE DZIMMER0(M, N, F, LDF, G, LDG, V, LDV, MAXCYC, TOL, H, K, SIGMA, NSWEEP, NROT, INFO)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: M, N, LDF, LDG, LDV, MAXCYC
  DOUBLE PRECISION, INTENT(IN) :: TOL
  DOUBLE PRECISION, INTENT(INOUT) :: F(LDF,N), G(LDG,N)
  DOUBLE PRECISION, INTENT(OUT) :: V(LDV,N), H(*), K(*), SIGMA(*)
  INTEGER, INTENT(OUT) :: NSWEEP, NROT(2), INFO

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
