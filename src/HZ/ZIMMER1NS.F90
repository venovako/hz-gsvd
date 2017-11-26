SUBROUTINE MY_DZIMMER1NS(FAST, M, N, F, LDF, G, LDG, V, LDV, MAXCYC, TOL, LDAC, NBMAXS, IPART,&
     H, K, SIGMA, NC, IFC, FB, GB, VB, NBSIZE, NSWEEP, NROT, INFO)

  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: FAST
  INTEGER, INTENT(IN) :: M, N, LDF, LDG, LDV, MAXCYC, LDAC, NBMAXS, IPART
  DOUBLE PRECISION, INTENT(IN) :: TOL
  DOUBLE PRECISION, INTENT(INOUT) :: F(LDF,N), G(LDG,N)
  DOUBLE PRECISION, INTENT(OUT) :: V(LDV,N), H(*), K(*), SIGMA(*)
  INTEGER, INTENT(OUT) :: NC(*), IFC(*), NBSIZE, NSWEEP, NROT(2), INFO(2)
  DOUBLE PRECISION, INTENT(OUT) :: FB(LDAC,*), GB(LDAC,*), VB(LDAC,*)

  INTEGER :: ITER, Q
  LOGICAL :: TRANSF

  INTEGER :: NBL, LDAC5
  INTEGER :: IBL, NC1, IFE1
  INTEGER :: JBL, NC2, IFE2
  INTEGER :: NCF, IFRS1, IFRS2
  INTEGER :: ITEMP
  INTEGER :: NROTIN(2)
  DOUBLE PRECISION :: D, FCT

  DOUBLE PRECISION, EXTERNAL :: DNRM2
  EXTERNAL :: DGEMM, DLACPY, DLASET, DPOTRF, DSCAL, DSYRK

  !DIR$ ASSUME_ALIGNED F:64,G:64,V:64, FB:64,GB:64,VB:64, H:64,K:64,SIGMA:64, NC:64,IFC:64
  !DIR$ ASSUME (MOD(LDF, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDG, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDV, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDAC, 8) .EQ. 0)

  ! Test the input arguments (assume already done).

  INFO = 0

  ! =====================================================================

  ! Prepare all parameters for inner HZ.

  LDAC5 = 5 * LDAC
  NBSIZE = NBMAXS
  NSWEEP = 0
  NROT = 0

  CALL DLASET('A', M, N, D_ZERO, D_ONE, V, LDV)

  ! =====================================================================

  ! Block HZ - Step 1.
  ! Make a partition of F (and G) into NBL blocks with target size NBMAXS.
  ! Note: NBMAXS is given an input parameter to Block HZ.
  ! It is copyied into NBSIZE, which is then used to make a partition.
  ! It may be modified to contain the actual maximal block size.

  ! For each block V( IBL ), IBL = 1, NBL, we need the following
  ! partitioning information:
  ! IFC( IBL ) = index of the first column of V( IBL )
  !              in global matrix V,
  ! NC( IBL )  = number of columns in V( IBL )
  !              (``size'' of the block).

  CALL PARTBL(IPART, N, NBSIZE, NBL, NC, IFC)

  ! =====================================================================

  ! HZ block iterations (block cycles or sweeps) loop.

  DO ITER = 1, MAXCYC

     ! Sweep initializations.

     TRANSF = .FALSE.

     ! Block HZ - Steps 2, 3.
     ! Block sweeps on blocks.

     ! Block sweep = row by row (top to bottom), row = left to right.

     DO IBL = 1, NBL

        ! Block HZ - step 2 - single diagonal block.

        ! Make a 1 x 1 (``half-sized'') diagonal block in FB (and GB).
        ! FB (and GB) is the Gram matrix of a single block F( IBL ) (or G( IBL )).

        NC1 = NC(IBL)

        ! Check for trivial diagonal block of size 1.

        IF (NC1 .LE. 0) CYCLE

        IFE1 = IFC(IBL)

        ! Nontrivial diagonal block.

        CALL DSYRK('U', 'T', NC1, M, D_ONE, F(1, IFE1), LDF, D_ZERO, FB, LDAC)
        CALL DSYRK('U', 'T', NC1, M, D_ONE, G(1, IFE1), LDG, D_ZERO, GB, LDAC)

        ! Call Cholesky on FB, GB.

        CALL DPOTRF('U', NC1, FB, LDAC, INFO(2))
        IF (INFO(2) .NE. 0) THEN
           INFO(1) = 1
           RETURN
        END IF

        CALL DPOTRF('U', NC1, GB, LDAC, INFO(2))
        IF (INFO(2) .NE. 0) THEN
           INFO(1) = 2
           RETURN
        END IF

        ! The strictly lower triangle of FB, GB has to be annihilated before HZ.

        CALL DLASET('L', NC1 - 1, NC1 - 1, D_ZERO, D_ZERO, FB(2, 1), LDAC)
        CALL DLASET('L', NC1 - 1, NC1 - 1, D_ZERO, D_ZERO, GB(2, 1), LDAC)

        ! Call inner FL on (FB, GB).
        ! VB is the transformation matrix accumulated by HZ.

        IFE2 = NC1 + 1
        CALL MY_DZIMMER0(.TRUE., NC1, NC1, 0, FB, LDAC, GB, LDAC, VB(1, IFE2), LDAC, 1, TOL, H, K, SIGMA, ITEMP, NROTIN, INFO(2))
        IF (INFO(2) .LT. 0) THEN
           INFO(1) = 3
           RETURN
        END IF

        ! If there were no ``inner'' rotations (NROTIN = 0),
        ! then skip the block transformation.

        NROT(1) = NROT(1) + NROTIN(1)
        NROT(2) = NROT(2) + NROTIN(2)

        IF (NROTIN(1) .EQ. 0) CYCLE

        ! Update block convergence tests TRANSF.           
        ! Note: Actually, the only ``normal'' convergence test is that
        ! all block are completed with a single inner sweep, NSWIN = 1.

        ! Slow convergence test: no rotations, or NROTIN = 0, for all blocks
        ! throughout the current sweep.

        TRANSF = .TRUE.

        ! The transformation matrix is VB, returned by HZ.
        ! Update the original ``tall'' columns in block F( IBL ) (and G( IBL )).

        DO ITEMP = 1, M, LDAC5
           Q = MIN(LDAC5, M - ITEMP + 1)

           CALL DGEMM('N', 'N', Q, NC1, NC1, D_ONE, F(ITEMP, IFE1), LDF, VB(1, IFE2), LDAC, D_ZERO, FB, LDAC5)
           CALL DLACPY('A', Q, NC1, FB, LDAC5, F(ITEMP, IFE1), LDF)

           CALL DGEMM('N', 'N', Q, NC1, NC1, D_ONE, G(ITEMP, IFE1), LDG, VB(1, IFE2), LDAC, D_ZERO, FB, LDAC5)
           CALL DLACPY('A', Q, NC1, FB, LDAC5, G(ITEMP, IFE1), LDG)

           CALL DGEMM('N', 'N', Q, NC1, NC1, D_ONE, V(ITEMP, IFE1), LDV, VB(1, IFE2), LDAC, D_ZERO, FB, LDAC5)
           CALL DLACPY('A', Q, NC1, FB, LDAC5, V(ITEMP, IFE1), LDV)
        END DO

        ! End of code for nontrivial diagonal block.

        ! End of code for a single diagonal block.

     END DO

     DO IBL = 1, NBL

        ! Block HZ - step 3 - inner loop, offdiagonal blocks.

        DO JBL = IBL + 1, NBL

           ! Make a 2 x 2 (``full-sized'') block in FB (and GB).
           ! FB (and GB) is the Gram matrix of a pair of blocks
           ! ( F( IBL ), F( JBL ) ) (or G( IBL ), G( JBL ) ).

           NC1 = NC(IBL)
           IFE1 = IFC(IBL)

           NC2 = NC(JBL)
           IFE2 = IFC(JBL)

           NCF = NC1 + NC2

           ! Full block muliply for diagonal blocks FB_11, FB_22.

           CALL DSYRK('U', 'T', NC1, M, D_ONE, F(1, IFE1), LDF, D_ZERO, FB(1, 1), LDAC)
           CALL DSYRK('U', 'T', NC2, M, D_ONE, F(1, IFE2), LDF, D_ZERO, FB(NC1 + 1, NC1 + 1), LDAC)

           ! Full block multiply for diagonal blocks GB_11, GB_22.

           CALL DSYRK('U', 'T', NC1, M, D_ONE, G(1, IFE1), LDG, D_ZERO, GB(1, 1), LDAC)
           CALL DSYRK('U', 'T', NC2, M, D_ONE, G(1, IFE2), LDG, D_ZERO, GB(NC1 + 1, NC1 + 1), LDAC)

           ! Nontrivial off-diagonal block.
           ! Upper offdiagonal block is FB_12 = FB( JBL )^T * FB( IBL ).

           CALL DGEMM('T', 'N', NC1, NC2, M, D_ONE, F(1, IFE1), LDF, F(1, IFE2), LDF, D_ZERO, FB(1, NC1 + 1), LDAC)
           CALL DGEMM('T', 'N', NC1, NC2, M, D_ONE, G(1, IFE1), LDG, G(1, IFE2), LDG, D_ZERO, GB(1, NC1 + 1), LDAC)

           ! Call Cholesky on FB, GB.

           CALL DPOTRF('U', NCF, FB, LDAC, INFO(2))
           IF (INFO(2) .NE. 0) THEN
              INFO(1) = 4
              RETURN
           END IF

           CALL DPOTRF('U', NCF, GB, LDAC, INFO(2))
           IF (INFO(2) .NE. 0) THEN
              INFO(1) = 5
              RETURN
           END IF

           CALL DLASET('L', NCF - 1, NCF - 1, D_ZERO, D_ZERO, FB(2, 1), LDAC)
           CALL DLASET('L', NCF - 1, NCF - 1, D_ZERO, D_ZERO, GB(2, 1), LDAC)

           ! Call inner HZ on ( FB, GB ).
           ! VB is the transformation matrix accumulated by HZ.

           CALL MY_DZIMMER0(.TRUE., NCF, NCF, NC1, FB, LDAC, GB, LDAC, VB, LDAC, 1, TOL, H, K, SIGMA, ITEMP, NROTIN, INFO(2))
           IF (INFO(2) .LT. 0) THEN
              INFO(1) = 6
              RETURN
           END IF

           ! Update statistics for inner and block transformations.
           ! If there were no ``inner'' rotations (NROTIN = 0),
           ! then skip the block transformation.

           NROT(1) = NROT(1) + NROTIN(1)
           NROT(2) = NROT(2) + NROTIN(2)

           IF (NROTIN(1) .EQ. 0) CYCLE

           ! Update block convergence tests TRANSF.
           ! Note: Actually, the only ``normal'' convergence test is that
           ! all block are completed with a single inner sweep, NSWIN = 1.

           ! Slow convergence test: no rotations, or NROTIN = 0, for all blocks
           ! throughout the current sweep.

           TRANSF = .TRUE.

           ! The transformation matrix is VB, returned by HZ.
           ! Update the original ``tall'' columns in blocks
           ! F( IBL ), F( JBL ) (and G( IBL ), G( JBL )).

           DO ITEMP = 1, M, LDAC
              Q = MIN(LDAC, M - ITEMP + 1)

              CALL DLACPY('A', Q, NC1, F(ITEMP, IFE1), LDF, FB(1, 1), LDAC)
              CALL DLACPY('A', Q, NC2, F(ITEMP, IFE2), LDF, FB(1, NC1 + 1), LDAC)
              CALL DGEMM('N', 'N', Q, NCF, NCF, D_ONE, FB, LDAC, VB, LDAC, D_ZERO, GB, LDAC)
              CALL DLACPY('A', Q, NC1, GB(1, 1), LDAC, F(ITEMP, IFE1), LDF)
              CALL DLACPY('A', Q, NC2, GB(1, NC1 + 1), LDAC, F(ITEMP, IFE2), LDF)

              CALL DLACPY('A', Q, NC1, G(ITEMP, IFE1), LDG, FB(1, 1), LDAC)
              CALL DLACPY('A', Q, NC2, G(ITEMP, IFE2), LDG, FB(1, NC1 + 1), LDAC)
              CALL DGEMM('N', 'N', Q, NCF, NCF, D_ONE, FB, LDAC, VB, LDAC, D_ZERO, GB, LDAC)
              CALL DLACPY('A', Q, NC1, GB(1, 1), LDAC, G(ITEMP, IFE1), LDG)
              CALL DLACPY('A', Q, NC2, GB(1, NC1 + 1), LDAC, G(ITEMP, IFE2), LDG)

              CALL DLACPY('A', Q, NC1, V(ITEMP, IFE1), LDV, FB(1, 1), LDAC)
              CALL DLACPY('A', Q, NC2, V(ITEMP, IFE2), LDV, FB(1, NC1 + 1), LDAC)
              CALL DGEMM('N', 'N', Q, NCF, NCF, D_ONE, FB, LDAC, VB, LDAC, D_ZERO, GB, LDAC)
              CALL DLACPY('A', Q, NC1, GB(1, 1), LDAC, V(ITEMP, IFE1), LDV)
              CALL DLACPY('A', Q, NC2, GB(1, NC1 + 1), LDAC, V(ITEMP, IFE2), LDV)
           END DO

           ! End of code for offdiagonal block.

        END DO
     END DO

     ! Convergence Tests.

     ! Slow convergence test - no block transformations per sweep.

     IF (.NOT. TRANSF) THEN
        NSWEEP = ITER
        EXIT
     END IF

     ! End of block sweeps.

  END DO

  ! No convergence.

  IF (NSWEEP .EQ. 0) THEN
     ! No convergence
     NSWEEP = MAXCYC
     INFO(2) = 1
  ELSE
     INFO(2) = 0
  END IF

  ! Postprocessing after block sweeps.

  ! =====================================================================

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

END SUBROUTINE MY_DZIMMER1NS

SUBROUTINE DZIMMER1NS(M, N, F, LDF, G, LDG, V, LDV, MAXCYC, TOL, NBMAXS, IPART,&
     H, K, SIGMA, WORK, LWORK, IWORK, LIWORK, NBSIZE, NSWEEP, NROT, INFO)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: M, N, LDF, LDG, LDV, MAXCYC, NBMAXS, IPART, LWORK, LIWORK
  DOUBLE PRECISION, INTENT(IN) :: TOL
  DOUBLE PRECISION, INTENT(INOUT) :: F(LDF,N), G(LDG,N)
  DOUBLE PRECISION, INTENT(OUT) :: V(LDV,N), H(*), K(*), SIGMA(*), WORK(*)
  INTEGER, INTENT(OUT) :: IWORK(*), NBSIZE, NSWEEP, NROT(2), INFO(2)

  DOUBLE PRECISION :: MYTOL
  INTEGER :: MYIPRT
  LOGICAL :: FAST

  INTEGER :: I, MCYCLE, MY_LWORK, MY_LIWORK, OPT_LWORK, OPT_LIWORK
  INTEGER :: LDAC, LDAC2, MAXNBL, A_MAXNBL
  INTEGER :: I_NC, I_IFC, I_FB, I_GB, I_VB
  LOGICAL :: Q_LWORK, Q_LIWORK

  !DIR$ ASSUME_ALIGNED F:64,G:64,V:64, H:64,K:64,SIGMA:64, WORK:64,IWORK:64
  !DIR$ ASSUME (MOD(LDF, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDG, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDV, 8) .EQ. 0)

  INFO(2) = 0

  IF (IPART .LT. 0) THEN
     FAST = .TRUE.
     MYIPRT = -IPART
     Q_LWORK = .FALSE.
     Q_LIWORK = .FALSE.
  ELSE
     FAST = .FALSE.
     IF (IPART .EQ. 0) THEN
        MYIPRT = 2
     ELSE
        MYIPRT = IPART
     END IF
     Q_LWORK = (LWORK .EQ. -1)
     Q_LIWORK = (LIWORK .EQ. -1)
  END IF

  IF (Q_LWORK) THEN
     MY_LWORK = 0
  ELSE
     MY_LWORK = LWORK
  END IF

  IF (Q_LIWORK) THEN
     MY_LIWORK = 0
  ELSE
     MY_LIWORK = LIWORK
  END IF

  IF (FAST) THEN
     MCYCLE = MAXCYC
     MYTOL = TOL
     INFO(1) = 0
  ELSE
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
        INFO(1) = -1
     ELSE IF (N .LT. 0) THEN
        INFO(1) = -2
     ELSE IF (N .GT. M) THEN
        INFO(1) = -2
     ELSE IF (LDF .LT. M) THEN
        INFO(1) = -4
     ELSE IF (LDG .LT. M) THEN
        INFO(1) = -6
     ELSE IF (LDV .LT. M) THEN
        INFO(1) = -8
     ELSE IF (MCYCLE .LT. 0) THEN
        INFO(1) = -9
     ELSE IF (MYTOL .LT. D_ZERO) THEN
        INFO(1) = -10
     ELSE IF (NBMAXS .LE. 0) THEN
        INFO(1) = -11
     ELSE IF (NBMAXS .GE. N) THEN
        INFO(1) = -11
     ELSE IF (MYIPRT .LT. 1) THEN
        INFO(1) = -12
     ELSE IF (MYIPRT .GT. 2) THEN
        INFO(1) = -12
     ELSE IF (MY_LWORK .LT. 0) THEN
        INFO(1) = -14
     ELSE IF (MY_LIWORK .LT. 0) THEN
        INFO(1) = -16
     ELSE
        INFO(1) = 0
     END IF

     IF (INFO(1) .NE. 0) RETURN
     ! Quick return, if possible.
     IF (M .EQ. 0) RETURN
     IF (N .EQ. 0) RETURN
  END IF

  LDAC = 2 * NBMAXS
  I = MOD(LDAC, 8)
  IF (I .GT. 0) LDAC = LDAC + (8 - I)
  LDAC2 = LDAC * 2 * NBMAXS

  MAXNBL = (N + NBMAXS - 1) / NBMAXS
  I = MOD(MAXNBL, 16)
  IF (I .EQ. 0) THEN
     A_MAXNBL = MAXNBL
  ELSE
     A_MAXNBL = MAXNBL + (16 - I)
  END IF

  ! Calculate optimal LWORK.
  OPT_LWORK = 3 * LDAC2
  ! Calculate optimal LIWORK.
  OPT_LIWORK = 2 * A_MAXNBL

  IF (.NOT. FAST) THEN
     IF (Q_LWORK) WORK(1) = DBLE(OPT_LWORK)
     IF (Q_LIWORK) IWORK(1) = OPT_LIWORK
     IF (Q_LWORK .OR. Q_LIWORK) RETURN

     IF (MY_LWORK .EQ. 0) THEN
        MY_LWORK = OPT_LWORK
     ELSE IF (MY_LWORK .LT. OPT_LWORK) THEN
        INFO(1) = -14
        INFO(2) = OPT_LWORK
     END IF
     IF (MY_LIWORK .EQ. 0) THEN
        MY_LIWORK = OPT_LIWORK
     ELSE IF (MY_LIWORK .LT. OPT_LIWORK) THEN
        INFO(1) = -16
        INFO(2) = OPT_LIWORK
     END IF
     IF (INFO(1) .NE. 0) RETURN
  END IF

  I_NC = 1
  I_IFC = A_MAXNBL + I_NC

  I_FB = 1
  I_GB = LDAC2 + I_FB
  I_VB = LDAC2 + I_GB

  CALL MY_DZIMMER1NS(FAST, M, N, F, LDF, G, LDG, V, LDV, MCYCLE, TOL, LDAC, NBMAXS, MYIPRT,&
       H, K, SIGMA, IWORK(I_NC), IWORK(I_IFC), WORK(I_FB), WORK(I_GB), WORK(I_VB), NBSIZE, NSWEEP, NROT, INFO)

END SUBROUTINE DZIMMER1NS
