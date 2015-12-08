SUBROUTINE MY_DZIMMER2BO(RANK, FAST, M, N, F, LDF, G, LDG, V, LDV, NBBL, NBMAXS, MAXCYC, TOL,&
     NBC, IFCSRC, IFCDST, H, K, SIGMA, LCOLV, FB, GB, VB, LDA, WORK, LWORK, IWORK, LIWORK, NBSWP, NBROT, INFO)

  IMPLICIT NONE

  INTEGER, PARAMETER :: IPART = -2

  LOGICAL, INTENT(IN) :: FAST
  INTEGER, INTENT(IN) :: RANK, M, N, LDF, LDG, LDV, NBBL, NBMAXS, MAXCYC(2), NBC(*), LCOLV, LDA, LWORK, LIWORK
  INTEGER, INTENT(INOUT) :: IFCSRC(*), IFCDST(*)
  DOUBLE PRECISION, INTENT(IN) :: TOL
  DOUBLE PRECISION, INTENT(INOUT) :: F(LDF,*), G(LDG,*), V(LDV,*)
  DOUBLE PRECISION, INTENT(OUT) :: FB(LDA,*), GB(LDA,*), VB(LDA,*), H(*), K(*), SIGMA(*), WORK(*)
  INTEGER, INTENT(OUT) :: IWORK(*), NBSWP, INFO(2)
  INTEGER(8), INTENT(OUT) :: NBROT(2)

  INTEGER :: ITER, ISTP, ITH(2)
  INTEGER :: NBCF, NBSIZE, NBSW
  INTEGER :: INC, JNC
  INTEGER :: IP, JP, IBLK, JBLK
  INTEGER :: ISRC, JSRC, IDST, JDST
  INTEGER(8) :: NROT(2), MYROT(2)
  DOUBLE PRECISION :: SCL
#ifndef USE_DIV
  DOUBLE PRECISION :: D
#endif

  DOUBLE PRECISION, EXTERNAL :: DNRM2
  EXTERNAL :: DCOPY, DGEMM, DLASET, DPOTRF, DSYRK
  
  !DIR$ ASSUME_ALIGNED F:64,G:64,V:64, NBC:64,IFCSRC:64,IFCDST:64, FB:64,GB:64,VB:64, WORK:64,IWORK:64, H:64,K:64,SIGMA:64
  !DIR$ ASSUME (MOD(LDF, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDG, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDV, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDA, 8) .EQ. 0)
  
  CALL STRINI(RANK, NBBL, IP, JP, IBLK, JBLK)

  INC = NBC(IBLK)
  JNC = NBC(JBLK)
  NBCF = INC + JNC
  ISRC = IFCSRC(IBLK)
  JSRC = IFCSRC(JBLK)
  IDST = IFCDST(IBLK)
  JDST = IFCDST(JBLK)

  ! V = I
  ITH(1) = ISRC - 1
  ITH(2) = LDV - ITH(1)
  CALL DLASET('A', ITH(1), INC, D_ZERO, D_ZERO, V(1, ISRC),    LDV)
  CALL DLASET('A', ITH(2), INC, D_ZERO, D_ONE,  V(ISRC, ISRC), LDV)
  ITH(1) = JSRC - 1
  ITH(2) = LDV - ITH(1)
  CALL DLASET('A', ITH(1), JNC, D_ZERO, D_ZERO, V(1, JSRC),    LDV)
  CALL DLASET('A', ITH(2), JNC, D_ZERO, D_ONE,  V(JSRC, JSRC), LDV)

  DO ITER = 1, MAXCYC(1)
     MYROT = 0_8

     DO ISTP = 1, NBBL
        !$OMP BARRIER

        CALL DSYRK('U', 'T', INC, M, D_ONE, F(1,ISRC), LDF, D_ZERO, FB(1,1),         LDA)
        CALL DSYRK('U', 'T', JNC, M, D_ONE, F(1,JSRC), LDF, D_ZERO, FB(INC+1,INC+1), LDA)

        CALL DSYRK('U', 'T', INC, M, D_ONE, G(1,ISRC), LDG, D_ZERO, GB(1,1),         LDA)
        CALL DSYRK('U', 'T', JNC, M, D_ONE, G(1,JSRC), LDG, D_ZERO, GB(INC+1,INC+1), LDA)

        CALL DGEMM('T', 'N', INC, JNC, M, D_ONE, F(1,ISRC), LDF, F(1,JSRC), LDF, D_ZERO, FB(1,INC+1), LDA)
        CALL DGEMM('T', 'N', INC, JNC, M, D_ONE, G(1,ISRC), LDG, G(1,JSRC), LDG, D_ZERO, GB(1,INC+1), LDA)

        CALL DPOTRF('U', NBCF, FB, LDA, ITH(1))
        IF (ITH(1) .NE. 0) THEN
           !$OMP ATOMIC UPDATE
           INFO(2) = IOR(INFO(2), ITH(1))
           !$OMP END ATOMIC
        ELSE
           CALL DLASET('L', NBCF - 1, NBCF - 1, D_ZERO, D_ZERO, FB(2, 1), LDA)
        END IF
        !$OMP BARRIER
        !$OMP ATOMIC READ
        ITH(2) = INFO(2)
        !$OMP END ATOMIC
        IF (ITH(2) .NE. 0) GOTO 3
        !$OMP BARRIER

        CALL DPOTRF('U', NBCF, GB, LDA, ITH(1))
        IF (ITH(1) .NE. 0) THEN
           !$OMP ATOMIC UPDATE
           INFO(2) = IOR(INFO(2), ITH(1))
           !$OMP END ATOMIC
        ELSE
           CALL DLASET('L', NBCF - 1, NBCF - 1, D_ZERO, D_ZERO, GB(2, 1), LDA)
        END IF
        !$OMP BARRIER
        !$OMP ATOMIC READ
        ITH(2) = INFO(2)
        !$OMP END ATOMIC
        IF (ITH(2) .NE. 0) GOTO 4
        !$OMP BARRIER

        CALL DZIMMER1BO(NBCF, NBCF, FB, LDA, GB, LDA, VB, LDA, MAXCYC(2), TOL, NBMAXS, IPART, H, K, SIGMA,&
             WORK, LWORK, IWORK, LIWORK, NBSIZE, NBSW, NROT, ITH)
        IF (ITH(1) .NE. 0) THEN
           !$OMP ATOMIC UPDATE
           INFO(2) = IOR(INFO(2), ITH(1))
           !$OMP END ATOMIC
        ELSE IF (NROT(1) .GT. 0_8) THEN
           MYROT(1) = MYROT(1) + NROT(1)
           MYROT(2) = MYROT(2) + NROT(2)

           CALL DGEMM('N', 'N', M, INC, JNC, D_ONE, F(1,JSRC), LDF, VB(INC+1,1),     LDA, D_ZERO, F(1,IDST), LDF)
           CALL DGEMM('N', 'N', M, INC, INC, D_ONE, F(1,ISRC), LDF, VB(1,1),         LDA, D_ONE,  F(1,IDST), LDF)
           CALL DGEMM('N', 'N', M, JNC, INC, D_ONE, F(1,ISRC), LDF, VB(1,INC+1),     LDA, D_ZERO, F(1,JDST), LDF)
           CALL DGEMM('N', 'N', M, JNC, JNC, D_ONE, F(1,JSRC), LDF, VB(INC+1,INC+1), LDA, D_ONE,  F(1,JDST), LDF)

           CALL DGEMM('N', 'N', M, INC, JNC, D_ONE, G(1,JSRC), LDG, VB(INC+1,1),     LDA, D_ZERO, G(1,IDST), LDG)
           CALL DGEMM('N', 'N', M, INC, INC, D_ONE, G(1,ISRC), LDG, VB(1,1),         LDA, D_ONE,  G(1,IDST), LDG)
           CALL DGEMM('N', 'N', M, JNC, INC, D_ONE, G(1,ISRC), LDG, VB(1,INC+1),     LDA, D_ZERO, G(1,JDST), LDG)
           CALL DGEMM('N', 'N', M, JNC, JNC, D_ONE, G(1,JSRC), LDG, VB(INC+1,INC+1), LDA, D_ONE,  G(1,JDST), LDG)

           CALL DGEMM('N', 'N', M, INC, JNC, D_ONE, V(1,JSRC), LDV, VB(INC+1,1),     LDA, D_ZERO, V(1,IDST), LDV)
           CALL DGEMM('N', 'N', M, INC, INC, D_ONE, V(1,ISRC), LDV, VB(1,1),         LDA, D_ONE,  V(1,IDST), LDV)
           CALL DGEMM('N', 'N', M, JNC, INC, D_ONE, V(1,ISRC), LDV, VB(1,INC+1),     LDA, D_ZERO, V(1,JDST), LDV)
           CALL DGEMM('N', 'N', M, JNC, JNC, D_ONE, V(1,JSRC), LDV, VB(INC+1,INC+1), LDA, D_ONE,  V(1,JDST), LDV)

           IFCSRC(IBLK) = IDST
           IFCSRC(JBLK) = JDST
           IFCDST(IBLK) = ISRC
           IFCDST(JBLK) = JSRC
        END IF
        !$OMP BARRIER
        !$OMP ATOMIC READ
        ITH(2) = INFO(2)
        !$OMP END ATOMIC
        IF (ITH(2) .NE. 0) GOTO 5
        !$OMP BARRIER

        CALL MMSTEP(NBBL, IP, JP, IBLK, JBLK)

        INC = NBC(IBLK)
        JNC = NBC(JBLK)
        NBCF = INC + JNC
        ISRC = IFCSRC(IBLK)
        JSRC = IFCSRC(JBLK)
        IDST = IFCDST(IBLK)
        JDST = IFCDST(JBLK)
     END DO

     !$OMP BARRIER
     !$OMP ATOMIC READ
     NROT(1) = NBROT(1)
     !$OMP END ATOMIC
     !$OMP ATOMIC READ
     NROT(2) = NBROT(2)
     !$OMP END ATOMIC
     !$OMP BARRIER
     !$OMP ATOMIC UPDATE
     NBROT(1) = NBROT(1) + MYROT(1)
     !$OMP END ATOMIC
     !$OMP ATOMIC UPDATE
     NBROT(2) = NBROT(2) + MYROT(2)
     !$OMP END ATOMIC
     !$OMP BARRIER
     !$OMP ATOMIC READ
     MYROT(1) = NBROT(1)
     !$OMP END ATOMIC
     !$OMP ATOMIC READ
     MYROT(2) = NBROT(2)
     !$OMP END ATOMIC
     !$OMP BARRIER

     MYROT(1) = MYROT(1) - NROT(1)
     MYROT(2) = MYROT(2) - NROT(2)

     !$OMP SINGLE
     NBSWP = ITER
     IF (.NOT. FAST) WRITE (*,'(A,I3,A,I13,A,I13)') 'ITER=', ITER, ', MYROT(1)=', MYROT(1), ', MYROT(2)=', MYROT(2)
     !$OMP END SINGLE

     IF (FAST) THEN
        IF (MYROT(1) .EQ. 0_8) GOTO 8
     ELSE
        IF (MYROT(2) .EQ. 0_8) GOTO 8
     END IF
  END DO

  !$OMP SINGLE
  INFO(2) = 1
  !$OMP END SINGLE
  GOTO 9

3 CONTINUE
  !$OMP SINGLE
  INFO(1) = 3
  !$OMP END SINGLE
  GOTO 9

4 CONTINUE
  !$OMP SINGLE
  INFO(1) = 4
  !$OMP END SINGLE
  GOTO 9

5 CONTINUE
  !$OMP SINGLE
  INFO(1) = 5
  !$OMP END SINGLE
  GOTO 9

8 DO IP = 1, INC
     ITH(2) = ISRC + IP - 1
     IF (ITH(2) .GT. N) THEN
        ITH(1) = ITH(2) - N
     ELSE
        ITH(1) = ITH(2)
     END IF

     IF (FAST) THEN
        ! Normalize V.
        SCL = HYPOT(DNRM2(M, F(1, ITH(2)), 1), DNRM2(M, G(1, ITH(2)), 1))
#ifdef USE_DIV
        IF (SCL .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_DIV_SCPY(M, V(1, ITH(2)), V(1, ITH(1)), SCL)
           ELSE
              CALL DARR_DIV_SCAL(M, V(1, ITH(1)), SCL)
           END IF
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, V(1, ITH(2)), 1, V(1, ITH(1)), 1)
        END IF
#else
        D = D_ONE / SCL
        IF (D .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_MUL_SCPY(M, V(1, ITH(2)), V(1, ITH(1)), D)
           ELSE
              CALL DARR_MUL_SCAL(M, V(1, ITH(1)), D)
           END IF
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, V(1, ITH(2)), 1, V(1, ITH(1)), 1)
        END IF
#endif
     ELSE
        H(ITH(1)) = DNRM2(M, F(1, ITH(2)), 1)
        K(ITH(1)) = DNRM2(M, G(1, ITH(2)), 1)

        SCL = K(ITH(1))
#ifdef USE_DIV
        IF (SCL .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_DIV_SCPY(M, G(1, ITH(2)), G(1, ITH(1)), SCL)
           ELSE
              CALL DARR_DIV_SCAL(M, G(1, ITH(1)), SCL)
           END IF
           SIGMA(ITH(1)) = H(ITH(1)) / SCL
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, G(1, ITH(2)), 1, G(1, ITH(1)), 1)
           SIGMA(ITH(1)) = H(ITH(1))
        END IF
#else
        D = D_ONE / SCL
        IF (D .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_MUL_SCPY(M, G(1, ITH(2)), G(1, ITH(1)), D)
           ELSE
              CALL DARR_MUL_SCAL(M, G(1, ITH(1)), D)
           END IF
           SIGMA(ITH(1)) = H(ITH(1)) * D
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, G(1, ITH(2)), 1, G(1, ITH(1)), 1)
           SIGMA(ITH(1)) = H(ITH(1))
        END IF
#endif

        SCL = H(ITH(1))
#ifdef USE_DIV
        IF (SCL .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_DIV_SCPY(M, F(1, ITH(2)), F(1, ITH(1)), SCL)
           ELSE
              CALL DARR_DIV_SCAL(M, F(1, ITH(1)), SCL)
           END IF
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, F(1, ITH(2)), 1, F(1, ITH(1)), 1)
        END IF
#else
        D = D_ONE / SCL
        IF (D .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_MUL_SCPY(M, F(1, ITH(2)), F(1, ITH(1)), D)
           ELSE
              CALL DARR_MUL_SCAL(M, F(1, ITH(1)), D)
           END IF
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, F(1, ITH(2)), 1, F(1, ITH(1)), 1)
        END IF
#endif

        SCL = HYPOT(H(ITH(1)), K(ITH(1)))
#ifdef USE_DIV
        IF (SCL .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_DIV_SCPY(M, V(1, ITH(2)), V(1, ITH(1)), SCL)
           ELSE
              CALL DARR_DIV_SCAL(M, V(1, ITH(1)), SCL)
           END IF
           H(ITH(1)) = H(ITH(1)) / SCL
           K(ITH(1)) = K(ITH(1)) / SCL
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, V(1, ITH(2)), 1, V(1, ITH(1)), 1)
        END IF
#else
        D = D_ONE / SCL
        IF (D .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_MUL_SCPY(M, V(1, ITH(2)), V(1, ITH(1)), D)
           ELSE
              CALL DARR_MUL_SCAL(M, V(1, ITH(1)), D)
           END IF
           H(ITH(1)) = H(ITH(1)) * D
           K(ITH(1)) = K(ITH(1)) * D
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, V(1, ITH(2)), 1, V(1, ITH(1)), 1)
        END IF
#endif
     END IF
  END DO

  DO JP = 1, JNC
     ITH(2) = JSRC + JP - 1
     IF (ITH(2) .GT. N) THEN
        ITH(1) = ITH(2) - N
     ELSE
        ITH(1) = ITH(2)
     END IF

     IF (FAST) THEN
        ! Normalize V.
        SCL = HYPOT(DNRM2(M, F(1, ITH(2)), 1), DNRM2(M, G(1, ITH(2)), 1))
#ifdef USE_DIV
        IF (SCL .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_DIV_SCPY(M, V(1, ITH(2)), V(1, ITH(1)), SCL)
           ELSE
              CALL DARR_DIV_SCAL(M, V(1, ITH(1)), SCL)
           END IF
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, V(1, ITH(2)), 1, V(1, ITH(1)), 1)
        END IF
#else
        D = D_ONE / SCL
        IF (D .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_MUL_SCPY(M, V(1, ITH(2)), V(1, ITH(1)), D)
           ELSE
              CALL DARR_MUL_SCAL(M, V(1, ITH(1)), D)
           END IF
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, V(1, ITH(2)), 1, V(1, ITH(1)), 1)
        END IF
#endif
     ELSE
        H(ITH(1)) = DNRM2(M, F(1, ITH(2)), 1)
        K(ITH(1)) = DNRM2(M, G(1, ITH(2)), 1)

        SCL = K(ITH(1))
#ifdef USE_DIV
        IF (SCL .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_DIV_SCPY(M, G(1, ITH(2)), G(1, ITH(1)), SCL)
           ELSE
              CALL DARR_DIV_SCAL(M, G(1, ITH(1)), SCL)
           END IF
           SIGMA(ITH(1)) = H(ITH(1)) / SCL
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, G(1, ITH(2)), 1, G(1, ITH(1)), 1)
           SIGMA(ITH(1)) = H(ITH(1))
        END IF
#else
        D = D_ONE / SCL
        IF (D .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_MUL_SCPY(M, G(1, ITH(2)), G(1, ITH(1)), D)
           ELSE
              CALL DARR_MUL_SCAL(M, G(1, ITH(1)), D)
           END IF
           SIGMA(ITH(1)) = H(ITH(1)) * D
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, G(1, ITH(2)), 1, G(1, ITH(1)), 1)
           SIGMA(ITH(1)) = H(ITH(1))
        END IF
#endif

        SCL = H(ITH(1))
#ifdef USE_DIV
        IF (SCL .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_DIV_SCPY(M, F(1, ITH(2)), F(1, ITH(1)), SCL)
           ELSE
              CALL DARR_DIV_SCAL(M, F(1, ITH(1)), SCL)
           END IF
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, F(1, ITH(2)), 1, F(1, ITH(1)), 1)
        END IF
#else
        D = D_ONE / SCL
        IF (D .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_MUL_SCPY(M, F(1, ITH(2)), F(1, ITH(1)), D)
           ELSE
              CALL DARR_MUL_SCAL(M, F(1, ITH(1)), D)
           END IF
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, F(1, ITH(2)), 1, F(1, ITH(1)), 1)
        END IF
#endif

        SCL = HYPOT(H(ITH(1)), K(ITH(1)))
#ifdef USE_DIV
        IF (SCL .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_DIV_SCPY(M, V(1, ITH(2)), V(1, ITH(1)), SCL)
           ELSE
              CALL DARR_DIV_SCAL(M, V(1, ITH(1)), SCL)
           END IF
           H(ITH(1)) = H(ITH(1)) / SCL
           K(ITH(1)) = K(ITH(1)) / SCL
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, V(1, ITH(2)), 1, V(1, ITH(1)), 1)
        END IF
#else
        D = D_ONE / SCL
        IF (D .NE. D_ONE) THEN
           IF (ITH(1) .NE. ITH(2)) THEN
              CALL DARR_MUL_SCPY(M, V(1, ITH(2)), V(1, ITH(1)), D)
           ELSE
              CALL DARR_MUL_SCAL(M, V(1, ITH(1)), D)
           END IF
           H(ITH(1)) = H(ITH(1)) * D
           K(ITH(1)) = K(ITH(1)) * D
        ELSE
           IF (ITH(1) .NE. ITH(2)) CALL DCOPY(M, V(1, ITH(2)), 1, V(1, ITH(1)), 1)
        END IF
#endif
     END IF
  END DO

9 CONTINUE

END SUBROUTINE MY_DZIMMER2BO

SUBROUTINE DZIMMER2BO(P, BP, M, N, F, LDF, G, LDG, V, LDV, MAXCYC, TOL, NBMAXS, H, K, SIGMA, WORK, LWORK, IWORK, LIWORK,&
     NSWEEP, NROT, INFO)

  USE OMP_LIB
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: P, BP, M, N, LDF, LDG, LDV, MAXCYC(2), NBMAXS, LWORK, LIWORK
  DOUBLE PRECISION, INTENT(IN) :: TOL
  DOUBLE PRECISION, INTENT(INOUT) :: F(LDF,*), G(LDG,*), V(LDV,*)
  DOUBLE PRECISION, INTENT(OUT) :: H(*), K(*), SIGMA(*), WORK(*)
  INTEGER, INTENT(OUT) :: IWORK(*), NSWEEP, INFO(2)
  INTEGER(8), INTENT(OUT) :: NROT(2)

  INTEGER :: MY_P, NBBL, A_NBBL, N_P, LDA, N_2P, NM2P, LCOLV, A_B, OPT_LWORK, OPT_LIWORK, LW, LIW, RANK, I, NT
  INTEGER :: MCYC(2), I_NBC, I_IFCSRC, I_IFCDST, I_INIWRK
  LOGICAL :: FAST
  DOUBLE PRECISION :: MTOL

  !DIR$ ASSUME_ALIGNED F:64,G:64,V:64, MAXCYC:64, H:64,K:64,SIGMA:64, WORK:64,IWORK:64
  !DIR$ ASSUME (MOD(LDF, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDG, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDV, 8) .EQ. 0)

  INFO(2) = 0
  NSWEEP = 0
  NROT = 0_8

  NT = GET_NTHR()
  IF (P .LT. 0) THEN
     FAST = .TRUE.
     MY_P = -P
  ELSE
     FAST = .FALSE.
     IF (P .EQ. 0) THEN
        MY_P = NT
     ELSE
        MY_P = P
     END IF
  END IF

  IF (FAST) THEN
     MCYC = MAXCYC
     MTOL = TOL
     INFO(1) = 0
  ELSE
     IF (MAXCYC(1) .EQ. -1) THEN
        MCYC(1) = HUGE(MAXCYC(1))
     ELSE
        MCYC(1) = MAXCYC(1)
     END IF
     IF (MAXCYC(2) .EQ. -1) THEN
        MCYC(2) = HUGE(MAXCYC(2))
     ELSE
        MCYC(2) = MAXCYC(2)
     END IF

     IF (.NOT. (TOL .EQ. TOL)) THEN ! QNaN
        MTOL = D_MONE
     ELSE IF (TOL .EQ. D_MONE) THEN ! Compute own TOL.
        MTOL = HUGE(TOL)
     ELSE IF (TOL .EQ. D_ZERO) THEN ! May be +0 or -0.
        MTOL = ABS(TOL)
     ELSE
        MTOL = TOL
     END IF

     IF ((MY_P .LT. 2) .OR. (MY_P .GT. NT)) THEN
        INFO(1) = -1
     ELSE IF ((BP .LT. 1) .OR. ((BP * MY_P) .GT. NT)) THEN
        INFO(1) = -2
     ELSE IF (N .LT. 4) THEN
        INFO(1) = -4
     ELSE IF (M .LT. N) THEN
        INFO(1) = -3
     ELSE IF (LDF .LT. M) THEN
        INFO(1) = -6
     ELSE IF (LDG .LT. M) THEN
        INFO(1) = -8
     ELSE IF (LDV .LT. M) THEN
        INFO(1) = -10
     ELSE IF (MCYC(1) .LT. 0) THEN
        INFO(1) = -11
     ELSE IF (MCYC(2) .LT. 0) THEN
        INFO(1) = -11
     ELSE IF (MTOL .LT. D_ZERO) THEN
        INFO(1) = -12
     ELSE IF (NBMAXS .LE. 0) THEN
        INFO(1) = -13
     ELSE
        INFO(1) = 0
     END IF
     IF (INFO(1) .NE. 0) RETURN
  END IF

  NBBL = 2 * MY_P
  I = MOD(NBBL, 16)
  IF (I .EQ. 0) THEN
     A_NBBL = NBBL
  ELSE
     A_NBBL = NBBL + (16 - I)
  END IF

  N_P = 2 * ((N + NBBL - 1) / NBBL)
  I = MOD(N_P, 8)
  IF (I .EQ. 0) THEN
     LDA = N_P
  ELSE
     LDA = N_P + (8 - I)
  END IF

  CALL DZIMMER1BO(LDA, N_P, F, LDF, G, LDG, V, LDV, MCYC(2), MTOL, NBMAXS, 0,&
     H, K, SIGMA, WORK, -1, IWORK, -1, I, NSWEEP, NROT, INFO)
  IF (INFO(1) .EQ. 0) THEN
     LCOLV = INFO(2)

     A_B = LDA * LCOLV
     I = CEILING(WORK(1))
     OPT_LWORK = (3 * A_B + I) * MY_P
     OPT_LIWORK = 3 * A_NBBL + IWORK(1) * MY_P

     IF (LWORK .EQ. -1) THEN
        WORK(1) = DBLE(OPT_LWORK)
     ELSE IF (LWORK .LT. OPT_LWORK) THEN
        INFO(1) = -17
        INFO(2) = OPT_LWORK
     END IF
     IF (LIWORK .EQ. -1) THEN
        IWORK(1) = OPT_LIWORK
     ELSE IF (LIWORK .LT. OPT_LIWORK) THEN
        INFO(1) = -19
        INFO(2) = OPT_LIWORK
     END IF
     IF ((LWORK .EQ. -1) .OR. (LIWORK .EQ. -1)) THEN
        INFO(2) = 2 * N
        RETURN
     ELSE
        INFO(2) = 0
     END IF
  ELSE
     RETURN
  END IF

  LW = I
  LIW = IWORK(1)

  I_NBC = 1
  I_IFCSRC = I_NBC + A_NBBL
  I_IFCDST = I_IFCSRC + A_NBBL
  I_INIWRK = I_IFCDST + A_NBBL

  N_2P = N / NBBL
  NM2P = MOD(N, NBBL)

  DO I = 0, NBBL - 1
     IWORK(I_NBC + I) = N_2P
  END DO
  DO I = 0, NM2P - 1
     IWORK(I_NBC + I) = IWORK(I_NBC + I) + 1
  END DO

  IWORK(I_IFCSRC) = 1
  IWORK(I_IFCDST) = N + 1
  DO I = 1, NBBL - 1
     IWORK(I_IFCSRC + I) = IWORK(I_IFCSRC + I - 1) + IWORK(I_NBC + I - 1)
     IWORK(I_IFCDST + I) = IWORK(I_IFCSRC + I) + N
  END DO

  !$OMP PARALLEL NUM_THREADS(MY_P) DEFAULT(SHARED) PRIVATE(RANK,I,NT)
  RANK = OMP_GET_THREAD_NUM()
  I = (3 * A_B + LW) * RANK + 1
  NT = BLAS_SET_NUM_THREADS(BP)
  CALL MY_DZIMMER2BO(RANK, FAST, M, N, F, LDF, G, LDG, V, LDV, NBBL, NBMAXS, MCYC, MTOL,&
       IWORK(I_NBC), IWORK(I_IFCSRC), IWORK(I_IFCDST), H, K, SIGMA,&
       LCOLV, WORK(I), WORK(A_B + I), WORK(2 * A_B + I), LDA,&
       WORK(3 * A_B + I), LW, IWORK(I_INIWRK + RANK * LIW), LIW,&
       NSWEEP, NROT, INFO)
  I = BLAS_SET_NUM_THREADS(NT)
  !$OMP END PARALLEL

END SUBROUTINE DZIMMER2BO
