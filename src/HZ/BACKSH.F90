SUBROUTINE BACKSH(NBSIZE, NBL, NC, IPL, IFCS, INVP, M, N, F, LDF, G, LDG, V, LDV)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NBSIZE, NBL, M, N, LDF, LDG, LDV
  INTEGER, INTENT(IN) :: NC(*), IFCS(*)
  INTEGER, INTENT(INOUT) :: IPL(*)
  INTEGER, INTENT(OUT) :: INVP(*)
  DOUBLE PRECISION, INTENT(INOUT) :: F(LDF,*), G(LDG,*), V(LDV,*)

  ! Purpose
  ! =======
  !
  ! Restores original columns from shuffled places for NBL blocks
  ! with 2 additional places.
  ! Needs:
  !    NBSIZE = target maximal block size for partition,
  !    NBL = the number of blocks NBL,
  ! and for each block IBL (IBL = 1, NBL),
  !    NC( IBL ) = the number of columns.
  !    IPL( IBL ) = ``block-index'' in which current block the original
  !                 block V( IBL ) actually resides,
  !    IFCS( IBL ) = index of the first ``spread'' column.
  ! Working array INVP holds the inverse of the permutation IPL.

  INTEGER :: IBL, I, IC0, IC1, IC2, IP, INP, NCMIN, NCINP

  EXTERNAL :: DCOPY, DSWAP

  !DIR$ ASSUME_ALIGNED F:64,G:64,V:64, NC:64,IFCS:64, IPL:64,INVP:64
  !DIR$ ASSUME (MOD(LDF, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDG, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDV, 8) .EQ. 0)

  ! Invert the permutation IPL and store it in INVP.

  DO IBL = 1, NBL + 2
     INVP(IPL(IBL)) = IBL
  END DO

  ! Back-shuffle and ``compress'' the blocks into their original positions.

  IC0 = 0
  IC1 = 0

  DO IBL = 1, NBL

     ! Bring the block V( IBL ) into its proper place, at columns
     ! with indices IC0 + 1 through IC0 + NC( IBL ).
     ! Note: IC0 = IFC( IBL ) - 1,
     !       IC1 = IFCS( IBL ) - 1 = (IBL - 1) * NBSIZE.
     ! Current state:
     ! The original V( IBL ) is at ``spread'' V( IPL( IBL ) ), and
     ! the original V( INVP( IBL ) ) is currently at ``spread'' V( IBL ).

     ! Check if we have to swap blocks V( IBL ) and V( IPL( IBL ) ).

     IF (IPL(IBL) .NE. IBL) THEN

        ! Swap blocks V( IBL ) and V( IPL( IBL ) ), and ``compress''
        ! the final V( IBL ), if necessary.

        IP = IPL(IBL)

        ! Note: IC2 = IFCS( IP ) - 1 = (IP - 1) * NBSIZE, and IC2 <> IC1.

        IC2 = IFCS(IP) - 1

        ! Keep track of the sizes NC( IBL ), NC( INP = INVP( IBL ) ).
        ! First check if INP points to a ``free block'', i.e., INP > NBL.

        INP = INVP(IBL)
        IF (INP .LE. NBL) THEN
           NCINP = NC(INP)
        ELSE
           NCINP = 0
        END IF

        NCMIN = MIN(NC(IBL), NCINP)

        ! Check if ``compression'' is necessary (IC0 < IC1).

        IF (IC0 .EQ. IC1) THEN

           ! No ``compression''. Just swap, up to the smaller block size.

           DO I = 1, NCMIN
              CALL DSWAP(M, F(1, IC0 + I), 1, F(1, IC2 + I), 1)
              CALL DSWAP(M, G(1, IC0 + I), 1, G(1, IC2 + I), 1)
              CALL DSWAP(M, V(1, IC0 + I), 1, V(1, IC2 + I), 1)
           END DO

        ELSE

           ! We have ``compression'' (IC0 < IC1). Instead of swap, and then
           ! ``compress'', we immediately copy each column into its proper
           ! place. First ``copy in'' (IC2 + I -> IC0 + I), then ``copy out''
           ! (IC1 + I -> IC2 + I), up to the smaller block size.

           DO I = 1, NCMIN
              CALL DCOPY(M, F(1, IC2 + I), 1, F(1, IC0 + I), 1)
              CALL DCOPY(M, G(1, IC2 + I), 1, G(1, IC0 + I), 1)
              CALL DCOPY(M, V(1, IC2 + I), 1, V(1, IC0 + I), 1)
              CALL DCOPY(M, F(1, IC1 + I), 1, F(1, IC2 + I), 1)
              CALL DCOPY(M, G(1, IC1 + I), 1, G(1, IC2 + I), 1)
              CALL DCOPY(M, V(1, IC1 + I), 1, V(1, IC2 + I), 1)
           END DO

        END IF

        ! If one of the blocks is larger in size, then copy the remaining
        ! columns of that block into the proper place to complete the job.
        ! At most one of the following two loops will be executed.

        DO I = NCMIN + 1, NC(IBL)
           CALL DCOPY(M, F(1, IC2 + I), 1, F(1, IC0 + I), 1)
           CALL DCOPY(M, G(1, IC2 + I), 1, G(1, IC0 + I), 1)
           CALL DCOPY(M, V(1, IC2 + I), 1, V(1, IC0 + I), 1)
        END DO
        DO I = NCMIN + 1, NCINP
           CALL DCOPY(M, F(1, IC1 + I), 1, F(1, IC2 + I), 1)
           CALL DCOPY(M, G(1, IC1 + I), 1, G(1, IC2 + I), 1)
           CALL DCOPY(M, V(1, IC1 + I), 1, V(1, IC2 + I), 1)
        END DO

        ! Record the swap in both IPL and INVP.
        ! New values: INVP( IPL( IBL ) ) = INVP( IBL ),
        !             IPL( INVP( IBL ) ) = IPL( IBL ).

        INVP(IP) = INP
        IPL(INP) = IP

        ! Not necessary to set: IPL( IBL ) = IBL, INVP( IBL ) = IBL.

     ELSE

        ! No swap (IBL = IPL( IBL ), or IC1 = IC2).
        ! Just ``compress'' the block V( IBL ), if necessary (IC0 < IC1).

        IF (IC0 .LT. IC1) THEN
           DO I = 1, NC(IBL)
              CALL DCOPY(M, F(1, IC1 + I), 1, F(1, IC0 + I), 1)
              CALL DCOPY(M, G(1, IC1 + I), 1, G(1, IC0 + I), 1)
              CALL DCOPY(M, V(1, IC1 + I), 1, V(1, IC0 + I), 1)
           END DO
        END IF

     END IF

     ! Update starting positions for the next block.

     IC0 = IC0 + NC(IBL)
     IC1 = IC1 + NBSIZE

  END DO

END SUBROUTINE BACKSH
