SUBROUTINE SINITSH(NBSIZE, NBL, NC, IPL, IFCS, M, N, F, LDF, G, LDG, V, LDV)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NBSIZE, NBL, M, N, LDF, LDG, LDV
  INTEGER, INTENT(IN) :: NC(*)
  INTEGER, INTENT(OUT) :: IPL(*), IFCS(*)
  REAL, INTENT(INOUT) :: F(LDF,*), G(LDG,*), V(LDV,*)

  ! Purpose
  ! =======
  !
  ! Initialize shuffler for NBL blocks with 2 additional places.
  ! Needs:
  !    NBSIZE = target maximal block size for partition,
  !    NBL = the number of blocks NBL,
  ! and for each block IBL (IBL = 1, NBL),
  !    NC( IBL ) = the number of columns.
  ! Returns:
  !    IPL( IBL ) = ``block-index'' in which current block the original
  !                 block V( IBL ) actually resides,
  !    IFCS( IBL ) = index of the first ``spread'' column.

  INTEGER :: IBL, I, IC1, IC2

  EXTERNAL :: SCOPY

  !DIR$ ASSUME_ALIGNED F:64,G:64,V:64, NC:64,IPL:64,IFCS:64
  !DIR$ ASSUME (MOD(LDF, 16) .EQ. 0)
  !DIR$ ASSUME (MOD(LDG, 16) .EQ. 0)
  !DIR$ ASSUME (MOD(LDV, 16) .EQ. 0)

  DO IBL = 1, NBL + 2
     IPL(IBL) = IBL
     IFCS(IBL) = 1 + (IBL - 1) * NBSIZE
  END DO

  ! For non-greedy partitioning, rearrange (``spread'') the columns
  ! so that all blocks have the same artificial size NBSIZE, as the
  ! largest one. After that, all ``blocks'' can accomodate the
  ! largest block, when shuffling.

  IC1 = N - NC(NBL)
  IC2 = (NBL - 1) * NBSIZE
  IBL = NBL

  DO WHILE (IC1 .LT. IC2)

     ! Spread the columns (backward loop).
     ! Note: IC1 = IFC( IBL ) - 1, IC2 = IFCS( IBL ) - 1
     ! and IC2 - IC1 = NSHIFT for copying.

     DO I = NC(IBL), 1, -1
        CALL SCOPY(M, F(1, IC1 + I), 1, F(1, IC2 + I), 1)
        CALL SCOPY(M, G(1, IC1 + I), 1, G(1, IC2 + I), 1)
        CALL SCOPY(M, V(1, IC1 + I), 1, V(1, IC2 + I), 1)
     END DO

     IBL = IBL - 1
     IF (IBL .LE. 0) EXIT
     IC1 = IC1 - NC(IBL)
     IC2 = IC2 - NBSIZE

  END DO

END SUBROUTINE SINITSH

SUBROUTINE DINITSH(NBSIZE, NBL, NC, IPL, IFCS, M, N, F, LDF, G, LDG, V, LDV)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NBSIZE, NBL, M, N, LDF, LDG, LDV
  INTEGER, INTENT(IN) :: NC(*)
  INTEGER, INTENT(OUT) :: IPL(*), IFCS(*)
  DOUBLE PRECISION, INTENT(INOUT) :: F(LDF,*), G(LDG,*), V(LDV,*)

  ! Purpose
  ! =======
  !
  ! Initialize shuffler for NBL blocks with 2 additional places.
  ! Needs:
  !    NBSIZE = target maximal block size for partition,
  !    NBL = the number of blocks NBL,
  ! and for each block IBL (IBL = 1, NBL),
  !    NC( IBL ) = the number of columns.
  ! Returns:
  !    IPL( IBL ) = ``block-index'' in which current block the original
  !                 block V( IBL ) actually resides,
  !    IFCS( IBL ) = index of the first ``spread'' column.

  INTEGER :: IBL, I, IC1, IC2

  EXTERNAL :: DCOPY

  !DIR$ ASSUME_ALIGNED F:64,G:64,V:64, NC:64,IPL:64,IFCS:64
  !DIR$ ASSUME (MOD(LDF, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDG, 8) .EQ. 0)
  !DIR$ ASSUME (MOD(LDV, 8) .EQ. 0)

  DO IBL = 1, NBL + 2
     IPL(IBL) = IBL
     IFCS(IBL) = 1 + (IBL - 1) * NBSIZE
  END DO

  ! For non-greedy partitioning, rearrange (``spread'') the columns
  ! so that all blocks have the same artificial size NBSIZE, as the
  ! largest one. After that, all ``blocks'' can accomodate the
  ! largest block, when shuffling.

  IC1 = N - NC(NBL)
  IC2 = (NBL - 1) * NBSIZE
  IBL = NBL

  DO WHILE (IC1 .LT. IC2)

     ! Spread the columns (backward loop).
     ! Note: IC1 = IFC( IBL ) - 1, IC2 = IFCS( IBL ) - 1
     ! and IC2 - IC1 = NSHIFT for copying.

     DO I = NC(IBL), 1, -1
        CALL DCOPY(M, F(1, IC1 + I), 1, F(1, IC2 + I), 1)
        CALL DCOPY(M, G(1, IC1 + I), 1, G(1, IC2 + I), 1)
        CALL DCOPY(M, V(1, IC1 + I), 1, V(1, IC2 + I), 1)
     END DO

     IBL = IBL - 1
     IF (IBL .LE. 0) EXIT
     IC1 = IC1 - NC(IBL)
     IC2 = IC2 - NBSIZE

  END DO

END SUBROUTINE DINITSH
