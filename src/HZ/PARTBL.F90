SUBROUTINE PARTBL(IPART, N, NBSIZE, NBL, NC, IFC)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: IPART, N
  INTEGER, INTENT(INOUT) :: NBSIZE
  INTEGER, INTENT(OUT) :: NBL
  INTEGER, INTENT(OUT) :: NC(*), IFC(*)

  ! Purpose
  ! =======
  !
  ! Partition N columns into blocks of given maximal size NBSIZE.
  ! IPART selects the type of partitioning:
  !    IPART <= 1  - ``greedy'' partitioning,
  !    IPART  = 2  - ``equal sized'' partitioning.
  ! Currently:
  !    IPART >= 2 signals ``equal sized'' partitioning.
  !
  ! Returns:
  !    NBL = the number of blocks NBL,
  !    NBSIZE = the actual maximal block size,
  ! and for each block IBL (IBL = 1, NBL),
  !    NC( IBL ) = the number of columns,
  !    IFC( IBL ) = index of the first column.

  INTEGER :: IBL, NR

  !DIR$ ASSUME_ALIGNED NC:64, IFC:64

  IF (N .GT. NBSIZE) THEN

     ! More than one block. Set the number of blocks.

     NBL = (N + NBSIZE - 1) / NBSIZE

     IF (IPART .LE. 1) THEN

        ! ``Greedy'' partition.
        ! N = (N / NBSIZE) * NBSIZE + NR, with NR = MOD( N, NBSIZE ) >= 0.

        NR = MOD(N, NBSIZE)

        ! Large blocks (length NBSIZE) first (N/NBSIZE of them).

        DO IBL = 1, N / NBSIZE
           NC(IBL) = NBSIZE
        END DO

        ! Small block (length NR) last (at most one, if NR > 0).

        IF (NR .GT. 0) NC(NBL) = NR

     ELSE

        ! ``Equal sized'' partition.
        ! N = (NBL - NR) * (N / NBL) + NR * (N / NBL + 1),
        ! with NR = MOD( N, NBL ) >= 0.

        NBSIZE = N / NBL
        NR = MOD(N, NBL)

        ! Smaller blocks (length N/NBL) last (NBL-NR of them).

        DO IBL = NR + 1, NBL
           NC(IBL) = NBSIZE
        END DO

        ! Larger blocks (length N/NBL + 1) first (NR of them, may be none).

        IF (NR .GT. 0) THEN
           NBSIZE = NBSIZE + 1
           DO IBL = 1, NR
              NC(IBL) = NBSIZE
           END DO
        END IF

     END IF

     ! Set IFC via NC. New IFC is previous IFC + previous block length.

     IFC(1) = 1
     DO IBL = 2, NBL
        IFC(IBL) = IFC(IBL - 1) + NC(IBL - 1)
     END DO

     ! End of partitioning.

  ELSE

     ! Trivial partition with a single block.

     NBSIZE = N
     NBL = 1
     NC(1) = N
     IFC(1) = 1

  END IF

END SUBROUTINE PARTBL
