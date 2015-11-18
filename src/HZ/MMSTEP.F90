PURE SUBROUTINE MMSTEP(NBL, IP, JP, IBLK, JBLK)

  ! Purpose
  ! =======
  !
  ! MMSTEP applies a single step of MODULAR-PACKED strategy in Jacobi
  ! for a given processor, and updates processor information.
  !
  !
  ! It computes new position (IP, JP) and new indices (IBLK, JBLK)
  ! of block columns for processor RANK. It also returns all necessary
  ! communication info for this step.
  !
  !
  ! This is for MODULAR PACKED strategy with EVEN number of blocks.
  !
  ! Note: the arguments IP, JP, IBLK, JBLK should NOT be altered
  ! between successive calls of this routine!
  !
  !
  ! Arguments
  ! =========
  !
  ! NBL     (input) INTEGER
  ! The total number of block columns in the Jacobi process.
  ! NBL > 0. NBL must be EVEN (not checked).
  !
  ! IP      (input/output) INTEGER
  ! Imaginary row index of the processor position (IP, JP)
  ! on the block matrix map. 0 < IP <= NBL.
  !
  ! JP      (input/output) INTEGER
  ! Imaginary column index of the processor position (IP, JP)
  ! on the block matrix map. 0 < IP <= JP <= NBL.
  !
  ! IBLK    (input/output) INTEGER
  ! Index of the first block column in the processor.
  ! 0 < IBLK <= NBL.
  !
  ! JBLK    (input/output) INTEGER
  ! Index of the second block column in the processor.
  ! 0 < IBLK < JBLK <= NBL.

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NBL
  INTEGER, INTENT(INOUT) :: IP, JP, IBLK, JBLK

  ! =====================================================================
  !
  ! This is for MODULAR PACKED strategy with EVEN number of blocks.

  IF ((IP + JP) .GT. NBL) THEN

     ! At or below the SW-NE antidiagonal.
     ! Move down <-> increment IP.
     ! If on diagonal after move, jump up to smaller diagonal position
     ! shifted by NBL / 2. As a consequence, there are no diagonal
     ! processors in the bottom triangle.

     ! Move down, increment IP.

     IP = IP + 1
     IF (IP .EQ. JP) THEN

        ! Moved to diagonal, already has its new IBLK in its old JBLK.
        ! Immediately jumps up by NBL / 2, and receives new IBLK.
        ! Keeps its old JBLK.

        IP = IP - NBL / 2
        JP = IP
     END IF

     IBLK = IP

  ELSE

     ! Here, IP + JP <= NBL. Strictly above the SW-NE antidiagonal.
     ! Move right <-> increment JP.
     !
     ! Diagonal processor keeps its IBLK, but sends its JBLK,
     ! as if it is still at bottom diagonal position.

     ! Now move right, increment JP.

     JP = JP + 1
     JBLK = JP

  END IF

END SUBROUTINE MMSTEP
