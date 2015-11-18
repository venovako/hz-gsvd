PURE SUBROUTINE STRINI(RANK, NBL, IP, JP, IBLK, JBLK)

  ! Purpose
  ! =======
  !
  ! STRINI performs the initialization step of MODULAR strategies in
  ! Jacobi for a given processor.
  !
  !
  ! It computes the initial position (IP, JP) and the initial indices
  ! (IBLK, JBLK) of block columns for processor RANK.
  !
  !
  ! This is for MODULAR and ROUND-ROBIN MODULAR strategies in Jacobi.
  ! The sum of indices is NBL + 1 (main antidiagonal), and the ordering
  ! is NE -> SW.
  !
  ! Note: the arguments IP, JP, IBLK, JBLK should NOT be altered
  ! before or between successive calls of the step routine!
  !
  !
  ! Arguments
  ! =========
  !
  ! RANK    (input) INTEGER
  !         The index of the processor. 0 <= RANK < NPROC (not checked).
  !
  ! NBL     (input) INTEGER
  !         The total number of block columns in the Jacobi process.
  !         NBL > 0. NBL must be EVEN (not checked).
  !
  ! IP      (output) INTEGER
  !         Imaginary row index of the processor position (IP, JP)
  !         on the block matrix map. 0 < IP <= NBL.
  !
  ! JP      (output) INTEGER
  !         Imaginary column index of the processor position (IP, JP)
  !         on the block matrix map. 0 < IP <= JP <= NBL.
  !
  ! IBLK    (output) INTEGER
  !         Index of the first block column in the processor.
  !         0 < IBLK <= NBL.
  !
  ! JBLK    (output) INTEGER
  !         Index of the second block column in the processor.
  !         0 < IBLK < JBLK <= NBL.

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: RANK, NBL
  INTEGER, INTENT(OUT) :: IP, JP, IBLK, JBLK

  ! =====================================================================
  !
  ! This is for MODULAR and ROUND-ROBIN MODULAR strategies in Jacobi.
  !
  ! The sum of indices is IP + JP = NBL + 1 (main antidiagonal).
  ! The ordering is NE -> SW, i.e., from the top-right corner.

  IP = RANK + 1
  JP = NBL - RANK
  IBLK = IP
  JBLK = JP

END SUBROUTINE STRINI
