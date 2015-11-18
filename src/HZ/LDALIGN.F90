! Computes LDA >= N such that LDA mod (A/S) = 0
! S: size of datatype in bytes
! A: alignment in bytes
PURE INTEGER FUNCTION LDALIGN(N, S, A)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N, S, A

  INTEGER :: AA, K, L, M

  LDALIGN = 0
  IF (N .EQ. 0) RETURN

  LDALIGN = -1
  IF (N .LT. 0) RETURN

  LDALIGN = -2
  IF (S .LE. 0) RETURN

  LDALIGN = -3
  IF (A .LE. 0) RETURN

  AA = LCM(A, S)

  K = AA / S
  L = N
  M = MOD(L, K)
  IF (M .GT. 0) L = L + (K - M)
  LDALIGN = L

END FUNCTION LDALIGN
