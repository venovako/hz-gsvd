PROGRAM S0

  USE H5LT
  USE HDF5

  USE HZ

  IMPLICIT NONE

  CHARACTER(LEN=255) :: H5FN, H5RN
  CHARACTER(LEN=5) :: H5GN
  LOGICAL :: FEXIST

  INTEGER(HID_T) :: FID, GID
  INTEGER(HSIZE_T) :: DIMS1(1), DIMS2(2)

  INTEGER :: IDADIM(4), LDF, LDG, LDV, M, N
  EQUIVALENCE (IDADIM(1), LDF), (IDADIM(1), LDG), (IDADIM(1), LDV), (IDADIM(1), M), (IDADIM(1), N)

  INTEGER :: BP, OLDBP, P, Q, LWORK, ISTATS(4), NSWEEP, INFO
  INTEGER :: NROT(2), WALLTM(3)
  DOUBLE PRECISION :: DSTATS(2)
  !DIR$ ATTRIBUTES ALIGN: 64:: ISTATS, DSTATS
  REAL :: TOL

  REAL, ALLOCATABLE :: F(:,:), G(:,:), V(:,:), VQ(:,:)
  !DIR$ ATTRIBUTES ALIGN: 64:: F, G, V, VQ
  REAL, ALLOCATABLE :: SIGMA(:), H(:), K(:), TAU(:), WORK(:)
  !DIR$ ATTRIBUTES ALIGN: 64:: SIGMA, H, K, TAU, WORK

  EXTERNAL :: SLACPY, SLASET, SGEQLF, SORGQL

  CALL READCL(H5FN, H5GN, H5RN, BP)
  IF (BP .LT. 0) STOP 'BTHRDS < 0'
  OLDBP = BLAS_SET_NUM_THREADS(BP)

  INQUIRE(FILE=H5FN, EXIST=FEXIST)
  IF (.NOT. FEXIST) STOP 'Input file nonexistent!'

  CALL h5open_f(INFO)
  IF (INFO .NE. 0) STOP 'Error initializing HDF5!'

  CALL h5fopen_f(H5FN, H5F_ACC_RDONLY_F, FID, INFO)
  IF (INFO .NE. 0) STOP 'Error opening the input file!'

  CALL h5gopen_f(FID, H5GN, GID, INFO)
  IF (INFO .NE. 0) STOP 'Error opening the input group!'

  DIMS1(1) = 4
  CALL h5ltread_dataset_int_f(GID, 'IDADIM', IDADIM, DIMS1, INFO)
  IF (INFO .NE. 0) STOP 'Error reading IDADIM!'

  ALLOCATE(F(LDF, N))
  ALLOCATE(G(LDG, N))
  ALLOCATE(V(LDV, N))
  ALLOCATE(VQ(LDV, N))

  ALLOCATE(SIGMA(N))
  ALLOCATE(H(N))
  ALLOCATE(K(N))
  ALLOCATE(TAU(N))

  LWORK = -1
  CALL SGEQLF(M, N, V, LDV, TAU, SIGMA, LWORK, INFO)
  IF (INFO .NE. 0) STOP 'SGEQLF workspace query'
  P = CEILING(SIGMA(1))
  CALL SORGQL(M, N, N, VQ, LDV, TAU, SIGMA, LWORK, INFO)
  IF (INFO .NE. 0) STOP 'SORGQL workspace query'
  Q = CEILING(SIGMA(1))
  LWORK = MAX(P, Q)
  ALLOCATE(WORK(LWORK))

  DIMS2(1) = LDF
  DIMS2(2) = N
  CALL h5ltread_dataset_float_f(GID, 'F', F, DIMS2, INFO)
  IF (INFO .NE. 0) STOP 'Error reading F!'

  DIMS2(1) = LDG
  DIMS2(2) = N
  CALL h5ltread_dataset_float_f(GID, 'G', G, DIMS2, INFO)
  IF (INFO .NE. 0) STOP 'Error reading G!'

  CALL h5gclose_f(GID, INFO)
  IF (INFO .NE. 0) STOP 'Error closing the input group!'

  CALL h5fclose_f(FID, INFO)
  IF (INFO .NE. 0) STOP 'Error closing the input file!'

  IF (LDF .GT. M) THEN
     P = LDV - M
     Q = M + 1
     CALL SLASET('A', P, N, S_ZERO, S_ZERO, F(Q, 1), LDF)
  END IF

  IF (LDG .GT. M) THEN
     P = LDV - M
     Q = M + 1
     CALL SLASET('A', P, N, S_ZERO, S_ZERO, G(Q, 1), LDG)
  END IF

  IF (LDV .GT. M) THEN
     P = LDV - M
     Q = M + 1
     CALL SLASET('A', P, N, S_ZERO, S_ZERO, V(Q, 1), LDV)
     CALL SLASET('A', P, N, S_ZERO, S_ZERO, VQ(Q, 1), LDV)
  END IF

  WRITE (*,*) N
  CALL TIMER_START(WALLTM)

  TOL = S_MONE
  NSWEEP = 0
  NROT = 0
  INFO = 0

  CALL SZIMMER0(M, N, F, LDF, G, LDG, V, LDV, -1, TOL, H, K, SIGMA, NSWEEP, NROT, INFO)

  IF (INFO .EQ. 0) THEN
     CALL SGEQLF(M, N, V, LDV, TAU, WORK, LWORK, INFO)
     IF (INFO .NE. 0) STOP 'SGEQLF'

     CALL SLACPY('A', M, N, V, LDV, VQ, LDV)
     CALL SLASET('U', M, N - 1, S_ZERO, S_ZERO, V(1, 2), LDV)

     CALL SORGQL(M, N, N, VQ, LDV, TAU, WORK, LWORK, INFO)
     IF (INFO .NE. 0) STOP 'SORGQL'
  ELSE
     SIGMA = S_ZERO
  END IF

  CALL TIMER_STOP(WALLTM)

  DSTATS(1) = TIMER2DBLE(WALLTM)
  DSTATS(2) = TOL

  WRITE (*,*) 'WALL_TIME ', DSTATS(1)
  WRITE (*,*) 'NSWEEP ', NSWEEP
  WRITE (*,*) 'NROT ', NROT(1), ' ', NROT(2)
  WRITE (*,*) 'INFO ', INFO

  ISTATS(1) = NROT(1)
  ISTATS(2) = NROT(2)
  ISTATS(3) = NSWEEP
  ISTATS(4) = INFO

  IF (TRIM(H5RN) .NE. '-') THEN
     INQUIRE(FILE=H5RN, EXIST=FEXIST)
     IF (FEXIST) THEN
        CALL h5fopen_f(H5RN, H5F_ACC_RDWR_F, FID, INFO)
        IF (INFO .NE. 0) STOP 'Error opening the output file!'
     ELSE
        CALL h5fcreate_f(H5RN, H5F_ACC_TRUNC_F, FID, INFO)
        IF (INFO .NE. 0) STOP 'Error creating the output file!'
     END IF

     CALL h5gcreate_f(FID, H5GN, GID, INFO)
     IF (INFO .NE. 0) STOP 'Error creating the output group!'

     DIMS1(1) = 4
     CALL h5ltmake_dataset_int_f(GID, 'IDADIM', 1, DIMS1, IDADIM, INFO)
     IF (INFO .NE. 0) STOP 'Error writing IDADIM!'

     DIMS1(1) = 4
     CALL h5ltmake_dataset_int_f(GID, 'ISTATS', 1, DIMS1, ISTATS, INFO)
     IF (INFO .NE. 0) STOP 'Error writing ISTATS!'

     DIMS1(1) = 2
     CALL h5ltmake_dataset_double_f(GID, 'DSTATS', 1, DIMS1, DSTATS, INFO)
     IF (INFO .NE. 0) STOP 'Error writing DSTATS!'

     DIMS1(1) = N
     CALL h5ltmake_dataset_float_f(GID, 'SIGMA', 1, DIMS1, SIGMA, INFO)
     IF (INFO .NE. 0) STOP 'Error writing SIGMA!'
     CALL h5ltmake_dataset_float_f(GID, 'H', 1, DIMS1, H, INFO)
     IF (INFO .NE. 0) STOP 'Error writing H!'
     CALL h5ltmake_dataset_float_f(GID, 'K', 1, DIMS1, K, INFO)
     IF (INFO .NE. 0) STOP 'Error writing K!'

#ifdef USE_MTXOUT
     DIMS2(1) = LDV
     DIMS2(2) = N
     CALL h5ltmake_dataset_float_f(GID, 'F', 2, DIMS2, F, INFO)
     IF (INFO .NE. 0) STOP 'Error writing F!'
     CALL h5ltmake_dataset_float_f(GID, 'G', 2, DIMS2, G, INFO)
     IF (INFO .NE. 0) STOP 'Error writing G!'
     CALL h5ltmake_dataset_float_f(GID, 'V', 2, DIMS2, V, INFO)
     IF (INFO .NE. 0) STOP 'Error writing V!'
     CALL h5ltmake_dataset_float_f(GID, 'VQ', 2, DIMS2, VQ, INFO)
     IF (INFO .NE. 0) STOP 'Error writing VQ!'
#endif

     CALL h5gclose_f(GID, INFO)
     IF (INFO .NE. 0) STOP 'Error closing the output group!'

     CALL h5fclose_f(FID, INFO)
     IF (INFO .NE. 0) STOP 'Error closing the output file!'
  END IF

  DEALLOCATE(WORK)
  DEALLOCATE(TAU)

  DEALLOCATE(K)
  DEALLOCATE(H)
  DEALLOCATE(SIGMA)

  DEALLOCATE(VQ)
  DEALLOCATE(V)
  DEALLOCATE(F)
  DEALLOCATE(G)

  CALL h5close_f(INFO)
  IF (INFO .NE. 0) STOP 'Error closing HDF5!'
  BP = BLAS_SET_NUM_THREADS(OLDBP)

CONTAINS

  SUBROUTINE READCL(H5FN, H5GN, H5RN, BP)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(OUT) :: H5FN, H5GN, H5RN
    INTEGER, INTENT(OUT) :: BP

    CHARACTER(LEN=8) :: CLA

    IF (COMMAND_ARGUMENT_COUNT() .GT. 4) STOP 's0.exe H5FN H5GN H5RN BTHRDS'

    IF (COMMAND_ARGUMENT_COUNT() .GE. 1) THEN
       CALL GET_COMMAND_ARGUMENT(1, H5FN)
    ELSE
       WRITE (*,'(A)',ADVANCE='NO') 'H5FN: '
       READ (*,'(A)') H5FN
    END IF
    IF (COMMAND_ARGUMENT_COUNT() .GE. 2) THEN
       CALL GET_COMMAND_ARGUMENT(2, H5GN)
    ELSE
       WRITE (*,'(A)',ADVANCE='NO') 'H5GN: '
       READ (*,'(A)') H5GN
    END IF
    IF (COMMAND_ARGUMENT_COUNT() .GE. 3) THEN
       CALL GET_COMMAND_ARGUMENT(3, H5RN)
    ELSE
       WRITE (*,'(A)',ADVANCE='NO') 'H5RN: '
       READ (*,'(A)') H5RN
    END IF
    IF (COMMAND_ARGUMENT_COUNT() .GE. 4) THEN
       CALL GET_COMMAND_ARGUMENT(4, CLA)
       READ (CLA,*) BP
    ELSE
       WRITE (*,'(A)',ADVANCE='NO') 'BTHRDS: '
       READ (*,*) BP
    END IF

  END SUBROUTINE READCL

END PROGRAM S0
