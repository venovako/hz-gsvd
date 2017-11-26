PROGRAM EXE2

  USE H5LT
  USE HDF5

  USE HZ

  IMPLICIT NONE

  CHARACTER(LEN=255) :: H5FN, H5RN
  CHARACTER(LEN=5) :: H5GN
  LOGICAL :: FEXIST

  INTEGER(HID_T) :: FID, GID
  INTEGER(HSIZE_T) :: DIMS1(1), DIMS2(2)

  INTEGER :: IDADIM(4)
  INTEGER :: M, N, LDA, P, BP, OLDBP, MAXCYC(2), NBMAXS, LCOLV
  EQUIVALENCE (IDADIM(1), LDA), (IDADIM(1), M), (IDADIM(1), N)

  INTEGER :: LWORK, LIWORK, NSWEEP, INFO(2), ISTATS(4)
  INTEGER :: NROT(2), CLK(3)
  DOUBLE PRECISION :: TOL, DSTATS(2)
  !DIR$ ATTRIBUTES ALIGN: 64:: ISTATS, DSTATS

  DOUBLE PRECISION, ALLOCATABLE :: F(:,:), G(:,:), V(:,:), VQ(:,:)
  !DIR$ ATTRIBUTES ALIGN: 64:: F, G, V, VQ
  DOUBLE PRECISION, ALLOCATABLE :: H(:), K(:), SIGMA(:), TAU(:), WORK(:)
  !DIR$ ATTRIBUTES ALIGN: 64:: H, K, SIGMA, TAU, WORK
  INTEGER, ALLOCATABLE :: IWORK(:)
  !DIR$ ATTRIBUTES ALIGN: 64:: IWORK

  EXTERNAL :: DLACPY, DLASET, DGEQLF, DORGQL

  CALL READCL(H5FN, H5GN, H5RN, NBMAXS, P, BP)
  IF (BP .LT. 0) STOP 'BTHRDS < 0'
  OLDBP = BLAS_SET_NUM_THREADS(P)

  INQUIRE(FILE=H5FN, EXIST=FEXIST)
  IF (.NOT. FEXIST) STOP 'Input file nonexistent!'

  CALL h5open_f(INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error initializing HDF5!'

  CALL h5fopen_f(H5FN, H5F_ACC_RDONLY_F, FID, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error opening the input file!'

  CALL h5gopen_f(FID, H5GN, GID, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error opening the input group!'

  DIMS1(1) = 4
  CALL h5ltread_dataset_int_f(GID, 'IDADIM', IDADIM, DIMS1, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error reading IDADIM!'

  LDA = IDADIM(1)
  M = IDADIM(1)
  N = IDADIM(1)

  IF (LDA .LT. M) STOP 'LDA < M'
  IF (M .LT. N) STOP 'M < N'

  LCOLV = 2 * N
  DIMS2(1) = LDA
  DIMS2(2) = N

  ALLOCATE(F(LDA, LCOLV))
  CALL h5ltread_dataset_double_f(GID, 'F', F, DIMS2, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error reading F!'

  ALLOCATE(G(LDA, LCOLV))
  CALL h5ltread_dataset_double_f(GID, 'G', G, DIMS2, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error reading G!'

  CALL h5gclose_f(GID, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error closing the input group!'

  CALL h5fclose_f(FID, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error closing the input file!'

  ALLOCATE(V(LDA, LCOLV))
  ALLOCATE(VQ(LDA, N))

  ALLOCATE(H(N))
  ALLOCATE(K(N))
  ALLOCATE(SIGMA(N))

  MAXCYC(1) = HUGE(MAXCYC(1))
  MAXCYC(2) = 1
  TOL = D_MONE
  LWORK = -1
  LIWORK = -1
  NSWEEP = 0
  NROT = 0

  CALL DZIMMER2BO(P, BP, M, N, F, LDA, G, LDA, V, LDA, MAXCYC, TOL, NBMAXS, H, K, SIGMA, DSTATS, LWORK, ISTATS, LIWORK,&
       NSWEEP, NROT, INFO)
  IF (INFO(1) .NE. 0) STOP 'Error querying the worksize for DZIMMER2BO!'
  LIWORK = ISTATS(1)
  ISTATS(2) = MAX(CEILING(DSTATS(1)), 1)
  IF (LCOLV .NE. INFO(2)) STOP 'LCOLV inconsistent!'

  ALLOCATE(TAU(N))

  DSTATS = D_ZERO
  LWORK = -1
  CALL DGEQLF(M, N, V, LDA, TAU, DSTATS(1), LWORK, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'DGEQLF workspace query'
  CALL DORGQL(M, N, N, VQ, LDA, TAU, DSTATS(2), LWORK, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'DORGQL workspace query'
  LWORK = MAX(CEILING(MAX(DSTATS(1), DSTATS(2))), ISTATS(2))

  ALLOCATE(WORK(LWORK))
  ALLOCATE(IWORK(LIWORK))

  CALL TIMER_START(CLK)

  CALL DZIMMER2BO(P, BP, M, N, F, LDA, G, LDA, V, LDA, MAXCYC, TOL, NBMAXS, H, K, SIGMA, WORK, LWORK, IWORK, LIWORK,&
       NSWEEP, NROT, INFO)
  IF (INFO(1) .EQ. 0) THEN
    CALL DGEQLF(M, N, V, LDA, TAU, WORK, LWORK, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'DGEQLF'

    CALL DLACPY('A', M, N, V, LDA, VQ, LDA)
    CALL DLASET('U', M, N - 1, D_ZERO, D_ZERO, V(1, 2), LDA)

    CALL DORGQL(M, N, N, VQ, LDA, TAU, WORK, LWORK, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'DORGQL'
  ELSE
    SIGMA = D_ZERO
  END IF

  CALL TIMER_STOP(CLK)

  DEALLOCATE(IWORK)
  DEALLOCATE(WORK)
  DEALLOCATE(TAU)

  DSTATS(1) = TIMER2DBLE(CLK)
  DSTATS(2) = TOL

  WRITE (*,*) 'WALL_TIME ', DSTATS(1)
  WRITE (*,*) 'NSWEEP ', NSWEEP, ' BP ', BP
  WRITE (*,*) 'NROT ', NROT(1), ' ', NROT(2)
  WRITE (*,*) 'INFO ', INFO(1), ' ', INFO(2)

  ISTATS(1) = NROT(1)
  ISTATS(2) = NROT(2)
  ISTATS(3) = NSWEEP
  ISTATS(4) = INFO(2)

  IF (TRIM(H5RN) .NE. '-') THEN
    INQUIRE(FILE=H5RN, EXIST=FEXIST)
    IF (FEXIST) THEN
      CALL h5fopen_f(H5RN, H5F_ACC_RDWR_F, FID, INFO(2))
      IF (INFO(2) .NE. 0) STOP 'Error opening the output file!'
    ELSE
      CALL h5fcreate_f(H5RN, H5F_ACC_TRUNC_F, FID, INFO(2))
      IF (INFO(2) .NE. 0) STOP 'Error creating the output file!'
    END IF

    CALL h5gcreate_f(FID, H5GN, GID, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error creating the output group!'

    DIMS1(1) = 4
    CALL h5ltmake_dataset_int_f(GID, 'IDADIM', 1, DIMS1, IDADIM, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error writing IDADIM!'

    DIMS1(1) = 4
    CALL h5ltmake_dataset_int_f(GID, 'ISTATS', 1, DIMS1, ISTATS, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error writing ISTATS!'

    DIMS1(1) = 2
    CALL h5ltmake_dataset_double_f(GID, 'DSTATS', 1, DIMS1, DSTATS, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error writing DSTATS!'

    DIMS1(1) = N
    CALL h5ltmake_dataset_double_f(GID, 'H', 1, DIMS1, H, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error writing H!'
    CALL h5ltmake_dataset_double_f(GID, 'K', 1, DIMS1, K, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error writing K!'
    CALL h5ltmake_dataset_double_f(GID, 'SIGMA', 1, DIMS1, SIGMA, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error writing SIGMA!'

#ifdef USE_MTXOUT
    DIMS2(1) = LDA
    DIMS2(2) = N
    CALL h5ltmake_dataset_double_f(GID, 'F', 2, DIMS2, F, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error writing F!'
    CALL h5ltmake_dataset_double_f(GID, 'G', 2, DIMS2, G, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error writing G!'
    CALL h5ltmake_dataset_double_f(GID, 'V', 2, DIMS2, V, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error writing V!'
    CALL h5ltmake_dataset_double_f(GID, 'VQ', 2, DIMS2, VQ, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error writing VQ!'
#endif

    CALL h5gclose_f(GID, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error closing the output group!'

    CALL h5fclose_f(FID, INFO(2))
    IF (INFO(2) .NE. 0) STOP 'Error closing the output file!'
  END IF

  DEALLOCATE(SIGMA)
  DEALLOCATE(K)
  DEALLOCATE(H)

  DEALLOCATE(VQ)
  DEALLOCATE(V)
  DEALLOCATE(F)
  DEALLOCATE(G)

  CALL h5close_f(INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error closing HDF5!'
  BP = BLAS_SET_NUM_THREADS(OLDBP)

CONTAINS

  SUBROUTINE READCL(H5FN, H5GN, H5RN, NBMAXS, P, BP)

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(OUT) :: H5FN, H5GN, H5RN
  INTEGER, INTENT(OUT) :: NBMAXS, P, BP

  CHARACTER(LEN=8) :: CLA

  IF (COMMAND_ARGUMENT_COUNT() .GT. 6) STOP 'exe2.exe H5FN H5GN H5RN NBMAXS NTHRDS BTHRDS'

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
    READ (CLA,*) NBMAXS
  ELSE
    WRITE (*,'(A)',ADVANCE='NO') 'NBMAXS: '
    READ (*,*) NBMAXS
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 5) THEN
    CALL GET_COMMAND_ARGUMENT(5, CLA)
    READ (CLA,*) P
  ELSE
    WRITE (*,'(A)',ADVANCE='NO') 'NTHRDS: '
    READ (*,*) P
  END IF
  IF (COMMAND_ARGUMENT_COUNT() .GE. 6) THEN
    CALL GET_COMMAND_ARGUMENT(6, CLA)
    READ (CLA,*) BP
  ELSE
    WRITE (*,'(A)',ADVANCE='NO') 'BTHRDS: '
    READ (*,*) BP
  END IF

  END SUBROUTINE READCL

END PROGRAM EXE2
