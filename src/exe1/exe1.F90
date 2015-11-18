PROGRAM EXE1

  USE H5LT
  USE HDF5

  USE HZ

  IMPLICIT NONE

  INTEGER, PARAMETER :: IPART = 2

  CHARACTER(LEN=255) :: H5FN, H5RN
  CHARACTER(LEN=5) :: H5GN
  LOGICAL :: FEXIST

  INTEGER(HID_T) :: FID, GID
  INTEGER(HSIZE_T) :: DIMS1(1), DIMS2(2)

  INTEGER :: IDADIM(4), LDF, LDG, LDV, M, N
  EQUIVALENCE (IDADIM(1), LDF), (IDADIM(1), LDG), (IDADIM(1), LDV), (IDADIM(1), M), (IDADIM(1), N)

  INTEGER :: NBMAXS
  INTEGER :: P, Q, LWORK, LIWORK, ISTATS(4), INFO(2)
  INTEGER :: NBSIZE, NBSW, LCOLV
  INTEGER(8) :: NROT(2), WALLTM(3)
  DOUBLE PRECISION :: TOL, DSTATS(2)
  !DIR$ ATTRIBUTES ALIGN: 64:: ISTATS, INFO, NROT, WALLTM, DSTATS

  DOUBLE PRECISION, ALLOCATABLE :: F(:,:), G(:,:), V(:,:), VQ(:,:)
  !DIR$ ATTRIBUTES ALIGN: 64:: F, G, V, VQ
  DOUBLE PRECISION, ALLOCATABLE :: SIGMA(:), H(:), K(:), TAU(:), WORK(:)
  !DIR$ ATTRIBUTES ALIGN: 64:: SIGMA, H, K, TAU, WORK
  INTEGER, ALLOCATABLE :: IWORK(:)
  !DIR$ ATTRIBUTES ALIGN: 64:: IWORK

  DOUBLE PRECISION, EXTERNAL :: DNRM2, DLAMCH
  EXTERNAL :: DLACPY, DLASET, DGEQLF, DORGQL

  CALL READCL(H5FN, H5GN, H5RN, NBMAXS)

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

  TOL = D_MONE
  LWORK = -1
  LIWORK = -1
  DSTATS(1) = D_ZERO
  ISTATS(1) = 0
  CALL DZIMMER1BO(M, N, F, LDF, G, LDG, V, LDV, -1, TOL, NBMAXS, IPART,&
     H, K, SIGMA, DSTATS(1), LWORK, ISTATS(1), LIWORK, NBSIZE, NBSW, NROT, INFO)
  IF (INFO(1) .NE. 0) STOP 'Error querying workspace size for DZIMMER1BO!'

  LWORK = CEILING(DSTATS(1))
  LIWORK = ISTATS(1)
  LCOLV = INFO(2)

  ALLOCATE(F(LDF, LCOLV))
  ALLOCATE(G(LDG, LCOLV))
  ALLOCATE(V(LDV, LCOLV))
  ALLOCATE(VQ(LDV, LCOLV))

  ALLOCATE(SIGMA(N))
  ALLOCATE(H(N))
  ALLOCATE(K(N))
  ALLOCATE(TAU(N))

  DSTATS = D_ZERO
  ISTATS(2) = -1
  INFO(2) = 0
  CALL DGEQLF(M, N, V, LDV, TAU, DSTATS(1), ISTATS(2), INFO(2))
  IF (INFO(2) .NE. 0) STOP 'DGEQLF workspace query'
  CALL DORGQL(M, N, N, VQ, LDV, TAU, DSTATS(2), ISTATS(2), INFO(2))
  IF (INFO(2) .NE. 0) STOP 'DORGQL workspace query'
  LWORK = MAX(LWORK, CEILING(MAX(DSTATS(1), DSTATS(2))))

  ALLOCATE(WORK(LWORK))
  ALLOCATE(IWORK(LIWORK))

  DIMS2(1) = LDF
  DIMS2(2) = N
  CALL h5ltread_dataset_double_f(GID, 'F', F, DIMS2, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error reading F!'

  DIMS2(1) = LDG
  DIMS2(2) = N
  CALL h5ltread_dataset_double_f(GID, 'G', G, DIMS2, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error reading G!'

  CALL h5gclose_f(GID, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error closing the input group!'

  CALL h5fclose_f(FID, INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error closing the input file!'

  IF (LDF .GT. M) THEN
     P = LDV - M
     Q = M + 1
     CALL DLASET('A', P, N, D_ZERO, D_ZERO, F(Q, 1), LDF)
  END IF

  IF (LDG .GT. M) THEN
     P = LDV - M
     Q = M + 1
     CALL DLASET('A', P, N, D_ZERO, D_ZERO, G(Q, 1), LDG)
  END IF

  IF (LDV .GT. M) THEN
     P = LDV - M
     Q = M + 1
     CALL DLASET('A', P, N, D_ZERO, D_ZERO, V(Q, 1), LDV)
     CALL DLASET('A', P, N, D_ZERO, D_ZERO, VQ(Q, 1), LDV)
  END IF

  WRITE (*,*) N, NBMAXS
  CALL TIMER_START(WALLTM)

  NBSIZE = 0
  NBSW = 0
  NROT = 0_8
  INFO = 0

  CALL DZIMMER1BO(M, N, F, LDF, G, LDG, V, LDV, -1, TOL, NBMAXS, IPART,&
     H, K, SIGMA, WORK, LWORK, IWORK, LIWORK, NBSIZE, NBSW, NROT, INFO)

  IF ((INFO(1) .EQ. 0) .AND. (INFO(2) .EQ. 0)) THEN
     CALL DGEQLF(M, N, V, LDV, TAU, WORK, LWORK, INFO(2))
     IF (INFO(2) .NE. 0) STOP 'DGEQLF'

     CALL DLACPY('A', M, N, V, LDV, VQ, LDV)
     CALL DLASET('U', M, N - 1, D_ZERO, D_ZERO, V(1, 2), LDV)

     CALL DORGQL(M, N, N, VQ, LDV, TAU, WORK, LWORK, INFO(2))
     IF (INFO(2) .NE. 0) STOP 'DORGQL'
  ELSE
     SIGMA = D_ZERO
  END IF

  CALL TIMER_STOP(WALLTM)

  DSTATS(1) = TIMER2DBLE(WALLTM)
  DSTATS(2) = TOL

  WRITE (*,*) 'WALL_TIME ', DSTATS(1)
  WRITE (*,*) 'NBSIZE ', NBSIZE
  WRITE (*,*) 'NBSW ', NBSW
  WRITE (*,*) 'NROT ', NROT(1), ' ', NROT(2)
  WRITE (*,*) 'INFO ', INFO(1), ' ', INFO(2)

  ISTATS(1) = INT(MOD(NROT(1), 2147483648_8))
  ISTATS(2) = INT(NROT(1) / 2147483648_8)
  ISTATS(3) = NBSW
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
     CALL h5ltmake_dataset_double_f(GID, 'SIGMA', 1, DIMS1, SIGMA, INFO(2))
     IF (INFO(2) .NE. 0) STOP 'Error writing SIGMA!'
     CALL h5ltmake_dataset_double_f(GID, 'H', 1, DIMS1, H, INFO(2))
     IF (INFO(2) .NE. 0) STOP 'Error writing H!'
     CALL h5ltmake_dataset_double_f(GID, 'K', 1, DIMS1, K, INFO(2))
     IF (INFO(2) .NE. 0) STOP 'Error writing K!'

     ! DIMS2(1) = LDV
     ! DIMS2(2) = N
     ! CALL h5ltmake_dataset_double_f(GID, 'F', 2, DIMS2, F, INFO(2))
     ! IF (INFO(2) .NE. 0) STOP 'Error writing F!'
     ! CALL h5ltmake_dataset_double_f(GID, 'G', 2, DIMS2, G, INFO(2))
     ! IF (INFO(2) .NE. 0) STOP 'Error writing G!'
     ! CALL h5ltmake_dataset_double_f(GID, 'V', 2, DIMS2, V, INFO(2))
     ! IF (INFO(2) .NE. 0) STOP 'Error writing V!'
     ! CALL h5ltmake_dataset_double_f(GID, 'VQ', 2, DIMS2, VQ, INFO(2))
     ! IF (INFO(2) .NE. 0) STOP 'Error writing VQ!'

     CALL h5gclose_f(GID, INFO(2))
     IF (INFO(2) .NE. 0) STOP 'Error closing the output group!'

     CALL h5fclose_f(FID, INFO(2))
     IF (INFO(2) .NE. 0) STOP 'Error closing the output file!'
  END IF

  DEALLOCATE(IWORK)
  DEALLOCATE(WORK)
  DEALLOCATE(TAU)

  DEALLOCATE(K)
  DEALLOCATE(H)
  DEALLOCATE(SIGMA)

  DEALLOCATE(VQ)
  DEALLOCATE(V)
  DEALLOCATE(F)
  DEALLOCATE(G)

  CALL h5close_f(INFO(2))
  IF (INFO(2) .NE. 0) STOP 'Error closing HDF5!'

CONTAINS

  SUBROUTINE READCL(H5FN, H5GN, H5RN, NBMAXS)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(OUT) :: H5FN, H5GN, H5RN
    INTEGER, INTENT(OUT) :: NBMAXS

    CHARACTER(LEN=8) :: CLA

    IF (COMMAND_ARGUMENT_COUNT() .GT. 4) STOP 'exe1.exe H5FN H5GN H5RN NBMAXS'

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

  END SUBROUTINE READCL

END PROGRAM EXE1
