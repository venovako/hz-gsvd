PROGRAM DGSVD

  USE H5LT
  USE HDF5

  USE HZ

  IMPLICIT NONE

  CHARACTER, PARAMETER :: JOBU = 'U', JOBV = 'V', JOBQ = 'Q'

  CHARACTER(LEN=255) :: H5FN, H5RN
  CHARACTER(LEN=5) :: H5GN
  LOGICAL :: FEXIST

  INTEGER(HID_T) :: FID, GID
  INTEGER(HSIZE_T) :: DIMS1(1), DIMS2(2)
  !DIR$ ATTRIBUTES ALIGN: 64:: DIMS1, DIMS2

  INTEGER :: IDADIM(4), LDA, LDB, LDU, LDV, LDQ, M, N, P, K, L
  !DIR$ ATTRIBUTES ALIGN: 64:: IDADIM
  EQUIVALENCE (IDADIM(1), LDA), (IDADIM(1), LDB), (IDADIM(1), LDU), (IDADIM(1), LDV), (IDADIM(1), LDQ)
  EQUIVALENCE (IDADIM(1), M), (IDADIM(1), N), (IDADIM(1), P)
  EQUIVALENCE (IDADIM(3), K)
  EQUIVALENCE (IDADIM(4), L)

  DOUBLE PRECISION :: TOL(2), SMAX, TEMP
  !DIR$ ATTRIBUTES ALIGN: 64:: TOL
  INTEGER :: NCYCLE, I, J, IBND, ISUB, INFO

  INTEGER :: ISTATS(2)
  !DIR$ ATTRIBUTES ALIGN: 64:: ISTATS
  INTEGER(8) :: WALLTM(3)
  !DIR$ ATTRIBUTES ALIGN: 64:: WALLTM
  DOUBLE PRECISION :: DSTATS(1)
  !DIR$ ATTRIBUTES ALIGN: 64:: DSTATS

  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), U(:,:), V(:,:), Q(:,:)
  !DIR$ ATTRIBUTES ALIGN: 64:: A, B, U, V, Q
  DOUBLE PRECISION, ALLOCATABLE :: ALPHA(:), BETA(:), WORK(:)
  !DIR$ ATTRIBUTES ALIGN: 64:: ALPHA, BETA, WORK
  INTEGER, ALLOCATABLE :: IWORK(:)
  !DIR$ ATTRIBUTES ALIGN: 64:: IWORK

  EXTERNAL :: DTGSJA, DCOPY, DLASET

  CALL READCL(H5FN, H5GN, H5RN)

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

  ALLOCATE(A(LDA,MAX(1,N)))
  DIMS2(1) = LDA
  DIMS2(2) = N
  CALL h5ltread_dataset_double_f(GID, 'F', A, DIMS2, INFO)
  IF (INFO .NE. 0) STOP 'Error reading F to A!'
  IF (LDA .GT. M) CALL DLASET('A', LDA - M, N, D_ZERO, D_ZERO, A(M+1, 1), LDA)

  ALLOCATE(B(LDB,MAX(1,N)))
  DIMS2(1) = LDB
  DIMS2(2) = N
  CALL h5ltread_dataset_double_f(GID, 'G', B, DIMS2, INFO)
  IF (INFO .NE. 0) STOP 'Error reading G to B!'
  IF (LDB .GT. P) CALL DLASET('A', LDB - P, N, D_ZERO, D_ZERO, B(P+1, 1), LDB)

  ALLOCATE(U(LDU,MAX(1,M)))
  DIMS2(1) = LDU
  DIMS2(2) = M
  CALL h5ltread_dataset_double_f(GID, 'U', U, DIMS2, INFO)
  IF (INFO .NE. 0) STOP 'Error reading U!'
  IF (LDU .GT. M) CALL DLASET('A', LDU - M, M, D_ZERO, D_ZERO, U(M+1, 1), LDU)

  ALLOCATE(V(LDV,MAX(1,P)))
  DIMS2(1) = LDV
  DIMS2(2) = P
  CALL h5ltread_dataset_double_f(GID, 'V', V, DIMS2, INFO)
  IF (INFO .NE. 0) STOP 'Error reading V!'
  IF (LDV .GT. P) CALL DLASET('A', LDV - P, P, D_ZERO, D_ZERO, V(P+1, 1), LDV)

  ALLOCATE(Q(LDQ,MAX(1,N)))
  DIMS2(1) = LDQ
  DIMS2(2) = N
  CALL h5ltread_dataset_double_f(GID, 'Q', Q, DIMS2, INFO)
  IF (INFO .NE. 0) STOP 'Error reading Q!'
  IF (LDQ .GT. N) CALL DLASET('A', LDQ - N, N, D_ZERO, D_ZERO, Q(N+1, 1), LDQ)

  DIMS1(1) = 2
  CALL h5ltread_dataset_double_f(GID, 'TOL', TOL, DIMS1, INFO)
  IF (INFO .NE. 0) STOP 'Error reading TOL!'

  CALL h5gclose_f(GID, INFO)
  IF (INFO .NE. 0) STOP 'Error closing the input group!'

  CALL h5fclose_f(FID, INFO)
  IF (INFO .NE. 0) STOP 'Error closing the input file!'

  ALLOCATE(ALPHA(MAX(1,N))); ALPHA = D_ZERO
  ALLOCATE(BETA(MAX(1,N))); BETA = D_ZERO

  ALLOCATE(IWORK(MAX(1,N))); IWORK = 0
  ALLOCATE(WORK(MAX(1,2*N))); WORK = D_ZERO

  NCYCLE = 0; INFO = 0
  CALL TIMER_START(WALLTM)

  CALL DTGSJA(JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB, TOL(1), TOL(2), ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, &
       WORK, NCYCLE, INFO)
  IF (INFO .EQ. 0) THEN
     CALL DCOPY(N, ALPHA, 1, WORK, 1)
     IBND = MIN(L, M-K)
     DO I = 1, IBND
        ISUB = I
        SMAX = WORK(K+I)
        DO J = I + 1, IBND
           TEMP = WORK(K+J)
           IF (TEMP .GT. SMAX) THEN
              ISUB = J
              SMAX = TEMP
           END IF
        END DO
        IF (ISUB .NE. I) THEN
           WORK(K+ISUB) = WORK(K+I)
           WORK(K+I) = SMAX
           IWORK(K+I) = K + ISUB
        ELSE
           IWORK(K+I) = K + I
        END IF
     END DO
  END IF

  CALL TIMER_STOP(WALLTM)
  DEALLOCATE(WORK)

  DSTATS(1) = TIMER2DBLE(WALLTM)
  WRITE (*,*) 'WALL_TIME    ', DSTATS(1), ' s'
  ISTATS(1) = INFO
  WRITE (*,*) 'INFO         ', ISTATS(1)
  ISTATS(2) = NCYCLE
  WRITE (*,*) 'NCYCLE       ', ISTATS(2)

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

     DIMS1(1) = 2
     CALL h5ltmake_dataset_int_f(GID, 'ISTATS', 1, DIMS1, ISTATS, INFO)
     IF (INFO .NE. 0) STOP 'Error writing ISTATS!'

     DIMS1(1) = 1
     CALL h5ltset_attribute_double_f(GID, 'ISTATS', 'WALL_TIME', DSTATS, DIMS1(1), INFO)
     IF (INFO .NE. 0) STOP 'Error setting attribute WALL_TIME!'

     DIMS1(1) = N
     CALL h5ltmake_dataset_double_f(GID, 'ALPHA', 1, DIMS1, ALPHA, INFO)
     IF (INFO .NE. 0) STOP 'Error writing ALPHA!'
     CALL h5ltmake_dataset_double_f(GID, 'BETA', 1, DIMS1, BETA, INFO)
     IF (INFO .NE. 0) STOP 'Error writing BETA!'
     CALL h5ltmake_dataset_int_f(GID, 'IWORK', 1, DIMS1, IWORK, INFO)
     IF (INFO .NE. 0) STOP 'Error writing IWORK!'

     ! DIMS2(1) = LDA
     ! DIMS2(2) = N
     ! CALL h5ltmake_dataset_double_f(GID, 'A', 2, DIMS2, A, INFO)
     ! IF (INFO .NE. 0) STOP 'Error writing A!'
     ! DIMS2(1) = LDB
     ! DIMS2(2) = N
     ! CALL h5ltmake_dataset_double_f(GID, 'B', 2, DIMS2, B, INFO)
     ! IF (INFO .NE. 0) STOP 'Error writing B!'
     ! DIMS2(1) = LDU
     ! DIMS2(2) = M
     ! CALL h5ltmake_dataset_double_f(GID, 'U', 2, DIMS2, U, INFO)
     ! IF (INFO .NE. 0) STOP 'Error writing U!'
     ! DIMS2(1) = LDV
     ! DIMS2(2) = P
     ! CALL h5ltmake_dataset_double_f(GID, 'V', 2, DIMS2, V, INFO)
     ! IF (INFO .NE. 0) STOP 'Error writing V!'
     ! DIMS2(1) = LDQ
     ! DIMS2(2) = N
     ! CALL h5ltmake_dataset_double_f(GID, 'Q', 2, DIMS2, Q, INFO)
     ! IF (INFO .NE. 0) STOP 'Error writing Q!'

     CALL h5gclose_f(GID, INFO)
     IF (INFO .NE. 0) STOP 'Error closing the output group!'

     CALL h5fclose_f(FID, INFO)
     IF (INFO .NE. 0) STOP 'Error closing the output file!'
  END IF

  DEALLOCATE(IWORK)
  DEALLOCATE(BETA)
  DEALLOCATE(ALPHA)

  DEALLOCATE(Q)
  DEALLOCATE(V)
  DEALLOCATE(U)
  DEALLOCATE(B)
  DEALLOCATE(A)

  CALL h5close_f(INFO)
  IF (INFO .NE. 0) STOP 'Error closing HDF5!'

CONTAINS

  SUBROUTINE READCL(H5FN, H5GN, H5RN)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(OUT) :: H5FN, H5GN, H5RN

    IF (COMMAND_ARGUMENT_COUNT() .GT. 3) STOP 'DGSVD.exe H5FN H5GN H5RN'

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

  END SUBROUTINE READCL

END PROGRAM DGSVD
