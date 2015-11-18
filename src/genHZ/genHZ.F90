program genHZ

  use HDF5
  use H5LT

  implicit none

  include 'WP.f'

  double precision, parameter :: zero = 0.0d0

  ! command-line parameters
  integer :: seed, n, idist_f, idist_g, idist_x, info
  double precision :: eps_f, scal_f, eps_g, scal_g, eps_x, scal_x
  character(len=256) :: sigma_f, sigma_g, lambda_x, fil, grp

  logical :: fex
  integer(hid_t) :: fid, gid

  integer :: oseed(4), iseed(4), lda, npos, k, l, m, p
  double precision :: tola, tolb, ulp, unfl
  real(WP) :: h

  integer, allocatable :: iwork(:)

  double precision, allocatable :: ds_f(:), ds_g(:), ds(:), dl_x(:), tau(:), dwork(:)
  real(WP), allocatable :: xs_f(:), xs_g(:), xl_x(:), xwork(:)

  double precision, allocatable :: df(:,:), dg(:,:), du(:,:), dv(:,:), dq(:,:)
  real(WP), allocatable :: xf(:,:), xg(:,:), xx(:,:)

  double precision, external :: dlamch, dlange
  external :: dggsvp

  call readcl(sigma_f, sigma_g, lambda_x, seed, n, fil, grp, &
       idist_f, eps_f, scal_f, idist_g, eps_g, scal_g, idist_x, eps_x, scal_x, info)
  if (info .ne. 0) then
    write (*,*) info
    stop 'readcl'
  end if

  call seedix(seed, oseed, info)
  if (info .ne. 0) then
    write (*,*) info
    stop 'seedix'
  end if
  iseed = oseed

  lda = n
  m = n
  p = mod(lda, 4)
  !if (p .ne. 0) lda = lda + (4 - p)
  if (p .ne. 0) stop '{LDA=N} mod 4 <> 0'

  allocate(ds_f(n))
  if (idist_f .ne. 0) then
    call genlam(n, iseed, idist_f, eps_f, scal_f, ds_f, npos, info)
  else
    call txtlam(sigma_f, n, ds_f, npos, info)
  end if
  if (info .ne. 0) then
    write (*,*) info
    stop 'sigma_f'
  end if
  if (npos .ne. n) then
    write (*,*) npos, '<', n
    stop 'sigma_f'
  end if

  allocate(ds_g(n))
  if (idist_g .ne. 0) then
    call genlam(n, iseed, idist_g, eps_g, scal_g, ds_g, npos, info)
  else
    call txtlam(sigma_g, n, ds_g, npos, info)
  end if
  if (info .ne. 0) then
    write (*,*) info
    stop 'sigma_g'
  end if
  if (npos .ne. n) then
    write (*,*) npos, '<', n
    stop 'sigma_g'
  end if

  allocate(ds(n))
  do p = 1, n
    ds(p) = ds_f(p) / ds_g(p)
  end do

  allocate(dl_x(n))
  if (idist_x .ne. 0) then
    call genlam(n, iseed, idist_x, eps_x, scal_x, dl_x, npos, info)
  else
    call txtlam(lambda_x, n, dl_x, npos, info)
  end if
  if (info .ne. 0) then
    write (*,*) info
    stop 'lambda_x'
  end if

  allocate(df(lda,n))
  allocate(dg(lda,n))
  allocate(du(lda,n))

  allocate(xf(n,n))
  allocate(xg(n,n))
  allocate(xx(n,n))

  allocate(xs_f(n))
  allocate(xs_g(n))
  allocate(xl_x(n))

  do p = 1, n
    xs_f(p) = real(ds_f(p),WP)
    xs_g(p) = real(ds_g(p),WP)
    xl_x(p) = real(dl_x(p),WP)

    h = hypot(xs_f(p), xs_g(p))
    xs_f(p) = xs_f(p) / h
    xs_g(p) = xs_g(p) / h

    ds_f(p) = dble(xs_f(p))
    ds_g(p) = dble(xs_g(p))
  end do

  p = 3 * n
  allocate(xwork(p))
  call gendat(n, iseed, xs_f, xs_g, xl_x, xf, xg, xx, df, dg, du, lda, xwork, info)
  if (info .ne. 0) then
    write (*,*) info
    stop 'gendat'
  end if
  deallocate(xwork)

  deallocate(xl_x)
  deallocate(xs_g)
  deallocate(xs_f)

  deallocate(xx)
  deallocate(xg)
  deallocate(xf)

  allocate(dv(lda,n))
  allocate(dq(lda,n))

  allocate(tau(n))
  allocate(dwork(p))
  allocate(iwork(n))

  p = n
  tola = dlange('1', m, n, df, lda, dwork)
  tolb = dlange('1', p, n, dg, lda, dwork)
  ulp = dlamch('Precision')
  unfl = dlamch('Safe Minimum')
  tola = max(m,n) * max(tola,unfl) * ulp
  tolb = max(p,n) * max(tolb,unfl) * ulp

  call dggsvp('U', 'V', 'Q', m, p, n, df, lda, dg, lda, tola, tolb, k, l, du, lda, dv, lda, dq, lda, iwork, tau, dwork, info)
  if (info .ne. 0) then
    write (*,*) info
    stop 'dggsvp'
  end if

  deallocate(iwork)
  deallocate(dwork)
  deallocate(tau)

  call h5open_f(info)
  if (info .ne. 0) then
    write (*,*) info
    stop 'h5open_f'
  end if

  inquire(file=fil, exist=fex)
  if (fex) then
    call h5fopen_f(fil, H5F_ACC_RDWR_F, fid, info)
  else
    call h5fcreate_f(fil, H5F_ACC_TRUNC_F, fid, info)
  end if
  if (info .ne. 0) then
    write (*,*) info
    if (fex) then
      stop 'h5open_f'
    else
      stop 'h5create_f'
    end if
  end if

  call h5gcreate_f(fid, grp, gid, info)
  if (info .ne. 0) then
    write (*,*) info
    stop 'h5gcreate_f'
  end if

  call h5wrds(gid, lda, n, npos, k, l, seed, oseed, ds_f, ds_g, ds, dl_x, df, dg, du, dv, dq, &
       idist_f, idist_g, idist_x, eps_f, eps_g, eps_x, scal_f, scal_g, scal_x, tola, tolb, info)
  if (info .ne. 0) then
    write (*,*) info
    stop 'h5wrds'
  end if

  deallocate(dq)
  deallocate(dv)
  deallocate(du)
  deallocate(dg)
  deallocate(df)

  deallocate(dl_x)
  deallocate(ds)
  deallocate(ds_g)
  deallocate(ds_f)

  call h5gclose_f(gid, info)
  if (info .ne. 0) then
    write (*,*) info
    stop 'h5gclose_f'
  end if

  call h5fclose_f(fid, info)
  if (info .ne. 0) then
    write (*,*) info
    stop 'h5fclose_f'
  end if

  call h5close_f(info)
  if (info .ne. 0) then
    write (*,*) info
    stop 'h5close_f'
  end if

contains

  subroutine readcl(sigma_f, sigma_g, lambda_x, seed, n, fil, grp, &
       idist_f, eps_f, scal_f, idist_g, eps_g, scal_g, idist_x, eps_x, scal_x, info)

    implicit none

    integer, intent(out) :: seed, n, idist_f, idist_g, idist_x, info
    double precision, intent(out) :: eps_f, scal_f, eps_g, scal_g, eps_x, scal_x
    character(len=*), intent(out) :: sigma_f, sigma_g, lambda_x, fil, grp

    integer, parameter :: nrqp = 7
    double precision, parameter :: zero = 0.0d0
    integer :: nxta
    character(len=256) :: cas

    seed = 0
    n = 0
    idist_f = 0
    idist_g = 0
    idist_x = 0
    info = 0

    eps_f = zero
    scal_f = zero
    eps_g = zero
    scal_g = zero
    eps_x = zero
    scal_x = zero

    sigma_f = ''
    sigma_g = ''
    lambda_x = ''
    fil = ''
    grp = ''

    nxta = nrqp
    cas = ''

    if (command_argument_count() .lt. nrqp) then
      write (*,*) 'genHZ.exe SIGMA_F SIGMA_G LAMBDA_X SEEDIX N FILE GROUP [ SIG|LAM_PARAMS ]'
      write (*,*) '>> COMMAND LINE (INPUT) ARGUMENTS <<'
      write (*,*) 'SIGMA_F : \Sigma(F); 1, 3, or FILENAME'
      write (*,*) 'SIGMA_G : \Sigma(G); 1, 3, or FILENAME'
      write (*,*) 'LAMBDA_X: \Lambda(X); 1, 2, 3, or FILENAME'
      write (*,*) 'IDIST123: 1 [uniform (0,1)], 2 [uniform(-1,1)], or 3 [normal(0,1)]'
      write (*,*) 'FILENAME: SIG|LAM.txt: max 256 chars, >= N lines [each line = one real value]'
      write (*,*) 'SEEDIX  : index of hard-coded pRNG seed (see seedix.F90); 1 or 2'
      write (*,*) 'N       : order of the output matrix: > 0'
      write (*,*) 'FILE.h5 : output HDF5 file (may exist): max 256 chars'
      write (*,*) 'GROUP   : output HDF5 group (must NOT exist): max 256 chars'
      write (*,*) 'SIG|LAM ; SIG|LAM_PARAMS if SIG|LAM is IDIST123'
      write (*,*) ' EPS_F  : \sigma''_i survives iff |\sigma''_i| > EPS_F'
      write (*,*) ' SCALE_F: final \sigma_i = \sigma''_i * SCALE_F'
      write (*,*) ' EPS_G  : \sigma''_i survives iff |\sigma''_i| > EPS_G'
      write (*,*) ' SCALE_G: final \sigma_i = \sigma''_i * SCALE_G'
      write (*,*) ' EPS_X  : \lambda''_i survives iff |\lambda''_i| > EPS_X'
      write (*,*) ' SCALE_X: final \lambda_i = \lambda''_i * SCALE_X'
      write (*,*) '<< OUTPUT DATASETS IN FILE.h5/GROUP >>'
      write (*,*) 'IDADIM  : integer(4) { {LDA=M=P=}N, NPOS_X, K, L }'
      write (*,*) 'ISEED   : integer(4); initial seed for (d|x)laran pRNG (see LaPACK dlaran.f)'
      write (*,*) 'SIGMA_F : double precision(N); normalized: \sigma_F^2 + \sigma_G^2 = 1'
      write (*,*) 'SIGMA_G : double precision(N); normalized: \sigma_F^2 + \sigma_G^2 = 1'
      write (*,*) 'SIGMA   : double precision(N); \sigma_F / \sigma_G'
      write (*,*) 'LAMBDA_X: double precision(N); as read/generated'
      write (*,*) 'F       : double precision(LDA,N) = U_F \SIGMA_F X'
      write (*,*) 'G       : double precision(LDA,N) = U_G \SIGMA_G X'
      write (*,*) 'U       : double precision(LDA,N); see LaPACK dggsvp.f'
      write (*,*) 'V       : double precision(LDA,N); see LaPACK dggsvp.f'
      write (*,*) 'Q       : double precision(LDA,N); see LaPACK dggsvp.f'
      write (*,*) 'TOL     : double precision(2) = (/ TOLA, TOLB /); see LaPACK dggsvd.f'
      info = 1
      return
    end if

    call get_command_argument(1, sigma_f, status=info)
    if (info .ne. 0) then
      info = -1
      return
    end if
    if (trim(sigma_f) .eq. '1') then
      idist_f = 1
    else if (trim(sigma_f) .eq. '3') then
      idist_f = 3
    end if

    call get_command_argument(2, sigma_g, status=info)
    if (info .ne. 0) then
      info = -2
      return
    end if
    if (trim(sigma_g) .eq. '1') then
      idist_g = 1
    else if (trim(sigma_g) .eq. '3') then
      idist_g = 3
    end if

    call get_command_argument(3, lambda_x, status=info)
    if (info .ne. 0) then
      info = -3
      return
    end if
    if (trim(lambda_x) .eq. '1') then
      idist_x = 1
    else if (trim(lambda_x) .eq. '2') then
      idist_x = 2
    else if (trim(lambda_x) .eq. '3') then
      idist_x = 3
    end if

    call get_command_argument(4, cas, status=info)
    if (info .ne. 0) then
      info = -4
      return
    end if
    read (cas,*) seed
    if (seed .le. 0) then
      info = -4
      return
    end if

    call get_command_argument(5, cas, status=info)
    if (info .ne. 0) then
      info = -5
      return
    end if
    read (cas,*) n
    if (n .le. 0) then
      info = -5
      return
    end if

    call get_command_argument(6, fil, status=info)
    if (info .ne. 0) then
      info = -6
      return
    end if

    call get_command_argument(7, grp, status=info)
    if (info .ne. 0) then
      info = -7
      return
    end if

    if (idist_f .ne. 0) then
      nxta = nxta + 1
      call get_command_argument(nxta, cas, status=info)
      if (info .ne. 0) then
        info = -nxta
        return
      end if
      read (cas,*) eps_f
      if (eps_f .lt. zero) then
        info = -nxta
        return
      end if
      nxta = nxta + 1
      call get_command_argument(nxta, cas, status=info)
      if (info .ne. 0) then
        info = -nxta
        return
      end if
      read (cas,*) scal_f
      if (scal_f .eq. zero) then
        info = -nxta
        return
      end if
    end if

    if (idist_g .ne. 0) then
      nxta = nxta + 1
      call get_command_argument(nxta, cas, status=info)
      if (info .ne. 0) then
        info = -nxta
        return
      end if
      read (cas,*) eps_g
      if (eps_g .lt. zero) then
        info = -nxta
        return
      end if
      nxta = nxta + 1
      call get_command_argument(nxta, cas, status=info)
      if (info .ne. 0) then
        info = -nxta
        return
      end if
      read (cas,*) scal_g
      if (scal_g .eq. zero) then
        info = -nxta
        return
      end if
    end if

    if (idist_x .ne. 0) then
      nxta = nxta + 1
      call get_command_argument(nxta, cas, status=info)
      if (info .ne. 0) then
        info = -nxta
        return
      end if
      read (cas,*) eps_x
      if (eps_x .lt. zero) then
        info = -nxta
        return
      end if
      nxta = nxta + 1
      call get_command_argument(nxta, cas, status=info)
      if (info .ne. 0) then
        info = -nxta
        return
      end if
      read (cas,*) scal_x
      if (scal_x .eq. zero) then
        info = -nxta
        return
      end if
    end if
  end subroutine readcl

  subroutine h5wrds(gid, lda, n, npos, k, l, six, iseed, ds_f, ds_g, ds, dl_x, df, dg, du, dv, dq, &
       idist_f, idist_g, idist_x, eps_f, eps_g, eps_x, scal_f, scal_g, scal_x, tola, tolb, info)

    implicit none

    integer(hid_t), intent(in) :: gid
    integer, intent(in) :: lda, n, npos, k, l, six, iseed(4), idist_f, idist_g, idist_x
    double precision, intent(in) :: ds_f(n), ds_g(n), ds(n), dl_x(n), df(lda,n), dg(lda,n), du(lda,n), dv(lda,n), dq(lda,n), &
         eps_f, eps_g, eps_x, scal_f, scal_g, scal_x, tola, tolb
    integer, intent(out) :: info

    integer(hsize_t) :: dims(2)
    integer :: idadim(4), idist(4)
    double precision :: eps(3), scal(3), tol(2)

    info = 0

    idadim = (/ n, npos, k, l /)
    dims(1) = 4
    call h5ltmake_dataset_int_f(gid, 'IDADIM', 1, dims, idadim, info)
    if (info .ne. 0) then
      info = 1
      return
    end if

    dims(1) = 4
    call h5ltmake_dataset_int_f(gid, 'ISEED', 1, dims, iseed, info)
    if (info .ne. 0) then
      info = 2
      return
    end if
    
    dims(1) = n
    call h5ltmake_dataset_double_f(gid, 'SIGMA_F', 1, dims, ds_f, info)
    if (info .ne. 0) then
      info = 3
      return
    end if

    dims(1) = n
    call h5ltmake_dataset_double_f(gid, 'SIGMA_G', 1, dims, ds_g, info)
    if (info .ne. 0) then
      info = 4
      return
    end if

    dims(1) = n
    call h5ltmake_dataset_double_f(gid, 'SIGMA', 1, dims, ds, info)
    if (info .ne. 0) then
      info = 5
      return
    end if

    dims(1) = n
    call h5ltmake_dataset_double_f(gid, 'LAMBDA_X', 1, dims, dl_x, info)
    if (info .ne. 0) then
      info = 6
      return
    end if

    dims(1) = lda
    dims(2) = n
    call h5ltmake_dataset_double_f(gid, 'F', 2, dims, df, info)
    if (info .ne. 0) then
      info = 7
      return
    end if

    dims(1) = lda
    dims(2) = n
    call h5ltmake_dataset_double_f(gid, 'G', 2, dims, dg, info)
    if (info .ne. 0) then
      info = 8
      return
    end if

    dims(1) = lda
    dims(2) = n
    call h5ltmake_dataset_double_f(gid, 'U', 2, dims, du, info)
    if (info .ne. 0) then
      info = 9
      return
    end if

    dims(1) = lda
    dims(2) = n
    call h5ltmake_dataset_double_f(gid, 'V', 2, dims, dv, info)
    if (info .ne. 0) then
      info = 10
      return
    end if

    dims(1) = lda
    dims(2) = n
    call h5ltmake_dataset_double_f(gid, 'Q', 2, dims, dq, info)
    if (info .ne. 0) then
      info = 11
      return
    end if

    idist = (/ idist_f, idist_g, idist_x, six /)
    dims(1) = 4
    call h5ltmake_dataset_int_f(gid, 'IDIST', 1, dims, idist, info)
    if (info .ne. 0) then
      info = 12
      return
    end if

    eps = (/ eps_f, eps_g, eps_x /)
    dims(1) = 3
    call h5ltmake_dataset_double_f(gid, 'EPS', 1, dims, eps, info)
    if (info .ne. 0) then
      info = 13
      return
    end if

    scal = (/ scal_f, scal_g, scal_x /)
    dims(1) = 3
    call h5ltmake_dataset_double_f(gid, 'SCALE', 1, dims, scal, info)
    if (info .ne. 0) then
      info = 14
      return
    end if

    tol = (/ tola, tolb /)
    dims(1) = 2
    call h5ltmake_dataset_double_f(gid, 'TOL', 1, dims, tol, info)
    if (info .ne. 0) then
      info = 15
      return
    end if

  end subroutine h5wrds

end program genHZ
