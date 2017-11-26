subroutine sgenlam(n, iseed, idist, eps, scal, lam, npos, info)

  implicit none

  integer, intent(in) :: n, idist
  integer, intent(inout) :: iseed(4)
  real, intent(in) :: eps, scal
  real, intent(out) :: lam(n)
  integer, intent(out) :: npos, info

  real, parameter :: zero = 0.0e0, one = 1.0e0

  integer :: i

  real, external :: slarnd
  external :: seedok

  call seedok(iseed, i)

  if (n .lt. 0) then
    info = -1
  else if (i .ne. 0) then
    info = -2
  else if ((idist .lt. 1) .or. (idist .gt. 3)) then
    info = -3
  else if (eps .lt. zero) then
    info = -4
  else if (scal .eq. zero) then
    info = -5
  else
    info = 0
  end if
  if (info .ne. 0) return

  npos = 0
  i = 1

  do while (i .le. n)
    lam(i) = slarnd(idist, iseed)
    if (abs(lam(i)) .gt. eps) then
      if (lam(i) .gt. zero) npos = npos + 1
      i = i + 1
    end if
  end do

  if (scal .ne. one) then
    do i = 1, n
      lam(i) = lam(i) * scal
    end do
  end if

end subroutine sgenlam

subroutine dgenlam(n, iseed, idist, eps, scal, lam, npos, info)

  implicit none

  integer, intent(in) :: n, idist
  integer, intent(inout) :: iseed(4)
  double precision, intent(in) :: eps, scal
  double precision, intent(out) :: lam(n)
  integer, intent(out) :: npos, info

  double precision, parameter :: zero = 0.0d0, one = 1.0d0

  integer :: i

  double precision, external :: dlarnd
  external :: seedok

  call seedok(iseed, i)

  if (n .lt. 0) then
    info = -1
  else if (i .ne. 0) then
    info = -2
  else if ((idist .lt. 1) .or. (idist .gt. 3)) then
    info = -3
  else if (eps .lt. zero) then
    info = -4
  else if (scal .eq. zero) then
    info = -5
  else
    info = 0
  end if
  if (info .ne. 0) return

  npos = 0
  i = 1

  do while (i .le. n)
    lam(i) = dlarnd(idist, iseed)
    if (abs(lam(i)) .gt. eps) then
      if (lam(i) .gt. zero) npos = npos + 1
      i = i + 1
    end if
  end do

  if (scal .ne. one) then
    do i = 1, n
      lam(i) = lam(i) * scal
    end do
  end if

end subroutine dgenlam
