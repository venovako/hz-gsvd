subroutine gendat(n, iseed, xs_f, xs_g, xl_x, xf, xg, xx, df, dg, dx, lda, xwork, info)

  implicit none

  include 'WP.f'

  double precision, parameter :: dzero = 0.0d0
  real(WP), parameter :: zero = 0.0e0_WP, one = 1.0e0_WP

  integer, intent(in) :: n, iseed(4), lda
  real(WP), intent(in) :: xs_f(n), xs_g(n), xl_x(n)
  double precision, intent(out) :: df(lda,n), dg(lda,n), dx(lda,n)
  real(WP), intent(out) :: xf(n,n), xg(n,n), xx(n,n), xwork(*)
  integer, intent(out) :: info

  integer :: p, q
  external :: seedok, xgemm, xlagsy, xlaror

  info = 0

  if (n .lt. 0) then
    info = -1
    return
  end if
  call seedok(iseed, info)
  if (info .ne. 0) then
    info = -2
    return
  end if

  call xlagsy(n, n-1, xl_x, xx, n, iseed, xwork, info)
  if (info .ne. 0) return
  do q = 1, n
    do p = 1, n
      dx(p,q) = dble(xx(p,q))
    end do
    do p = n+1, lda
      dx(p,q) = dzero
    end do
  end do

  do q = 1, n
    do p = 1, q-1
      xf(p,q) = zero
    end do
    xf(q,q) = xs_f(q)
    do p = q+1, n
      xf(p,q) = zero
    end do
  end do

  call xlaror('L', 'N', n, n, xf, n, iseed, xwork, info)
  if (info .ne. 0) return

  call xgemm('N', 'N', n, n, n, one, xf, n, xx, n, zero, xg, n)
  do q = 1, n
    do p = 1, n
      df(p,q) = dble(xg(p,q))
    end do
    do p = n+1, lda
      df(p,q) = dzero
    end do
  end do

  do q = 1, n
    do p = 1, q-1
      xg(p,q) = zero
    end do
    xg(q,q) = xs_g(q)
    do p = q+1, n
      xg(p,q) = zero
    end do
  end do

  call xlaror('L', 'N', n, n, xg, n, iseed, xwork, info)
  if (info .ne. 0) return

  call xgemm('N', 'N', n, n, n, one, xg, n, xx, n, zero, xf, n)
  do q = 1, n
    do p = 1, n
      dg(p,q) = dble(xf(p,q))
    end do
    do p = n+1, lda
      dg(p,q) = dzero
    end do
  end do

end subroutine gendat
