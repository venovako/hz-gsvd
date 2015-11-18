subroutine txtlam(fn, n, lam, npos, info)

  implicit none

  include 'WP.f'

  character(len=*), intent(in) :: fn
  integer, intent(in) :: n
  double precision, intent(out) :: lam(n)
  integer, intent(out) :: npos, info

  double precision, parameter :: zero = 0.0d0

  integer :: i
  logical :: fex

  inquire(file=trim(fn), exist=fex)
  if (.not. fex) then
    info = -1
  else if (n .lt. 0) then
    info = -2
  else
    info = 0
  end if
  if (info .ne. 0) return

  open(unit=1,file=trim(fn),action='READ',status='OLD')

  npos = 0
  do i = 1, n
    read (1,*) lam(i)
    if (lam(i) .gt. zero) npos = npos + 1
  end do

  close(unit=1)

end subroutine txtlam
