pure subroutine seedok(iseed, info)

  implicit none

      INTEGER WP
#ifdef USE_GNU
      PARAMETER (WP=10)
#else
      PARAMETER (WP=16)
#endif

  integer, intent(in) :: iseed(4)
  integer, intent(out) :: info

  if ((iseed(1) .lt. 0) .or. (iseed(1) .gt. 4095)) then
     info = -1
  else if ((iseed(2) .lt. 0) .or. (iseed(2) .gt. 4095)) then
     info = -2
  else if ((iseed(3) .lt. 0) .or. (iseed(3) .gt. 4095)) then
     info = -3
  else if ((iseed(4) .lt. 0) .or. (iseed(4) .gt. 4095) .or. (mod(iseed(4),2) .eq. 0)) then
     info = -4
  else
     info = 0
  end if

end subroutine seedok
