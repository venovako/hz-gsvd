pure subroutine seedix(ix, iseed, info)

  implicit none

  integer, intent(in) :: ix
  integer, intent(out) :: iseed(4), info

  integer, parameter :: iseeds(4,2) = &
       reshape( (/ &
       2741, 1223, 1987, 541, &
       3281,  123,  456, 333  &
       /), (/ 4,2 /) )

  if ((ix .lt. 1) .or. (ix .gt. 2)) then
     info = -1
  else
     info = 0
     iseed = iseeds(:,ix)
  end if

end subroutine seedix
