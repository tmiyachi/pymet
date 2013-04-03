subroutine linear_interp(in,zin,zout,xnin,znin,znout,out)
  implicit none
  integer, intent(in) :: xnin, znin, znout
  real(4), intent(in) :: in(xnin,znin), zin(xnin,znin), zout(znout)
  real(4), intent(out) :: out(xnin,znout)
  integer idx(znout)
  real(4) a, x(znin), y(znin)
  integer i, j, jo

    do i = 1, xnin
       x = zin(i,:)
       y = in(i,:)
       call idxsearch(x, zout, znin, znout, idx)
       do j = 1, znout
          jo = idx(j)
          if (jo == 0) then ! outside of bottom level is constant
             out(i,j) = y(1)
          else if (jo == znin) then !outside of top level is constant
             out(i,j) = y(znin)
          else !linear interpolation
             a = (y(jo+1) - y(jo))/(x(jo+1) - x(jo))
             out(i,j) = a*(zout(j) - x(jo)) + y(jo)
          end if
       end do
    end do

  end subroutine linear_interp
subroutine idxsearch(ar1, ar2, k1, k2, idx)
  !ar1 and ar2  assumed to be monotonically decreasing or increasing
  !idx(k) means ar2(k) is between ar1(k) and ar1(k+1)
  !idx(k) = 0 or size(ar1) indicates ar2(k) is outside the range of ar1
  integer,intent(in) :: k1, k2
  real(kind=4), intent(in)      :: ar1(k1) !sequence value to search
  real(kind=4), intent(in)      :: ar2(k2) !set of value to serach for
  integer, intent(out) :: idx(k2) !interval locations 
  
  integer :: i, j
  real(kind=4) :: val
  
  if (ar1(1) < ar1(k1)) then
     !monotonically increasing case
     do i = 1, k2
        val = ar2(i)
        if (val <= ar1(1)) then
           idx(i) = 0
           cycle
        else if ( val >= ar1(k1)) then
           idx(i) = k1
           cycle
        end if
        do j = 1, k1-1
           if (val <= ar1(j+1)) then
              idx(i) = j
              exit
           end if
        end do
     end do
  else
     !monotonically decreasing case
     do i = 1, k2
        val = ar2(i)
        if (val >= ar1(1)) then
           idx(i) = 0
           cycle
        else if ( val <= ar1(k1)) then
           idx(i) = k1
           cycle
        end if
        do j = 1, k1-1
           if (val >= ar1(j+1)) then
              idx(i) = j
              exit
           end if
        end do
     end do
  end if
  
end subroutine idxsearch


