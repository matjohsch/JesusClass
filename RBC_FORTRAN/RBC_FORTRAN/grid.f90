module grid

implicit none

contains

! --------------------------------------------------------------------------------------
function func_makegrid(x1,x2,n,c)

! builds grids, linear if scaling factor c=1, triple exponential for c=3

implicit none
real(8)::grid(n),func_makegrid(n)
real(8),intent(in)::x1,x2,c
integer,intent(in)::n

integer::i
real(8)::scale

scale = x2-x1
grid(1) = x1
grid(n) = x2
do i = 2,n-1
	grid(i) = x1+scale*((i-1.0)/(n-1.0))**c
end do

func_makegrid = grid

end function func_makegrid
! --------------------------------------------------------------------------------------

end module grid