module RBC_functions
    
contains 
subroutine sub_makegrid(x1,x2,n,c,grid)

! builds grids, linear if scaling factor c=1, triple exponential for c=3

implicit none
real(8),intent(in)::x1,x2,c
integer,intent(in)::n
real(8),dimension(n)::grid!,intent(out)
integer::i
real(8)::scale

scale = x2-x1
grid(1) = x1
grid(n) = x2
do i = 2,n-1
	grid(i) = x1+scale*((i-1.0)/(n-1.0))**c
end do

end subroutine sub_makegrid


!---------------------------------------------------------------------------------------
subroutine sub_interpolation(gridold,valuegrid,np,point,value)

! Interpolation for EXGM

implicit none
integer::np
real(8),dimension(np),intent(in)::gridold,valuegrid
real(8),intent(in)::point
real(8),intent(out)::value
real(8)::vals(2)
integer inds(2)

call sub_basefun(gridold,np,point,vals,inds)

value = vals(1)*valuegrid(inds(1))+vals(2)*valuegrid(inds(2))

end subroutine sub_interpolation
 ! --------------------------------------------------------------------------------------
subroutine sub_basefun (grid,np,point,vals,inds) 

! interpolation: basis functions 

! this subroutine returns the values and the indices of the two basis
! functions that are positive on a given x in the grid_x

implicit none

real(8),intent(in) :: point
integer , intent(in):: np
real(8), intent(in) :: grid (np)
real(8), intent(out) ::vals(2)
integer ,intent(out) ::inds(2)
integer :: i

call sub_lookup(i,point,grid,np)

vals(2) = ( point-grid(i-1) )/(grid(i)-grid(i-1))
vals(1) = ( grid(i)-point )/(grid(i)-grid(i-1))
inds(2) = i
inds(1) = i-1

end subroutine sub_basefun
! --------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------
subroutine sub_lookup(i,point,grid,np)

! lookup position of x in gridx

integer::ju,jl,jm
integer,intent(in)::np
real(8),intent(in)::point,grid(np)
integer, intent(out)::i

jl = 1      
ju = np

do
	if (ju-jl .le. 1) exit
	jm = (ju+jl)/2
	if (point .ge. grid(jm)) then
		jl=jm
	else
		ju=jm
	endif
end do

i = jl+1

end subroutine sub_lookup
! --------------------------------------------------------------------------------------
end module RBC_functions