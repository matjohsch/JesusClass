module RBC_Functions
    
contains 
    
subroutine sub_makegrid(x1,x2,n,c,grid)

! builds grids, linear if scaling factor c=1, triple exponential for c=3

implicit none
real(8),intent(in)::x1,x2,c
integer,intent(in)::n
real(8),dimension(n)::grid
integer::i
real(8)::scale

scale = x2-x1
grid(1) = x1
grid(n) = x2
do i = 2,n-1
	grid(i) = x1+scale*((i-1.0)/(n-1.0))**c
end do

end subroutine sub_makegrid

subroutine sub_makerandomgrid(x1,x2,n,grid)

implicit none
real(8),intent(in)::x1,x2
integer,intent(in)::n
real(8),dimension(n)::grid
integer::i
real(8)::scale, num

call RANDOM_SEED

scale = x2-x1
grid(1) = x1
grid(n) = x2
do i = 2,n-1
    call RANDOM_NUMBER(num)
	grid(i) = x1+scale*num
end do

call sub_Bubble_Sort(grid)

end subroutine sub_makerandomgrid

!---------------------------------------------------------------------------------------
subroutine sub_interpolation(GridOld,valuegrid,np,point,value)

! Interpolation

implicit none
integer::np
real(8),dimension(np),intent(in)::GridOld,valuegrid
real(8),intent(in)::point
real(8),intent(out)::value
real(8)::vals(2)
integer inds(2)

call sub_basefun(GridOld,np,point,vals,inds)

value = vals(1)*valuegrid(inds(1))+vals(2)*valuegrid(inds(2))

if (maxval(abs(vals(:)))>1.0) then
    print* ,'extrapolation interpolation'
end if


end subroutine sub_interpolation

! --------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
subroutine sub_derivative(GridOld,valuegrid,np,point,derivative)

! Interpolation

implicit none
integer::np
real(8),dimension(np),intent(in)::GridOld,valuegrid
real(8),intent(in)::point
real(8),intent(out)::derivative
real(8)::vals(2)
real(8)::x,y

integer inds(2)

call sub_basefun(GridOld,np,point,vals,inds)

x=gridold(inds(2))-gridold(inds(1))
y=valuegrid(inds(2))-valuegrid(inds(1))

derivative=y/x

!if (maxval(abs(vals(:)))>1.0) then
!    print* ,'extrapolation derivative'
!end if


end subroutine sub_derivative

! --------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------
subroutine sub_derivative2(GridOld,valuegrid,np,point,derivative)

! Interpolation

implicit none
integer::np
real(8),dimension(np),intent(in)::GridOld,valuegrid
integer,intent(in)::point
real(8),intent(out)::derivative
real(8)::vals(2)
real(8)::x,y,derivative1,derivative2

integer inds(2)

if (point==1) then
    x=gridold(2)-gridold(1)
    y=valuegrid(2)-valuegrid(1)

    derivative=y/x
else if (point==np) then
    x=gridold(np)-gridold(np-1)
    y=valuegrid(np)-valuegrid(np-1)

    derivative=y/x
else
    x=gridold(point)-gridold(point-1)
    y=valuegrid(point)-valuegrid(point-1)
    
    derivative1=y/x
    
    x=gridold(point+1)-gridold(point)
    y=valuegrid(point+1)-valuegrid(point)
    derivative2=y/x
    
    derivative=(derivative1+derivative2)/2.0
end if


if (maxval(abs(vals(:)))>1.0) then
    print* ,'extrapolation'
end if


end subroutine sub_derivative2

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

SUBROUTINE sub_Bubble_Sort(a)
  REAL(8), INTENT(in out), DIMENSION(:) :: a
  REAL(8) :: temp
  INTEGER :: i, j
  LOGICAL :: swapped = .TRUE.
 
  DO j = SIZE(a)-1, 1, -1
    swapped = .FALSE.
    DO i = 1, j
      IF (a(i) > a(i+1)) THEN
        temp = a(i)
        a(i) = a(i+1)
        a(i+1) = temp
        swapped = .TRUE.
      END IF
    END DO
    IF (.NOT. swapped) EXIT
  END DO
END SUBROUTINE sub_Bubble_Sort

end module RBC_Functions