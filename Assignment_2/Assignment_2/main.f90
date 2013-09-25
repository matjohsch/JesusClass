!============================================================================
! Name : RBC_F90.f90
! Description : Basic RBC model with labour and partial depreciation
! Date : September 23, 2013
!============================================================================
    
program main_RBC_Labour

use RBC_var
use RBC_param
use RBC_functions

implicit none

integer i,nold,nnew
real(8), dimension(:),allocatable :: aaa,vGridCapital
real(8), dimension(:,:),allocatable :: bbb, mValueFunction, mPolicyCapital, mPolicyLabour!,expectedValueFunction
integer,dimension(3)::scheme
scheme(1)=100
scheme(2)=1000
scheme(3)=25000

!----------------------------------------------------------------
! 0. PreSettings
!----------------------------------------------------------------

! Steady state Labour Supply (Calibration Target)
! Either choose utility parameter phi or steady state labour supply
labourSteadyState=0.3333333333      
x=(aalpha*bbeta)/(1.0-bbeta+bbeta*ddelta)
pphi = ((1.0-aalpha)/(labourSteadyState**2.0))/(1.0-ddelta*x)

!pphi=0.0 !6.72  

! Define the initial guess for the Value function
! 1:= constant function equal zeros; 2:=constant function equal to Steady State Value;
! 3:= value function obtained from the deterministic version of the problem
initial_guess = 2

! Accelarator
! To use the accelarator Howard:=true
howard = 1

!----------------------------------------------------------------
! 1. Steady State
!----------------------------------------------------------------
call sub_steady_state

!----------------------------------------------------------------
! 2. Technology Morkov Chain and Initial Guess 
!----------------------------------------------------------------

! Productivity value
vProductivity = (/0.9792, 0.9896, 1.0000, 1.0106, 1.0212/)
! Transition matrix
mTransition = reshape( (/0.9727, 0.0273, 0.0, 0.0, 0.0, &
0.0041, 0.9806, 0.0153, 0.0, 0.0, &
0.0, 0.0082, 0.9837, 0.0082, 0.0, &
0.0, 0.0, 0.0153, 0.9806, 0.0041, &
0.0, 0.0, 0.0, 0.0273, 0.9727 /), (/5,5/))

mTransition = transpose(mTransition)

nold=scheme(1)
allocate (vGridCapital(nold))
call sub_makegrid(lowcapital,highcapital,nold,curv,vGridCapital)
allocate (mValueFunction(nold,nGridProductivity))
allocate (mPolicyCapital(nold,nGridProductivity))
allocate (mPolicyLabour(nold,nGridProductivity))
! Initial guess for value function
if (initial_guess .eq. 1) then

    mValueFunction = 0.0
else if (initial_guess .eq. 2) then    
    mValueFunction = consumptionSteadyState/(1.0-bbeta)
else 
    mValueFunction = 0.0
    ! Deterministic version
    vProductivity = (/1.0000, 1.0000, 1.0000, 1.0000, 1.0000/)     
    
    nnew=scheme(1)
    nold=nnew
    call sub_makegrid(lowcapital,highcapital,nnew,curv,vGridCapital)
    call sub_RBC_labour_solver(nold,nnew,vGridCapital,mValueFunction,vGridCapital,mValueFunction,mPolicyCapital,mPolicyLabour)

    vProductivity = (/0.9792, 0.9896, 1.0000, 1.0106, 1.0212/)
    call sub_makegrid(lowcapital,highcapital,nnew,curv,vGridCapital)
end if    

!----------------------------------------------------------------
! Main Iteration - RBC
!----------------------------------------------------------------

! solve RBC either with brute force or with a numerical solver for the labour supply given capital choice
! call sub_RBC_labour ! Bruce force - second grid for labour 
! using zbrent from numerical recipes

do i=1,3
    nnew=scheme(i)
       
    allocate (aaa(nnew))
    allocate (bbb(nnew,nGridProductivity))
    deallocate (mPolicyCapital)
    deallocate (mPolicyLabour)
    allocate (mPolicyCapital(nnew,nGridProductivity))
    allocate (mPolicyLabour(nnew,nGridProductivity))
    
    call sub_RBC_labour_Solver(nold,nnew,vGridCapital,mValueFunction,aaa,bbb,mPolicyCapital,mPolicyLabour) 
    
    deallocate(vGridCapital)
    allocate (vGridCapital(nnew))
    vGridCapital=aaa
    deallocate(mValueFunction)
    allocate (mValueFunction(nnew,nGridProductivity))
    mValueFunction=bbb
    deallocate (aaa)
    deallocate (bbb)
end do

!----------------------------------------------------------------
! 5. PRODUCE OUTPUT FILE
!----------------------------------------------------------------

!call sub_out

!----------------------------------------------------------------
! 5. PRINT MAIN RESULTS
!----------------------------------------------------------------
  
print *, 'Final Iteration:', iteration, 'Sup Diff:', MaxDifference
print *, ' '
timeInSeconds=times(3)*60.0+times(4)
print *, 'Computing Time', timeInSeconds
print *, ' '
!
ss=nnew/2
print *, 'My check Capital:', mPolicyCapital(ss,3)
print *, ' '
print *, 'My check Labour:', mPolicyLabour(ss,3)
print *, ' '

print *, 'Press enter to continue'
pause
    
end program main_RBC_Labour