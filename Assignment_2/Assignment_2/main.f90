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
initial_guess = 1

! Accelarator
! To use the accelarator Howard:=true
howard = 1

!----------------------------------------------------------------
! 1. Steady State
!----------------------------------------------------------------
call sub_steady_state

!----------------------------------------------------------------
! 2. Initial Guess and Grids 
!----------------------------------------------------------------

nGridCapital=178
allocate (vGridCapital(nGridCapital))
allocate (mValueFunction(nGridCapital,nGridProductivity))
allocate (mValueFunctionNew(nGridCapital,nGridProductivity))
allocate (mPolicyCapital(nGridCapital,nGridProductivity))
allocate (mPolicyLabour(nGridCapital,nGridProductivity))
allocate (mOutput(nGridCapital,nGridLabour,nGridProductivity))
allocate (expectedValueFunction(nGridCapital,nGridProductivity))
allocate (xxx(nGridLabour,nGridCapital))


! Productivity value
vProductivity = (/0.9792, 0.9896, 1.0000, 1.0106, 1.0212/)
call sub_grids

! Initial guess for value function
if (initial_guess .eq. 1) then
     mValueFunction = 0.0
else if (initial_guess .eq. 2) then    
    mValueFunction = consumptionSteadyState/(1.0-bbeta)
else 
    mValueFunction = 0.0
    ! Deterministic version
    vProductivity = (/1.0000, 1.0000, 1.0000, 1.0000, 1.0000/)     
    call sub_grids
    call sub_RBC_labour

    vProductivity = (/0.9792, 0.9896, 1.0000, 1.0106, 1.0212/)
    call sub_grids
end if    

!----------------------------------------------------------------
! Main Iteration - RBC
!----------------------------------------------------------------

! solve RBC either with brute force or with a numerical solver for the labour supply given capital choice
! call sub_RBC_labour ! Bruce force - second grid for labour 

call sub_RBC_labour_Solver ! using zbrent from numerical recipes
  
!----------------------------------------------------------------
! 5. PRODUCE OUTPUT FILE
!----------------------------------------------------------------

call sub_out

!----------------------------------------------------------------
! 5. PRINT MAIN RESULTS
!----------------------------------------------------------------
  
print *, 'Final Iteration:', iteration, 'Sup Diff:', MaxDifference
print *, ' '
timeInSeconds=times(3)*60.0+times(4)
print *, 'Computing Time', timeInSeconds
print *, ' '

ss=nGridCapital/2
print *, 'My check Capital:', mPolicyCapital(ss,3)
print *, ' '
print *, 'My check Labour:', mPolicyLabour(ss,3)
print *, ' '

print *, 'Press enter to continue'
pause
    
end program main_RBC_Labour