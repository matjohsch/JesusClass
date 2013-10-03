!============================================================================
! Name : RBC_F90.f90
! Description : Basic RBC model with labour and partial depreciation
! Date : October 02, 2013
!============================================================================
    
program RBC_Labour_Main

use RBC_Parameter
use RBC_Variable
use RBC_Functions
use clock

implicit none

integer,dimension(:), allocatable:: nMultiGrid
integer:: nGridOld, nGridNew, cMultigrid, nGrids

real(8)::capital, labour, timeInSeconds
real(8), dimension(:), allocatable   :: vGridCapitalSave, vGridCapital
real(8), dimension(:,:), allocatable :: mValueFuctionSave, mValueFunction, mPolicyCapital, mPolicyLabour, mPolicyConsumption

integer:: SteadyState

!----------------------------------------------------------------
! 0. PreSettings
!----------------------------------------------------------------

! MULTIGRID SCHEME
! nGrids:= How many grids implemented; nMultiGrid:= Define grid size in each step
nGrids = 3
allocate (nMultiGrid(nGrids))
nMultiGrid = (/ 100, 1000, 200000 /)

! INITIAL GUESS for Value function
! 1:= constant function equal zeros
! 2:=constant function equal to Steady State Value
! 3:= value function obtained from the deterministic version of the problem
INITIAL_GUESS = 2

! ACCELERATOR
! To use the accelarator HOWARD:=1
HOWARD = 1

! STOCHASTIC GRID
! STOChASTIC:= 0 - Equally spaced grid
! STOCHASTIC:= 1 - randomly distributed gridpoints 
STOCHASTIC = 0

! ENDOGENOUS GRID
! ENDOGENOUS:= 0 - VFI
! ENDOGENOUS:= 1 - Endogenous Grid Method
ENDOGENOUS = 0

Single=.false.
LabourSteady=.false.

!----------------------------------------------------------------
! 1. Calibration
!----------------------------------------------------------------
! Steady state Labour Supply (Calibration Target), for other parameter see module RBC_PARAMETER
! Either choose utility parameter pphi or steady state labour supply
labourSteadyState = 0.3333333333      
x = (aalpha*bbeta)/(1.0-bbeta+bbeta*ddelta)
pphi = ((1.0-aalpha)/(labourSteadyState**2.0))/(1.0-ddelta*x)

!pphi=7.6175  

!----------------------------------------------------------------
! 2. Steady State
!----------------------------------------------------------------
call sub_SteadyState

!----------------------------------------------------------------
! 3. Productivity Markov Chain and Initial Guess 
!----------------------------------------------------------------

! Productivity States
vProductivity = (/0.9792, 0.9896, 1.0000, 1.0106, 1.0212/)
! Transition matrix
mTransition = reshape( (/0.9727, 0.0273, 0.0,    0.0,    0.0, &
                         0.0041, 0.9806, 0.0153, 0.0,    0.0, &
                         0.0,    0.0082, 0.9837, 0.0082, 0.0, &
                         0.0,    0.0,    0.0153, 0.9806, 0.0041, &
                         0.0,    0.0,    0.0,    0.0273, 0.9727 /), (/5,5/))

mTransition = transpose(mTransition)

! Define Initial Grid for Capital
nGridOld=nMultiGrid(1)
allocate (vGridCapital(nGridOld))

if (STOCHASTIC==1) then
    call sub_makerandomgrid(lowcapital,highcapital,nGridOld,vGridCapital)
else
    call sub_makegrid(lowcapital,highcapital,nGridOld,curv,vGridCapital)
end if

allocate (mValueFunction(nGridOld,nGridProductivity))
allocate (mPolicyCapital(nGridOld,nGridProductivity))
allocate (mPolicyLabour(nGridOld,nGridProductivity))
allocate (mPolicyConsumption(nGridOld,nGridProductivity))

! Initial Guess for Value Function
if (INITIAL_GUESS == 1) then   
    mValueFunction = 0.0   
else if (INITIAL_GUESS == 2) then        
    mValueFunction = (log(consumptionSteadyState)-pphi*0.5*labourSteadyState**2)/(1.0-bbeta)
else     
    ! Deterministic Version   
    vProductivity = (/1.0000, 1.0000, 1.0000, 1.0000, 1.0000/)     
    mValueFunction = 0.0  
    nGridOld = nMultiGrid(1)    
    if (STOCHASTIC == 1) then
        call sub_makerandomgrid(lowcapital,highcapital,nGridOld,vGridCapital)
    else
        call sub_makegrid(lowcapital,highcapital,nGridOld,curv,vGridCapital)
    end if
    call sub_RBC_labour_solver(nGridOld,nGridOld,vGridCapital,mValueFunction,vGridCapital,mValueFunction,mPolicyCapital,mPolicyLabour, mPolicyConsumption)

    vProductivity = (/0.9792, 0.9896, 1.0000, 1.0106, 1.0212/)
    if (STOCHASTIC == 1) then
        call sub_makerandomgrid(lowcapital,highcapital,nGridOld,vGridCapital)
    else
        call sub_makegrid(lowcapital,highcapital,nGridOld,curv,vGridCapital)
    end if
end if    

!----------------------------------------------------------------
! 4. Main Iteration 
!----------------------------------------------------------------

call tic
if (.not. ENDOGENOUS) then
    !----------------------------------------------------------------
    ! Value Function Iteration with a numerical solver (zbrent from numerical recipes) for the labour supply given capital choice
    !----------------------------------------------------------------
    do cMultigrid=1,size(nMultiGrid)
        nGridNew=nMultiGrid(cMultigrid)
       
        allocate (vGridCapitalSave(nGridNew))
        allocate (mValueFuctionSave(nGridNew,nGridProductivity))
    
        deallocate (mPolicyCapital)
        deallocate (mPolicyLabour)
        deallocate (mPolicyConsumption)
        allocate (mPolicyConsumption(nGridNew,nGridProductivity))
        allocate (mPolicyCapital(nGridNew,nGridProductivity))
        allocate (mPolicyLabour(nGridNew,nGridProductivity))
    
        call sub_RBC_labour_Solver(nGridOld,nGridNew,vGridCapital,mValueFunction,vGridCapitalSave,mValueFuctionSave,mPolicyCapital,mPolicyLabour,mPolicyConsumption) 
    
        deallocate(vGridCapital)
        allocate (vGridCapital(nGridNew))
        vGridCapital=vGridCapitalSave
    
        deallocate(mValueFunction)
        allocate (mValueFunction(nGridNew,nGridProductivity))
        mValueFunction=mValueFuctionSave
    
        deallocate (vGridCapitalSave)
        deallocate (mValueFuctionSave)
        nGridOld=nGridNew
    end do
else
    !----------------------------------------------------------------
    ! Endogenous Grid Method
    !----------------------------------------------------------------
    nGridNew=nMultigrid(3)
    deallocate(vGridCapital)
    deallocate (mPolicyCapital)
    deallocate (mPolicyLabour)
    deallocate (mPolicyConsumption)
    deallocate (mValueFunction)  
    
    allocate (vGridCapital(nMultigrid(3)))
    allocate (mPolicyConsumption(nMultigrid(3),nGridProductivity))
    allocate (mPolicyCapital(nMultigrid(3),nGridProductivity))
    allocate (mPolicyLabour(nMultigrid(3),nGridProductivity))
    allocate (mValueFunction(nMultigrid(3),nGridProductivity))
    call sub_RBC_EGM_labour(nMultigrid(3), vGridCapital, mValueFunction, mPolicyConsumption, mPolicyCapital, mPolicyLabour)
end if
call toc

!----------------------------------------------------------------
! 5. MEASURING ACCURACY BY EULER EQUATION ERROR
!----------------------------------------------------------------

call sub_EulerError(size(vGridCapital), vGridCapital, mPolicyLabour, mPolicyConsumption)

!----------------------------------------------------------------
! 6. PRODUCE OUTPUT FILE
!----------------------------------------------------------------

call sub_Output(size(vGridCapital),vGridCapital, mValueFunction, mPolicyConsumption, mPolicyCapital, mPolicyLabour)

!----------------------------------------------------------------
! 7. PRINT MAIN RESULTS
!----------------------------------------------------------------

timeInSeconds = times(3)*60.0+times(4)
print *, 'Computing Time', timeInSeconds
print *, ' '
!
if (STOCHASTIC == 1 .OR. ENDOGENOUS == 1 ) then
    call sub_interpolation(vGridCapital,mPolicyCapital(:,3),nGridNew,capitalSteadyState,capital)
    call sub_interpolation(vGridCapital,mPolicyLabour(:,3),nGridNew,capitalSteadyState,labour)
else
    SteadyState = nGridNew/2
    capital = mPolicyCapital(SteadyState,3)
    labour = mPolicyLabour(SteadyState,3)
end if

print *, 'My check Capital:', capital
print *, ' '
print *, 'My check Labour:', labour
print *, ' '

pause
    
end program RBC_Labour_Main