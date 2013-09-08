!============================================================================
! Name : RBC_F90.f90
! Description : Basic RBC model with labour and partial depreciation
! Date : July 31, 2013
!============================================================================

program RBC_F90
  
!----------------------------------------------------------------
! 0. variables definition
!----------------------------------------------------------------
  
use  grid

implicit none

integer, parameter :: nGridCapital = 100
integer, parameter :: nGridLabour =50
integer, parameter :: nGridProductivity = 5
real (8), parameter :: curv=1.0
real (8), parameter :: tolerance = 0.0000001
  
integer :: nCapital, nCapitalNextPeriod, gridCapitalNextPeriod, nProductivity, nProductivityNextPeriod
integer :: nLabour
integer :: iteration
        
real :: elapsed(2), total

real(8) :: aalpha, bbeta, capitalSteadyState, labourSteadyState, outputSteadyState, consumptionSteadyState
real(8) :: lowcapital,highcapital,lowlabour,highlabour
real(8) :: valueHighSoFar1, valueProvisional1,valueHighSoFar2, valueProvisional2, consumption, capitalChoice	, labourChoice
real(8) :: maxDifference,diff,diffHighSoFar

real(8) :: pphi, ddelta,x

real(8), dimension(nGridProductivity) :: vProductivity
real(8), dimension(nGridProductivity,nGridProductivity) :: mTransition
real(8), dimension(nGridCapital) :: vGridCapital
real(8), dimension(nGridLabour) :: vGridLabour
real(8), dimension(nGridCapital,nGridProductivity) :: mValueFunction, mValueFunctionNew, mPolicyFunctionCapital, mPolicyFunctionLabour
real(8), dimension(nGridCapital,nGridLabour,nGridProductivity) :: mOutput
real(8), dimension(nGridCapital,nGridProductivity) :: expectedValueFunction
real(8), dimension(nGridLabour,nGridcapital) :: xxx
integer,dimension(2)::res 
!----------------------------------------------------------------
! 1. Calibration
!----------------------------------------------------------------
        
aalpha = 0.33333333333; ! Elasticity of output w.r.t
bbeta = 0.95 ! Discount factor
ddelta = 0.09
labourSteadyState=0.33333333333
x=(aalpha*bbeta)/(1.0-bbeta+bbeta*ddelta)
pphi = ((1.0-aalpha)/(labourSteadyState**2))/(1.0-ddelta*x) !6.72

! Productivity value
vProductivity = (/0.9792, 0.9896, 1.0000, 1.0106, 1.0212/)

! Transition matrix
mTransition = reshape( (/0.9727, 0.0273, 0., 0., 0., &
0.0041, 0.9806, 0.0153, 0.0, 0.0, &
0.0, 0.0082, 0.9837, 0.0082, 0.0, &
0.0, 0.0, 0.0153, 0.9806, 0.0041, &
0.0, 0.0, 0.0, 0.0273, 0.9727 /), (/5,5/))
   
!! Transition matrix
!mTransition = reshape( (/0.0, 0.0, 1.0, 0.0, 0.0, &
!0.0, 0.0, 1.0, 0.0, 0.0, &
!0.0, 0.0, 1.0, 0.0, 0.0, &
!0.0, 0.0, 1.0, 0.0, 0.0, &
!0.0, 0.0, 1.0, 0.0, 0.0 /), (/5,5/))
!mTransition = transpose(mTransition)

!----------------------------------------------------------------
! 2. Steady state
!----------------------------------------------------------------

capitalSteadyState = (x)**(1.0/(1.0-aalpha))*((1.0-aalpha)/pphi*1.0/(1.0-ddelta*x))**(-1/2)
outputSteadyState = capitalSteadyState**(aalpha)*labourSteadyState**(1.0-aalpha)
consumptionSteadyState = outputSteadyState-ddelta*capitalSteadyState
  
print *, 'Steady State values'
print *, 'Output: ', outputSteadyState, 'Capital: ', capitalSteadyState, 'Consumption: ', consumptionSteadyState
      
! Grid for capital
lowcapital=0.8*capitalSteadyState
highcapital =1.2*capitalSteadyState
vGridCapital(1:nGridCapital) = func_makegrid(lowcapital,highcapital,nGridCapital,curv) 
 
lowlabour=0.5*labourSteadyState
highlabour =1.5*labourSteadyState
vGridLabour(1:nGridLabour) = func_makegrid(lowlabour,highlabour,nGridLabour,curv) 

!----------------------------------------------------------------
! 3. Pre-build Output for each point in the grid
!----------------------------------------------------------------
  
do nProductivity = 1, nGridProductivity
    do nCapital = 1, nGridCapital
        do nLabour =1,nGridLabour
            mOutput(nCapital, nLabour, nProductivity) = vProductivity(nProductivity)*(vGridCapital(nCapital)**aalpha)*(vGridLabour(nLabour)**(1.0-aalpha))
        end do
    end do
end do

!----------------------------------------------------------------
! 4. Main Iteration
!----------------------------------------------------------------

maxDifference = 10.0
iteration = 0
  
do while (maxDifference>tolerance)
     
    expectedValueFunction = matmul(mValueFunction,transpose(mTransition));
   
    do nProductivity = 1,nGridProductivity
               

                
        do nCapital = 1,nGridCapital
                        
            valueHighSoFar1 = -1000000.0

            ! We start from previous choice (monotonicity of policy function)
            gridCapitalNextPeriod = 1  
            xxx=0.0
            do nCapitalNextPeriod = 1,nGridCapital !gridCapitalNextPeriod,nGridCapital
                valueHighSoFar2 = -100000.0
                do nLabour = 1,nGridLabour
                    
                    consumption = (1-ddelta)*vGridCapital(nCapital)+mOutput(nCapital,nLabour,nProductivity)-vGridCapital(nCapitalNextPeriod)
                    
                    xxx(nLabour,nCapitalNextperiod)  = log(consumption)-pphi*0.5*vGridLabour(nLabour)**2+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity) !(1.0-bbeta)*(
                    !valueProvisional2
                    if (consumption .lt. 0.0) then
                        xxx(nLabour,nCapitalNextPeriod) = -100000.0
                        valueProvisional2 = -100000.0
                    end if
                    !if (valueProvisional2>=valueHighSoFar2) then ! we break when we have achieved the max
                    !    valueHighSoFar2 = valueProvisional2
                    !else
                    !    x=2
                    !    exit
                    !end if
                end do
                !if (valueProvisional2>=valueHighSoFar1) then ! we break when we have achieved the max
                !    valueHighSoFar1 = valueProvisional2
                !    labourChoice = vGridLabour(nLabour)
                !    capitalChoice = vGridCapital(nCapitalNextPeriod)
                !    gridCapitalNextPeriod = nCapitalNextPeriod
                !else
                !    x=2
                !    exit
                !end if
            end do
            res=maxloc(xxx)
            
            mValueFunctionNew(nCapital,nProductivity) = xxx(res(1),res(2))
            mPolicyFunctionCapital(nCapital,nProductivity) = vGridCapital(res(2))
            mPolicyFunctionLabour(nCapital,nProductivity) = vGridLabour(res(1))
            
            !mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar1
            !mPolicyFunctionCapital(nCapital,nProductivity) = capitalChoice
            !mPolicyFunctionLabour(nCapital,nProductivity) = labourChoice
            x=2
        end do
    end do

    maxDifference = maxval((abs(mValueFunctionNew-mValueFunction)))
    mValueFunction = mValueFunctionNew
           
    iteration = iteration+1
    if (mod(iteration,10)==0 .OR. iteration==1) then
        print *, 'Iteration:', iteration, 'Sup Diff:', MaxDifference
    end if
end do
  
!----------------------------------------------------------------
! 5. PRINT RESULTS
!----------------------------------------------------------------
  
print *, 'Iteration:', iteration, 'Sup Diff:', MaxDifference
print *, ' '
print *, 'My check:', mPolicyFunctionCapital(50,3)
print *, ' '
print *, 'My check:', mPolicyFunctionLabour(50,3)
print *, ' '
!total = etime(elapsed)
!
!print *, 'Elapsed time is ', elapsed(1)
!pause
end program RBC_F90