!============================================================================
! Name : RBC_F90.f90
! Description : Basic RBC model with labour and partial depreciation
! Date : September 09, 2013
!============================================================================

program RBC_F90
  
!----------------------------------------------------------------
! 0. variables definition
!----------------------------------------------------------------
  
!use  grid

implicit none

integer, parameter :: nGridCapital = 100
integer, parameter :: nGridLabour =75
integer, parameter :: nGridProductivity = 5
real (8), parameter :: curv=1.0
real (8), parameter :: tolerance = 0.000001
  
integer :: nCapital, nLabour, nCapitalNextPeriod, gridCapitalNextPeriod, nProductivity, nProductivityNextPeriod
integer :: iteration
        
!real :: elapsed(2), total

real(8) :: aalpha, bbeta, capitalSteadyState, labourSteadyState, outputSteadyState, consumptionSteadyState
real(8) :: lowcapital,highcapital,lowlabour,highlabour
real(8) :: consumption
real(8) :: maxDifference,diff
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
        
aalpha = 0.3333333333          ! Elasticity of output w.r.t Capital
bbeta = 0.95       ! Discount factor
ddelta = 0.09      ! Capital depreciation rate

labourSteadyState=0.3333333333 ! Steady state Labour Supply (Calibration Target)
x=(aalpha*bbeta)/(1.0-bbeta+bbeta*ddelta)
pphi = ((1.0-aalpha)/(labourSteadyState**2.0))/(1.0-ddelta*x) !6.72    !Leisure Utility parameter

! Productivity value
vProductivity = (/0.9792, 0.9896, 1.0000, 1.0106, 1.0212/)

 !Transition matrix
mTransition = reshape( (/0.9727, 0.0273, 0.0, 0.0, 0.0, &
0.0041, 0.9806, 0.0153, 0.0, 0.0, &
0.0, 0.0082, 0.9837, 0.0082, 0.0, &
0.0, 0.0, 0.0153, 0.9806, 0.0041, &
0.0, 0.0, 0.0, 0.0273, 0.9727 /), (/5,5/))

mTransition = transpose(mTransition)

!----------------------------------------------------------------
! 2. Steady state
!----------------------------------------------------------------
labourSteadyState=(((1.0-aalpha)/pphi)*(1.0-ddelta*x)**(-1.0))**(1.0/2.0)
capitalSteadyState = (x)**(1.0/(1.0-aalpha))*labourSteadyState
outputSteadyState = capitalSteadyState**(aalpha)*labourSteadyState**(1.0-aalpha)
consumptionSteadyState = outputSteadyState-ddelta*capitalSteadyState
  
print *, 'Steady State values'
print *, 'Output: ', outputSteadyState, 'Capital: ', capitalSteadyState, 'Consumption: ', consumptionSteadyState
      
! Grid for capital
lowcapital=0.8*capitalSteadyState
highcapital =1.2*capitalSteadyState
vGridCapital(1:nGridCapital) = func_makegrid(lowcapital,highcapital,nGridCapital,curv) 
 
lowlabour=0.8*labourSteadyState
highlabour =1.2*labourSteadyState
vGridLabour(1:nGridLabour) = func_makegrid(lowlabour,highlabour,nGridLabour,curv) 

!----------------------------------------------------------------
! 3. Pre-build Output for each point in the grid
!----------------------------------------------------------------
  
do nProductivity = 1, nGridProductivity
    do nCapital = 1, nGridCapital
        do nLabour =1, nGridLabour
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

            xxx=0.0
            do nCapitalNextPeriod =1,nGridCapital
                do nLabour = 1,nGridLabour                   
                    
                    consumption = (1-ddelta)*vGridCapital(nCapital)+mOutput(nCapital,nLabour,nProductivity)-vGridCapital(nCapitalNextPeriod)                                     
                    
                    if (consumption .lt. 0.0) then
                        xxx(nLabour,nCapitalNextPeriod) = -100000.0
                    else
xxx(nLabour,nCapitalNextperiod) = log(consumption)-pphi*0.50*vGridLabour(nLabour)**2.0+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity)
                    end if
                end do
            end do    
            
            res=maxloc(xxx)
            
            mValueFunctionNew(nCapital,nProductivity) = xxx(res(1),res(2))
            mPolicyFunctionCapital(nCapital,nProductivity) = vGridCapital(res(2))
            mPolicyFunctionLabour(nCapital,nProductivity) = vGridLabour(res(1))           
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

pause

    contains 
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
end program RBC_F90