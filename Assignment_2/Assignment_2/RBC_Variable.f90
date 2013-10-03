module RBC_Variable
  
use RBC_Parameter

implicit none

logical :: HOWARD
logical :: STOCHASTIC
logical :: ENDOGENOUS
integer :: INITIAL_GUESS

logical :: Single           
logical :: LabourSteady  

real (8) :: pphi, x, labourSteadyState
real (8) :: capitalSteadyState, outputSteadyState, consumptionSteadyState
real (8) :: lowcapital, highcapital

real (8), dimension(nGridProductivity) :: vProductivity
real (8), dimension(nGridProductivity,nGridProductivity) :: mTransition

real (8), dimension(4):: times

! Numerical Solver
real (8) :: fsige,xa,xb
logical  :: succes

end module RBC_Variable