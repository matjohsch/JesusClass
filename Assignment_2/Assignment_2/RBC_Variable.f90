module RBC_Variable
  
use RBC_Parameter

implicit none

integer :: iteration

real(8), dimension(nGridProductivity) :: vProductivity
real(8), dimension(nGridProductivity,nGridProductivity) :: mTransition

real(8) :: pphi, x, labourSteadyState
real(8) :: capitalSteadyState, outputSteadyState, consumptionSteadyState
real(8) :: lowcapital, highcapital

integer :: INITIAL_GUESS
logical :: HOWARD
logical :: STOCHASTIC
logical :: Single ! just one iteration in VFI
logical :: LabourSteady
real(8), dimension(4):: times

! Numerical Solver
! Zbrent/brac

real(8):: fsige,xa,xb
logical:: succes

end module RBC_Variable