module RBC_var
  
use RBC_param

implicit none

integer :: nCapital, nLabour, nCapitalNextPeriod, gridCapitalNextPeriod, nProductivity, nProductivityNextPeriod
integer :: iteration, k
        
!real :: elapsed(2), total

real(8) :: pphi, x, labourSteadyState
real(8) :: capitalSteadyState, outputSteadyState, consumptionSteadyState
real(8) :: lowcapital, highcapital, lowlabour, highlabour
real(8) :: consumption
real(8) :: maxDifference, diff, ss


real(8), dimension(nGridProductivity) :: vProductivity
real(8), dimension(nGridProductivity,nGridProductivity) :: mTransition
real(8), dimension(nGridCapital) :: vGridCapital
real(8), dimension(nGridLabour) :: vGridLabour
real(8), dimension(nGridCapital,nGridProductivity) :: mValueFunction, mValueFunctionNew, mPolicyCapital, mPolicyLabour
real(8), dimension(nGridCapital,nGridLabour,nGridProductivity) :: mOutput
real(8), dimension(nGridCapital,nGridProductivity) :: expectedValueFunction
real(8), dimension(nGridLabour,nGridcapital) :: xxx
integer,dimension(2)::resi

integer initial_guess

real(8), dimension(4):: times
real(8)::timeInSeconds

end module RBC_var