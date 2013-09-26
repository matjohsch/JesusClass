subroutine sub_SteadyState
    
use RBC_Variable
use RBC_Parameter

implicit none

!----------------------------------------------------------------
! 2. Steady state
!----------------------------------------------------------------

labourSteadyState=(((1.0-aalpha)/pphi)*(1.0-ddelta*x)**(-1.0))**(1.0/2.0)
capitalSteadyState = (x)**(1.0/(1.0-aalpha))*labourSteadyState
outputSteadyState = capitalSteadyState**(aalpha)*labourSteadyState**(1.0-aalpha)
consumptionSteadyState = outputSteadyState-ddelta*capitalSteadyState

! for later grid generation 
lowcapital=0.5*capitalSteadyState
highcapital =1.5*capitalSteadyState

print *, 'Deterministi Steady State Values'
print *, 'Output: ', outputSteadyState, 'Capital: ', capitalSteadyState
print *, 'Consumption: ', consumptionSteadyState, 'Labour: ', labourSteadyState
      
end subroutine sub_SteadyState