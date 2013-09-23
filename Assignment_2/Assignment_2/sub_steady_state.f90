subroutine sub_steady_state
    
use RBC_var
use RBC_param

implicit none

!----------------------------------------------------------------
! 2. Steady state
!----------------------------------------------------------------

labourSteadyState=(((1.0-aalpha)/pphi)*(1.0-ddelta*x)**(-1.0))**(1.0/2.0)
capitalSteadyState = (x)**(1.0/(1.0-aalpha))*labourSteadyState
outputSteadyState = capitalSteadyState**(aalpha)*labourSteadyState**(1.0-aalpha)
consumptionSteadyState = outputSteadyState-ddelta*capitalSteadyState
  
print *, 'Steady State values'
print *, 'Output: ', outputSteadyState, 'Capital: ', capitalSteadyState
print *, 'Consumption: ', consumptionSteadyState, 'Labour: ', labourSteadyState
      
end subroutine sub_steady_state