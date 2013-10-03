subroutine sub_RBC_SteadyState
    
use RBC_Variable
use RBC_Parameter

implicit none

labourSteadyState  = (((1.0-aalpha) / pphi) * (1.0 - ddelta * x)**(-1.0))**(1.0/2.0)
capitalSteadyState = (x)**(1.0 / (1.0-aalpha)) * labourSteadyState
outputSteadyState  = capitalSteadyState**(aalpha) * labourSteadyState**(1.0-aalpha)
consumptionSteadyState = outputSteadyState - ddelta * capitalSteadyState

! VFI grid generation 
lowcapital = 0.5 * capitalSteadyState
highcapital = 1.5 * capitalSteadyState

print *, 'Deterministic Steady State Values'
print *, 'Output:      ', outputSteadyState,      'Capital: ', capitalSteadyState
print *, 'Consumption: ', consumptionSteadyState, 'Labour:  ', labourSteadyState
print *, ''      
end subroutine sub_RBC_SteadyState