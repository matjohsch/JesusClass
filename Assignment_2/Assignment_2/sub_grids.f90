subroutine sub_grids
    
use RBC_var
use RBC_param
use RBC_functions

implicit none

 !Transition matrix
mTransition = reshape( (/0.9727, 0.0273, 0.0, 0.0, 0.0, &
0.0041, 0.9806, 0.0153, 0.0, 0.0, &
0.0, 0.0082, 0.9837, 0.0082, 0.0, &
0.0, 0.0, 0.0153, 0.9806, 0.0041, &
0.0, 0.0, 0.0, 0.0273, 0.9727 /), (/5,5/))

mTransition = transpose(mTransition)

! Grid for capital
lowcapital=0.5*capitalSteadyState
highcapital =1.5*capitalSteadyState
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

end subroutine sub_grids