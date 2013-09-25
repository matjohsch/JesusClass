!subroutine sub_grids(n,vGridCapital)
!    
!use RBC_var
!use RBC_param
!use RBC_functions
!
!implicit none
!
!integer,intent(in)::n
!real(8),dimension(n),intent(out)::vGridCapital
!
!
!
!
!
!!call sub_makegrid(n,curv,vGridCapital) 
! 
!!lowlabour=0.8*labourSteadyState
!!highlabour =1.2*labourSteadyState
!!vGridLabour(1:nGridLabour) = func_makegrid(lowlabour,highlabour,nGridLabour,curv) 
!
!
!!----------------------------------------------------------------
!! 3. Pre-build Output for each point in the grid
!!----------------------------------------------------------------
!  
!!do nProductivity = 1, nGridProductivity
!!    do nCapital = 1, nGridCapital
!!        do nLabour =1, nGridLabour
!!            mOutput(nCapital, nLabour, nProductivity) = vProductivity(nProductivity)*(vGridCapital(nCapital)**aalpha)*(vGridLabour(nLabour)**(1.0-aalpha))
!!        end do
!!    end do
!!end do
!
!end subroutine sub_grids