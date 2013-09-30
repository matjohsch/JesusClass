module RBC_Parameter
    
implicit none
    
integer, parameter  :: nGridProductivity = 5    ! number of technology states
real (8), parameter :: curv = 1.0               ! curvature of grids
real (8), parameter :: tolerance = 0.000001     ! convergence tolerance level 
    
real (8), parameter :: aalpha = 0.3333333333    ! Elasticity of output w.r.t Capital
real (8), parameter :: bbeta  = 0.95            ! Discount factor
real (8), parameter :: ddelta = 1.00            ! Capital depreciation rate

end module RBC_Parameter