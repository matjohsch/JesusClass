module RBC_Parameter
    
implicit none
    
integer, parameter  :: nGridProductivity = 5    ! number of technology states
real (8), parameter :: curv = 1.0               ! curvature of grids
real (8), parameter :: tolerance = 1.0e-8     ! convergence tolerance level 
    
real (8), parameter :: aalpha = 0.3333333333    ! Elasticity of output w.r.t Capital
real (8), parameter :: bbeta  = 0.95            ! Discount factor
real (8), parameter :: ddelta = 0.09           ! Capital depreciation rate
real(8),parameter   :: epsi = 1.0e-12         ! absolute error level 

end module RBC_Parameter