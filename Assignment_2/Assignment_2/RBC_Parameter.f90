module RBC_Parameter
    
implicit none
    
integer, parameter  :: nGridProductivity = 5  ! number of technology states
real (8), parameter :: curv = 1.0             ! curvature of grids
real (8), parameter :: tolerance = 1.0e-5     ! convergence tolerance level 
    
real (8), parameter :: aalpha = 0.3333333333  ! Elasticity of output w.r.t Capital
real (8), parameter :: bbeta  = 0.95          ! Discount factor
real (8), parameter :: ddelta = 0.09          ! Capital depreciation rate

! Numerical Solver
real (8), parameter::tolbre = 1.0e-12         ! tolerance level of zbrent
real (8), parameter::ddx    = 1.0e-4          ! range of first bracket in zbrac
real (8), parameter::errabs = 1.0e-8          ! absolute error level 
real (8), parameter :: epsi = 1.0e-12         ! absolute error level 

end module RBC_Parameter