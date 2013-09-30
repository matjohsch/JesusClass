!subroutine sub_RBC_EGM_solver(nGrid,mValueFunction,mPolicyCapital,mPolicyLabour)
!    
!use RBC_Parameter    
!use RBC_Variable
!use RBC_Functions
!use Solver
!use clock
!    
!implicit none
!
!integer, intent(in)::nGrid
!
!real(8),dimension(nGrid)::vGridCapitalNext
!real(8),dimension(nGrid,nGridProductivity)::mValueFunction,GridY,mPolicyConsumption,mPolicyCapital,mPolicyLabour,ValueDerivative,ValueEndogenous,ValueDerivativeold
!real(8):: rhsc,cons, labour,maxDifference,q,value, capital,expectedvalue
!
!integer :: cCapital, cProductivity, cCapitalNextPeriod,cProductivitytomorrow
!
!! Numerical Solver
!! Zbrent/brac
!real(8),parameter::tolbre = 1.0e-12        ! tolerance level of zbrent
!real(8),parameter::ddx    = 1.0e-4         ! range of first bracket in zbrac
!real(8),parameter::errabs = 1.0e-6         ! absolute error level 
!real(8),parameter::epsi   = 1.0e-10         ! absolute error level 
!
!real(8):: fsige,xa,xb
!logical:: succes
!
!
!
!
!lowcapital=1.2*lowcapital
!highcapital=0.8*highcapital
!call sub_makegrid(lowcapital,highcapital,nGrid,curv,vGridCapitalNext)
!
!do cProductivity = 1,nGridProductivity 
!   do cCapitalNext =1, nGrid
!      GridY 
!   end do
!   
!end do
!
!!ValueDerivativeold=0.0
!!ValueDerivative=((1.0-ddelta)+aalpha* capitalSteadyState**aalpha*labourSteadyState**(1.0-aalpha))/consumptionSteadyState       
!!mValueFunction=(log(consumptionSteadyState)-pphi*0.5*labourSteadyState**2)/(1.0-bbeta)
!
!
!maxDifference = 10.0
!iteration = 0
!
!! Starting Value Function Iteration
!call tic
!do while (maxDifference>tolerance)
!    iteration = iteration+1
!    do cProductivity = 1,nGridProductivity        
!        do cCapital = 1,nGrid
!            rhsc=0.0
!            do cProductivitytomorrow = 1,nGridProductivity 
!            
!                call sub_derivative(GridEndogenous(:,cProductivitytomorrow),mValuefunction(:,cProductivitytomorrow),nGrid,vGridCapitalNext(cCapital),Value)
!
!                rhsc = rhsc + mtransition(cProductivity,cProductivitytomorrow)*Value
!                
!            end do 
!            mPolicyConsumption(cCapital,cProductivity) = 1.0/(bbeta**rhsc)
!            q =pphi*mPolicyConsumption(cCapital,cProductivity)/(vProductivity(cProductivity)*(1.0-aalpha))
!            xa = log(1.0)
!            xb = xa+log(ddx)
!
!            call  zbrac(f_labour,xa,xb,succes)
!            xb = zbrent(f_labour,xa,xb,tolbre)
!
!            fsige = f_labour(xb)
!            if (abs(fsige) .gt. 0.000001) then
!                print*, 'failed to compute labour supply'
!                pause
!            end if
!
!            mPolicyLabour(cCapital,cProductivity) = exp(xb)
!            
!            
!        end do
!    end do
!    
!    ! upadate guess
!    do cProductivity = 1,nGridProductivity        
!        do cCapital  = 1,nGrid
!            
!            GridEndogenous(cCapital,cProductivity) = q**(1.0/aalpha)*mPolicyLabour(cCapital,cProductivity)**((1.0+aalpha)/aalpha)
!
!            ValueDerivative(cCapital,cProductivity)=((1.0-ddelta)+aalpha*vProductivity(cProductivity)* GridEndogenous(cCapital,cProductivity)**aalpha*mPolicyLabour(cCapital,cProductivity)**(1.0-aalpha))/mPolicyConsumption(cCapital,cProductivity)        
!            expectedvalue=0.0
!            
!            do cProductivitytomorrow = 1,nGridProductivity 
!            
!                call sub_interpolation(GridEndogenous,mValueFunction(:,cProductivitytomorrow),nGrid,vGridCapitalNext(cCapital),Value)
!            
!                expectedvalue = expectedvalue + mtransition(cProductivity,cProductivitytomorrow)*Value
!                
!            end do 
!            
!            mValueFunction(cCapital,cProductivity)=(log(mPolicyConsumption(cCapital,cProductivity))-pphi*0.5*mPolicyLabour(cCapital,cProductivity)**2)+bbeta*expectedValue
!        end do
!    end do
!    maxDifference = maxval((abs(ValueDerivative-ValueDerivativeold)))
!    if (mod(iteration,10) == 0 .OR. iteration == 1) then
!        
!
!        print *, 'Iteration:', iteration, 'Sup Diff:', MaxDifference
!
!    end if
!    
!    ValueDerivativeold = ValueDerivative     
!end do
!
!pause
!
!contains
!
!!---------------------------------------------------------------------------------------
!function f_labour(x) 
!
!real(8),intent(in)::x
!real(8)::xtransformed
!real(8)::f_labour
!
!xtransformed = exp(x)
!
!f_labour = mPolicyConsumption(cCapital,cProductivity)+vGridCapitalNext(cCapital)-(1.0-ddelta)*q**(1.0/aalpha)*xtransformed**((1.0+aalpha)/aalpha)-vProductivity(cProductivity)*q*xtransformed**2.0
!
!end function f_labour
!!---------------------------------------------------------------------------------------
!
!end subroutine sub_RBC_EGM_solver