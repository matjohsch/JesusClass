subroutine sub_RBC_EGM(nGrid,mValueFunction,mPolicyCapital,mPolicyLabour)
    
use RBC_Parameter    
use RBC_Variable
use RBC_Functions
use Solver
use Clock
    
implicit none

integer, intent(in)::nGrid

real(8),dimension(nGrid)::vGridCapitalNext
real(8),dimension(nGrid,nGridProductivity)::mValueFunction,derivative,GridY,mPolicyConsumption,mPolicyCapital,mPolicyLabour,Vtilde,gridyendo,vtildeold,valueVendo,valuev,kendogenous
real(8):: maxDifference,y

integer :: cCapitalNext, cProductivity

 !Numerical Solver
 !Zbrent/brac
real(8),parameter::tolbre = 1.0e-12        ! tolerance level of zbrent
real(8),parameter::ddx    = 1.0e-4         ! range of first bracket in zbrac
real(8),parameter::errabs = 1.0e-6         ! absolute error level 
real(8),parameter::epsi   = 1.0e-10         ! absolute error level 

real(8):: fsige,xa,xb
logical:: succes


vProductivity = (/1.0000, 1.0000, 1.0000, 1.0000, 1.0000/) 

lowcapital = 0.05
highcapital = 0.03
call sub_makegrid(lowcapital,highcapital,nGrid,curv,vGridCapitalNext)

do cProductivity = 1,nGridProductivity
   do cCapitalNext =1, nGrid
      GridY(cCapitalNext,cProductivity)=vProductivity(cProductivity)* vGridCapitalNext(cCapitalNext)**aalpha+(1.0-ddelta)*vGridCapitalNext(cCapitalNext)
   end do   
end do
mValueFunction=GridY-1.0
Vtilde=bbeta*matmul(mValueFunction,transpose(mTransition))

!lowcapital =0.9*lowcapital
!highcapital = 1.1*highcapital
!do cProductivity = 1,nGridProductivity
!    call sub_makegrid(lowcapital,highcapital,nGrid,curv,GridY(:,cproductivity))
!end do
!mValueFunction=(log(consumptionSteadyState)-pphi*0.5*labourSteadyState**2)/(1.0-bbeta)

maxDifference = 10.0
iteration = 0

mPolicyConsumption=0.0
GridYendo=0.0
ValueVendo=0.0
ValueV=0.0
Vtildeold=0.0

! Starting Value Function Iteration
call tic
do while (maxDifference>tolerance)
    iteration = iteration + 1
    do cProductivity = 1,nGridProductivity        
        do cCapitalNext = 1,nGrid
            
           !Y=vProductivity(cProductivity)* vGridCapitalNext(cCapitalNext)**aalpha+(1.0-ddelta)*vGridCapitalNext(cCapitalNext)
           !
           !call sub_derivative(GridY(:,cProductivity),Vtilde(:,cProductivity),nGrid,Y,derivative(cCapitalNext,cProductivity))
           call sub_derivative2(GridY(:,cProductivity),Vtilde(:,cProductivity),nGrid,CCapitalNext,derivative(cCapitalNext,cProductivity))
           
           mPolicyConsumption(cCapitalNext,cProductivity) = 1.0/derivative(cCapitalNext,cProductivity)
            
           GridYendo(cCapitalNext,cProductivity) = mPolicyConsumption(cCapitalNext,cProductivity) + vGridCapitalNext(cCapitalNext)

           ValueVendo(cCapitalNext,cProductivity)= log(mPolicyConsumption(cCapitalNext,cProductivity)) + Vtilde(cCapitalNext,cProductivity)
        end do
    end do
    do cProductivity = 1,nGridProductivity
        do cCapitalNext = 1,nGrid
           call sub_interpolation(GridYendo(:,cProductivity),ValueVendo(:,cProductivity),nGrid,GridY(cCapitalNext,cProductivity),ValueV(cCapitalNext,cProductivity))           
        end do
    end do
    Vtilde=bbeta*matmul(ValueV,transpose(mTransition)) 
    
    maxDifference = maxval((abs(Vtildeold-Vtilde)))
    if (mod(iteration,10) == 0 .OR. iteration == 1) then        
        print *, 'Iteration:', iteration, 'Sup Diff:', MaxDifference
    end if
    
    Vtildeold = Vtilde     
end do

do cProductivity = 1,nGridProductivity   
    do cCapitalNext = 1,nGrid
        xa = log(1.0)
        xb = xa+log(ddx)

        call  zbrac(f_endo,xa,xb,succes)
        xb = zbrent(f_endo,xa,xb,tolbre)

        fsige = f_endo(xb)
        if (abs(fsige) .gt. 0.000001) then
            print*, 'failed to compute labour supply'
            pause
        end if

        kendogenous(cCapitalNext,cProductivity) = exp(xb)
    end do
end do    
pause

contains

!---------------------------------------------------------------------------------------
function f_endo(x) 

real(8),intent(in)::x
real(8)::xtransformed,value
real(8)::f_endo

xtransformed = exp(x)
call sub_interpolation(vGridCapitalNext,GridYendo(:,cProductivity),nGrid,xtransformed,value)
f_endo = value - vProductivity(cProductivity)*xtransformed**aalpha-(1.0-ddelta)*xtransformed

end function f_endo
!---------------------------------------------------------------------------------------

end subroutine sub_RBC_EGM

   