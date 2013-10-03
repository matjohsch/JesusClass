subroutine sub_RBC_EGM_Labour(nGridCapital, vGridCapitalNext, mValueFunction, mPolicyConsumption, mPolicyCapital, mPolicyLabour)
    
use RBC_Parameter    
use RBC_Variable
use RBC_Functions
use Solver
    
implicit none

integer, intent(in)::nGridCapital
real(8),dimension(nGridCapital),intent(out) :: vGridCapitalNext
real(8),dimension(nGridCapital,nGridProductivity),intent(out) :: mValueFunction, mPolicyConsumption, mPolicyCapital, mPolicyLabour

integer :: cCapital, cCapitalNext, cProd, iteration, iteration1

real(8) :: maxDifference, maxDifference1, derivative

real(8),dimension(nGridCapital,nGridProductivity) :: CapitalEndo, vGridCapitalEndo, LabourEndo, mValueFunctionold, Vtilde, Vtildeold, VEndogenous, VY


! Define Grids for Capital and prepare Guess for Valuefunction
lowcapital = 0.1
highcapital = 5.00
call sub_makegrid(lowcapital,highcapital,nGridCapital,curv,vGridCapitalNext)

!Step 1 Compute Value Function with EGM and SteadyState Labour Supply
labourSteady = .true.
call sub_RBC_EGM(nGridCapital, vGridCapitalEndo, mValueFunction)

! Step 2 - Recover PolicyFunction by VFI using the guess from EGM-SteadyState
VEndogenous = mValueFunction
do cProd = 1,nGridProductivity        
    do cCapitalNext = 1,nGridCapital
        call sub_interpolation(vGridCapitalEndo(:,cProd),VEndogenous(:,cProd),nGridCapital,vGridCapitalNext(cCapitalNext),mValueFunction(cCapitalNext,cProd))
    end do
end do
        
single = .true.
call sub_RBC_Labour_Solver(nGridCapital,nGridCapital,vGridCapitalNext,mValueFunction,vGridCapitalNext,mValueFunction,mPolicyCapital,mPolicyLabour,mPolicyConsumption)

! Step 3 - Using LabourPolicy Function from previous VFI to compute EGM 
maxDifference1 = 10.0
iteration1 = 0
do while (maxDifference1>tolerance)
    iteration1 = iteration1 + 1
    ! 3(a) -Find according CapitalEndo that leads to CapitalNextPeriod 
    do cProd = 1,nGridProductivity        
        do cCapitalNext = 1,nGridCapital
          call sub_interpolation(mPolicyCapital(:,cProd),vGridCapitalNext,nGridCapital,vGridCapitalNext(cCapitalNext),CapitalEndo(cCapitalNext,Cprod))
        end do
    end do
    
    !3(b) - Find the Labour Policy Function that belongs to CapitalEndo
    do cProd = 1,nGridProductivity        
        do cCapitalNext = 1,nGridCapital
          call sub_interpolation(vGridCapitalNext,mPolicyLabour(:,cProd),nGridCapital,CapitalEndo(cCapitalNext,Cprod),LabourEndo(cCapitalNext,Cprod))
        end do
    end do

    maxDifference = 10.0
    iteration = 0
    
    ! 3(c) Take Expectations of Value Function  
    Vtilde = bbeta*matmul(mValueFunction,transpose(mTransition))
    
    do while (maxDifference>tolerance .and. iteration.lt.100)
        iteration = iteration + 1
        
        ! 3(d) Compute given CapitalNext and LabourEndo optimal c and according CapitalToday
        do cProd = 1,nGridProductivity        
            do cCapitalNext = 1,nGridCapital
                call sub_derivative(vGridCapitalNext,Vtilde(:,cProd),nGridCapital,cCapitalNext,derivative)
                mPolicyConsumption(cCapitalNext,cProd) = 1.0/derivative
                
                ! Solve nonlinear function for CapitalToday
                xa = log(1.0)
                xb = xa+log(ddx)           
                call  zbrac(f_endolabour,xa,xb,succes)
                xb = zbrent(f_endolabour,xa,xb,tolbre)
                fsige = f_endolabour(xb)
                if (abs(fsige) .gt. 0.000001) then
                    print*, 'failed to compute capital out of cash on hand'
                    pause
                end if
                CapitalEndo(cCapitalNext,cProd) = exp(xb)
                
                ! 3(e) - Update ValueFunction ...
                VEndogenous(cCapitalNext,cProd) = log(mPolicyConsumption(cCapitalNext,cProd)) - pphi/2.0*LabourEndo(cCapitalNext,cProd)**2.0 + Vtilde(cCapitalNext,cProd)
            end do
        end do
        
        ! 3(e) ... and interpolate it on GridCapitalNext 
        do cProd = 1,nGridProductivity        
            do cCapitalNext = 1,nGridCapital
               call sub_interpolation(CapitalEndo(:,cProd),VEndogenous(:,cProd),nGridCapital, vGridCapitalNext(cCapitalNext),VY(cCapitalNext,cProd)) 
            end do 
        end do
        
        ! 3(f) - Find again the Labour Policy Function that belongs to the now new CapitalEndo
        do cProd = 1,nGridProductivity        
          do cCapitalNext = 1,nGridCapital
             call sub_interpolation(vGridCapitalNext,mPolicyLabour(:,cProd),nGridCapital,CapitalEndo(cCapitalNext,Cprod),LabourEndo(cCapitalNext,Cprod))
          end do
        end do
        
        ! 3(g) - Take expectations of Value Function
        Vtilde = bbeta*matmul(VY,transpose(mTransition))

        ! 3(h)- Compute convergence criterion
        maxDifference = maxval((abs(Vtildeold-Vtilde)))
        if (mod(iteration,10) == 0 .OR. iteration == 1) then        
            print *, 'Iteration:', iteration, 'Sup Diff:', MaxDifference
        end if
    
        Vtildeold = Vtilde     
    end do   
    
    ! back to step 2
    mValueFunction = VY
    call sub_RBC_Labour_Solver(nGridCapital,nGridCapital,vGridCapitalNext,mValueFunction,vGridCapitalNext,mValueFunction,mPolicyCapital,mPolicyLabour,mPolicyConsumption)
    
    maxDifference1 = maxval((abs(mValueFunctionold-mValueFunction)))
    print *, 'Iteration on Value Function:', iteration1, 'Sup Diff:', MaxDifference1
    print *,''
    mValueFunctionold=mValueFunction
end do

contains

!---------------------------------------------------------------------------------------
function f_endolabour(x) 

real(8),intent(in)::x
real(8)::xtransformed
real(8)::f_endolabour

xtransformed = exp(x)

f_endolabour = mPolicyConsumption(cCapitalNext,Cprod)+VgridCapitalNext(cCapitalNext)-vProductivity(cProd)*LabourEndo(cCapitalNext,Cprod)**(1.0-aalpha)*xtransformed**aalpha-(1.0-ddelta)*xtransformed

end function f_endolabour
!---------------------------------------------------------------------------------------

end subroutine sub_RBC_EGM_labour