subroutine sub_RBC_EGM(nGrid, vGridCapitalEndo, mValueFunction)
    
use RBC_Parameter    
use RBC_Variable
use RBC_Functions
use Solver
    
implicit none

integer, intent(in)::nGrid

integer :: iteration, cCapitalNext, cProd, gc

real(8) :: maxDifference,derivative
real(8),dimension(nGrid) :: vGridCapitalNext
real(8),dimension(nGrid,nGridProductivity) :: mValueFunction,VEndogenous,VY,Vtilde,Vtildeold,GridY,mPolicyConsumption,vGridCapitalEndo,GridYendo

call sub_makegrid(lowcapital,highcapital,nGrid,curv,vGridCapitalNext)

if (labourSteady) then 
    do cProd = 1,nGridProductivity
       do cCapitalNext =1, nGrid
          GridY(cCapitalNext,cProd)=vProductivity(cProd)* vGridCapitalNext(cCapitalNext)**aalpha*labourSteadyState**(1.0-aalpha)+(1.0-ddelta)*vGridCapitalNext(cCapitalNext)
          Vtilde(cCapitalNext,cProd)=log(vGridCapitalNext(cCapitalNext))-pphi/2.0*labourSteadyState**2.0
       end do   
    end do
else
    do cProd = 1,nGridProductivity
       do cCapitalNext =1, nGrid
          GridY(cCapitalNext,cProd)=vProductivity(cProd)* vGridCapitalNext(cCapitalNext)**aalpha+(1.0-ddelta)*vGridCapitalNext(cCapitalNext)
          Vtilde(cCapitalNext,cProd)=log(vGridCapitalNext(cCapitalNext))
       end do   
    end do
end if

maxDifference = 10.0
iteration = 0

do while (maxDifference>tolerance)
    iteration = iteration + 1
    do cProd = 1,nGridProductivity        
        do cCapitalNext = 1,nGrid
            call sub_derivative(vGridCapitalNext,Vtilde(:,cProd),ngrid,cCapitalNext,derivative)
            mPolicyConsumption(cCapitalNext,cProd)=1.0/derivative
            GridYendo(cCapitalNext,cProd)=mPolicyConsumption(cCapitalNext,cProd)+vGridCapitalNext(cCapitalNext)
            if (labourSteady) then
                VEndogenous(cCapitalNext,cProd)=log(mPolicyConsumption(cCapitalNext,cProd))-pphi/2.0*labourSteadyState**2.0+Vtilde(cCapitalNext,cProd)
            else
                VEndogenous(cCapitalNext,cProd)=log(mPolicyConsumption(cCapitalNext,cProd))+Vtilde(cCapitalNext,cProd)              
             end if            
        
        end do
    end do
    do cProd = 1,nGridProductivity        
        do cCapitalNext = 1,nGrid
           call sub_interpolation(GridYendo(:,cProd),VEndogenous(:,cProd),nGrid, GridY(cCapitalNext,cProd),VY(cCapitalNext,cProd)) 
        end do 
    end do
    Vtilde=bbeta*matmul(VY,transpose(mTransition))

    maxDifference = maxval((abs(Vtildeold-Vtilde)))
    if (mod(iteration,10) == 0 .OR. iteration == 1) then        
        print *, 'Iteration:', iteration, 'Sup Diff:', MaxDifference
    end if
    
    Vtildeold = Vtilde     
end do    

mValueFunction=VEndogenous

do cProd = 1,nGridProductivity   
    do cCapitalNext = 1,nGrid
        xa = log(1.0)
        xb = xa+log(ddx)

        call  zbrac(f_endo,xa,xb,succes)
        xb = zbrent(f_endo,xa,xb,tolbre)

        fsige = f_endo(xb)
        if (abs(fsige) .gt. 0.000001) then
            print*, 'failed to compute capital out of cash on hand'
            pause
        end if

        vGridCapitalendo(cCapitalNext,cProd) = exp(xb)
    end do
end do    

contains

!---------------------------------------------------------------------------------------
function f_endo(x) 

real(8),intent(in)::x
real(8)::xtransformed
real(8)::f_endo

xtransformed = exp(x)
if (labourSteady) then
    f_endo = gridYendo(cCapitalNext,cProd) - vProductivity(cProd)*xtransformed**aalpha*labourSteadyState**(1.0-aalpha)-(1.0-ddelta)*xtransformed
else
    f_endo = gridYendo(cCapitalNext,cProd) - vProductivity(cProd)*xtransformed**aalpha-(1.0-ddelta)*xtransformed
end if

end function f_endo
!---------------------------------------------------------------------------------------

end subroutine sub_RBC_EGM