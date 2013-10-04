subroutine sub_RBC_EGM(nGridCapital, mGridCapitalEndo, mValueFunction)
    
use RBC_Parameter    
use RBC_Variable
use RBC_Functions
use Solver
    
implicit none

integer, intent(in) :: nGridCapital

integer :: iteration, cCapitalNext, cProd, cOut

real (8) :: maxDifference, derivative
real (8), dimension(nGridCapital) :: vGridCapitalNext
real (8), dimension(nGridCapital,nGridProductivity) :: mValueFunction, mVEndogenous, mValueY, mValueTilde, mValueTildeold, mGridY, mPolicyConsumption, mGridCapitalEndo, mGridYendo

call sub_makegrid(lowcapital,highcapital,nGridCapital,curv,vGridCapitalNext)

if (labourSteady) then 
    do cProd = 1,nGridProductivity
       do cCapitalNext = 1,nGridCapital
          mGridY(cCapitalNext,cProd) = vProductivity(cProd) * vGridCapitalNext(cCapitalNext)**aalpha * labourSteadyState**(1.0 - aalpha) + (1.0 - ddelta) * vGridCapitalNext(cCapitalNext)
          mValueTilde(cCapitalNext,cProd) = log(vGridCapitalNext(cCapitalNext)) - pphi / 2.0 * labourSteadyState**2.0
       end do   
    end do
else
    do cProd = 1,nGridProductivity
       do cCapitalNext = 1,nGridCapital
          mGridY(cCapitalNext,cProd) = vProductivity(cProd) * vGridCapitalNext(cCapitalNext)**aalpha + (1.0 - ddelta) * vGridCapitalNext(cCapitalNext)
          mValueTilde(cCapitalNext,cProd) = log(vGridCapitalNext(cCapitalNext))
       end do   
    end do
end if

maxDifference = 10.0
iteration = 0

do while (maxDifference > tolerance)
    iteration = iteration + 1
    do cProd = 1,nGridProductivity        
        do cCapitalNext = 1,nGridCapital
            call sub_derivative(vGridCapitalNext,mValueTilde(:,cProd),nGridCapital,cCapitalNext,derivative)
            
            mPolicyConsumption(cCapitalNext,cProd) = 1.0 / derivative
            mGridYendo(cCapitalNext,cProd) = mPolicyConsumption(cCapitalNext,cProd) + vGridCapitalNext(cCapitalNext)
            
            if (labourSteady) then
                mVEndogenous(cCapitalNext,cProd) = log(mPolicyConsumption(cCapitalNext,cProd)) - pphi / 2.0 * labourSteadyState**2.0 + mValueTilde(cCapitalNext,cProd)
            else
                mVEndogenous(cCapitalNext,cProd) = log(mPolicyConsumption(cCapitalNext,cProd)) + mValueTilde(cCapitalNext,cProd)              
             end if                   
        end do
    end do
    do cProd = 1,nGridProductivity        
        do cCapitalNext = 1,nGridCapital
           call sub_interpolation(mGridYendo(:,cProd),mVEndogenous(:,cProd),nGridCapital, mGridY(cCapitalNext,cProd),mValueY(cCapitalNext,cProd)) 
        end do 
    end do
    mValueTilde = bbeta * matmul(mValueY,transpose(mTransition))

    maxDifference = maxval((abs(mValueTildeold-mValueTilde)))
    if (mod(iteration,10) == 0 .OR. iteration == 1) then        
        print *, 'Iteration:', iteration, 'Sup Diff:', MaxDifference
    end if
    
    mValueTildeold = mValueTilde     
end do    

mValueFunction = mVEndogenous

do cProd = 1,nGridProductivity   
    do cCapitalNext = 1,nGridCapital
        xa = log(1.0)
        xb = xa+log(ddx)

        call  zbrac(f_endo,xa,xb,succes)
        xb = zbrent(f_endo,xa,xb,tolbre)

        fsige = f_endo(xb)
        if (abs(fsige) .gt. 0.000001) then
            print*, 'failed to compute capital out of cash on hand'
            pause
        end if

        mGridCapitalEndo(cCapitalNext,cProd) = exp(xb)
    end do
end do    

contains

!---------------------------------------------------------------------------------------
function f_endo(x) 

real (8), intent(in) :: x
real (8) :: xtransformed
real (8) :: f_endo

xtransformed = exp(x)
if (labourSteady) then
    f_endo = mGridYendo(cCapitalNext,cProd) - vProductivity(cProd) * xtransformed**aalpha * labourSteadyState**(1.0-aalpha) - (1.0-ddelta) * xtransformed
else
    f_endo = mGridYendo(cCapitalNext,cProd) - vProductivity(cProd) * xtransformed**aalpha - (1.0-ddelta) * xtransformed
end if

end function f_endo
!---------------------------------------------------------------------------------------

end subroutine sub_RBC_EGM