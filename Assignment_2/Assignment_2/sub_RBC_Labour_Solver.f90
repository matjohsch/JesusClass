subroutine sub_RBC_labour_solver
    
use RBC_var
use RBC_param
use RBC_functions
use sing_dim
use clock

implicit none

real(8), dimension(nGridCapital,nGridProductivity)::mCapitalProd
real(8) :: valueHighSoFar, valueProvisional, capitalChoice, labourChoice

! Numerical Solver
! Zbrent/brac
real(8),parameter::tolbre = 1.0e-12        ! tolerance level of zbrent
real(8),parameter::ddx    = 1.0e-4         ! range of first bracket in zbrac
real(8),parameter::errabs = 1.0e-6         ! absolute error level 
real(8),parameter::epsi = 1.0e-10         ! absolute error level 

real(8):: fsige,xa,xb
logical:: succes

do nProductivity = 1, nGridProductivity
    do nCapital = 1, nGridCapital
            mCapitalProd(nCapital, nProductivity) = vProductivity(nProductivity)*vGridCapital(nCapital)**aalpha
    end do
end do

maxDifference = 10.0
iteration = 0

call tic
do while (maxDifference>tolerance)
     
    expectedValueFunction = matmul(mValueFunction,transpose(mTransition));
   
    do nProductivity = 1,nGridProductivity               
         
        gridCapitalNextPeriod = 1
        do nCapital = 1,nGridCapital
            
            valueHighSoFar = -100000.0
            
            if (mod(iteration,10)==0 .OR. iteration==1) then 
                do nCapitalNextPeriod = gridCapitalNextPeriod , nGridCapital
                
                    xa = log(epsi)
                    xb = xa+log(ddx)

                    call  zbrac(f_labour,xa,xb,succes)
                    xb = zbrent(f_labour,xa,xb,tolbre)

                    fsige = f_labour(xb)
                    if (abs(fsige) .gt. 0.0001) then
                       print*, 'labour failed to compute supply'
                       pause
                    end if

                    labourChoice = exp(xb)
                    consumption = (1-ddelta)*vGridCapital(nCapital) + mCapitalProd(nCapital,nProductivity) *  labourChoice**(1.0-aalpha) - vGridCapital(nCapitalNextPeriod)                                     
                    
                    if (consumption .lt. 0.0) then
                        print*, 'negative consumption'
                    else
                        valueProvisional = log(consumption)-pphi*0.50* labourChoice**2.0+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity)
                    end if
                    if (valueProvisional>valueHighSoFar) then ! we break when we have achieved the max
                        valueHighSoFar = valueProvisional
                        !wer=nCapitalNextPeriod
                        !capitalChoice = vGridCapital(nCapitalNextPeriod)
                        gridCapitalNextPeriod = nCapitalNextPeriod
                    else
                        exit
                    end if
                end do    
                
                mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar
                
                mPolicyCapital(nCapital,nProductivity) = vGridCapital(wer)
                mPolicyLabour(nCapital,nProductivity) = labourChoice                 
            else 
                consumption = (1-ddelta)*vGridCapital(nCapital) + mCapitalProd(nCapital,nProductivity) *  mPolicyLabour(nCapital,nProductivity)**(1.0-aalpha) - mPolicyCapital(nCapital,nProductivity)                                     
                    
                mValueFunctionNew(nCapital,nProductivity) = log(consumption)-pphi*0.50* mPolicyLabour(nCapital,nProductivity)**2.0+bbeta*expectedValueFunction(mPolicyCapital(nCapital,nProductivity),nProductivity) 
            end if
     
        end do
    end do

    maxDifference = maxval((abs(mValueFunctionNew-mValueFunction)))
    mValueFunction = mValueFunctionNew
           
    iteration = iteration+1
    if (mod(iteration,10)==0 .OR. iteration==1) then
        print *, 'Iteration:', iteration, 'Sup Diff:', MaxDifference
    end if
end do
call toc

contains

!---------------------------------------------------------------------------------------
function f_labour(x) 

real(8),intent(in)::x
real(8)::xtransformed
real(8)::f_labour

xtransformed = exp(x)

f_labour = mCapitalProd(nCapital,nProductivity)*(xtransformed**(1.0-aalpha)-(1.0-aalpha)/pphi*xtransformed**(-aalpha-1.0))- vGridCapital(nCapitalNextPeriod) + (1.0-ddelta)*vgridCapital(nCapital)

end function f_labour
!---------------------------------------------------------------------------------------
    
end subroutine sub_RBC_labour_solver