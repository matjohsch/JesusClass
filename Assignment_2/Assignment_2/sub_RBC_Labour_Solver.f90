subroutine sub_RBC_labour_solver(nold,nnew,gridold,valueold,vGridCapital,mValueFunction,mPolicyCapital,mPolicyLabour)
    
use RBC_var
use RBC_param
use RBC_functions
use sing_dim
use clock

implicit none

integer,intent(in)::nnew,nold
real(8), dimension(nold),intent(in)::gridold
real(8),dimension(nold,nGridProductivity),intent(in)::valueold

real(8), dimension(nnew),intent(out)::vGridCapital
real(8),dimension(nnew,nGridProductivity),intent(out)::mValueFunction

real(8), dimension(nnew,nGridProductivity) :: mCapitalProd, mexpectedValueFunction
real(8), dimension(nnew,nGridProductivity) :: mValueFunctionNew, mPolicyCapital, mPolicyLabour,expectedValueFunction

real(8) :: value, valueHighSoFar, valueProvisional, labourChoice
integer:: capitalChoice, counterProd, counterCap

! Numerical Solver
! Zbrent/brac
real(8),parameter::tolbre = 1.0e-12        ! tolerance level of zbrent
real(8),parameter::ddx    = 1.0e-4         ! range of first bracket in zbrac
real(8),parameter::errabs = 1.0e-6         ! absolute error level 
real(8),parameter::epsi = 1.0e-10         ! absolute error level 

real(8):: fsige,xa,xb
logical:: succes

call sub_makegrid(lowcapital,highcapital,nnew,curv,vGridCapital)

! interpolation
do counterProd=1,nGridProductivity 
    do counterCap=1,nnew
        call sub_interpolation(gridold,valueold(:,counterProd),nold,vGridCapital(counterCap),value)
        mValueFunction(counterCap,counterProd)=value
    end do
end do

do nProductivity = 1, nGridProductivity
    do nCapital = 1, nnew
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
        do nCapital = 1,nnew
            
            valueHighSoFar = -100000.0
            
            if (howard==0 .OR. mod(iteration,10)==0 .OR. iteration==1 .or. iteration .gt.500 ) then 
                do nCapitalNextPeriod = gridCapitalNextPeriod , nnew
                
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
                        capitalChoice = nCapitalNextPeriod
                        gridCapitalNextPeriod = nCapitalNextPeriod
                    else
                        exit
                    end if
                end do    
                
                mValueFunctionNew(nCapital,nProductivity) = valueHighSoFar  
                mexpectedValueFunction(nCapital,nProductivity) = expectedValueFunction(capitalchoice,nProductivity)
                mPolicyCapital(nCapital,nProductivity) = vGridCapital(capitalchoice)
                mPolicyLabour(nCapital,nProductivity) = labourChoice                 
            else 
                consumption = (1-ddelta)*vGridCapital(nCapital) + mCapitalProd(nCapital,nProductivity) *  mPolicyLabour(nCapital,nProductivity)**(1.0-aalpha) - mPolicyCapital(nCapital,nProductivity)                                     
                    
                mValueFunctionNew(nCapital,nProductivity) = log(consumption)-pphi*0.50* mPolicyLabour(nCapital,nProductivity)**2.0+bbeta*mexpectedValueFunction(nCapital,nProductivity) 
            end if
     
        end do
    end do
           
   
    if (howard==0 .OR. mod(iteration,10)==0 .OR. iteration==1) then
        maxDifference = maxval((abs(mValueFunctionNew-mValueFunction)))
        if (mod(iteration,10)==0 .OR. iteration==1) then
            print *, 'Iteration:', iteration, 'Sup Diff:', MaxDifference
        end if
    end if
    iteration = iteration+1
    mValueFunction = mValueFunctionNew

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