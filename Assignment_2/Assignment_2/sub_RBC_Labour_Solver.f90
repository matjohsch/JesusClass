subroutine sub_RBC_Labour_Solver(nGridOld,nGridNew,GridOld,ValueOld,vGridCapital,mValueFunction,mPolicyCapital,mPolicyLabour, mPolicyConsumption)

use RBC_Parameter    
use RBC_Variable
use RBC_Functions
use Solver

implicit none

integer, intent(in):: nGridNew, nGridOld
real(8), dimension(nGridOld), intent(in) :: GridOld
real(8), dimension(nGridOld,nGridProductivity), intent(in) :: ValueOld

real(8), dimension(nGridNew), intent(out) :: vGridCapital
real(8), dimension(nGridNew,nGridProductivity), intent(out) :: mValueFunction

integer :: iteration, cCapital, cProductivity, cCapitalNextPeriod
integer :: gridCapitalNextPeriod
integer :: capitalChoice

real(8) :: valueHighSoFar, valueProvisional
real(8) :: labour, cons, labourChoice
real(8) :: maxDifference

real(8), dimension(nGridNew,nGridProductivity) :: mCapitalProd, mValueFunctionNew, expectedValueFunction
real(8), dimension(nGridNew,nGridProductivity) :: mPolicyCapital, mPolicyCapitalIndex, mPolicyLabour, mPolicyConsumption


! Generation of New Grid 
if (STOCHASTIC == 1) then
    call sub_makerandomgrid(lowcapital,highcapital,nGridNew,vGridCapital)
else
    call sub_makegrid(lowcapital,highcapital,nGridNew,curv,vGridCapital)
end if

! Interpolation of new Gridpoints at old Grid
do cProductivity = 1,nGridProductivity 
    do cCapital = 1,nGridNew
        call sub_interpolation(GridOld,ValueOld(:,cProductivity),nGridOld,vGridCapital(cCapital),mValueFunction(cCapital,cProductivity))
    end do
end do

! Capital State times Productivity State
do cProductivity = 1, nGridProductivity
    do cCapital = 1, nGridNew
        mCapitalProd(cCapital, cProductivity) = vProductivity(cProductivity)*vGridCapital(cCapital)**aalpha
    end do
end do

maxDifference = 10.0
iteration = 0

! Starting Value Function Iteration
do while (maxDifference>tolerance)
    
    iteration = iteration+1
    
    expectedValueFunction = matmul(mValueFunction,transpose(mTransition));

    do cProductivity = 1,nGridProductivity               
         
        gridCapitalNextPeriod = 1
        do cCapital = 1,nGridNew
                 
            if (HOWARD == 0 .OR. mod(iteration,10) == 0 .OR. iteration <= 3) then 
                
                valueHighSoFar = -100000.0
               
                do cCapitalNextPeriod = gridCapitalNextPeriod , nGridNew
                
                    xa = log(1.0)
                    xb = xa + log(ddx)

                    call zbrac(f_LabourOpt,xa,xb,succes)
                    xb = zbrent(f_LabourOpt,xa,xb,tolbre)

                    fsige = f_LabourOpt(xb)
                    if (abs(fsige) .gt. 0.000001) then
                       print *,'failed to compute labour supply'
                       pause
                    end if

                    labour = exp(xb)
                    cons = (1-ddelta) * vGridCapital(cCapital) + mCapitalProd(cCapital,cProductivity) *  labour**(1.0-aalpha) - vGridCapital(cCapitalNextPeriod)                                     
                    
                    if (cons < 0.0) then
                        print *,'negative consumption'
                    else
                        valueProvisional = log(cons) - pphi / 2.0 * labour**2.0 + bbeta * expectedValueFunction(cCapitalNextPeriod,cProductivity)
                    end if
                    if (valueProvisional > valueHighSoFar) then ! we break when we have achieved the max
                        valueHighSoFar = valueProvisional
                        capitalChoice = cCapitalNextPeriod
                        gridCapitalNextPeriod = cCapitalNextPeriod
                        labourchoice = labour
                    else
                        exit
                    end if
                end do    
                
                mValueFunctionNew(cCapital,cProductivity) = valueHighSoFar  

                mPolicyCapitalIndex(cCapital,cProductivity) = capitalChoice
                mPolicyCapital(cCapital,cProductivity) = vGridCapital(capitalchoice)
                
                mPolicyLabour(cCapital,cProductivity) = labourchoice     
                
                mPolicyConsumption(cCapital,cProductivity) = (1-ddelta)*vGridCapital(cCapital) + mCapitalProd(cCapital,cProductivity) *  mPolicyLabour(cCapital,cProductivity)**(1.0-aalpha) - mPolicyCapital(cCapital,cProductivity)                                     

            else 
                                    
                mValueFunctionNew(cCapital,cProductivity) = log(mPolicyConsumption(cCapital,cProductivity)) - pphi /2.0 * mPolicyLabour(cCapital,cProductivity)**2.0 + bbeta * expectedValueFunction(mPolicyCapitalIndex(cCapital,cProductivity),cProductivity)
            
            end if
     
        end do
    end do      
   
    if (HOWARD == 0 .OR. mod(iteration,10) == 0 .OR. iteration <= 3) then
        maxDifference = maxval((abs(mValueFunctionNew-mValueFunction)))
        if (mod(iteration,10) == 0 .OR. iteration == 1) then
            print *, 'Iteration:', iteration, 'Sup Diff:', MaxDifference
        end if
    end if

    mValueFunction = mValueFunctionNew 
    if (single) exit
end do

print *, 'Final Iteration:', iteration, 'Sup Diff:', MaxDifference
print *, ' '

contains

!---------------------------------------------------------------------------------------
function f_LabourOpt(x) 

real(8), intent(in) :: x
real(8) :: xtransformed
real(8) :: f_LabourOpt

xtransformed = exp(x)

f_LabourOpt = mCapitalProd(cCapital,cProductivity) * (xtransformed**(1.0-aalpha) - (1.0-aalpha) / pphi * xtransformed**(-aalpha-1.0)) - vGridCapital(cCapitalNextPeriod) + (1.0-ddelta) * vgridCapital(cCapital)

end function f_LabourOpt
!---------------------------------------------------------------------------------------
    
end subroutine sub_RBC_Labour_Solver