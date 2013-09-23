subroutine sub_RBC_Labour
    
use RBC_var
use RBC_param
use RBC_functions
use clock

implicit none

maxDifference = 10.0
iteration = 0

call tic
do while (maxDifference>tolerance)
     
    expectedValueFunction = matmul(mValueFunction,transpose(mTransition));
   
    do nProductivity = 1,nGridProductivity               
                
        do nCapital = 1,nGridCapital

            xxx=-1000000.0
            do nCapitalNextPeriod =1,nGridCapital
                do nLabour = 1,nGridLabour                   
                    
                    consumption = (1-ddelta)*vGridCapital(nCapital)+mOutput(nCapital,nLabour,nProductivity)-vGridCapital(nCapitalNextPeriod)                                     
                    
                    if (consumption .lt. 0.0) then
                        xxx(nLabour,nCapitalNextPeriod) = -100000.0
                    else
                        xxx(nLabour,nCapitalNextperiod) = log(consumption)-pphi*0.50*vGridLabour(nLabour)**2.0+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity)
                    end if
                    if (nlabour.gt.1) then
                        if (xxx(nLabour,nCapitalNextperiod) .lt. xxx(nLabour-1,nCapitalNextperiod)) then
                           exit
                        end if   
                    end if
                end do
            end do    
            
            resi=maxloc(xxx)
            
            mValueFunctionNew(nCapital,nProductivity) = xxx(resi(1),resi(2))
            mPolicyCapital(nCapital,nProductivity) = vGridCapital(resi(2))
            mPolicyLabour(nCapital,nProductivity) = vGridLabour(resi(1))           
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

end subroutine sub_RBC_Labour