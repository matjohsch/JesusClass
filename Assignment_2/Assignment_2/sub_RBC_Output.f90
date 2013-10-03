subroutine sub_Output(n,vGridCapital,mValueFunction,mPolicyConsumption,mPolicyCapital,mPolicyLabour)

use RBC_Variable

real(8), dimension(n) :: vGridCapital
real(8), dimension(n,nGridProductivity) :: mValueFunction, mPolicyCapital, mPolicyLabour

open (unit = 1, file = 'policyfunction.txt', status = 'old')
do cOut = 1,n
   write (1,'(16f30.15)') vGridcapital(cOut),&
                          mPolicyCapital(cOut,1),mPolicyCapital(cOut,2),mPolicyCapital(cOut,3),mPolicyCapital(cOut,4),mPolicyCapital(cOut,5), &
                          mPolicyLabour(cOut,1), mPolicyLabour(cOut,2), mPolicyLabour(cOut,3), mPolicyLabour(cOut,4), mPolicyLabour(cOut,5),  &
                          mValueFunction(cOut,1),mValueFunction(cOut,2),mValueFunction(cOut,3),mValueFunction(cOut,4),mValueFunction(cOut,5)                          
end do
close (1)
    
end subroutine sub_Output
