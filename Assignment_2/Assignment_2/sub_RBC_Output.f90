subroutine sub_Output(n,vGridCapital,mValueFunction,mPolicyCapital,mPolicyLabour)

use RBC_Variable

real(8), dimension(n) :: vGridCapital
real(8), dimension(n,nGridProductivity) :: mValueFunction, mPolicyCapital, mPolicyLabour

open (unit = 1, file = 'policyfunction.txt', status = 'old')
do gc = 1,n
   write (1,'(16f30.15)') vGridcapital(gc),&
                          mPolicyCapital(gc,1),mPolicyCapital(gc,2),mPolicyCapital(gc,3),mPolicyCapital(gc,4),mPolicyCapital(gc,5), &
                          mPolicyLabour(gc,1), mPolicyLabour(gc,2), mPolicyLabour(gc,3), mPolicyLabour(gc,4), mPolicyLabour(gc,5),  &
                          mValueFunction(gc,1),mValueFunction(gc,2),mValueFunction(gc,3),mValueFunction(gc,4),mValueFunction(gc,5)                          
end do
close (1)
    
end subroutine sub_Output
