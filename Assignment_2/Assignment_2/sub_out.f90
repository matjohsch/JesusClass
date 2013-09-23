subroutine sub_out

use RBC_var
    
open (unit = 1, file = 'policyfunction.txt', status = 'old')
do gc = 1,nGridcapital
   write (1,'(16f30.15)') mpolicyCapital(gc,1),mpolicyCapital(gc,2),mpolicyCapital(gc,3),mpolicyCapital(gc,4),mpolicyCapital(gc,5), &
                          mpolicyLabour(gc,1), mpolicyLabour(gc,2), mpolicyLabour(gc,3), mpolicyLabour(gc,4), mpolicyLabour(gc,5),  &
                          mValueFunction(gc,1),mValueFunction(gc,2),mValueFunction(gc,3),mValueFunction(gc,4),mValueFunction(gc,5), &
                          vGridcapital(gc)
end do
close (1)
    
    
end subroutine sub_out
