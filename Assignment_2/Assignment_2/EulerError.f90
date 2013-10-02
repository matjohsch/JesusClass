subroutine sub_EulerError(nGridCapital, vGridCapital, mPolicyLabour, mPolicyConsumption)
    
use RBC_Parameter
use RBC_Variable
use RBC_Functions

implicit none

integer, intent(in):: nGridCapital
real(8),dimension(nGridCapital), intent(in):: vGridCapital
real(8),dimension(nGridCapital,nGridProductivity), intent(in):: mPolicyLabour, mPolicyConsumption

integer:: cCapital, cProd, gc
real(8):: cprime, kprime, rhs, cc, ll, auxiliar,average_eee,max_eee

real(8),dimension(nGridCapital)::eee,eee1

eee1=0.0000
do cCapital=1,nGridCapital
    kprime=vProductivity(3)*vGridCapital(cCapital)**aalpha*mPolicyLabour(cCapital,3)**(1.000-aalpha)+(1.000-ddelta)*vGridCapital(cCapital)-mPolicyConsumption(cCapital,3)
    rhs=0.00000

    do cProd=1,nGridProductivity
        call sub_interpolation (vGridCapital,mPolicyConsumption(:,3),nGridCapital,kprime,cc)
        call sub_interpolation (vGridCapital,mPolicyLabour(:,3),nGridCapital,kprime,ll)
        auxiliar=(1.0000/cc)*((1.0000-ddelta)+vProductivity(cProd)*aalpha*kprime**(aalpha-1.0000)*ll**(1.0000-aalpha))
        rhs=rhs+auxiliar*mtransition(3,cProd)
    end do    

    eee1(cCapital)=1.0000-mPolicyConsumption(cCapital,3)*bbeta*rhs!/cprime
    eee(cCapital)=log(abs(eee1(cCapital)))/log(10.000)
end do

average_eee=sum(eee)/real(nGridCapital)
max_eee=maxval(eee)
print *, ' '
print *, 'average Euler Error:', average_eee
print *, 'maximum Euler Error:', max_eee
print *, ' '

open (unit = 6, file = 'Error.txt', status = 'old')
do gc = 1,nGridCapital
   write (6,'(2f30.15)') vGridCapital(gc), eee(gc)                         
end do
close (6)
    

end subroutine sub_EulerError
    