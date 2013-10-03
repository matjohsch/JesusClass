subroutine sub_EulerError(nGridCapital, vGridCapital, mPolicyLabour, mPolicyConsumption)
    
use RBC_Parameter
use RBC_Variable
use RBC_Functions

implicit none

integer, intent(in) :: nGridCapital
real (8), dimension(nGridCapital), intent(in) :: vGridCapital
real (8), dimension(nGridCapital,nGridProductivity), intent(in) :: mPolicyLabour, mPolicyConsumption

integer  :: cCapital, cProd, cOut
real (8) :: cprime, kprime, rhs, cc, ll, auxiliar, average_eee, max_eee

real(8),dimension(nGridCapital)::eee, eee1

eee1 = 0.0
do cCapital = 1,nGridCapital
    kprime = vProductivity(3) * vGridCapital(cCapital)**aalpha * mPolicyLabour(cCapital,3)**(1.000-aalpha) + (1.000-ddelta) * vGridCapital(cCapital) - mPolicyConsumption(cCapital,3)

    rhs = 0.0
    do cProd = 1,nGridProductivity
        call sub_interpolation (vGridCapital,mPolicyConsumption(:,3),nGridCapital,kprime,cc)
        call sub_interpolation (vGridCapital,mPolicyLabour(:,3),nGridCapital,kprime,ll)
        auxiliar = (1.0 / cc) * ((1.0-ddelta) + vProductivity(cProd) * aalpha * kprime**(aalpha-1.0) * ll**(1.0-aalpha))
        rhs = rhs + auxiliar * mtransition(3,cProd)
    end do    

    eee1(cCapital) = 1.0 - mPolicyConsumption(cCapital,3) * bbeta * rhs
    eee(cCapital) = log(abs(eee1(cCapital))) / log(10.000)
end do

average_eee = sum(eee) / real(nGridCapital)
max_eee = maxval(eee)
print *, ' '
print *, 'Average Euler Error:', average_eee
print *, 'Maximum Euler Error:', max_eee
print *, ' '

open (unit = 2, file = 'Error.txt', status = 'old')
do cOut = 1,nGridCapital
   write (2,'(2f30.15)') vGridCapital(cOut), eee(cOut)                         
end do
close (2)

end subroutine sub_EulerError
    