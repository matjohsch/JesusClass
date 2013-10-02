subroutine sub_RBC_EGM(nGrid)
    
use RBC_Parameter    
use RBC_Variable
use RBC_Functions
use Solver
use Clock
    
implicit none

integer, intent(in)::nGrid

integer :: cCapitalNext, cProductivity , gc

real(8) :: maxDifference,derivative
real(8),dimension(nGrid) :: vGridCapitalNext
real(8),dimension(nGrid,nGridProductivity) :: mValueFunction,Vendo,VY,Vtilde,Vtildeold,GridY,mPolicyConsumption,vGridCapitalEndo,GridYendo

! Numerical Solver
! Zbrent/brac
real(8),parameter::tolbre = 1.0e-12        ! tolerance level of zbrent
real(8),parameter::ddx    = 1.0e-4         ! range of first bracket in zbrac
real(8),parameter::errabs = 1.0e-8         ! absolute error level 

real(8):: fsige,xa,xb
logical:: succes

integer:: cCapital, cProd
real(8):: cprime, kprime, rhs, cc, ll, auxiliar,average_eee,max_eee

real(8),dimension(nGrid)::eee,eee1

!vProductivity = (/1.0000, 1.0000, 1.0000, 1.0000, 1.0000/) 
vProductivity = (/0.9792, 0.9896, 1.0000, 1.0106, 1.0212/)

lowcapital = 0.1
highcapital = 6.00
call sub_makegrid(lowcapital,highcapital,nGrid,curv,vGridCapitalNext)

do cProd = 1,nGridProductivity
   do cCapitalNext =1, nGrid
      GridY(cCapitalNext,cProd)=vProductivity(cProd)* vGridCapitalNext(cCapitalNext)**aalpha+(1.0-ddelta)*vGridCapitalNext(cCapitalNext)
      Vtilde(cCapitalNext,cProd)=log(vGridCapitalNext(cCapitalNext))
   end do   
end do

maxDifference = 10.0
iteration = 0

call tic
do while (maxDifference>tolerance)
    iteration = iteration + 1
    do cProd = 1,nGridProductivity        
        do cCapitalNext = 1,nGrid
            call sub_derivative2(vGridCapitalNext,Vtilde(:,cProd),ngrid,cCapitalNext,derivative)
            mPolicyConsumption(cCapitalNext,cProd)=1.0/derivative
            GridYendo(cCapitalNext,cProd)=mPolicyConsumption(cCapitalNext,cProd)+vGridCapitalNext(cCapitalNext)
            Vendo(cCapitalNext,cProd)=log(mPolicyConsumption(cCapitalNext,cProd))+Vtilde(cCapitalNext,cProd)
        end do
    end do
    do cProd = 1,nGridProductivity        
        do cCapitalNext = 1,nGrid
           call sub_interpolation(GridYendo(:,cProd),Vendo(:,cProd),nGrid, GridY(cCapitalNext,cProd),VY(cCapitalNext,cProd)) 
        end do 
    end do
    Vtilde=bbeta*matmul(VY,transpose(mTransition))

    maxDifference = maxval((abs(Vtildeold-Vtilde)))
    if (mod(iteration,10) == 0 .OR. iteration == 1) then        
        print *, 'Iteration:', iteration, 'Sup Diff:', MaxDifference
    end if
    
    Vtildeold = Vtilde     
end do    

mValueFunction=Vendo

do cProductivity = 1,nGridProductivity   
    do cCapitalNext = 1,nGrid
        xa = log(1.0)
        xb = xa+log(ddx)

        call  zbrac(f_endo,xa,xb,succes)
        xb = zbrent(f_endo,xa,xb,tolbre)

        fsige = f_endo(xb)
        if (abs(fsige) .gt. 0.000001) then
            print*, 'failed to compute capital out of cash on hand'
            pause
        end if

        vGridCapitalendo(cCapitalNext,cProductivity) = exp(xb)
    end do
end do    
call toc

open (unit = 2, file = 'EGM.txt', status = 'old')
do gc = 1,nGrid
   write (2,'(3f30.15)') vGridCapitalEndo(gc,3),mValueFunction(gc,3),mPolicyConsumption(gc,3)                          
end do
close (2)
    
eee1=0.0000
do cCapital=1,nGrid
    kprime=vProductivity(3)*vGridCapitalEndo(cCapital,3)**aalpha+(1.000-ddelta)*vGridCapitalendo(cCapital,3)-mPolicyConsumption(cCapital,3) !*mPolicyLabour(cCapital,3)**(1.000-aalpha)
    rhs=0.00000

    do cProd=1,nGridProductivity
        call sub_interpolation (vGridCapitalendo,mPolicyConsumption(:,3),nGrid,kprime,cc)
        !call sub_interpolation (vGridCapitalendo,mPolicyLabour(:,3),nGridCapital,kprime,ll)
        auxiliar=(1.0000/cc)*((1.0000-ddelta)+vProductivity(cProd)*aalpha*kprime**(aalpha-1.0000))!*ll**(1.0000-aalpha))
        rhs=rhs+auxiliar*mtransition(3,cProd)
    end do    

    eee1(cCapital)=1.0000-mPolicyConsumption(cCapital,3)*bbeta*rhs!/cprime
    eee(cCapital)=log(abs(eee1(cCapital)))/log(10.000)
end do

average_eee=sum(eee)/real(nGrid)
max_eee=maxval(eee)
print *, ' '
print *, 'average Euler Error:', average_eee
print *, 'maximum Euler Error:', max_eee
print *, ' '

pause

contains

!---------------------------------------------------------------------------------------
function f_endo(x) 

real(8),intent(in)::x
real(8)::xtransformed
real(8)::f_endo

xtransformed = exp(x)
f_endo = gridYendo(cCapitalNext,cProductivity) - vProductivity(cProductivity)*xtransformed**aalpha-(1.0-ddelta)*xtransformed

end function f_endo
!---------------------------------------------------------------------------------------

end subroutine sub_RBC_EGM