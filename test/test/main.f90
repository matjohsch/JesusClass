program test
    
    implicit none
    
    real(8)::a,b,c
    integer:: i
    b=3.0
    c=4.0
    
    a=b+c
    i=3
    call sub_subroutine(i)
    
    
    end program test