subroutine sub_subroutine(i)
    
    implicit none
    integer,intent(in)::i
    real (8):: q,w,e
    real(8),dimension(i)::test
    
    test(1)=4
    test(2)=5
    test(3)=6
    e=10
    q=19
    w=23
    q=w*e
    
    end subroutine sub_subroutine