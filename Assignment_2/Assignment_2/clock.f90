!##############################################################################
! MODULE clock
! 
! Module for gmeasuring elapsed time.
!
! copyright: Fabian Kindermann
!            University of Wuerzburg
!            kindermann.fabian@uni-wuerzburg.de
!##############################################################################

module clock



implicit none

save

!##############################################################################
! Declaration of Variables
!##############################################################################

! starting time for cpu timer
real*8, private :: starttime_cpu

! starting time real time timer
integer, private :: starttime_real(8)


!##############################################################################
! Subroutines and functions                                               
!##############################################################################

contains


!##############################################################################
! SUBROUTINE tic
! 
! Starts cpu timer.
!##############################################################################
subroutine tic()

    
    !##### ROUTINE CODE #######################################################
    
    ! get cpu time
    call cpu_time(starttime_cpu)

end subroutine tic


!##############################################################################
! SUBROUTINE toc
! 
! Stops cpu timer.
!##############################################################################
subroutine toc(file)

    use RBC_var

    !##### INPUT/OUTPUT VARIABLES #############################################
    
    ! optional file identifier
    integer, intent(in), optional :: file


    !##### OTHER VARIABLES ####################################################
    
    real*8 :: time
    integer :: outfile

    
    
    
    !##### ROUTINE CODE #######################################################
    
    ! get output file identifier
    if(present(file))then
        outfile = file
    else
        outfile = 0
    endif
    
    ! get cpu time
    call cpu_time(time)            
    
    ! calculate time difference
    time = time - starttime_cpu    
    
    ! get number of days
    times(1) = floor(time/(24d0*60d0*60d0))
    time = time - times(1)*24d0*60d0*60d0    
        
    ! get number of hours
    times(2) = floor(time/(60d0*60d0))
    time = time - times(2)*60d0*60d0
        
    ! get number of minutes
    times(3) = floor(time/60d0)
    time = time - times(3)*60d0
        
    ! get number of seconds
    times(4) = time
    
!    call outTime(times, outfile)

end subroutine toc


!##############################################################################
! SUBROUTINE outTime
! 
! Writes time to file.
!##############################################################################
subroutine outTime(times, file)

    !##### INPUT/OUTPUT VARIABLES #############################################
    
    ! time as integer array
    real*8, intent(in) :: times(4)
    
    ! the output file identifier
    integer, intent(in) :: file
    
    !##### OTHER VARIABLES ####################################################
    
    character(len=1000) :: output
    
    
    !##### ROUTINE CODE #######################################################
    
    ! set up output
    write(output, '(a)')'Time elapsed: '
    
    ! write time values
    if(times(1) > 0d0) write(output, '(a,1x,i3,a)') &
        trim(output), int(times(1)), ' d  '
    if(times(2) > 0d0) write(output, '(a,1x,i3,a)') &
        trim(output), int(times(2)), ' h  '
    if(times(3) > 0d0) write(output, '(a,1x,i3,a)') &
        trim(output), int(times(3)), ' min  '
    write(output, '(a,1x,f7.3,a)') &
        trim(output), times(4), ' s  '
    
    if(file > 0) then
        write(file, '(/a/)')trim(output)
    else 
        write(*, '(/a/)')trim(output)
    endif
    
end subroutine outTime

end module clock