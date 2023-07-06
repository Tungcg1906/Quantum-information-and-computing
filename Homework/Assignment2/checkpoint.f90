! Quantum Information and Computing 2022-2023
! Nguyen Xuan Tung, ID: 2005491
! Week 2, Exercise 1: Checkpoint
!    Write a subroutine to be used as a checkpoint for debugging

!   (a) Include a control on a logical variable (Debug=.TRUE. or .FALSE.).
!   (b) Include an additional (optional) string to be printed.
!   (c) Include additional (optional) variables to be printed.
! --------------------------------------------------------------------------------------


subroutine check_point(DEBUG, realarg, line)
    logical        :: DEBUG   ! input, if DEBUG_ == TRUE we are in 'debug mode'
    real           :: realarg ! optional generic real argument
    integer        :: line

    if (DEBUG .eqv. .TRUE.) then
        print *, 'LINE:',line   ! print file and line
        print*,'arg:', realarg
    end if
end subroutine check_point

#define check_real_(realarg) check_point(DEBUG, realarg,__LINE__)

program test
    logical :: DEBUG = .TRUE.
    real    :: x = 3.14159265359, y = 9.53562951413

! DEBUG_ is true, this value should be printed
    call check_real_(x)

    DEBUG = .FALSE.
! DEBUG_ is now false, this value should NOT be printed
    call check_real_(y)

end program test

! Error appears when compiling, in order to compile sucessful, a flag -cpp is used
! gfortran -o checkpoint checkpoint.f90 -cpp
! ./checkpoint
