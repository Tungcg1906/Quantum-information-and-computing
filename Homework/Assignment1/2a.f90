! Quantum Information and Computing 2022-2023
! Nguyen Xuan Tung, ID: 2005491
! Week 1, Exercise 2: Number precision
!   Integer and real numbers have a finite precision. Explore the limits of INTEGER
!   and REAL in Fortran.

!   (a) Sum 2.000.000 and 1 with INTEGER*2 and INTEGER*4
!   (b) Sum 𝜋 · 10^32 and √2 · 10^21 with single and double precision.
! --------------------------------------------------------------------------------------
program num_prec_a
    implicit none

!   (a) Sum 2.000.000 and 1 with INTEGER*2 and INTEGER*4

    integer*2 :: x_2, y_2   
    integer*4 :: x_4, y_4

    x_2   = 2000000
    y_2   = 1

    x_4   = 2000000
    y_4   = 1

    print *, "Sum of 2.000.000 and 1 with INTEGER*2 using 2 bytes:", x_2 + y_2
    print *, "The result of Sum of 2.000.000 and 1 with INTEGER*2 using 2 bytes  is wrong!!!"
    print *, "Sum of 2.000.000 and 1 with INTEGER*2 using 4 bytes:", x_4 + y_4
    print *, "The result of Sum of 2.000.000 and 1 with INTEGER*2 using 4 bytes  is correct!!!"

end program num_prec_a
! Error appears when compiling, in order to compile sucessful, a flag -fno-range-check must be added
! gfortran 2a.f90 -o 2a -fno-range-check
! ./2a