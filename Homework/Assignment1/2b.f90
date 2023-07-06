! Quantum Information and Computing 2022-2023
! Nguyen Xuan Tung, ID: 2005491
! Week 1, Exercise 2: Number precision
!   Integer and real numbers have a finite precision. Explore the limits of INTEGER
!   and REAL in Fortran.

!   (a) Sum 2.000.000 and 1 with INTEGER*2 and INTEGER*4
!   (b) Sum 𝜋 · 10^32 and √2 · 10^21 with single and double precision.
! --------------------------------------------------------------------------------------
program num_prec_b
    implicit none

!   (b) Sum 𝜋 · 10^32 and √2 · 10^21 with single and double precision.

! PI = 4.D0*DATAN(1.0d0)
! This style ensures that the maximum precision available on ANY architecture
! is used when assigning a value to PI. - StackOverflow
    integer*2 :: int2
    integer*4 :: int4
    real*4     :: real4
    real*8     :: real8

    integer*2, parameter :: int2_number   = 2000000
    integer*4, parameter :: int4_number   = 2000000
    real*4, parameter    :: single_pi     = 4.D0*DATAN(1.0d0)
    real*8, parameter    :: double_pi     = 4.D0*DATAN(1.0d0)
    real*4, parameter    :: single_number = sqrt(2.0e0)*1.e21
    real*8, parameter    :: double_number = sqrt(2.0d0)*1.d21

    print *, "Sum pi * 10^32 and sqtr(2) * 10^21 with single precision:", (single_pi*1.e32) + single_number
    print *, "Sum pi * 10^32 and sqrt(2) * 10^21 with double precision:", (double_pi*1.d32) + double_number



end program num_prec_b
! Error appears when compiling, in order to compile sucessful, a flag -fno-range-check must be added
! gfortran 2b.f90 -o 2b -fno-range-check
! ./2b
