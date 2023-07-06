! Quantum Information and Computing 2022-2023
! Nguyen Xuan Tung, ID: 2005491
! Week 1, Exercise 1: Setup
!   (a) Create a working directory. 
!   (b) Open a code editor - i.e. emacs, Vim, Visual Studio editor, Atom - and write your first program in FORTRAN.
!   (c) Submit a test job.
!   (d) (Optional) Connect to CloudVeneto cluster via ssh and repeat the execution.

program Setup 
   implicit none

   character*20 :: none
   character (len=20) :: f_name, l_name
   print *, "What's your name?"

   read *, f_name, l_name
   print *, "Hello ", trim(f_name), " ", trim(l_name)


end program Setup



