program test_limits
  implicit none

  integer :: N, d
  double complex, dimension(:,:), allocatable :: rho
  integer*16     :: nbytes

  print*, "--------------------------------------------"
  print*, "+              ALLOCATE LIMITS             +"
  print*, "+------------------------------------------+"
  write(*,"(A)",advance='no') " + N, d: "
  read (*,*) N, d

  nbytes = 16*(d**(2*N))
  print*, "+ Trying to allocate", nbytes, "bytes"

  allocate(rho(d**N,d**N))
end program test_limits
