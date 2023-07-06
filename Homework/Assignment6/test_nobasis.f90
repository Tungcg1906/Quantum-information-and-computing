! AFTER GENERATING D^N NUMBERS (FOR THE UNSEPARABLE CASE) 
! A WAY TO CONVERT THE INDEX 1 <= i <= d^N ONE CAN USE 
! A BASE CHANGE (FROM 10 TO D)

program changebase
  implicit none

  integer                            :: N, b_to, ii
  integer, dimension(:), allocatable :: vec
  integer                            :: number, inv_result

  N = 10

  print*, "--------------------------------------------"
  print*, "+               BASIS CHANGE               +"
  print*, "+------------------------------------------+"
  write(*,"(A)",advance='no') " + Bto, number: "
  read (*,*) b_to, number

  allocate(vec(N))

  do ii = 1, N, 1
    vec(N - ii + 1) = modulo(number,b_to)
    number = number/b_to
  end do

  print*, "+ Output"
  do ii = 1, N, 1
    print*, "+", vec(ii)
  end do
  print*, "+------------------------------------------+"
  print*, "+ Inverse:                                 +"

  inv_result = 0
  do ii = 1, N, 1
    inv_result = inv_result + vec(N - ii + 1)*b_to**(ii-1)
  end do

  write(*,"(A)",advance='no') " + Your number was:"
  print*, inv_result
  print*, "+------------------------------------------+"
  print*, "+------------------------------------------+"
end program changebase
