! Quantum Information and Computing 2022-2023
! Nguyen Xuan Tung, ID: 2005491
! Week 7, Exercise: Ising model Consider 𝑁 spin-1/2 particles on a one-dimensional lattice, described by the Hamiltonian:

! 𝐻 = 𝜆𝑁∑︁𝑖𝜎𝑧+𝑁−1∑︁𝑖𝜎𝑥𝑖 𝜎𝑥𝑖+1

! where 𝜎𝑥 and 𝜎𝑧 are the Pauli matrices and 𝜆 is the interaction strength

! (a) Write a program to compute the 2𝑁 x 2𝑁 matrix representation of the Hamiltonian ˆ𝐻 for different 𝑁.
! (b) Diagonalize H for different 𝑁 = 1, ..., 𝑁𝑚𝑎𝑥 and 𝜆 ∈ [0 − 3]. What is the largest is 𝑁𝑚𝑎𝑥 you can reach?
! (c) Plot the first 𝑘 levels as a function of 𝜆 for different 𝑁 and comment on the spectrum.
 ! --------------------------------------------------------------------------------------
module isingmodel
  implicit none
  contains

  function ising_init_H(N,lambda) result(H)
    integer   :: N
    double precision :: lambda
    double complex, dimension(:,:), allocatable :: H, int_A, int_B

    integer :: ii,jj,kk,ll

    allocate(H(2**N,2**N))
    H = 0.0 * H

    ! External field part: \lambda \sum_i^N \sigma_z^i
    do ii = 1, N, 1
      do jj = 1, 2**N, 1
        H(jj,jj) = H(jj,jj) + -2*(modulo( (jj-1)/int(2**(N-ii)),2) ) +1
      end do
    end do
    H = lambda * H ! Adding the magnetization field factor

    ! Interaction part -\sum_i^{N-1}\sigma_x^{i+1}\sigma_x^i
    do ii = 1, N-1, 1
      allocate(int_A(2**N,2**N))
      allocate(int_B(2**N,2**N))
      int_A = int_A * 0.0
      int_B = int_B * 0.0
      do kk = 0,2**(ii-1)-1,1
        do jj=1,2**(N-ii),1
          int_A(kk*(2**(N-ii+1))  + 2**(N-ii)+jj, kk*(2**(N-ii+1))  + jj) = 1
          int_A(kk*(2**(N-ii+1))  + jj, kk*(2**(N-ii+1))  + 2**(N-ii)+jj) = 1
        end do
      end do
      do kk = 0,2**(ii)-1,1
        do jj=1,2**(N-ii-1),1
          int_B(kk*(2**(N-ii)) + 2**(N-ii-1)+jj, kk*(2**(N-ii)) + jj) = 1
          int_B(kk*(2**(N-ii)) + jj, kk*(2**(N-ii)) + 2**(N-ii-1)+jj) = 1
        end do
      end do
      if(.False. .eqv. .True.) then
      print*, "mata"
      do jj = 1, ubound(int_A, 1)
        print*, "|", real(int_A(jj, :)), "|"
      end do
      print*, "matB"
      do jj = 1, ubound(int_B, 1)
        print*, "|", real(int_B(jj, :)), "|"
      end do
      end if
      H = H - matmul(int_B,int_A)
      deallocate(int_A,int_B)
    end do
  end function
end module


module debugmod

  contains

  subroutine debugging(condition, msg, content)
    logical, intent(IN)                 :: condition
    character(*), intent(IN), optional  :: msg
    class(*), intent(IN), optional      :: content

    if(condition) then
      if (present(content)) then
        select type(content)
          type is (integer(1))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (integer(2))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (integer(4))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (integer(8))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (real(4))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (real(8))
            print*, "+ ", msg, " => [OK], Variable = ", content
          type is (logical)
            print*, "+ ", msg, " => [OK], Variable = ", content
        end select
      else
        print*, "+ ", msg, " => [OK]"
      end if
    end if
    end subroutine
end module debugmod



program ising
  use isingmodel
  use debugmod

  implicit none

  integer          :: N, d, SEP, nlams,ll
  double precision :: lambda
  character(20)    :: folder
  double complex, dimension(:,:), allocatable :: H

  ! Debugging
  logical :: DEBUG
  integer :: iostat,ii

  ! CPU times
  real*8 :: cpu_start, cpu_finish ! for the CPU times

  ! LAPACK variables
  double precision, dimension(:), allocatable   :: RWORK
  integer                                       :: INFO, LWORK
  integer, parameter                            :: LWMAX = 100000
  complex*16                                    :: WORK(LWMAX)
  complex(kind=8), dimension(:,:), allocatable  :: VR
  double precision, dimension(:), allocatable   :: eigenv

  ! PARAMETERS
  d   = 2 ! Spin 1/2 paticles
  SEP = 0
  DEBUG = .FALSE.

  print*, "+                ISING MODEL               +"
  print*, "+ This program investigates quantities for +"
  print*, "+ the Ising Model                          +"
  print*, "+------------------------------------------+"
  print*, "+ N: Number of spin 1/2 particles          +"
  print*, "+ nlams: Number of lambdas (range is [0,3])+"
  print*, "+------------------------------------------+"
  write(*,"(A)",advance='no') " + Type N, nlams, folder: "
  read (*,*, iostat=iostat) N, nlams, folder

  call system("mkdir -p ./data/lambdas/"//trim(folder))
  open(1, file = "./data/lambdas/"//trim(folder)//'/cputime')
  open(2, file = "./data/lambdas/"//trim(folder)//'/eigenvalue')

  ! Check if parameters are the right types
  if(iostat /= 0) then
    print*, "+ !!! INVALID PARAMETER TYPES !!!          +"
    print*, "+ N = integer                              +"
    print*, "+ nlams = integer                          +"
    print*, "+ Finish                                   +"

    stop
  end if

  ! Check if parameters are in the right ranges
  if(N < 2) then
    print*, "+ !!! INVALID PARAMETER RANGES !!!         +"
    print*, "+ (N > 1)                                  +"
    print*, "+ Finish                                   +"
    stop
  end if

  ! Check on READ
  if(DEBUG) then
    write(*,"(A)",advance='no') " + N:"
    print*, N
    write(*,"(A)",advance='no') " + d:"
    print*, d
    write(*,"(A)",advance='no') " + nlams:"
    print*, nlams
    if(SEP==1) then
      print*, "+ System is separable"
    else
      print*, "+ System is NOT separable"
    end if
  end if

  do ll = 0, nlams, 1
    lambda = (ll)*3/float(nlams)
    ! Preparation

    ! Building the Hamiltonian
    call cpu_time(cpu_start)  ! start time

    ! Diagonalizing and stuff
    H = ising_init_H(N,lambda)
    if(DEBUG) then
      print*, "+ lambda:", lambda
      do ii = 1, ubound(H, 1)
        print*, "|", real(H(ii, :)), "|"
      end do
    end if

    allocate(RWORK(3*(d**N)-2))
    allocate(VR(d**N,d**N))
    allocate(eigenv(d**N))
    ! Compute optimal size of workspace
    LWORK = -1

    call ZHEEV('N', 'U', d**N, H, d**N, eigenv, WORK,LWORK,RWORK,INFO)
    LWORK = min(LWMAX, int(WORK(1)))

    ! Compute eigenvalues
    call ZHEEV('N', 'U', d**N, H, d**N, eigenv, WORK,LWORK,RWORK,INFO)

    if(DEBUG) then
      write(*,'(A)',advance='no') " +   "
      do ii = 1, 3, 1
        write(*,'(A,ES14.5,A)', advance='no') "  ", eigenv(d**N - ii + 1), ", "
      end do
      print*, ""
    end if

    deallocate(RWORK, VR, H)

    call cpu_time(cpu_finish) ! end time

    write(1,*) cpu_finish - cpu_start
    write(2,*) eigenv

    deallocate(eigenv)
  end do

  print*, "+ Saving results in:                       +"
  print*, "+   ./data/lambdas/"//trim(folder)
  print*, "+------------------------------------------+"

end program ising


! To run the program: gfortran ising_mod.f90 -o ising_mod -llapack
! ./ising_mod
