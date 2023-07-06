! Quantum Information and Computing 2022-2023
! Nguyen Xuan Tung, ID: 2005491
! Week 1, Exercise 3: Performance testing
!   Matrix-matrix multiplication is many times the bottleneck of linear algebra computations.

!   (a) Write explicitly the matrix-matrix multiplication loop in two different orders.
!   (b) Use the FORTRAN intrinsic function.
!   (c) Increase the matrix size and track the code performance using FORTRAN basic date and time
!   routines (i.e. CPU TIME).
!   (d)Use different optimization flags available with the compiler and compare the performance with
!   increasing matrix size.
! --------------------------------------------------------------------------------------

!Module matrix multiplication
module matrix_product
    implicit none
    contains
!create a matrix returns random values NxM
    function fill_no_matrix(N, M, rand_range) result(matrixA)
        integer :: N ,M, rand_range
        real*4, dimension(N,M) :: matrixA
        call random_number(matrixA)
        matrixA = rand_range * matrixA
        
    end function fill_no_matrix

!Assume A is lhs matrix, B rhs second, C is the result matrix
!first index "i" runs slower in c_ij 
    function matrix_multiplication(matrixA, matrixB) result(matrixC)
        integer :: ii, jj, kk
        logical :: check
        real*4, dimension(:,:)   :: matrixA, matrixB
        real*4, dimension(size(matrixA,2),size(matrixB,1)) :: matrixC

! Check if multiplication is possible (shapes)
        if (size(matrixA,2) .eq. size(matrixB,1)) then
            check = .TRUE.
        else
            print*, "Input matrices cannot be multiplied"
            check = .FALSE.
        end if

! Matrix multiplication
        if (check .eqv. .TRUE.) then
            do ii = 1, size(matrixA, 1), 1
                do jj = 1, size(matrixB, 2), 1
                    do kk = 1, size(matrixB, 1), 1
                        if (kk == 0) then
                            matrixC(ii,jj) = 0
                        end if
                        matrixC(ii,jj) = matrixC(ii,jj) + matrixA(ii,kk)*matrixB(kk,jj)
                    enddo
                enddo
            enddo
        end if
    end function matrix_multiplication


!The second orders (with transposed matrix)
    function matrix_multiplication_transposed (matrixA, matrixB) result(matrixC)
        integer :: ii, jj, kk
        logical :: check
        real*4, dimension(:,:)   :: matrixA, matrixB
        real*4, dimension(size(matrixA,2),size(matrixB,1)) :: matrixC

! Check if multiplication is possible (shapes)
        if (size(matrixA,2) .eq. size(matrixB,1)) then
            check = .TRUE.
        else
            print*, "Input matrices cannot be multiplied"
            check = .FALSE.
        end if

! Matrix multiplication
        if (check .eqv. .TRUE.) then
            do kk = 1, size(matrixA,1), 1
                do jj = 1, size(matrixB,2), 1
                    do ii = 1, size(matrixB,1), 1
                        if (kk == 1) then
                            matrixC(ii,jj) = 0
                        end if
                        matrixC(ii,jj) = matrixC(ii,jj) + matrixA(ii,kk)*matrixB(kk,jj)
                    enddo
                enddo
            enddo
        end if
    end function matrix_multiplication_transposed 


end module matrix_product

! Printing matrix
module matrix_print
    implicit none
    contains

    subroutine printmatrix(matrixA)
    integer                 :: ii
    real*4, dimension(:,:)  :: matrixA

    do ii = 1, ubound(matrixA, 1)
        print*, "|", matrixA(ii, :), "|"
    end do

    end subroutine printmatrix

end module matrix_print


program matrix_performance
    use matrix_product
    use matrix_print


    real*4, dimension(:,:), allocatable :: matrixA, matrixB, matrixC
    real*8 :: start_time, finish_time
    real*8 :: time_delta1, time_delta2, time_delta3

! Perform a NxN NxN matrix multiplication
! Generate two N by N random matrices and
    n = 1000
    allocate(matrixA(n, n))
    allocate(matrixB(n, n))

! Fill the matrix with radom values
    matrixA = fill_no_matrix(n,n,10)
    matrixB = fill_no_matrix(n,n,10)

! Performance of matrix multiplication

    call cpu_time(start_time)  ! start time
    matrixC = matrix_multiplication(matrixA,matrixB)
    call cpu_time(finish_time) ! end time
    
    time_delta1 = finish_time - start_time
    print*, "Time passed using matrix multiplication", time_delta1 , "s" 

    
! Performance of transposed matrix multiplication

    call cpu_time(start_time)  ! start time
    matrixC = matrix_multiplication_transposed(matrixA,matrixB)
    call cpu_time(finish_time) ! end time

    time_delta2 = finish_time - start_time
    print*, "Time passed using transposed matrix multiplication", time_delta2 , "s"  

! Performance of matmul

    call cpu_time(start_time)  ! start time
    matrixC = matmul(matrixA,matrixB)
    call cpu_time(finish_time) ! end time

    time_delta3 = finish_time - start_time
    print*, "Time passed using transposed matrix multiplication", time_delta3 , "s"  


    deallocate(matrixA)
    deallocate(matrixB)


end program matrix_performance
