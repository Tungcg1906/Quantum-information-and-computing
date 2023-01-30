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

program matrix_mult_internal
    use matrix_product

    real*4, dimension(:,:), allocatable :: mat_A, mat_B, mat_C
    real*8 :: start, finish ! for the CPU times

! Perform a NxN NxN matrix multiplication
! Generate two N by N random matrices
    open(1, file = './matrix_mul.csv', status = 'new')
    print*, "t | Size"
! We save the times to perform a n by n matrix multiplication
    do n = 5, 100, 1
        allocate(mat_A(n, n))
        allocate(mat_B(n, n))

        mat_A = fill_no_matrix(n,n,10)
        mat_B = fill_no_matrix(n,n,10)
    ! TEST call graphics_printmatrix(mat_A)
        call cpu_time(start)  ! start time
        mat_C = matrix_multiplication(mat_A,mat_B)
        call cpu_time(finish) ! end time

    ! finish - start is Delta T
        print*, finish - start , "|", n      ! print on terminal
        write(1,*) finish - start , ",", n   ! save on file

        deallocate(mat_A)
        deallocate(mat_B)

    end do

end program matrix_mult_internal