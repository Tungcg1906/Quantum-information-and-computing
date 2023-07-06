! Quantum Information and Computing 2022-2023
! Nguyen Xuan Tung, ID: 2005491
! Week 3, Exercise 1: Scaling of the matrix-matrix multiplication
!    Consider the code developed in the Exercise 3 from
!    Assignment 1 (matrix-matrix multiplication):


!   (a) Write a python script that changes N between two values N𝑚𝑖𝑛 and N𝑚𝑎𝑥, and launches the program.
!   (b) Store the results of the execution time in different files depending on the multiplication method
!   used.
!   (c) Fit the scaling of the execution time for different methods as a function of the input size. Consider
!   the largest possible difference between N𝑚𝑖𝑛 and N𝑚𝑎𝑥.
!   (d) Plot results for different multiplication methods.
! --------------------------------------------------------------------------------------
!Module matrix multiplication
module matrix_product
    implicit none
    contains
!> @brief function creates a matrix returns random values NxM
!!
!! @param[in] N : integer, number of columns in matrix 
!! @param[in] M : integer, number of rows in matrix
!! @param[in] rand_range : integer, range of random numbers
!! 
!! @return a matrixA = real*4, dimension(N,M), real random matrix
!!---------------------------------------------------------------------
    function fill_no_matrix(N, M, rand_range) result(matrixA)
        integer :: N ,M, rand_range
        real*4, dimension(N,M) :: matrixA
        call random_number(matrixA)
        matrixA = rand_range * matrixA
        
    end function fill_no_matrix

!> @brief function perfrom matrix multiplication through a loop method
!! @param[in] matrixA = real*4, dimension(N,M), real random matrix
!! @param[in] matrixB = real*4, dimension(N,M), real random matrix
!!
!! @return a matrixC = real*4, dimension(N,M), matrixA*matrixB
!! ----------------------------------------------------------------------
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

!> @brief function perfrom matrix multiplication through a 2nd order loop method
!! @param[in] matrixA = real*4, dimension(N,M), real random matrix
!! @param[in] matrixB = real*4, dimension(N,M), real random matrix
!!
!! @return a matrixC = real*4, dimension(N,M), matrixA*matrixB
!! --------------------------------------------------------------------------------
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
!> @brief print matrix in the terminal
!! @param[in] matrixA = real*4, dimension(N,M), real random matrix 
!!-----------------------------------------------------------------------
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
    real*8 :: start_time, finish_time ! for the CPU times



    open(1, file = './matrix_mul.csv', status = 'old') 
    print*, "t | Size"
    do n = 5, 1000, 1  

        allocate(matrixA(N, N))
        allocate(matrixB(N, N))  

! Creat random matrices
        matrixA = fill_no_matrix(N,N,10)
        matrixB = fill_no_matrix(N,N,10)

! TEST call graphics_printmatrix(mat_A)
        call cpu_time(start_time)  ! start time
        matrixC = mat_mul(matrixA,matrixB)
        call cpu_time(finish_time) ! end time

! finish - start is Delta T
!print*, finish - start , "|", n      ! print on terminal
        write(1,*) finish_time - start_time , ",", N   ! save on file

        deallocate(matrixA)
        deallocate(matrixB)
    end do

end program matrix_performance
