! Quantum Information and Computing 2022-2023
! Nguyen Xuan Tung, ID: 2005491
! Week 2, Exercise 2: Documentation
!   Rewrite Exercise 3 from Assignment 1 including

!   (a) Documentation.
!   (b) Comments.
!   (c) Pre- and post- conditions.
!   (d) Error handling.
!   (e) Checkpoints. 
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

! Module for debugging function
! Initialize DEBUG_ and set it .TRUE. to be in debug mode
module debug
    use matrix_product
    use matrix_print
    contains

!> @brief print a real value and the line   
!! @param[in] realarg = real
!!---------------------------------------------------------------------------------
    subroutine check_real(DEBUG, realarg, line)
        logical        :: DEBUG   ! input, if DEBUG_ == TRUE we are in 'debug mode'
        real           :: realarg ! optional generic real argument
        integer        :: line

        if (DEBUG .eqv. .TRUE.) then
            print *, 'LINE:',line   ! print file and line
            print*,'arg:', realarg
        end if
    end subroutine check_real

!> @brief print a matrix and its lines
!! @param[in] matrix = real*4, dimension(:,:)
!!--------------------------------------------------------------------------------------
    subroutine check_matrix(DEBUG, mat, line)
        logical        :: DEBUG   ! input, if DEBUG_ == TRUE we are in 'debug mode'
        integer        :: line
        real*4, dimension(:,:), allocatable :: mat

        if (DEBUG .eqv. .TRUE.) then
            print *, 'LINE:',line   ! print file and line
            call printmatrix(mat)
        end if
    end subroutine check_matrix

!> @brief check for custom implemented matrix multiplication using matmul 
!! @param[in] matrixA = real*4, dimension(:,:)
!! @param[in] matrixB = real*4, dimension(:,:)
!!--------------------------------------------------------------------------------------
    subroutine check_custom_matrix_multiplication(DEBUG, mat_A, mat_B)
        logical        :: DEBUG   ! input, if DEBUG_ == TRUE we are in 'debug mode'
        real*4, dimension(:,:), allocatable :: mat_A, mat_B, mat_C, mat_D, mat_E

! Initialize matrix according to our function
        mat_C = matrix_multiplication(mat_A,mat_B)

! Initialize matrix using fortran function
        mat_D = matmul(mat_A,mat_B)

        mat_E = mat_C - mat_D
        call printmatrix(mat_E)

    end subroutine check_custom_matrix_multiplication

end module debug


! Macro for passing __LINE__ and DEBUG_ automatically
#define check_real(realarg) check_real(DEBUG, realarg,__LINE__)
#define check_matrix(mat) check_matrix(DEBUG, mat,__LINE__)
#define check_custom_matrix_multiplication(mat_A, mat_B) check_custom_matrix_multiplication_(DEBUG, mat_A, mat_B)



program matrix_performance
    use matrix_product
    use matrix_print
    implicit none

    integer :: row1, row2, col1, col2
    real*4, dimension(:,:), allocatable :: matrixA, matrixB, matrixC
    real*8 :: start_time, finish_time 

! Testing matrix multiplication
! Input from the terminal for the dimension of the matrices
    print*, "Matrix multiplication A*B"
    print*, "Enter the size of matrix A:"
    read*, row1, col1
    print*, "Enter the size of matrix B:"
    read*, row2, col2

! Print the sizes of matrices
    print*, "size A :", "(", row1, ",", col1, ")"
    print*, "size B :", "(", row2, ",", col2, ")"

! Creat random matrices
    matrixA = fill_no_matrix(row1,col1,10)
    matrixB = fill_no_matrix(row2,col2,10)

! Print matrices for testing!
    print*,""
    print*,"A:"
    call printmatrix(matrixA)

    print*,""
    print*,"B:"
    call printmatrix(matrixB)

! Multiplication matrix
    matrixC = matrix_multiplication(matrixA,matrixB)
    print*,""
    print*,"C:"
    call printmatrix(matrixC)

! Multiplication transposed matrix
    matrixC = matrix_multiplication_transposed (matrixA,matrixB)
    print*,""
    print*,"C transposed:"
    call printmatrix(matrixC)

! Multiplication matmul
    matrixC = matmul (matrixA,matrixB)
    print*,""
    print*,"C matmul"
    call printmatrix(matrixC)
    
end program matrix_performance

! Error appears when compiling, in order to compile sucessful, a flag -cpp is used
! gfortran -o documentation documentation.f90 -cpp
! ./documentation