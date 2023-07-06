! Quantum Information and Computing 2022-2023
! Nguyen Xuan Tung, ID: 2005491
! Week 2, Exercise 2: Derived types
!   In Fortran90 write a MODULE containing a double complex matrix derived TYPE that
!   includes the components: Matrix elements, Matrix Dimensions, Matrix Trace, and Matrix Adjoint.


!   (a) Define the correspondent TYPE.
!   (b) Define a function/subroutine that initializes this new TYPE.
!   (c) Define the functions/subroutines Trace and Adjoint.
!   (d) Define the correspondent Interfaces of the previous points.
!   (e) Define a subroutine that writes on file the Matrix TYPE in a readable form.
!   (f) Include everything in a test program.
! --------------------------------------------------------------------------------------

module matrices
    implicit none

!> @brief structure for a 2D matrix
!! dim: dimension of the matrix                 
!! element: actual array of elements               
!! trace: trace of the matrix               
!! det: determinant                        
!!---------------------------------------------------------------------------------------------------------
    type cmatrix
        integer, dimension(2)                   :: dim        ! dimension of the matrix
        complex*16, dimension(:,:), allocatable :: element
        complex*16                              :: trace, det
    end type cmatrix

    interface operator(.Adj.)
        module procedure cmatrix_adjoint
    end interface

    interface operator(.Trace.)
        module procedure cmatrix_trace
    end interface

    contains

!> @brief computes the trace 
!! @param[in] cmat = type(cmatrix), input complex matrix    
!! @return trace = complex*16, trace of thecomplex matrix  
!!------------------------------------------------------------------------------------------------
    function cmatrix_trace(cmat) result(trace)
        integer                   :: ii
        complex*16                :: trace
        type(cmatrix), intent(IN) :: cmat

        if(cmat%dim(1) == cmat%dim(2)) then ! iff the matrix is square
            trace = 0
            do ii = 1, size(cmat%element,1), 1
                trace = trace + cmat%element(ii,ii)
            end do
        else
            trace = 0
            trace = trace/trace
        end if
    end function cmatrix_trace

!> @brief initializes complex matrix randomly     
!! @param[in] nrow = integer, number of rows of matrix 
!! @param[in] ncol = integer, number of columns of matrix  
!! @param[in] range = real, range of numbers: 0<=x<=range
!! @return cmat = type(cmatrix), complex matrix  
!!-----------------------------------------------------------------------------------------
    function cmatrix_randinit(nrow,ncol,range) result(cmat)
        integer                 :: nrow, ncol
        real                    :: range
        real*16, dimension(:,:), allocatable :: elem_real,elem_imag

        type(cmatrix)          :: cmat

        cmat%dim = (/ nrow, ncol /) ! Assign the dimension

        allocate(elem_real(ncol,nrow))
        allocate(elem_imag(ncol,nrow))


        call random_number(elem_real)
        call random_number(elem_imag)

        elem_real = range*elem_real
        elem_imag = range*elem_imag

! Combining the matrices
        cmat%element = cmplx(elem_real,elem_imag)
        cmat%trace   = cmatrix_trace(cmat)
        cmat%det     = 0
    end function cmatrix_randinit

!> @brief initializes complex matrix type given the 2d array      
!! @param[in] array2d = real*16, dimension(:,:), input 2d array     
!! @return cmat = type(cmatrix), complex matrix   
!!------------------------------------------------------------------------------
    function cmatrix_init(array2d) result(cmat)
        integer                 :: nrow, ncol
        real                    :: range
        real*16, dimension(:,:), allocatable :: elem_real,elem_imag
        complex*16, dimension(:,:)           :: array2d

        type(cmatrix)          :: cmat

        cmat%dim = (/ size(array2d,1), size(array2d,2) /) ! Assign the dimension
        cmat%element = array2d
        cmat%trace   = cmatrix_trace(cmat)
        cmat%det   = 0
    end function cmatrix_init

!> @brief outputs the adjoints of the input complex matrix 
!! @param[in] cmat = type(cmatrix), input complex matrix  
!! @return cmat = type(cmatrix), adjoint matrix 
!!-----------------------------------------------------------------------------
    function cmatrix_adjoint(cmat) result(cmat_adj)
        type(cmatrix), intent(IN) :: cmat
        type(cmatrix)             :: cmat_adj

        cmat_adj = cmatrix_init( conjg(transpose(cmat%element)) )

    end function cmatrix_adjoint

!> @brief prints matrix on terminal  
!! @param[in] cmat = type(cmatrix), matrix to print               
!!-------------------------------------------------------------------------------
    subroutine cmatrix_print(cmat)
        type(cmatrix) :: cmat
        integer       :: ii

        do ii = 1, ubound(cmat%element, 1)
            print*, "|", cmat%element(ii, :), "|"
        end do

    end subroutine cmatrix_print

!> @brief writes matrix to file
!! @param[in] cmat = type(cmatrix), matrix to store              
!!------------------------------------------------------------------------------
    subroutine cmatrix_write(filename, cmat)
        character(len = 8) :: filename
        type(cmatrix)       :: cmat
        integer             :: ii

        open(1, file = filename, status = 'old')

        do ii = 1, ubound(cmat%element, 1)
            write(1, '(*(F0.8,SP,F0.8,"i",", "))')  cmat%element(ii, :)
        end do
    end subroutine cmatrix_write

end module matrices

program test
    use matrices

    type(cmatrix)               :: cmat, cmat_adj
    complex*16, dimension(2,2)  :: array

    array = transpose(reshape((/ 1, 2, 3, 4 /), shape(array)))
    cmat  = cmatrix_init(array)

    call cmatrix_print(cmat)
    cmat_adj = .Adj.cmat
    print*, "Dimension: ", cmat%dim
    print*, "Trace: ",     .Trace.cmat
    print*, "Adjoint: "
    call cmatrix_print(cmat_adj)

! write the matrix in file
    call cmatrix_write('matrices.csv', cmat)
end program test





