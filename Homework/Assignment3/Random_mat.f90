! Quantum Information and Computing 2022-2023
! Nguyen Xuan Tung, ID: 2005491
! Week 3, Exercise 3: Random matrix theory
!    Study 𝑃(𝑠), the distribution of normalized spacing 𝑠 defined in the previous
!    exercise, accumulating values from different random matrices of size at least N = 1000.

!   (a) Compute 𝑃(𝑠) for a random hermitian matrix.
!   (b) Compute 𝑃(𝑠) for a diagonal matrix with real random entries.
!   (c) Fit the corresponding distributions with the function: P(s) = as α exp(−bs β ) and report and report α, β, a, b.
! --------------------------------------------------------------------------------------

module cmatrices
implicit none

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

    interface operator(.Eigens.)
        module procedure cmatrix_heigenvalues
    end interface

    interface operator(.Eigenspacing.)
        module procedure cmatrix_heigenspacing
    end interface

    contains

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
! Weird way to assign trace as NaN
        trace = 0
        trace = trace/trace
        end if
    end function cmatrix_trace

    function cmatrix_randinit(nrow,ncol,range) result(cmat)
        integer                 :: nrow, ncol
        real                    :: range
        real*16, dimension(:,:), allocatable :: elem_real,elem_imag

        type(cmatrix)          :: cmat

        cmat%dim = (/ nrow, ncol /) ! Assign the dimension

! To define the complex matrix, we initialize 2 real matrices for the
! real and imaginary components and then we combine them together
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

    function cmatrix_randinit_hermitian(n,range) result(cmat)
        integer                 :: n
        real                    :: range
        real*16, dimension(:,:), allocatable :: elem_real,elem_imag
        integer                 :: ii,jj,kk

        type(cmatrix)          :: cmat

        cmat%dim = (/ n, n /) ! Assign the dimension

! To define the complex matrix, we initialize 2 real matrices for the
! real and imaginary components and then we combine them together
        allocate(elem_real(n,n))
        allocate(elem_imag(n,n))

! Imaginary part diagonal
        do ii = 1, n, 1
            call random_number(elem_real(ii, ii))
            call random_number(elem_imag(ii, ii))
        end do

! Upper right part
        do ii = 1, n, 1
            do jj = ii + 1, n, 1
                elem_real(jj, ii) = (RAND(0) - 0.5 )*range*2
                elem_imag(jj, ii) = (RAND(0) - 0.5 )*range*2
            end do
        end do

! Flip negatively imaginary values
        do ii = 2, n, 1
            do jj = 1,ii-1, 1
                elem_real(jj, ii) = + elem_real(ii,jj)
                elem_imag(jj, ii) = - elem_imag(ii,jj)
            end do
        end do

! Combining the matrices
        cmat%element = cmplx(elem_real,elem_imag)
        cmat%trace   = cmatrix_trace(cmat)
        cmat%det     = 0
    end function cmatrix_randinit_hermitian

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

    function cmatrix_adjoint(cmat) result(cmat_adj)
        type(cmatrix), intent(IN) :: cmat
        type(cmatrix)             :: cmat_adj

        cmat_adj = cmatrix_init( conjg(transpose(cmat%element)) )
    end function cmatrix_adjoint

    function cmatrix_heigenvalues(cmat) result(eigens)
        type(cmatrix), intent(IN)                   :: cmat
        double precision, dimension(:), allocatable :: eigens, RWORK
        integer                                     :: INFO, LWORK
        integer                                     :: N
        integer, parameter                          :: LWMAX = 100000
        complex*16                                  :: WORK(LWMAX)

! Check if matrix is squared
        if(cmat%dim(1) == cmat%dim(2)) then
            N = cmat%dim(1)

            allocate(eigens(N))
            allocate(RWORK(3*N-2))

! Compute optimal size of workspace
            LWORK = -1
            call ZHEEV('N', 'U', N, cmat%element, N, eigens, WORK, LWORK, RWORK, INFO)
            LWORK = min(LWMAX, int(WORK(1)))

! Compute eigenvalues
            call ZHEEV('N', 'U', N, cmat%element, N, eigens, WORK, LWORK, RWORK, INFO)
        end if
    end function cmatrix_heigenvalues

    function cmatrix_heigenspacing(cmat) result(spacing)
        type(cmatrix), intent(IN)                   :: cmat
        double precision, dimension(:), allocatable :: eigens, spacing
        double precision                            :: inverseaveragelambda
        integer                                     :: ii, N

! Check if matrix is squared
        if(cmat%dim(1) == cmat%dim(2)) then
            N = cmat%dim(1)
            allocate(eigens(N))
            allocate(spacing(N-1))

            eigens = cmatrix_heigenvalues(cmat)
            inverseaveragelambda = (N-1)/ (eigens(N)-eigens(1))

            do ii = 1, N-1, 1
                spacing(ii) = inverseaveragelambda * (eigens(ii+1) - eigens(ii))
            end do

        end if
    end function cmatrix_heigenspacing

    subroutine cmatrix_print(cmat)
        type(cmatrix) :: cmat
        integer       :: ii

        do ii = 1, ubound(cmat%element, 1)
            print*, "|", cmat%element(ii, :), "|"
        end do

    end subroutine cmatrix_print

    subroutine cmatrix_write(filename, cmat)
        character(len = 25) :: filename
        type(cmatrix)       :: cmat
        integer             :: ii

        open(1, file = filename, status = 'new')

        do ii = 1, ubound(cmat%element, 1)
            write(1, '(*(F0.4, 1x, SP, F0.4, "i", 2x))')  cmat%element(ii, :)
        end do
    end subroutine cmatrix_write
end module cmatrices

program randommat
    use cmatrices

    type(cmatrix), allocatable                  :: cmat, cmat_diag
    double precision, dimension(:), allocatable :: spacing,spacing_diag
    integer, parameter                          :: N = 1000, NMATS = 100
    integer                                     :: ii, jj
    complex*16, dimension(:,:), allocatable     :: diag
    real, parameter                             :: range = 1

    open(101, file = './s_1.csv', status = 'old')
    open(102, file = './s_2.csv', status = 'old')

    call RANDOM_SEED()
    do ii = 1, NMATS, 1

    !------------------------------------------!
    !------ 3.a RANDOM HERMITIAN MATRIX -------!
    !------------------------------------------!
        cmat = cmatrix_randinit_hermitian(N,range)
        spacing = .Eigenspacing.cmat
        write(101,*) spacing   ! save on file

    !------------------------------------------!
    !------ 3.b RANDOM DIAGONAL MATRIX --------!
    !------------------------------------------!
        allocate(diag(N,N))
        diag = diag * 0.0
        do jj = 1, N, 1
            diag(jj,jj) = (RAND(0) - .5)*range*2
        end do
        cmat_diag = cmatrix_init(diag)
        spacing_diag = .Eigenspacing.cmat_diag
        write(102,*) spacing_diag   ! save on file

        deallocate(cmat)
        deallocate(diag)
    end do
end program randommat