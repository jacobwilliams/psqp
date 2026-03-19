!*******************************************************************************
!> author: Jacob Williams
!
!  Unit test for the sparse matrix module.

program test_sparse_matrix

    use psqp_sparse_module
    use psqp_kind_module, only: wp => psqp_wp

    implicit none

    logical :: all_tests_passed

    write(*,'(A)') ''
    write(*,'(A)') '========================================='
    write(*,'(A)') 'Testing psqp_sparse_module'
    write(*,'(A)') '========================================='
    write(*,'(A)') ''

    all_tests_passed = .true.

    ! Run tests
    all_tests_passed = test_csr_basic() .and. all_tests_passed
    all_tests_passed = test_matvec() .and. all_tests_passed
    all_tests_passed = test_row_operations() .and. all_tests_passed
    all_tests_passed = test_pattern() .and. all_tests_passed

    write(*,'(A)') ''
    write(*,'(A)') '========================================='
    if (all_tests_passed) then
        write(*,'(A)') 'All tests PASSED!'
    else
        write(*,'(A)') 'Some tests FAILED!'
        error stop 1
    end if
    write(*,'(A)') '========================================='
    write(*,'(A)') ''

contains

!*******************************************************************************
!>
!  Test basic CSR initialization and destruction.

    function test_csr_basic() result(passed)

        logical :: passed
        type(sparse_matrix_csr) :: A
        integer, parameter :: nrows = 3, ncols = 4, nnz = 5
        integer :: row_ptr(nrows+1), col_ind(nnz)

        write(*,'(A)') 'Test: CSR basic initialization...'

        ! Create a simple sparse matrix:
        ! [1.0  0.0  2.0  0.0]
        ! [0.0  3.0  0.0  4.0]
        ! [5.0  0.0  0.0  0.0]

        row_ptr = [1, 3, 5, 6]  ! row 1 has 2 nz, row 2 has 2 nz, row 3 has 1 nz
        col_ind = [1, 3, 2, 4, 1]  ! column indices

        call A%initialize(nrows, ncols, row_ptr, col_ind)

        passed = A%is_initialized()
        passed = passed .and. (A%nrows == nrows)
        passed = passed .and. (A%ncols == ncols)
        passed = passed .and. (A%nnz == nnz)

        ! Fill in values
        A%values = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]

        call A%destroy()

        passed = passed .and. (.not. A%is_initialized())

        write(*,'(A,L1)') '  Result: ', passed

    end function test_csr_basic

!*******************************************************************************
!>
!  Test sparse matrix-vector multiplication.

    function test_matvec() result(passed)

        logical :: passed
        type(sparse_matrix_csr) :: A
        integer, parameter :: nrows = 3, ncols = 4
        integer :: row_ptr(nrows+1), col_ind(5)
        real(wp) :: x(ncols), y(nrows), y_expected(nrows)
        real(wp), parameter :: tol = 1.0e-14_wp

        write(*,'(A)') 'Test: Sparse matvec (A*x)...'

        ! Matrix A:
        ! [1.0  0.0  2.0  0.0]
        ! [0.0  3.0  0.0  4.0]
        ! [5.0  0.0  0.0  0.0]

        row_ptr = [1, 3, 5, 6]
        col_ind = [1, 3, 2, 4, 1]

        call A%initialize(nrows, ncols, row_ptr, col_ind)
        A%values = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]

        ! Test vector
        x = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]

        ! Compute y = A*x
        call A%matvec(x, y)

        ! Expected: [1*1 + 2*3, 3*2 + 4*4, 5*1] = [7, 22, 5]
        y_expected = [7.0_wp, 22.0_wp, 5.0_wp]

        passed = all(abs(y - y_expected) < tol)

        write(*,'(A,L1)') '  Result: ', passed
        if (.not. passed) then
            write(*,'(A,3F10.4)') '  Expected: ', y_expected
            write(*,'(A,3F10.4)') '  Got:      ', y
        end if

        call A%destroy()

    end function test_matvec

!*******************************************************************************
!>
!  Test row operations.

    function test_row_operations() result(passed)

        logical :: passed
        type(sparse_matrix_csr) :: A
        integer, parameter :: nrows = 3, ncols = 4
        integer :: row_ptr(nrows+1), col_ind(5)
        real(wp) :: row_vec(ncols), expected(ncols)
        real(wp) :: dot_result, expected_dot
        real(wp), parameter :: tol = 1.0e-14_wp
        logical :: test1, test2, test3

        write(*,'(A)') 'Test: Row operations...'

        ! Matrix A:
        ! [1.0  0.0  2.0  0.0]
        ! [0.0  3.0  0.0  4.0]
        ! [5.0  0.0  0.0  0.0]

        row_ptr = [1, 3, 5, 6]
        col_ind = [1, 3, 2, 4, 1]

        call A%initialize(nrows, ncols, row_ptr, col_ind)
        A%values = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]

        ! Test 1: get_row
        call A%get_row(2, row_vec)
        expected = [0.0_wp, 3.0_wp, 0.0_wp, 4.0_wp]
        test1 = all(abs(row_vec - expected) < tol)

        ! Test 2: row_dot
        dot_result = A%row_dot(1, [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp])
        expected_dot = 1.0_wp + 2.0_wp  ! row 1: [1, 0, 2, 0] dot [1,1,1,1]
        test2 = abs(dot_result - expected_dot) < tol

        ! Test 3: get_nnz_in_row
        test3 = (A%get_nnz_in_row(1) == 2) .and. &
                (A%get_nnz_in_row(2) == 2) .and. &
                (A%get_nnz_in_row(3) == 1)

        passed = test1 .and. test2 .and. test3

        write(*,'(A,L1)') '  Result: ', passed
        if (.not. passed) then
            write(*,'(A,L1)') '    get_row test:      ', test1
            write(*,'(A,L1)') '    row_dot test:      ', test2
            write(*,'(A,L1)') '    get_nnz_in_row test: ', test3
        end if

        call A%destroy()

    end function test_row_operations

!*******************************************************************************
!>
!  Test sparse pattern.

    function test_pattern() result(passed)

        logical :: passed
        type(sparse_pattern) :: pattern
        integer, parameter :: nrows = 2, ncols = 3
        integer :: row_ptr(nrows+1), col_ind(4)

        write(*,'(A)') 'Test: Sparse pattern...'

        row_ptr = [1, 3, 5]  ! row 1 has 2 nz, row 2 has 2 nz
        col_ind = [1, 2, 2, 3]

        call pattern%initialize(nrows, ncols, row_ptr, col_ind)

        passed = pattern%is_initialized()
        passed = passed .and. (pattern%nrows == nrows)
        passed = passed .and. (pattern%ncols == ncols)
        passed = passed .and. (pattern%nnz == 4)

        call pattern%destroy()

        passed = passed .and. (.not. pattern%is_initialized())

        write(*,'(A,L1)') '  Result: ', passed

    end function test_pattern

end program test_sparse_matrix
