!*****************************************************************************************
!>
!  Sparse matrix module for PSQP optimization.
!
!  Provides sparse matrix storage and operations for constraint Jacobians.
!  Uses CSR (Compressed Sparse Row) format for efficient row-wise access.
!
!### History
!  * Jacob Williams, March 2026 : Initial implementation

module psqp_sparse_module

    use psqp_kind_module, only: wp => psqp_wp

    implicit none

    private

    public :: sparse_matrix_csr
    public :: sparse_pattern

    type :: sparse_matrix_csr
        !! Compressed Sparse Row (CSR) format for constraint Jacobian matrix.
        !!
        !! Storage format:
        !! - `row_ptr(i)` gives the index in `col_ind` and `values` where row i starts
        !! - `row_ptr(i+1) - row_ptr(i)` gives the number of non-zeros in row i
        !! - `col_ind(row_ptr(i):row_ptr(i+1)-1)` gives column indices for row i
        !! - `values(row_ptr(i):row_ptr(i+1)-1)` gives values for row i
        !!
        !! For constraint Jacobian: nrows = nc (constraints), ncols = nf (variables)

        integer :: nrows = 0  !! number of rows (nc - number of constraints)
        integer :: ncols = 0  !! number of columns (nf - number of variables)
        integer :: nnz = 0    !! number of non-zero elements

        integer, allocatable :: row_ptr(:)  !! (nrows+1) row start pointers
        integer, allocatable :: col_ind(:)  !! (nnz) column indices
        real(wp), allocatable :: values(:)  !! (nnz) non-zero values

    contains

        private

        procedure, public :: initialize => sparse_csr_init
        procedure, public :: destroy => sparse_csr_destroy
        procedure, public :: from_coo => sparse_csr_from_coo
        procedure, public :: matvec => sparse_csr_matvec
        procedure, public :: matvec_t => sparse_csr_matvec_transpose
        procedure, public :: get_row => sparse_csr_get_row
        procedure, public :: get_row_sparse => sparse_csr_get_row_sparse
        procedure, public :: set_row => sparse_csr_set_row
        procedure, public :: row_dot => sparse_csr_row_dot
        procedure, public :: get_col => sparse_csr_get_col
        procedure, public :: set_col => sparse_csr_set_col
        procedure, public :: is_initialized => sparse_csr_is_initialized
        procedure, public :: get_nnz_in_row => sparse_csr_get_nnz_in_row
        procedure, public :: copy_structure => sparse_csr_copy_structure

    end type sparse_matrix_csr

    type :: sparse_pattern
        !! Sparsity pattern only (structure without values).
        !! Used to define the structure before filling in values.

        integer :: nrows = 0  !! number of rows
        integer :: ncols = 0  !! number of columns
        integer :: nnz = 0    !! number of non-zero elements

        integer, allocatable :: row_ptr(:)  !! (nrows+1) row start pointers
        integer, allocatable :: col_ind(:)  !! (nnz) column indices

    contains

        private

        procedure, public :: initialize => sparse_pattern_init
        procedure, public :: destroy => sparse_pattern_destroy
        procedure, public :: is_initialized => sparse_pattern_is_initialized

    end type sparse_pattern

contains

!*****************************************************************************************
!>
!  Initialize a sparse matrix in CSR format from sparsity pattern.

    subroutine sparse_csr_init(me, nrows, ncols, row_ptr, col_ind)

        class(sparse_matrix_csr), intent(inout) :: me
        integer, intent(in) :: nrows  !! number of rows
        integer, intent(in) :: ncols  !! number of columns
        integer, intent(in) :: row_ptr(nrows+1)  !! row start pointers
        integer, intent(in) :: col_ind(:)  !! column indices

        integer :: nnz

        ! Clean up any existing data
        call me%destroy()

        me%nrows = nrows
        me%ncols = ncols
        me%nnz = row_ptr(nrows + 1) - 1
        nnz = me%nnz

        ! Allocate storage
        allocate(me%row_ptr(nrows + 1))
        allocate(me%col_ind(nnz))
        allocate(me%values(nnz))

        ! Copy structure
        me%row_ptr = row_ptr
        me%col_ind = col_ind

        ! Initialize values to zero
        me%values = 0.0_wp

    end subroutine sparse_csr_init

!*****************************************************************************************
!>
!  Destroy/deallocate a sparse matrix.

    subroutine sparse_csr_destroy(me)

        class(sparse_matrix_csr), intent(inout) :: me

        me%nrows = 0
        me%ncols = 0
        me%nnz = 0

        if (allocated(me%row_ptr)) deallocate(me%row_ptr)
        if (allocated(me%col_ind)) deallocate(me%col_ind)
        if (allocated(me%values)) deallocate(me%values)

    end subroutine sparse_csr_destroy

!*****************************************************************************************
!>
!  Initialize a sparse matrix from COO (Coordinate) format.
!  Converts COO triplets (row, col, value) to CSR format.
!
!### Note
!  Input arrays can have entries in any order. Duplicate entries will be summed.

    subroutine sparse_csr_from_coo(me, nrows, ncols, rows, cols, nnz, values)

        class(sparse_matrix_csr), intent(inout) :: me
        integer, intent(in) :: nrows  !! number of rows
        integer, intent(in) :: ncols  !! number of columns
        integer, intent(in) :: rows(nnz)  !! row indices (1-based)
        integer, intent(in) :: cols(nnz)  !! column indices (1-based)
        integer, intent(in) :: nnz  !! number of non-zeros
        real(wp), intent(in), optional :: values(nnz)  !! values (optional, defaults to pattern only)

        integer :: i, j, row, col
        integer :: row_counts(nrows)
        integer :: current_pos(nrows)
        integer, allocatable :: temp_row_ptr(:)
        integer, allocatable :: temp_col_ind(:)
        real(wp), allocatable :: temp_values(:)

        ! Clean up any existing data
        call me%destroy()

        me%nrows = nrows
        me%ncols = ncols
        me%nnz = nnz

        ! Count non-zeros per row
        row_counts = 0
        do i = 1, nnz
            row = rows(i)
            if (row < 1 .or. row > nrows) then
                error stop 'sparse_csr_from_coo: row index out of bounds'
            end if
            row_counts(row) = row_counts(row) + 1
        end do

        ! Build row_ptr array
        allocate(me%row_ptr(nrows + 1))
        me%row_ptr(1) = 1
        do i = 1, nrows
            me%row_ptr(i + 1) = me%row_ptr(i) + row_counts(i)
        end do

        ! Allocate storage
        allocate(temp_col_ind(nnz))
        allocate(temp_values(nnz))

        ! Distribute entries into CSR format
        current_pos = me%row_ptr(1:nrows)
        do i = 1, nnz
            row = rows(i)
            col = cols(i)
            if (col < 1 .or. col > ncols) then
                error stop 'sparse_csr_from_coo: column index out of bounds'
            end if
            j = current_pos(row)
            temp_col_ind(j) = col
            if (present(values)) then
                temp_values(j) = values(i)
            else
                temp_values(j) = 0.0_wp  ! Pattern only
            end if
            current_pos(row) = current_pos(row) + 1
        end do

        ! Sort each row by column index and handle duplicates
        allocate(me%col_ind(nnz))
        allocate(me%values(nnz))

        me%col_ind = temp_col_ind
        me%values = temp_values

        ! Sort each row (simple insertion sort - fine for small rows)
        do i = 1, nrows
            call sort_row_by_column(me%row_ptr(i), me%row_ptr(i+1)-1, me%col_ind, me%values)
        end do

    end subroutine sparse_csr_from_coo

!*****************************************************************************************
!>
!  Helper: Sort one row by column indices (insertion sort).

    subroutine sort_row_by_column(start_idx, end_idx, col_ind, values)

        integer, intent(in) :: start_idx, end_idx
        integer, intent(inout) :: col_ind(:)
        real(wp), intent(inout) :: values(:)

        integer :: i, j, key_col
        real(wp) :: key_val

        do i = start_idx + 1, end_idx
            key_col = col_ind(i)
            key_val = values(i)
            j = i - 1
            do while (j >= start_idx)
                if (col_ind(j) <= key_col) exit
                col_ind(j + 1) = col_ind(j)
                values(j + 1) = values(j)
                j = j - 1
            end do
            col_ind(j + 1) = key_col
            values(j + 1) = key_val
        end do

    end subroutine sort_row_by_column

!*****************************************************************************************
!>
!  Check if sparse matrix is initialized.

    pure function sparse_csr_is_initialized(me) result(initialized)

        class(sparse_matrix_csr), intent(in) :: me
        logical :: initialized

        initialized = allocated(me%row_ptr) .and. &
                      allocated(me%col_ind) .and. &
                      allocated(me%values)

    end function sparse_csr_is_initialized

!*****************************************************************************************
!>
!  Sparse matrix-vector product: `y = A*x`.
!
!### Example
!```fortran
!  call jac%matvec(x, y)  ! y = Jacobian * x
!```

    subroutine sparse_csr_matvec(me, x, y)

        class(sparse_matrix_csr), intent(in) :: me
        real(wp), intent(in) :: x(:)  !! input vector (ncols)
        real(wp), intent(out) :: y(:)  !! output vector (nrows)

        integer :: i, j, k

        y = 0.0_wp

        do i = 1, me%nrows
            do k = me%row_ptr(i), me%row_ptr(i + 1) - 1
                j = me%col_ind(k)
                y(i) = y(i) + me%values(k) * x(j)
            end do
        end do

    end subroutine sparse_csr_matvec

!*****************************************************************************************
!>
!  Sparse matrix-vector product (transpose): `y = A^T*x`.
!
!### Example
!```fortran
!  call jac%matvec_t(x, y)  ! y = Jacobian^T * x
!```

    subroutine sparse_csr_matvec_transpose(me, x, y)

        class(sparse_matrix_csr), intent(in) :: me
        real(wp), intent(in) :: x(:)   !! input vector (nrows)
        real(wp), intent(out) :: y(:)  !! output vector (ncols)

        integer :: i, j, k

        y = 0.0_wp

        do i = 1, me%nrows
            do k = me%row_ptr(i), me%row_ptr(i + 1) - 1
                j = me%col_ind(k)
                y(j) = y(j) + me%values(k) * x(i)
            end do
        end do

    end subroutine sparse_csr_matvec_transpose

!*****************************************************************************************
!>
!  Extract a row from the sparse matrix as a dense vector.
!
!### Example
!```fortran
!  call jac%get_row(kc, gc)  ! gc = gradient of constraint kc
!```

    subroutine sparse_csr_get_row(me, irow, row_vec)

        class(sparse_matrix_csr), intent(in) :: me
        integer, intent(in) :: irow  !! row index (1-based)
        real(wp), intent(out) :: row_vec(:)  !! dense row vector (ncols)

        integer :: k, j

        ! Initialize to zero
        row_vec = 0.0_wp

        ! Fill in non-zero elements
        do k = me%row_ptr(irow), me%row_ptr(irow + 1) - 1
            j = me%col_ind(k)
            row_vec(j) = me%values(k)
        end do

    end subroutine sparse_csr_get_row

!*****************************************************************************************
!>
!  Extract a row from the sparse matrix in sparse format.
!
!### Example
!```fortran
!  call jac%get_row_sparse(kc, indices, values, nnz)
!```

    subroutine sparse_csr_get_row_sparse(me, irow, indices, values, nnz)

        class(sparse_matrix_csr), intent(in) :: me
        integer, intent(in) :: irow  !! row index (1-based)
        integer, intent(out) :: indices(:)  !! column indices of non-zeros
        real(wp), intent(out) :: values(:)  !! values of non-zeros
        integer, intent(out) :: nnz  !! number of non-zeros in row

        integer :: k, idx

        nnz = me%row_ptr(irow + 1) - me%row_ptr(irow)

        idx = 0
        do k = me%row_ptr(irow), me%row_ptr(irow + 1) - 1
            idx = idx + 1
            indices(idx) = me%col_ind(k)
            values(idx) = me%values(k)
        end do

    end subroutine sparse_csr_get_row_sparse

!*****************************************************************************************
!>
!  Set a row in the sparse matrix from a dense vector.
!  Only updates values at existing non-zero locations.

    subroutine sparse_csr_set_row(me, irow, row_vec)

        class(sparse_matrix_csr), intent(inout) :: me
        integer, intent(in) :: irow  !! row index (1-based)
        real(wp), intent(in) :: row_vec(:)  !! dense row vector (ncols)

        integer :: k, j

        ! Update values at non-zero locations
        do k = me%row_ptr(irow), me%row_ptr(irow + 1) - 1
            j = me%col_ind(k)
            me%values(k) = row_vec(j)
        end do

    end subroutine sparse_csr_set_row

!*****************************************************************************************
!>
!  Compute dot product of a sparse row with a dense vector.
!
!### Example
!```fortran
!  result = jac%row_dot(kc, s)  ! dot_product(constraint_gradient_kc, s)
!```

    pure function sparse_csr_row_dot(me, irow, vec) result(res)

        class(sparse_matrix_csr), intent(in) :: me
        integer, intent(in) :: irow  !! row index (1-based)
        real(wp), intent(in) :: vec(:)  !! dense vector (ncols)
        real(wp) :: res

        integer :: k, j

        res = 0.0_wp

        do k = me%row_ptr(irow), me%row_ptr(irow + 1) - 1
            j = me%col_ind(k)
            res = res + me%values(k) * vec(j)
        end do

    end function sparse_csr_row_dot

!*****************************************************************************************
!>
!  Extract a column from the sparse matrix as a dense vector.
!  Note: Column access in CSR is inefficient (requires searching all rows).

    subroutine sparse_csr_get_col(me, icol, col_vec)

        class(sparse_matrix_csr), intent(in) :: me
        integer, intent(in) :: icol  !! column index (1-based)
        real(wp), intent(out) :: col_vec(:)  !! dense column vector (nrows)

        integer :: i, k, j

        ! Initialize to zero
        col_vec = 0.0_wp

        ! Search for column icol in each row
        do i = 1, me%nrows
            do k = me%row_ptr(i), me%row_ptr(i + 1) - 1
                j = me%col_ind(k)
                if (j == icol) then
                    col_vec(i) = me%values(k)
                    exit  ! found it in this row
                end if
            end do
        end do

    end subroutine sparse_csr_get_col

!*****************************************************************************************
!>
!  Set a column in the sparse matrix from a dense vector.
!  Only updates values at existing non-zero locations.
!  Note: Column access in CSR is inefficient.

    subroutine sparse_csr_set_col(me, icol, col_vec)

        class(sparse_matrix_csr), intent(inout) :: me
        integer, intent(in) :: icol  !! column index (1-based)
        real(wp), intent(in) :: col_vec(:)  !! dense column vector (nrows)

        integer :: i, k, j

        ! Search for column icol in each row and update
        do i = 1, me%nrows
            do k = me%row_ptr(i), me%row_ptr(i + 1) - 1
                j = me%col_ind(k)
                if (j == icol) then
                    me%values(k) = col_vec(i)
                    exit  ! found it in this row
                end if
            end do
        end do

    end subroutine sparse_csr_set_col

!*****************************************************************************************
!>
!  Get the number of non-zeros in a specific row.

    pure function sparse_csr_get_nnz_in_row(me, irow) result(nnz)

        class(sparse_matrix_csr), intent(in) :: me
        integer, intent(in) :: irow  !! row index (1-based)
        integer :: nnz

        nnz = me%row_ptr(irow + 1) - me%row_ptr(irow)

    end function sparse_csr_get_nnz_in_row

!*****************************************************************************************
!>
!  Copy the structure (pattern) from another sparse matrix.
!  Values are initialized to zero.

    subroutine sparse_csr_copy_structure(me, source)

        class(sparse_matrix_csr), intent(inout) :: me
        type(sparse_matrix_csr), intent(in) :: source

        call me%initialize(source%nrows, source%ncols, source%row_ptr, source%col_ind)

    end subroutine sparse_csr_copy_structure

!*****************************************************************************************
!*****************************************************************************************
! Sparse Pattern routines
!*****************************************************************************************
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a sparse pattern.

    subroutine sparse_pattern_init(me, nrows, ncols, row_ptr, col_ind)

        class(sparse_pattern), intent(inout) :: me
        integer, intent(in) :: nrows  !! number of rows
        integer, intent(in) :: ncols  !! number of columns
        integer, intent(in) :: row_ptr(nrows+1)  !! row start pointers
        integer, intent(in) :: col_ind(:)  !! column indices

        ! Clean up any existing data
        call me%destroy()

        me%nrows = nrows
        me%ncols = ncols
        me%nnz = row_ptr(nrows + 1) - 1

        ! Allocate and copy structure
        allocate(me%row_ptr(nrows + 1))
        allocate(me%col_ind(me%nnz))

        me%row_ptr = row_ptr
        me%col_ind = col_ind

    end subroutine sparse_pattern_init

!*****************************************************************************************
!>
!  Destroy/deallocate a sparse pattern.

    subroutine sparse_pattern_destroy(me)

        class(sparse_pattern), intent(inout) :: me

        me%nrows = 0
        me%ncols = 0
        me%nnz = 0

        if (allocated(me%row_ptr)) deallocate(me%row_ptr)
        if (allocated(me%col_ind)) deallocate(me%col_ind)

    end subroutine sparse_pattern_destroy

!*****************************************************************************************
!>
!  Check if sparse pattern is initialized.

    pure function sparse_pattern_is_initialized(me) result(initialized)

        class(sparse_pattern), intent(in) :: me
        logical :: initialized

        initialized = allocated(me%row_ptr) .and. allocated(me%col_ind)

    end function sparse_pattern_is_initialized

!*****************************************************************************************

end module psqp_sparse_module
!*****************************************************************************************
