!*******************************************************************************
!> author: Jacob Williams
!>
!> Test case for sparse Jacobian with nonlinear constraints.
!>
!> Problem: Minimize a quadratic objective with sparse nonlinear constraints.
!>
!> Minimize: f(x) = sum((x(i) - i)^2)  for i=1..6
!>
!> Subject to:
!>   x(1)^2 + x(2)   = 3    (constraint 1)
!>   x(2)^2 + x(3)   = 6    (constraint 2)
!>   x(3)^2 + x(4)   = 11   (constraint 3)
!>   x(4)^2 + x(5)   = 18   (constraint 4)
!>   x(5)^2 + x(6)   = 27   (constraint 5)
!>   x(i) >= 0  for all i   (simple bounds)
!>
!> The constraint Jacobian is sparse (bidiagonal structure):
!>   Row 1: [2*x(1),  1,      0,      0,      0,      0]
!>   Row 2: [0,       2*x(2), 1,      0,      0,      0]
!>   Row 3: [0,       0,      2*x(3), 1,      0,      0]
!>   Row 4: [0,       0,      0,      2*x(4), 1,      0]
!>   Row 5: [0,       0,      0,      0,      2*x(5), 1]
!>
!> Sparsity: 10 nonzeros out of 30 total (33% sparse)
!>
!> The solution minimizes the objective while satisfying the nonlinear constraints.

program test_sparse_nonlinear

    use psqp_module, only: psqp_class, wp => psqp_wp
    use psqp_sparse_module

    implicit none

    integer, parameter :: nf = 6   !! number of variables
    integer, parameter :: nc = 5   !! number of nonlinear constraints
    integer, parameter :: nnz = 10 !! number of nonzeros in Jacobian
    integer, parameter :: nb = 1   !! use simple bounds

    type(psqp_class) :: solver
    integer :: iterm

    ! Variables and bounds
    real(wp), dimension(nf) :: x, xl, xu
    integer, dimension(nf) :: ix

    ! Constraints and bounds
    real(wp), dimension(nc+1) :: cf
    real(wp), dimension(nc) :: cl, cu
    integer, dimension(nc) :: ic

    ! Parameters for psqpn
    integer, dimension(6) :: ipar
    real(wp), dimension(5) :: rpar
    integer, parameter :: iprnt = 1

    ! Sparse Jacobian in COO format
    integer, dimension(nnz) :: irow, jcol
    real(wp), dimension(nnz) :: vals

    ! Results
    real(wp) :: f, cmax, gmax

    ! Expected solution
    real(wp), dimension(nf) :: x_expected
    real(wp) :: tol
    integer :: i
    logical :: test_passed

    write(*,'(A)') ''
    write(*,'(A)') '======================================='
    write(*,'(A)') 'Sparse Nonlinear Jacobian Test'
    write(*,'(A)') '======================================='
    write(*,'(A)') ''

    ! Initialize variables
    x = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]  ! start at unconstrained optimum

    ! Simple bounds: x(i) >= 0
    xl = 0.0_wp
    xu = 20.0_wp
    ix = 3  ! two-sided bounds

    ! Nonlinear constraint bounds (all equality constraints)
    cl = [3.0_wp, 6.0_wp, 11.0_wp, 18.0_wp, 27.0_wp]
    cu = cl  ! equality constraints
    ic = 5   ! equality constraints

    ! Set parameters
    ipar = [1000, &  ! maximum number of iterations
            1000, &  ! maximum number of function evaluations
            0, &     ! not used
            0, &     ! not used
            1, &     ! variable metric update
            1]       ! correction of variable metric update

    rpar = [0.0_wp, &                ! maximum stepsize
            2*epsilon(1.0_wp), &     ! tolerance for change of variables
            2*epsilon(1.0_wp), &     ! tolerance for constraint violations
            0.0_wp, &                ! tolerance for gradient of Lagrangian
            0.0_wp]                  ! penalty coefficient

    ! Define sparse Jacobian sparsity pattern in COO (Coordinate) format
    ! Each constraint i has nonzeros in columns i and i+1
    ! Constraint 1: grad = [2*x(1), 1, 0, 0, 0, 0]
    irow(1) = 1; jcol(1) = 1
    irow(2) = 1; jcol(2) = 2

    ! Constraint 2: grad = [0, 2*x(2), 1, 0, 0, 0]
    irow(3) = 2; jcol(3) = 2
    irow(4) = 2; jcol(4) = 3

    ! Constraint 3: grad = [0, 0, 2*x(3), 1, 0, 0]
    irow(5) = 3; jcol(5) = 3
    irow(6) = 3; jcol(6) = 4

    ! Constraint 4: grad = [0, 0, 0, 2*x(4), 1, 0]
    irow(7) = 4; jcol(7) = 4
    irow(8) = 4; jcol(8) = 5

    ! Constraint 5: grad = [0, 0, 0, 0, 2*x(5), 1]
    irow(9) = 5; jcol(9) = 5
    irow(10) = 5; jcol(10) = 6

    write(*,'(A)') 'Problem setup:'
    write(*,'(A,I0)') '  Number of variables: ', nf
    write(*,'(A,I0)') '  Number of constraints: ', nc
    write(*,'(A,I0)') '  Jacobian nonzeros: ', nnz
    write(*,'(A,F6.2,A)') '  Sparsity: ', real(nnz, wp) / real(nf*nc, wp) * 100.0_wp, '%'
    write(*,'(A)') ''
    write(*,'(A)') 'Constraint structure:'
    write(*,'(A)') '  c(i) = x(i)^2 + x(i+1) = rhs(i)'
    write(*,'(A)') '  This creates a sparse bidiagonal Jacobian pattern'
    write(*,'(A)') ''

    ! Set up sparse mode with Jacobian pattern
    call solver%set_jacobian_pattern(nc, nf, irow, jcol, nnz)
    call solver%set_sparse_jacobian_callback(sparse_jac_func)

    write(*,'(A)') 'Solving with sparse nonlinear Jacobian mode...'

    ! Solve the problem
    call solver%psqpn(nf, nb, nc, x, ix, xl, xu, cf, ic, cl, cu, &
                      ipar, rpar, f, cmax, gmax, iprnt, iterm, &
                      obj_func, dobj_func, con_func, dcon_func)

    write(*,'(A)') ''
    write(*,'(A,I0)') 'Optimization complete. iterm: ', iterm
    write(*,'(A)') ''

    ! Display results
    write(*,'(A)') 'Solution:'
    do i = 1, nf
        write(*,'(A,I0,A,F10.6)') '  x(', i, ') = ', x(i)
    end do
    write(*,'(A)') ''
    write(*,'(A,F10.6)') 'Objective value: ', f
    write(*,'(A,E10.3)') 'Max constraint violation: ', cmax
    write(*,'(A,E10.3)') 'Max Lagrangian gradient: ', gmax
    write(*,'(A)') ''

    ! Verify constraints
    write(*,'(A)') 'Constraint verification:'
    write(*,'(A,F10.6,A,F10.6)') '  x(1)^2 + x(2) = ', x(1)**2 + x(2), '  (target: 3.0)'
    write(*,'(A,F10.6,A,F10.6)') '  x(2)^2 + x(3) = ', x(2)**2 + x(3), '  (target: 6.0)'
    write(*,'(A,F10.6,A,F10.6)') '  x(3)^2 + x(4) = ', x(3)**2 + x(4), '  (target: 11.0)'
    write(*,'(A,F10.6,A,F10.6)') '  x(4)^2 + x(5) = ', x(4)**2 + x(5), '  (target: 18.0)'
    write(*,'(A,F10.6,A,F10.6)') '  x(5)^2 + x(6) = ', x(5)**2 + x(6), '  (target: 27.0)'
    write(*,'(A)') ''

    ! Expected solution (approximately)
    ! The optimizer finds the solution that satisfies the nonlinear constraints
    ! while minimizing the distance to [1,2,3,4,5,6]
    x_expected = [1.089_wp, 1.814_wp, 2.709_wp, 3.662_wp, 4.587_wp, 5.959_wp]

    ! Check solution accuracy
    tol = 1.0e-2_wp  ! tolerance for nonlinear problem
    test_passed = .true.

    write(*,'(A)') 'Solution accuracy check:'
    do i = 1, nf
        write(*,'(A,I0,A,F10.6,A,F10.6,A,E10.3)') &
            '  x(', i, ') = ', x(i), '  expected: ', x_expected(i), &
            '  error: ', abs(x(i) - x_expected(i))
        if (abs(x(i) - x_expected(i)) > tol) test_passed = .false.
    end do
    write(*,'(A)') ''

    if (test_passed) then
        write(*,'(A)') '✓ Test PASSED - Solution matches expected values'
    else
        write(*,'(A)') '✗ Test FAILED - Solution does not match expected values'
        error stop 1
    end if

    write(*,'(A)') ''
    write(*,'(A)') '======================================='
    write(*,'(A)') ''

contains

    !***************************************************************************
    subroutine obj_func(me, nf, x, f)
        !! Objective function: f(x) = sum((x(i) - i)^2)

        implicit none
        class(psqp_class), intent(inout) :: me
        integer :: nf
        real(wp) :: x(nf)
        real(wp) :: f

        integer :: i

        f = 0.0_wp
        do i = 1, nf
            f = f + (x(i) - real(i, wp))**2
        end do

    end subroutine obj_func
    !***************************************************************************

    !***************************************************************************
    subroutine dobj_func(me, nf, x, g)
        !! Objective gradient: g(i) = 2*(x(i) - i)

        implicit none
        class(psqp_class), intent(inout) :: me
        integer :: nf
        real(wp) :: x(nf)
        real(wp) :: g(nf)

        integer :: i

        do i = 1, nf
            g(i) = 2.0_wp * (x(i) - real(i, wp))
        end do

    end subroutine dobj_func
    !***************************************************************************

    !***************************************************************************
    subroutine con_func(me, nf, kc, x, fc)
        !! Constraint function: c(i) = x(i)^2 + x(i+1)
        !! Computes constraint kc at point x

        implicit none
        class(psqp_class), intent(inout) :: me
        integer :: nf, kc
        real(wp) :: x(nf), fc

        if (kc >= 1 .and. kc <= nf-1) then
            fc = x(kc)**2 + x(kc+1)
        else
            error stop 'invalid constraint index'
        end if

    end subroutine con_func
    !***************************************************************************

    !***************************************************************************
    subroutine dcon_func(me, nf, kc, x, gc)
        !! Gradient of constraint kc (used when sparse Jacobian not provided)
        !! Not actually used in sparse mode, but required by interface

        implicit none
        class(psqp_class), intent(inout) :: me
        integer :: nf, kc
        real(wp) :: x(nf), gc(nf)

        error stop 'dcon_func should not be called in sparse mode'

    end subroutine dcon_func
    !***************************************************************************

    !***************************************************************************
    subroutine sparse_jac_func(me, nf, nc1, x, values)
        !! Sparse Jacobian callback - returns nonzero values
        !! The sparsity pattern was already provided via set_jacobian_pattern
        !!
        !! For constraint i: c(i) = x(i)^2 + x(i+1)
        !! Gradient: dc/dx(i) = 2*x(i), dc/dx(i+1) = 1
        !!
        !! The values array follows the same order as the COO input:
        !! values(1:2)   = gradient of constraint 1
        !! values(3:4)   = gradient of constraint 2
        !! values(5:6)   = gradient of constraint 3
        !! values(7:8)   = gradient of constraint 4
        !! values(9:10)  = gradient of constraint 5

        implicit none
        class(psqp_class), intent(inout) :: me
        integer, intent(in) :: nf, nc1
        real(wp), intent(in) :: x(nf)
        real(wp), intent(out) :: values(:)

        integer :: i, idx

        ! Fill in Jacobian values for each constraint
        idx = 1
        do i = 1, nc1
            values(idx) = 2.0_wp * x(i)  ! d/dx(i) of x(i)^2 + x(i+1)
            values(idx+1) = 1.0_wp       ! d/dx(i+1) of x(i)^2 + x(i+1)
            idx = idx + 2
        end do

    end subroutine sparse_jac_func
    !***************************************************************************

end program test_sparse_nonlinear
