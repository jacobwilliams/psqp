!*******************************************************************************
!> author: Jacob Williams
!>
!> Test case for sparse Jacobian optimization.
!>
!> Problem: Minimize a quadratic objective with sparse linear constraints.
!>
!> Minimize: f(x) = 0.5 * sum((x(i) - i)^2)  for i=1..n
!>
!> Subject to:
!>   x(1) + x(2) = 3           (constraint 1)
!>   x(2) + x(3) = 5           (constraint 2)
!>   x(3) + x(4) = 7           (constraint 3)
!>   x(4) + x(5) = 9           (constraint 4)
!>   x(i) >= 0  for all i      (simple bounds)
!>
!> The constraint Jacobian is sparse (tridiagonal band structure):
!>   [1 1 0 0 0]
!>   [0 1 1 0 0]
!>   [0 0 1 1 0]
!>   [0 0 0 1 1]
!>
!> Known solution: The constraints form a chain, and the objective function
!> wants x(i) = i. This is compatible with all constraints, giving:
!>   x* = [1, 2, 3, 4, 5]  with f* = 0

program test_sparse_optimization

    use psqp_module, only: psqp_class, wp => psqp_wp
    use psqp_sparse_module

    implicit none

    integer, parameter :: nf = 5   !! number of variables
    integer, parameter :: nc = 4   !! number of linear constraints
    integer, parameter :: nnz = 8  !! number of nonzeros in Jacobian
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
    write(*,'(A)') 'Sparse Jacobian Optimization Test'
    write(*,'(A)') '======================================='
    write(*,'(A)') ''

    ! Initialize variables at origin
    x = 0.0_wp

    ! Simple bounds: x(i) >= 0
    xl = 0.0_wp
    xu = 10.0_wp  ! upper bounds (not binding)
    ix = 3  ! two-sided bounds

    ! Linear constraint bounds (all equality constraints)
    cl = [3.0_wp, 5.0_wp, 7.0_wp, 9.0_wp]
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

    ! Define sparse Jacobian in COO (Coordinate) format
    ! Constraint 1: x(1) + x(2) = 3
    irow(1) = 1; jcol(1) = 1; vals(1) = 1.0_wp
    irow(2) = 1; jcol(2) = 2; vals(2) = 1.0_wp

    ! Constraint 2: x(2) + x(3) = 5
    irow(3) = 2; jcol(3) = 2; vals(3) = 1.0_wp
    irow(4) = 2; jcol(4) = 3; vals(4) = 1.0_wp

    ! Constraint 3: x(3) + x(4) = 7
    irow(5) = 3; jcol(5) = 3; vals(5) = 1.0_wp
    irow(6) = 3; jcol(6) = 4; vals(6) = 1.0_wp

    ! Constraint 4: x(4) + x(5) = 9
    irow(7) = 4; jcol(7) = 4; vals(7) = 1.0_wp
    irow(8) = 4; jcol(8) = 5; vals(8) = 1.0_wp

    write(*,'(A)') 'Problem setup:'
    write(*,'(A,I0)') '  Number of variables: ', nf
    write(*,'(A,I0)') '  Number of constraints: ', nc
    write(*,'(A,I0)') '  Jacobian nonzeros: ', nnz
    write(*,'(A,F6.2,A)') '  Sparsity: ', real(nnz, wp) / real(nf*nc, wp) * 100.0_wp, '%'
    write(*,'(A)') ''

    ! Set up sparse mode with Jacobian pattern
    call solver%set_jacobian_pattern(nc, nf, irow, jcol, nnz)
    call solver%set_sparse_jacobian_callback(sparse_jac_func)

    write(*,'(A)') 'Solving with sparse Jacobian mode...'

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
    write(*,'(A,F10.6,A,F10.6)') '  x(1) + x(2) = ', x(1) + x(2), '  (target: 3.0)'
    write(*,'(A,F10.6,A,F10.6)') '  x(2) + x(3) = ', x(2) + x(3), '  (target: 5.0)'
    write(*,'(A,F10.6,A,F10.6)') '  x(3) + x(4) = ', x(3) + x(4), '  (target: 7.0)'
    write(*,'(A,F10.6,A,F10.6)') '  x(4) + x(5) = ', x(4) + x(5), '  (target: 9.0)'
    write(*,'(A)') ''

    ! Expected solution (computed analytically from KKT conditions)
    ! The solution satisfies the constraints and minimizes the objective
    ! Since f(x) = 0.5 * sum((x(i) - i)^2), the minimum is at x(i) = i
    ! And x = [1,2,3,4,5] satisfies all constraints:
    !   1 + 2 = 3 ✓
    !   2 + 3 = 5 ✓
    !   3 + 4 = 7 ✓
    !   4 + 5 = 9 ✓
    x_expected = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]

    ! Check solution accuracy
    tol = 1.0e-4_wp
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
        !! Objective function: f(x) = 0.5 * sum((x(i) - i)^2)

        implicit none
        class(psqp_class), intent(inout) :: me
        integer :: nf
        real(wp) :: x(nf)
        real(wp) :: f

        integer :: i

        f = 0.0_wp
        do i = 1, nf
            f = f + 0.5_wp * (x(i) - real(i, wp))**2
        end do

    end subroutine obj_func
    !***************************************************************************

    !***************************************************************************
    subroutine dobj_func(me, nf, x, g)
        !! Objective gradient: g(i) = x(i) - i

        implicit none
        class(psqp_class), intent(inout) :: me
        integer :: nf
        real(wp) :: x(nf)
        real(wp) :: g(nf)

        integer :: i

        do i = 1, nf
            g(i) = x(i) - real(i, wp)
        end do

    end subroutine dobj_func
    !***************************************************************************

    !***************************************************************************
    subroutine con_func(me, nf, kc, x, fc)
        !! Constraint function (linear constraints)
        !! Computes constraint kc at point x

        implicit none
        class(psqp_class), intent(inout) :: me
        integer :: nf, kc
        real(wp) :: x(nf), fc

        select case (kc)
        case(1); fc = x(1) + x(2)
        case(2); fc = x(2) + x(3)
        case(3); fc = x(3) + x(4)
        case(4); fc = x(4) + x(5)
        case default
            error stop 'invalid constraint index'
        end select

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

        error stop 'not used in sparse mode'

    end subroutine dcon_func
    !***************************************************************************

    !***************************************************************************
    subroutine sparse_jac_func(me, nf, nc1, x, values)
        !! Sparse Jacobian callback - returns nonzero values
        !! The sparsity pattern was already provided via set_jacobian_pattern

        implicit none
        class(psqp_class), intent(inout) :: me
        integer, intent(in) :: nf, nc1
        real(wp), intent(in) :: x(nf)
        real(wp), intent(out) :: values(:)

        ! For this linear problem, the Jacobian is constant
        ! All elements are 1.0
        values = 1.0_wp

    end subroutine sparse_jac_func
    !***************************************************************************

end program test_sparse_optimization
