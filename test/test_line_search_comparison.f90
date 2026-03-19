!***********************************************************************
!>
!  Test case comparing Extended vs Wolfe line search methods.
!
!  This test solves the same optimization problem using both line search
!  methods to demonstrate the potential efficiency gains from using
!  directional derivatives (Wolfe strong conditions).
!
!  Problem: Rosenbrock function
!    Minimize: (1-x1)^2 + 100*(x2-x1^2)^2
!    Subject to: x1^2 + x2^2 <= 2
!
!  The Wolfe line search should typically require fewer function
!  evaluations due to better step length selection.

program test_line_search_comparison

   use psqp_module, only: psqp_class, wp => psqp_wp

   implicit none

   type(psqp_class) :: solver_extended, solver_wolfe
   integer, parameter :: nf = 2
   integer, parameter :: nc = 1
   real(wp), dimension(nf) :: x_ext, x_wolfe, xl, xu
   real(wp), dimension(nc) :: cl, cu
   real(wp), dimension(nc+1) :: cf_ext, cf_wolfe
   integer, dimension(nf) :: ix
   integer, dimension(nc) :: ic
   integer, dimension(6) :: ipar
   real(wp), dimension(5) :: rpar
   integer, parameter :: iprnt = 1  ! Print final results only
   integer :: iterm_ext, iterm_wolfe
   real(wp) :: f_ext, f_wolfe, cmax_ext, cmax_wolfe, gmax_ext, gmax_wolfe
   integer :: i

   write(*,'(A)') ''
   write(*,'(A)') '======================================='
   write(*,'(A)') 'Line Search Method Comparison Test'
   write(*,'(A)') '======================================='
   write(*,'(A)') ''

   write(*,'(A)') 'Problem: Rosenbrock function'
   write(*,'(A)') '  Minimize: (1-x1)^2 + 100*(x2-x1^2)^2'
   write(*,'(A)') '  Subject to: x1^2 + x2^2 <= 2'
   write(*,'(A)') ''

   ! Initial guess
   x_ext(1) = -1.0_wp
   x_ext(2) = 0.5_wp
   x_wolfe = x_ext  ! Same starting point for both

   write(*,'(A,2F10.3)') 'Initial guess: x = ', x_ext
   write(*,'(A)') ''

   ! No bounds
   ix(:) = 0
   xl(:) = -10.0_wp
   xu(:) = 10.0_wp

   ! Circle constraint: x1^2 + x2^2 <= 2
   ic(1) = 2  ! upper bound constraint
   cl(1) = -1.0e10_wp
   cu(1) = 2.0_wp

   ! Optimization parameters
   ipar = [100, &    ! maximum iterations
           1000, &   ! maximum function evaluations
           0, 0, &   ! unused
           1, &      ! BFGS
           2]        ! Powell correction

   rpar = [1.0e3_wp, &    ! max stepsize
           1.0e-8_wp, &   ! tol x
           1.0e-8_wp, &   ! tol c
           1.0e-8_wp, &   ! tol g
           1.0e-4_wp]     ! penalty coefficient

   ! ==========================================
   ! Test 1: Extended line search (original)
   ! ==========================================
   write(*,'(A)') 'Test 1: Extended Line Search (no directional derivatives)'
   write(*,'(A)') '---------------------------------------------------------------'

   call solver_extended%psqpn(nf, 0, nc, x_ext, ix, xl, xu, cf_ext, ic, cl, cu, &
                              ipar, rpar, f_ext, gmax_ext, cmax_ext, iprnt, iterm_ext, &
                              obj_rosenbrock, dobj_rosenbrock, con_circle, dcon_circle, &
                              line_search_method=1)

   write(*,'(A,I4)') 'Termination code: ', iterm_ext
   write(*,'(A,2F12.6)') 'Solution: x = ', x_ext
   write(*,'(A,F12.6)') 'Objective value: ', f_ext
   write(*,'(A,1PE12.3)') 'Constraint violation: ', cmax_ext
   write(*,'(A,I6)') 'Function evaluations: ', solver_extended%nfv
   write(*,'(A,I6)') 'Gradient evaluations: ', solver_extended%nfg
   write(*,'(A,I6)') 'Iterations: ', solver_extended%nit
   write(*,'(A)') ''

   ! ==========================================
   ! Test 2: Wolfe line search (with derivatives)
   ! ==========================================
   write(*,'(A)') 'Test 2: Wolfe Line Search (strong Wolfe conditions)'
   write(*,'(A)') '---------------------------------------------------------------'

   call solver_wolfe%psqpn(nf, 0, nc, x_wolfe, ix, xl, xu, cf_wolfe, ic, cl, cu, &
                           ipar, rpar, f_wolfe, gmax_wolfe, cmax_wolfe, iprnt, iterm_wolfe, &
                           obj_rosenbrock, dobj_rosenbrock, con_circle, dcon_circle, &
                           line_search_method=2)  ! Use default c1=1e-4, c2=0.1, max_iter=20

   write(*,'(A,I4)') 'Termination code: ', iterm_wolfe
   write(*,'(A,2F12.6)') 'Solution: x = ', x_wolfe
   write(*,'(A,F12.6)') 'Objective value: ', f_wolfe
   write(*,'(A,1PE12.3)') 'Constraint violation: ', cmax_wolfe
   write(*,'(A,I6)') 'Function evaluations: ', solver_wolfe%nfv
   write(*,'(A,I6)') 'Gradient evaluations: ', solver_wolfe%nfg
   write(*,'(A,I6)') 'Iterations: ', solver_wolfe%nit
   write(*,'(A)') ''

   ! ==========================================
   ! Comparison
   ! ==========================================
   write(*,'(A)') '======================================='
   write(*,'(A)') 'Comparison Summary'
   write(*,'(A)') '======================================='
   write(*,'(A)') ''
   write(*,'(A,I6,A,I6)') 'Function evals:  Extended =', solver_extended%nfv, &
          '  Wolfe =', solver_wolfe%nfv
   write(*,'(A,I6,A,I6)') 'Gradient evals:  Extended =', solver_extended%nfg, &
          '  Wolfe =', solver_wolfe%nfg
   write(*,'(A,I6,A,I6)') 'Iterations:      Extended =', solver_extended%nit, &
          '  Wolfe =', solver_wolfe%nit
   write(*,'(A)') ''

   if (solver_wolfe%nfv < solver_extended%nfv) then
      write(*,'(A,F6.1,A)') '✓ Wolfe line search used ', &
         100.0_wp*(1.0_wp - real(solver_wolfe%nfv,wp)/real(solver_extended%nfv,wp)), &
         '% fewer function evaluations'
   else if (solver_wolfe%nfv == solver_extended%nfv) then
      write(*,'(A)') '≈ Both methods used same number of function evaluations'
   else
      write(*,'(A)') '  Extended actually used fewer evaluations for this problem'
   end if

   ! Check both converged to similar solutions
   if (abs(f_ext - f_wolfe) < 1.0e-4_wp .and. &
       iterm_ext > 0 .and. iterm_wolfe > 0) then
      write(*,'(A)') '✓ Both methods converged to similar solutions'
   end if

   write(*,'(A)') ''
   write(*,'(A)') '======================================='
   write(*,'(A)') ''

contains

   subroutine obj_rosenbrock(me, nf, x, ff)
      class(psqp_class), intent(inout) :: me
      integer :: nf
      real(wp) :: x(nf)
      real(wp) :: ff
      ! Rosenbrock function: (1-x1)^2 + 100*(x2-x1^2)^2
      ff = (1.0_wp - x(1))**2 + 100.0_wp*(x(2) - x(1)**2)**2
   end subroutine obj_rosenbrock

   subroutine dobj_rosenbrock(me, nf, x, gf)
      class(psqp_class), intent(inout) :: me
      integer :: nf
      real(wp) :: x(nf)
      real(wp) :: gf(nf)
      ! Gradient of Rosenbrock
      gf(1) = -2.0_wp*(1.0_wp - x(1)) - 400.0_wp*x(1)*(x(2) - x(1)**2)
      gf(2) = 200.0_wp*(x(2) - x(1)**2)
   end subroutine dobj_rosenbrock

   subroutine con_circle(me, nf, kc, x, fc)
      class(psqp_class), intent(inout) :: me
      integer :: nf, kc
      real(wp) :: x(nf)
      real(wp) :: fc
      ! Constraint: x1^2 + x2^2 <= 2
      fc = x(1)**2 + x(2)**2
   end subroutine con_circle

   subroutine dcon_circle(me, nf, kc, x, gc)
      class(psqp_class), intent(inout) :: me
      integer :: nf, kc
      real(wp) :: x(nf)
      real(wp) :: gc(nf)
      ! Gradient: [2*x1, 2*x2]
      gc(1) = 2.0_wp * x(1)
      gc(2) = 2.0_wp * x(2)
   end subroutine dcon_circle

end program test_line_search_comparison
