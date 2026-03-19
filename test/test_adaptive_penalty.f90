!***********************************************************************
!>
!  Test case demonstrating adaptive penalty parameter updates.
!
!  This test uses a deliberately poor initial penalty parameter (rpf)
!  to show how the adaptive mechanism improves convergence.
!
!  Problem: Minimize sum(x(i)^2) subject to:
!    c1: x(1)^2 + x(2)^2 = 25
!    c2: x(3)^2 + x(4)^2 = 16
!    c3: x(1) + x(2) + x(3) + x(4) = 0
!
!  Start with poor initial guess and very small rpf to trigger adaptation.

program test_adaptive_penalty

   use psqp_module, only: psqp_class, wp => psqp_wp

   implicit none

   type(psqp_class) :: solver
   integer, parameter :: nf = 4
   integer, parameter :: nc = 3
   real(wp), dimension(nf) :: x, xl, xu
   real(wp), dimension(nc) :: cl, cu
   real(wp), dimension(nc+1) :: cf  ! Constraint function values
   integer, dimension(nf) :: ix
   integer, dimension(nc) :: ic
   integer, dimension(6) :: ipar
   real(wp), dimension(5) :: rpar
   integer, parameter :: iprnt = 2  ! Print iterations
   integer :: iterm
   real(wp) :: f, cmax, gmax
   integer :: i
   logical :: test_passed
   real(wp) :: rpf_initial

   write(*,'(A)') ''
   write(*,'(A)') '======================================='
   write(*,'(A)') 'Adaptive Penalty Parameter Test'
   write(*,'(A)') '======================================='
   write(*,'(A)') ''

   ! Problem setup
   write(*,'(A)') 'Problem:'
   write(*,'(A)') '  Minimize: (x(1)-3)^2 + (x(2)-2)^2 + (x(3)+1)^2 + (x(4)-1)^2'
   write(*,'(A)') '  Subject to:'
   write(*,'(A)') '    x(1) + 2*x(2) = 5'
   write(*,'(A)') '    x(3) + x(4) = 2'
   write(*,'(A)') '    x(1)^2 + x(3)^2 <= 10'
   write(*,'(A)') ''

   ! Initial guess
   x(1) = 0.0_wp
   x(2) = 0.0_wp
   x(3) = 0.0_wp
   x(4) = 0.0_wp

   write(*,'(A)') 'Initial guess:'
   write(*,'(A,4F10.3)') '  x = ', x
   write(*,'(A)') ''

   ! No bounds
   ix(:) = 0
   xl(:) = -1.0e10_wp
   xu(:) = 1.0e10_wp

   ! Constraints
   ic(1) = 5  ! x(1) + 2*x(2) = 5 (equality)
   ic(2) = 5  ! x(3) + x(4) = 2 (equality)
   ic(3) = 2  ! x(1)^2 + x(3)^2 <= 10 (inequality)
   cl(1) = 5.0_wp
   cl(2) = 2.0_wp
   cl(3) = -1.0e10_wp
   cu(1) = 5.0_wp
   cu(2) = 2.0_wp
   cu(3) = 10.0_wp

   ! Deliberately poor initial penalty parameter (way too small)
   rpf_initial = 1.0e-6_wp
   write(*,'(A,1PE12.3)') 'Initial penalty parameter (deliberately small): rpf =', rpf_initial
   write(*,'(A)') 'Adaptive mode enabled - rpf will increase as needed'
   write(*,'(A)') ''

   ! Set parameters for psqpn
   ipar = [100, &    ! maximum number of iterations
           1000, &   ! maximum number of function evaluations
           0, &      ! not used
           0, &      ! not used
           1, &      ! variable metric update (BFGS)
           2]        ! correction of variable metric update (Powell)

   rpar = [1.0e3_wp, &             ! maximum stepsize
           1.0e-6_wp, &            ! tolerance for change of variables
           1.0e-6_wp, &            ! tolerance for constraint violations
           1.0e-6_wp, &            ! tolerance for Lagrangian gradient
           rpf_initial]            ! penalty coefficient (poor initial value)

   ! Configure adaptive penalty parameters
   solver%rpf_adaptive = .true.
   solver%rpf_min = 1.0e-8_wp
   solver%rpf_max = 1.0e6_wp
   solver%rpf_increase_factor = 10.0_wp
   solver%rpf_stagnation_threshold = 0.95_wp
   solver%rpf_stagnation_limit = 2

   write(*,'(A)') 'Solving...'
   write(*,'(A)') ''

   ! Solve (callbacks passed as arguments)
   call solver%psqpn(nf, 0, nc, x, ix, xl, xu, cf, ic, cl, cu, ipar, rpar, &
                     f, gmax, cmax, iprnt, iterm, obj_func, dobj_func, con_func, dcon_func)

   write(*,'(A)') ''
   write(*,'(A,I4)') 'Optimization complete. iterm:', iterm
   write(*,'(A)') ''

   ! Report final penalty parameter
   write(*,'(A,1PE12.3)') 'Final penalty parameter: rpf =', rpar(5)
   write(*,'(A,1PE12.3)') 'RPF increase factor:', rpar(5) / rpf_initial
   write(*,'(A)') ''

   ! Display solution
   write(*,'(A)') 'Solution:'
   do i = 1, nf
      write(*,'(A,I0,A,F12.6)') '  x(', i, ') = ', x(i)
   end do
   write(*,'(A)') ''

   write(*,'(A,1PE12.3)') 'Objective value:', f
   write(*,'(A,1PE12.3)') 'Max constraint violation:', cmax
   write(*,'(A,1PE12.3)') 'Max Lagrangian gradient:', gmax
   write(*,'(A)') ''

   ! Verify constraints
   write(*,'(A)') 'Constraint verification:'
   write(*,'(A,F12.6,A)') '  x(1) + 2*x(2) =', x(1) + 2.0_wp*x(2), '  (target: 5.0)'
   write(*,'(A,F12.6,A)') '  x(3) + x(4) =', x(3) + x(4), '  (target: 2.0)'
   write(*,'(A,F12.6,A)') '  x(1)^2 + x(3)^2 =', x(1)**2 + x(3)**2, '  (must be <= 10.0)'
   write(*,'(A)') ''

   ! Check convergence
   test_passed = .true.
   if (abs(x(1) + 2.0_wp*x(2) - 5.0_wp) > 1.0e-3_wp) test_passed = .false.
   if (abs(x(3) + x(4) - 2.0_wp) > 1.0e-3_wp) test_passed = .false.
   if (x(1)**2 + x(3)**2 > 10.0_wp + 1.0e-3_wp) test_passed = .false.
   if (cmax > 1.0e-4_wp) test_passed = .false.
   if (iterm < 0 .and. iterm /= -6) test_passed = .false.

   if (test_passed) then
      write(*,'(A)') '✓ Test PASSED - Adaptive penalty enabled successful convergence'
   else
      write(*,'(A)') '✗ Test FAILED - Did not converge to acceptable solution'
   end if
   write(*,'(A)') ''
   write(*,'(A)') '======================================='
   write(*,'(A)') ''

   if (.not. test_passed) error stop 1

contains

   subroutine obj_func(me, nf, x, ff)
      class(psqp_class), intent(inout) :: me
      integer :: nf
      real(wp) :: x(nf)
      real(wp) :: ff
      ! Minimize: (x(1)-3)^2 + (x(2)-2)^2 + (x(3)+1)^2 + (x(4)-1)^2
      ff = (x(1) - 3.0_wp)**2 + (x(2) - 2.0_wp)**2 + (x(3) + 1.0_wp)**2 + (x(4) - 1.0_wp)**2
   end subroutine obj_func

   subroutine dobj_func(me, nf, x, gf)
      class(psqp_class), intent(inout) :: me
      integer :: nf
      real(wp) :: x(nf)
      real(wp) :: gf(nf)
      gf(1) = 2.0_wp * (x(1) - 3.0_wp)
      gf(2) = 2.0_wp * (x(2) - 2.0_wp)
      gf(3) = 2.0_wp * (x(3) + 1.0_wp)
      gf(4) = 2.0_wp * (x(4) - 1.0_wp)
   end subroutine dobj_func

   subroutine con_func(me, nf, kc, x, fc)
      class(psqp_class), intent(inout) :: me
      integer :: nf, kc
      real(wp) :: x(nf)
      real(wp) :: fc
      select case(kc)
      case(1)
         ! x(1) + 2*x(2) = 5
         fc = x(1) + 2.0_wp*x(2)
      case(2)
         ! x(3) + x(4) = 2
         fc = x(3) + x(4)
      case(3)
         ! x(1)^2 + x(3)^2 <= 10
         fc = x(1)**2 + x(3)**2
      end select
   end subroutine con_func

   subroutine dcon_func(me, nf, kc, x, gc)
      class(psqp_class), intent(inout) :: me
      integer :: nf, kc
      real(wp) :: x(nf)
      real(wp) :: gc(nf)
      gc = 0.0_wp
      select case(kc)
      case(1)
         ! grad(x(1) + 2*x(2))
         gc(1) = 1.0_wp
         gc(2) = 2.0_wp
      case(2)
         ! grad(x(3) + x(4))
         gc(3) = 1.0_wp
         gc(4) = 1.0_wp
      case(3)
         ! grad(x(1)^2 + x(3)^2)
         gc(1) = 2.0_wp * x(1)
         gc(3) = 2.0_wp * x(3)
      end select
   end subroutine dcon_func

end program test_adaptive_penalty
