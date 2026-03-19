!***********************************************************************
!>
!  PSQP: SQP variable metric method for general
!  nonlinear programming problems.
!
!### History
!  * [Original Fortran 77 code](http://www.cs.cas.cz/~luksan/subroutines.html)
!    by Ladislav Luksan.
!  * Jacob Williams, Aug 2017,
!    Significant refactoring to modern Fortran.

module psqp_module

   use psqp_matrix_module
   use psqp_sparse_module
   use psqp_kind_module, only: wp => psqp_wp

   implicit none

   private

   integer,parameter,public :: psqp_wp = wp !! export the working precision

   type, public :: psqp_class

      !! The main PSQP class to use.

      private

      ! these were formerly in the `stat` common block:
      integer, public :: nres = 0 !! number of restarts.
      integer, public :: ndec = 0 !! number of matrix decomposition.
      integer, public :: nrem = 0 !! number of constraint deletions.
      integer, public :: nadd = 0 !! number of constraint additions.
      integer, public :: nit = 0 !! number of iterations.
      integer, public :: nfv = 0 !! number of function evaluations.
      integer, public :: nfg = 0 !! number of gradient evaluations.
      integer, public :: nfh = 0 !! number of hessian evaluations.

      ! formerly saved variables in extended_line_search:
      integer  :: mtyp = 0
      integer  :: mode = 0
      integer  :: mes1 = 0
      integer  :: mes2 = 0
      real(wp) :: rl = 0.0_wp
      real(wp) :: fl = 0.0_wp
      real(wp) :: ru = 0.0_wp
      real(wp) :: fu = 0.0_wp
      real(wp) :: ri = 0.0_wp
      real(wp) :: fi = 0.0_wp

      procedure(obj_func), pointer  :: obj  => null() !! objective function
      procedure(dobj_func), pointer :: dobj => null() !! gradient of the objective function
      procedure(con_func), pointer  :: con  => null() !! constraint function
      procedure(dcon_func), pointer :: dcon => null() !! gradient of the constraint function

      procedure(report_f),pointer :: report => null() !! iteration report function

      ! Sparse Jacobian support
      logical :: use_sparse = .false. !! Use sparse Jacobian mode
      type(sparse_matrix_csr), allocatable :: jac_sparse !! Sparse constraint Jacobian matrix
      procedure(sparse_jac_func), pointer :: sparse_jac => null() !! Sparse Jacobian callback

      ! Adaptive penalty parameter control (set via psqpn optional arguments)
      logical :: rpf_adaptive = .true. !! Enable adaptive penalty parameter updates
      real(wp) :: rpf_min = 1.0e-6_wp !! Minimum penalty parameter value
      real(wp) :: rpf_max = 1.0e6_wp !! Maximum penalty parameter value
      real(wp) :: rpf_increase_factor = 10.0_wp !! Factor to increase rpf when stagnating
      real(wp) :: rpf_stagnation_threshold = 0.9_wp !! Threshold for detecting constraint stagnation
      integer :: rpf_stagnation_limit = 3 !! Iterations of stagnation before increasing rpf

      ! Internal state for adaptive penalty updates
      real(wp) :: cmax_previous = huge(1.0_wp) !! Previous major iteration's constraint violation
      integer :: rpf_stagnation_count = 0 !! Counter for consecutive stagnations

   contains

      private

      procedure, public :: psqpn
      procedure, public :: psqp

      procedure :: compute_obj_and_dobj
      procedure :: dual_range_space_quad_prog
      procedure :: ops_after_constr_deletion
      procedure :: compute_con_and_dcon
      procedure :: extended_line_search

      ! Sparse Jacobian helper methods
      procedure, public :: set_jacobian_pattern
      procedure, public :: set_sparse_jacobian_callback
      procedure :: get_constraint_gradient
      procedure :: constraint_gradient_dot
      procedure :: set_constraint_gradient

   end type psqp_class

   abstract interface

      subroutine obj_func(me, nf, x, ff)
         !! objective function interface
         import :: wp, psqp_class
         implicit none
         class(psqp_class), intent(inout) :: me
         integer  :: nf    !! the number of variables
         real(wp) :: x(nf) !! a vector of variables
         real(wp) :: ff    !! the value of the objective function
      end subroutine obj_func

      subroutine dobj_func(me, nf, x, gf)
         !! gradient of the objective function interface
         import :: wp, psqp_class
         implicit none
         class(psqp_class), intent(inout) :: me
         integer  :: nf     !! the number of variables
         real(wp) :: x(nf)  !! a vector of variables
         real(wp) :: gf(nf) !! the gradient of the objective function
      end subroutine dobj_func

      subroutine con_func(me, nf, kc, x, fc)
         !! constraint function interface
         import :: wp, psqp_class
         implicit none
         class(psqp_class), intent(inout) :: me
         integer  :: nf      !! the number of variables
         integer  :: kc      !! the index of the constraint function
         real(wp) :: x(nf)   !! a vector of variables
         real(wp) :: fc      !! the value of the constraint function
      end subroutine con_func

      subroutine dcon_func(me, nf, kc, x, gc)
         !! gradient of the constraint function interface
         import :: wp, psqp_class
         implicit none
         class(psqp_class), intent(inout) :: me
         integer  :: nf     !! the number of variables
         integer  :: kc     !! the index of the constraint function
         real(wp) :: x(nf)  !! a vector of variables and
         real(wp) :: gc(nf) !! the gradient of the constraint function
      end subroutine dcon_func

      subroutine report_f(me,iter,x,j,f)
         !! Report function to call once per iteration to report the solution.
         import :: wp, psqp_class
         implicit none
         class(psqp_class), intent(inout) :: me
         integer, intent(in) :: iter           !! Iteration number
         real(wp),dimension(:),intent(in) :: x !! optimization variables
         real(wp),intent(in) :: j              !! Objective function value
         real(wp),dimension(:),intent(in) :: f !! Constraint functions
      end subroutine report_f

      subroutine sparse_jac_func(me, nf, nc, x, values)
         !! Sparse Jacobian evaluation interface.
         !! User fills in the non-zero values of the constraint Jacobian.
         !! The sparsity pattern must be set beforehand via set_jacobian_pattern.
         import :: wp, psqp_class
         implicit none
         class(psqp_class), intent(inout) :: me
         integer, intent(in) :: nf        !! number of variables
         integer, intent(in) :: nc        !! number of constraints
         real(wp), intent(in) :: x(nf)    !! current variable values
         real(wp), intent(out) :: values(:) !! non-zero Jacobian values (size = nnz)
      end subroutine sparse_jac_func

   end interface

contains
!***********************************************************************

!***********************************************************************
!> date: 97/01/22
!
! easy to use subroutine for general nonlinear programming problems.

   subroutine psqpn(me, nf, nb, nc, x, bound_type, xl, xu, cf, constraint_type, &
                    cl, cu, ipar, rpar, f, gmax, &
                    cmax, iprnt, iterm, obj, dobj, con, dcon, report, &
                    rpf_adaptive, rpf_min, rpf_max, rpf_increase_factor, &
                    rpf_stagnation_threshold, rpf_stagnation_limit)

      class(psqp_class), intent(inout) :: me
      integer, intent(in) :: nf  !! number of variables
      integer, intent(in) :: nb  !! choice of simple bounds.
                                 !!
                                 !! * `nb=0` -- simple bounds suppressed.
                                 !! * `nb>0` -- simple bounds accepted.
      integer, intent(in) :: nc  !! number of general nonlinear constraints.
      real(wp),dimension(nf),intent(inout) :: x !! x(nf) vector of variables.
      integer, dimension(nf), intent(in) :: bound_type !! ix(nf) vector containing types of bounds.
                                                       !!
                                                       !! * `ix(i) = 0` -- variable x(i) is unbounded.
                                                       !! * `ix(i) = 1` -- lower bound xl(i) <= x(i).
                                                       !! * `ix(i) = 2` -- upper bound x(i) <= xu(i).
                                                       !! * `ix(i) = 3` -- two side bound xl(i) <= x(i) <= xu(i).
                                                       !! * `ix(i) = 5` -- variable x(i) is fixed.
      real(wp),dimension(nf),intent(in) :: xl   !! xl(nf) vector containing lower bounds for variables.
      real(wp),dimension(nf),intent(in) :: xu   !! xu(nf) vector containing upper bounds for variables.
      real(wp),dimension(nc+1),intent(out) :: cf !! cf(nc+1) vector containing values of the constraint functions.
      integer,dimension(nc),intent(in) :: constraint_type  !! ic(nc) vector containing types of constraints:
                                                           !!
                                                           !! * `ic(kc) = 0` -- constraint cf(kc) is not used.
                                                           !! * `ic(kc) = 1` -- lower constraint cl(kc) <= cf(kc).
                                                           !! * `ic(kc) = 2` -- upper constraint cf(kc) <= cu(kc).
                                                           !! * `ic(kc) = 3` -- two side constraint cl(kc) <= cf(kc) <= cu(kc).
                                                           !! * `ic(kc) = 5` -- equality constraint cf(kc) == cl(kc).
      real(wp),dimension(nc),intent(in) :: cl  !! cl(nc) vector containing lower bounds for constraint functions.
      real(wp),dimension(nc),intent(in) :: cu  !! cu(nc) vector containing upper bounds for constraint functions.
      integer,dimension(6),intent(in) :: ipar  !! integer paremeters:
                                               !!
                                               !! * `ipar(1)`  maximum number of iterations.
                                               !! * `ipar(2)`  maximum number of function evaluations.
                                               !! * `ipar(3)`  this parameter is not used in the subroutine psqp.
                                               !! * `ipar(4)`  this parameter is not used in the subroutine psqp.
                                               !! * `ipar(5)`  variable metric update used.
                                               !!   `ipar(5)=1` - the bfgs update.
                                               !!   `ipar(5)=2` - the hoshino update.
                                               !! * `ipar(6)`  correction of the variable metric update if a negative
                                               !!   curvature occurs.
                                               !!   `ipar(6)=1` - no correction.
                                               !!   `ipar(6)=2` - powell's correction.
      real(wp),dimension(5),intent(in) :: rpar   !! real parameters:
                                                 !!
                                                 !! * `rpar(1)` -- maximum stepsize.
                                                 !! * `rpar(2)` -- tolerance for change of variables.
                                                 !! * `rpar(3)` -- tolerance for constraint violations.
                                                 !! * `rpar(4)` -- tolerance for the gradient of the lagrangian function.
                                                 !! * `rpar(5)` -- penalty coefficient.
      real(wp),intent(out) :: f     !! value of the objective function.
      real(wp),intent(out) :: gmax  !! maximum partial derivative of the lagrangian function.
      real(wp),intent(out) :: cmax  !! maximum constraint violation.
      integer,intent(in) :: iprnt  !! print specification:
                                   !!
                                   !! * `iprnt=0`      -- no print.
                                   !! * `abs(iprnt)=1` -- print of final results.
                                   !! * `abs(iprnt)=2` -- print of final results and iterations.
                                   !! * `iprnt>0`      -- basic final results.
                                   !! * `iprnt<0`      -- extended final results.
      integer, intent(out) :: iterm !! variable that indicates the cause of termination.
                                    !!
                                    !! * `iterm=1`  -- if abs(x-xo) was less than or equal to tolx in mtesx (usually two) subsequent iterations.
                                    !! * `iterm=2`  -- if abs(f-fo) was less than or equal to tolf in mtesf (usually two) subsequent iterations.
                                    !! * `iterm=3`  -- if f is less than or equal to tolb.
                                    !! * `iterm=4`  -- if gmax is less than or equal to tolg.
                                    !! * `iterm=11` -- if nit exceeded mit. iterm=12 - if nfv exceeded mfv.
                                    !! * `iterm=13` -- if nfg exceeded mfg. iterm<0 - if the method failed.
                                    !! * `iterm=-6` -- then the termination criterion has not been satisfied, but the point obtained if usually acceptable.
      procedure(obj_func)  :: obj  !! computation of the value of the objective function
      procedure(dobj_func) :: dobj !! computation of the gradient of the objective function
      procedure(con_func)  :: con  !! computation of the value of the constraint function
      procedure(dcon_func) :: dcon !! computation of the gradient of the constraint function
      procedure(report_f),optional :: report !! iteration report function. Note: this is independent of `iprnt`.
                                             !! If this function is associated, each iteration will be reported
      logical,optional,intent(in) :: rpf_adaptive !! Enable adaptive penalty parameter updates (default: .true.)
      real(wp),optional,intent(in) :: rpf_min !! Minimum penalty parameter value (default: 1.0e-6)
      real(wp),optional,intent(in) :: rpf_max !! Maximum penalty parameter value (default: 1.0e6)
      real(wp),optional,intent(in) :: rpf_increase_factor !! Factor to increase rpf when stagnating (default: 10.0)
      real(wp),optional,intent(in) :: rpf_stagnation_threshold !! Threshold for detecting constraint stagnation (default: 0.9)
      integer,optional,intent(in) :: rpf_stagnation_limit !! Iterations of stagnation before increasing rpf (default: 3)

      integer :: lcfd, lcfo, lcg, lcp, lcr, lcz, lg, lgc, lgf, lgo, lh, lia, ls, lxo
      integer, dimension(:), allocatable :: ia
      real(wp), dimension(:), allocatable :: ra
      integer,dimension(:),allocatable :: ic !! local copy of `constraint_type` since it is modified
      integer,dimension(:),allocatable :: ix !! local copy of `bound_type` since it is modified

      ! set the functions:
      me%obj  => obj
      me%dobj => dobj
      me%con  => con
      me%dcon => dcon

      if (present(report)) me%report => report

      ! Set adaptive penalty parameters from optional arguments
      if (present(rpf_adaptive)) me%rpf_adaptive = rpf_adaptive
      if (present(rpf_min)) me%rpf_min = rpf_min
      if (present(rpf_max)) me%rpf_max = rpf_max
      if (present(rpf_increase_factor)) me%rpf_increase_factor = rpf_increase_factor
      if (present(rpf_stagnation_threshold)) me%rpf_stagnation_threshold = rpf_stagnation_threshold
      if (present(rpf_stagnation_limit)) me%rpf_stagnation_limit = rpf_stagnation_limit

      ! Reset adaptive penalty state
      me%cmax_previous = huge(1.0_wp)
      me%rpf_stagnation_count = 0

      ! Conditional allocation based on sparse mode
      if (me%use_sparse) then
         ! Sparse mode: no need to allocate dense Jacobian storage
         ! Still need space for: cfo(nc+1), cfd(nc), gc(nf), cr(nf*(nf+1)/2),
         ! cz(nf), cp(nc), gf(nf), g(nf), h(nf*(nf+1)/2), s(nf), xo(nf), go(nf)
         allocate (ia(nf), ra(nf*(nf + 8) + 3*nc + 1))
      else
         ! Dense mode: allocate space for cg(nf*nc)
         allocate (ia(nf), ra((nf + nc + 8)*nf + 3*nc + 1))
      end if
      allocate (ic(nc))
      allocate (ix(nf))

      ic = constraint_type ! make a copy of input
      ix = bound_type ! make a copy of input

      if (me%use_sparse) then
         ! Sparse mode: lcg dummy, reduced memory layout
         lcg = 1  ! dummy value, cg array not used in sparse mode
         lcfo = 1
         lcfd = lcfo + nc + 1
         lgc = lcfd + nc
         lcr = lgc + nf
         lcz = lcr + nf*(nf + 1)/2
         lcp = lcz + nf
         lgf = lcp + nc
         lg = lgf + nf
         lh = lg + nf
         ls = lh + nf*(nf + 1)/2
         lxo = ls + nf
         lgo = lxo + nf
      else
         ! Dense mode: normal memory layout with cg storage
         lcg = 1  ! dense Jacobian storage starts here
         lcfo = lcg + nf*nc
         lcfd = lcfo + nc + 1
         lgc = lcfd + nc
         lcr = lgc + nf
         lcz = lcr + nf*(nf + 1)/2
         lcp = lcz + nf
         lgf = lcp + nc
         lg = lgf + nf
         lh = lg + nf
         ls = lh + nf*(nf + 1)/2
         lxo = ls + nf
         lgo = lxo + nf
      end if
      lia = 1
      call me%psqp(nf, nb, nc, x, ix, xl, xu, cf, ic, cl, cu, ra, ra(lcfo), ra(lcfd), &
                   ra(lgc), ia, ra(lcr), ra(lcz), ra(lcp), ra(lgf), ra(lg), ra(lh), &
                   ra(ls), ra(lxo), ra(lgo), rpar(1), rpar(2), rpar(3), rpar(4), &
                   rpar(5), cmax, gmax, f, ipar(1), ipar(2), ipar(5), ipar(6), &
                   iprnt, iterm)

      deallocate (ia, ra, ic, ix)

   end subroutine psqpn
!***********************************************************************

!***********************************************************************
!> date: 97/01/22
!
!  recursive quadratic programming method with the bfgs variable metric
!  update for general nonlinear programming problems.
!
!### Method
!  recursive quadratic programming method with the bfgs variable metric
!  update.

   subroutine psqp(me, nf, nb, nc, x, ix, xl, xu, cf, ic, cl, cu, cg, cfo, cfd, gc, ica, &
                   cr, cz, cp, gf, g, h, s, xo, go, xmax, tolx, tolc, tolg, rpf, &
                   cmax, gmax, f, mit, mfv, met, mec, iprnt, iterm)

      class(psqp_class), intent(inout) :: me
      real(wp) :: f         !! value of the objective function.
      real(wp) :: cmax      !! maximum constraint violation.
      real(wp) :: gmax      !! maximum partial derivative of the lagrangian function.
      real(wp) :: rpf       !! penalty coefficient.
      real(wp) :: tolx      !! tolerance for change of variables.
      real(wp) :: tolc      !! tolerance for constraint violations.
      real(wp) :: tolg      !! tolerance for the gradient of the lagrangian function.
      real(wp) :: told      !!
      real(wp) :: tols      !!
      real(wp) :: xmax      !! maximum stepsize.
      integer :: iprnt      !! print specification.
                            !!
                            !! * iprnt=0      - no print.
                            !! * abs(iprnt)=1 - print of final results.
                            !! * abs(iprnt)=2 - print of final results and iterations.
                            !! * iprnt>0      - basic final results.
                            !! * iprnt<0      - extended final results.
      integer :: iterm      !! variable that indicates the cause of termination.
                            !!
                            !! * iterm=1-if abs(x-xo) was less than or equal to tolx in
                            !!   mtesx (usually two) subsequent iterations.
                            !! * iterm=2-if abs(f-fo) was less than or equal to tolf in
                            !!   mtesf (usually two) subsequent iterations.
                            !! * iterm=3-if f is less than or equal to tolb.
                            !! * iterm=4-if gmax is less than or equal to tolg.
                            !! * iterm=11-if nit exceeded mit.
                            !! * iterm=12-if nfv exceeded mfv.
                            !! * iterm=13-if nfg exceeded mfg.
                            !! * iterm<0-if the method failed.
                            !! * if iterm=-6, then the termination criterion has not been
                            !!   satisfied, but the point obtained if usually acceptable.
      integer :: met        !! variable metric update used.
                            !!
                            !! * met=1 - the bfgs update.
                            !! * met=2 - the hoshino update.
      integer :: met1       !!
      integer :: mec        !! correction if the negative curvature occurs.
                            !!
                            !! * mec=1 - correction suppressed.
                            !! * mec=2 - powell's correction.
      integer :: mes        !!
      integer :: mfv        !! maximum number of function evaluations.
      integer :: mit        !! maximum number of iterations.
      integer :: nb         !! choice of simple bounds.
                            !!
                            !! * nb=0 - simple bounds suppressed.
                            !! * nb>0 - simple bounds accepted.
      integer :: nc         !! number of linear constraints.
      integer :: nf         !! number of variables.
      real(wp) :: cf(*)     !! cf(nc+1)  vector containing values of the constraint functions.
      real(wp) :: cg(*)     !! cg(nf*nc)  matrix whose columns are normals of the linear constraints.
      real(wp) :: cfo(*)    !! cfo(nc)  vector containing saved values of the constraint functions.
      real(wp) :: cfd(*)    !! cfd(nc)  vector containing increments of the constraint functions.
      real(wp) :: cl(*)     !! cl(nc)  vector containing lower bounds for constraint functions.
      real(wp) :: cu(*)     !! cu(nc)  vector containing upper bounds for constraint functions.
      real(wp) :: cp(*)     !!
      real(wp) :: cr(*)     !! cr(nf*(nf+1)/2)  triangular decomposition of kernel
                            !! of the orthogonal projection.
      real(wp) :: cz(*)     !! cz(nf)  vector of lagrange multipliers.
      real(wp) :: g(*)      !! g(nf)  gradient of the objective function.
      real(wp) :: gc(*)     !! gc(nf)  gradient of the selected constraint function.
      real(wp) :: gf(*)     !! gf(nf)  gradient of the model function.
      real(wp) :: go(*)     !! go(nf)  gradients difference.
      real(wp) :: h(*)      !! h(nf*(nf+1)/2)  triangular decomposition or inversion of
                            !! the hessian matrix approximation.
      real(wp) :: s(*)      !! s(nf)  direction vector.
      real(wp) :: x(*)      !! x(nf)  vector of variables.
      real(wp) :: xl(*)     !! xl(nf)  vector containing lower bounds for variables.
      real(wp) :: xu(*)     !! xu(nf)  vector containing upper bounds for variables.
      real(wp) :: xo(*)     !! xo(nf)  vectors of variables difference.
      integer,intent(inout) :: ic(*) !! ic(nc)  vector containing types of constraints.
                                     !!
                                     !! * ic(kc)=0 - constraint cf(kc) is not used.
                                     !! * ic(kc)=1 - lower constraint cl(kc)<=cf(kc).
                                     !! * ic(kc)=2 - upper constraint cf(kc)<=cu(kc).
                                     !! * ic(kc)=3 - two side constraint cl(kc)<=cf(kc)<=cu(kc).
                                     !! * ic(kc)=5 - equality constraint cf(kc)==cl(kc).
      integer :: ica(*)     !! ica(nf)  vector containing indices of active constraints.
      integer,intent(inout) :: ix(*) !! ix(nf)  vector containing types of bounds.
                                     !!
                                     !! * ix(i)=0 - variable x(i) is unbounded.
                                     !! * ix(i)=1 - lover bound xl(i)<=x(i).
                                     !! * ix(i)=2 - upper bound x(i)<=xu(i).
                                     !! * ix(i)=3 - two side bound xl(i)<=x(i)<=xu(i).
                                     !! * ix(i)=5 - variable x(i) is fixed.

      real(wp) :: alf1, alf2, cmaxo, dmax, eps7, eps9, eta0, &
                  eta2, eta9, fmax, fmin, fo, gnorm, p, po, &
                  r, rmax, rmin, ro, snorm, tolb, tolf, &
                  umax, rp, fp, pp, ff, fc
      integer :: i, idecf, iext, irest, iterd, iterh, iterq, &
                 iters, kbc, kbf, kc, kd, kit, ld, mred, mtesf, &
                 mtesx, n, k, ntesx, iest, inits, kters, maxst, &
                 isys, mfp, nred, ipom, lds

      if (abs(iprnt) > 1) write (6, '(1x,"entry to psqp :")')

      ! initiation
      kbf = 0
      kbc = 0
      if (nb > 0) kbf = 2
      if (nc > 0) kbc = 2
      me%nres = 0
      me%ndec = 0
      me%nrem = 0
      me%nadd = 0
      me%nit = 0
      me%nfv = 0
      me%nfg = 0
      me%nfh = 0
      isys = 0
      iest = 0
      iext = 0
      mtesx = 2
      mtesf = 2
      inits = 1
      iterm = 0
      iters = 0
      iterd = 0
      iterq = 0
      mred = 20
      irest = 1
      iters = 2
      kters = 5
      idecf = 1
      eta0 = 1.0e-15_wp
      eta2 = 1.0e-15_wp
      eta9 = huge(1.0_wp) !1.0e60_wp
      eps7 = 1.0e-15_wp
      eps9 = 1.0e-8_wp
      alf1 = 1.0e-10_wp
      alf2 = 1.0e10_wp
      fmax = huge(1.0_wp) !1.0e60_wp
      fmin = -fmax
      tolb = -fmax
      dmax = eta9
      tolf = 1.0e-16_wp
      if (xmax <= 0.0_wp) xmax = 1.0e+16_wp
      if (tolx <= 0.0_wp) tolx = 1.0e-16_wp
      if (tolg <= 0.0_wp) tolg = 1.0e-6_wp
      if (tolc <= 0.0_wp) tolc = 1.0e-6_wp
      told = 1.0e-8_wp
      tols = 1.0e-4_wp
      if (rpf <= 0.0_wp) rpf = 1.0e-4_wp
      if (met <= 0) met = 1
      met1 = 2
      if (mec <= 0) mec = 2
      mes = 1
      if (mit <= 0) mit = 1000
      if (mfv <= 0) mfv = 2000
      kd = 1
      ld = -1
      kit = 0
      call mxvset(nc, 0.0_wp, cp)
      ! initial operations with simple bounds
      if (kbf > 0) then
         do i = 1, nf
            if ((ix(i) == 3 .or. ix(i) == 4) .and. xu(i) <= xl(i)) then
               xu(i) = xl(i)
               ix(i) = 5
            elseif (ix(i) == 5 .or. ix(i) == 6) then
               xl(i) = x(i)
               xu(i) = x(i)
               ix(i) = 5
            end if
            if (ix(i) == 1 .or. ix(i) == 3) x(i) = max(x(i), xl(i))
            if (ix(i) == 2 .or. ix(i) == 3) x(i) = min(x(i), xu(i))
         end do
      end if
      ! initial operations with general constraints
      if (kbc > 0) then
         k = 0
         do kc = 1, nc
            if ((ic(kc) == 3 .or. ic(kc) == 4) .and. cu(kc) <= cl(kc)) then
               cu(kc) = cl(kc)
               ic(kc) = 5
            elseif (ic(kc) == 5 .or. ic(kc) == 6) then
               cu(kc) = cl(kc)
               ic(kc) = 5
            end if
            k = k + nf
         end do
      end if
      if (kbf > 0) then
         do i = 1, nf
            if (ix(i) >= 5) ix(i) = -ix(i)
            if (ix(i) <= 0) then
            elseif ((ix(i) == 1 .or. ix(i) == 3) .and. x(i) <= xl(i)) then
               x(i) = xl(i)
            elseif ((ix(i) == 2 .or. ix(i) == 3) .and. x(i) >= xu(i)) then
               x(i) = xu(i)
            end if
            call test_simple_bound(x, ix, xl, xu, eps9, i)
            if (ix(i) > 10) ix(i) = 10 - ix(i)
         end do
      end if
      fo = fmin
      gmax = eta9
      dmax = eta9

      main: do

         lds = ld
         call me%compute_obj_and_dobj(nf, x, gf, gf, ff, f, kd, ld, iext)
         ld = lds
         call me%compute_con_and_dcon(nf, nc, x, fc, cf, cl, cu, ic, gc, cg, cmax, kd, ld)
         cf(nc + 1) = f
         ! JW : seems to start with iter 0, so add 1 to iterations for report output
         if (associated(me%report)) call me%report(me%nit+1,x(1:nf),j=f,f=cf(1:nc))
         if (abs(iprnt) > 1) &
            write (6, '(1x,"nit=",i9,2x,"nfv=",i9,2x,"nfg=",i9,2x,"f=",g13.6,2x,"c=",e8.1,2x,"g=",e8.1)') &
            me%nit, me%nfv, me%nfg, f, cmax, gmax
         ! start of the iteration with tests for termination.
         if (iterm < 0) exit main
         if (iters /= 0) then
            if (f <= tolb) then
               iterm = 3
               exit main
            end if
            if (dmax <= tolx) then
               iterm = 1
               ntesx = ntesx + 1
               if (ntesx >= mtesx) exit main
            else
               ntesx = 0
            end if
         end if
         if (me%nit >= mit) then
            iterm = 11
            exit main
         end if
         if (me%nfv >= mfv) then
            iterm = 12
            exit main
         end if
         iterm = 0
         me%nit = me%nit + 1

         restart: do

            ! restart
            n = nf
            if (irest > 0) then
               call mxdsmi(n, h)
               ld = min(ld, 1)
               idecf = 1
               if (kit < me%nit) then
                  me%nres = me%nres + 1
                  kit = me%nit
               else
                  iterm = -10
                  if (iters < 0) iterm = iters - 5
                  exit main
               end if
            end if
            ! direction determination using a quadratic programming procedure
            call mxvcop(nc + 1, cf, cfo)
            mfp = 2
            ipom = 0
            dir_loop: do
               call me%dual_range_space_quad_prog(nf, nc, x, ix, xl, xu, cf, cfd, ic, ica, &
                                                  cl, cu, cg, cr, cz, g, gf, h, s, mfp, kbf, &
                                                  kbc, idecf, eta2, eta9, eps7, &
                                                  eps9, umax, gmax, n, iterq)
               if (iterq < 0) then
                  if (ipom < 10) then
                     ipom = ipom + 1
                     call transform_incompatible_qp_subproblem(nc, cf, ic, cl, cu, kbc)
                     cycle dir_loop
                  end if
                  iterd = iterq - 10
               else
                  ipom = 0
                  iterd = 1
                  gmax = mxvmax(nf, g)
                  gnorm = sqrt(mxvdot(nf, g, g))
                  snorm = sqrt(mxvdot(nf, s, s))
               end if
               exit dir_loop
            end do dir_loop
            if (iterd < 0) iterm = iterd
            if (iterm == 0) then
               call mxvcop(nc + 1, cfo, cf)
               ! test for sufficient descent
               p = mxvdot(nf, g, s)
               irest = 1
               if (snorm <= 0.0_wp) then
               elseif (p + told*gnorm*snorm <= 0.0_wp) then
                  irest = 0
               end if
               if (irest /= 0) cycle restart
               nred = 0
               rmin = alf1*gnorm/snorm
               rmax = min(alf2*gnorm/snorm, xmax/snorm)
               if (gmax <= tolg .and. cmax <= tolc) then
                  iterm = 4
                  exit main
               end if
               call compute_new_penalty_parameters(nf, n, nc, ica, cz, cp)

               ! Adaptive penalty parameter update
               ! Increase rpf if constraints are not improving sufficiently
               if (me%rpf_adaptive .and. me%nit > 0) then
                  ! Check if constraint violation is stagnating (not improving enough)
                  if (cmax > me%rpf_stagnation_threshold * me%cmax_previous) then
                     me%rpf_stagnation_count = me%rpf_stagnation_count + 1
                     if (me%rpf_stagnation_count >= me%rpf_stagnation_limit) then
                        ! Increase penalty parameter
                        rpf = min(rpf * me%rpf_increase_factor, me%rpf_max)
                        me%rpf_stagnation_count = 0
                        if (abs(iprnt) > 1) then
                           write (6, '(1x,"Increased penalty parameter rpf to",1p,e11.4," (stagnation)")') rpf
                        end if
                     end if
                  else
                     ! Constraints improving - reset stagnation counter
                     me%rpf_stagnation_count = 0
                  end if

                  ! Also check if Lagrange multipliers are large compared to rpf
                  ! This indicates the penalty parameter is too small
                  if (nf - n > 0) then
                     ! Compute max magnitude of active Lagrange multipliers
                     rp = 0.0_wp
                     do i = 1, nf - n
                        rp = max(rp, abs(cz(i)))
                     end do
                     ! If max multiplier >> rpf, increase rpf
                     if (rp > 100.0_wp * rpf) then
                        rpf = min(rp / 10.0_wp, me%rpf_max)
                        me%rpf_stagnation_count = 0
                        if (abs(iprnt) > 1) then
                           write (6, '(1x,"Increased penalty parameter rpf to",1p,e11.4," (large multipliers)")') rpf
                        end if
                     end if
                  end if
               end if
               me%cmax_previous = cmax

               call mxvina(nc, ic)
               call compute_augmented_lagrangian(nf, n, nc, cf, ic, ica, cl, cu, cz, rpf, fc, f)
               ! preparation of line search
               ro = 0.0_wp
               fo = f
               po = p
               cmaxo = cmax
               call mxvcop(nf, x, xo)
               call mxvcop(nf, g, go)
               call mxvcop(nf, gf, cr)
               call mxvcop(nc + 1, cf, cfo)

               line_search: do

                  ! line search without directional derivatives
                  call me%extended_line_search(r, ro, rp, f, fo, fp, po, pp, fmin, fmax, &
                                               rmin, rmax, tols, kd, ld, me%nit, kit, nred, &
                                               mred, maxst, iest, inits, iters, kters, mes, isys)
                  if (isys == 0) then
                     kd = 1
                     ! decision after unsuccessful line search
                     if (iters <= 0) then
                        r = 0.0_wp
                        f = fo
                        p = po
                        call mxvcop(nf, xo, x)
                        call mxvcop(nf, cr, gf)
                        call mxvcop(nc + 1, cfo, cf)
                        irest = 1
                        ld = kd
                        cycle restart
                     end if
                     ! computation of the value and the gradient of the objective
                     ! function together with the values and the gradients of the
                     ! approximated functions
                     if (kd > ld) then
                        lds = ld
                        call me%compute_obj_and_dobj(nf, x, gf, gf, ff, f, kd, ld, iext)
                        ld = lds
                        call me%compute_con_and_dcon(nf, nc, x, fc, cf, cl, cu, ic, gc, &
                                                     cg, cmax, kd, ld)
                     end if
                     ! preparation of variable metric update
                     call mxvcop(nf, gf, g)
                     call dual_range_space_qp(me, nf, n, x, xo, ica, cg, cz, g, go, r, f, fo, p, po, &
                                              cmax, cmaxo, dmax, kd, ld, iters)
                     ! variable metric update
                     call bfgs_variable_metric_update(n, h, g, s, xo, go, r, po, me%nit, &
                                                      kit, iterh, met, met1, mec)
                     ! if (mer>0.and.iterh>0) irest=1
                     cycle main   ! end of the iteration
                  else
                     ! go to (11174,11172) isys+1
                     call mxvdir(nf, r, s, xo, x)
                     lds = ld
                     call me%compute_obj_and_dobj(nf, x, gf, g, ff, f, kd, ld, iext)
                     ld = lds
                     call me%compute_con_and_dcon(nf, nc, x, fc, cf, cl, cu, ic, gc, cg, cmax, kd, ld)
                     cf(nc + 1) = f
                     call compute_augmented_lagrangian(nf, n, nc, cf, ic, ica, cl, cu, cz, rpf, fc, f)
                  end if

               end do line_search

            end if

            exit restart
         end do restart

         exit main
      end do main

      if (iprnt > 1 .or. iprnt < 0) write (6, '(1x,"exit from psqp :")')
      if (iprnt /= 0) &
         write (6, '(1x,"nit=",i4,2x,"nfv=",i4,2x,"nfg=",i4,2x,"f=",g13.6,2x,"c=",e8.1,2x,"g=",e8.1,2x,"iterm=",i3)') &
                  me%nit, me%nfv, me%nfg, f, cmax, gmax, iterm
      if (iprnt < 0) write (6, '(1x,"x=",5(g14.7,1x):/(3x,5(g14.7,1x)))') (x(i), i=1, nf)

   end subroutine psqp

!***********************************************************************
!> date: 97/12/01
!
! computation of the value and the gradient of the constraint function.
!
!@note This routine was formerly called `pc1f01`.

   subroutine compute_con_and_dcon(me, nf, nc, x, fc, cf, cl, cu, ic, gc, cg, cmax, kd, ld)

      class(psqp_class), intent(inout) :: me
      real(wp) :: fc      !! value of the selected constraint function.
      real(wp) :: cmax    !! maximum constraint violation.
      integer  :: kd      !! degree of required derivatives.
      integer  :: ld      !! degree of previously computed derivatives.
      integer  :: nc      !! number of constraints.
      integer  :: nf      !! number of variables.
      real(wp) :: cf(*)   !! cf(nc) vector containing values of constraint functions.
      real(wp) :: cl(*)   !! cl(nc) vector containing lower bounds for constraint functions.
      real(wp) :: cu(*)   !! cu(nc) vector containing upper bounds for constraint functions.
      integer  :: ic(*)   !! ic(nc) vector containing types of constraints.
      real(wp) :: gc(nf)  !! gc(nf) gradient of the selected constraint function.
      real(wp) :: cg(*)   !! cg(nf*nc) matrix whose columns are gradients of constraint functions.
      real(wp) :: x(nf)   !! x(nf) vector of variables.

      real(wp) :: pom, temp
      integer :: kc

      if (kd <= ld) return
      if (ld < 0) cmax = 0.0_wp

      ! Optimization: Compute sparse Jacobian once before loop if needed
      if (kd >= 1 .and. ld < 1 .and. me%use_sparse .and. associated(me%sparse_jac)) then
         ! Call sparse Jacobian to get all constraint gradients at once
         call me%sparse_jac(nf, nc, x, me%jac_sparse%values)
      end if

      do kc = 1, nc
         if (kd >= 0) then
            if (ld < 0) then
               call me%con(nf, kc, x, fc)
               cf(kc) = fc
               if (ic(kc) > 0) then
                  pom = 0.0_wp
                  temp = cf(kc)
                  if (ic(kc) == 1 .or. ic(kc) >= 3) pom = min(pom, temp - cl(kc))
                  if (ic(kc) == 2 .or. ic(kc) >= 3) pom = min(pom, cu(kc) - temp)
                  if (pom < 0.0_wp) cmax = max(cmax, -pom)
               end if
            else
               fc = cf(kc)
            end if
            if (kd >= 1) then
               if (ld >= 1) then
                  ! Retrieve previously computed gradient
                  call me%get_constraint_gradient(nf, nc, kc, cg, gc)
               else
                  ! Compute gradient
                  if (me%use_sparse .and. associated(me%sparse_jac)) then
                     ! Sparse mode: extract gradient from pre-computed Jacobian
                     call me%get_constraint_gradient(nf, nc, kc, cg, gc)
                  else
                     ! Dense mode: compute individual gradient via dcon
                     call me%dcon(nf, kc, x, gc)
                  end if
                  ! Store gradient in cg array
                  call me%set_constraint_gradient(nf, nc, kc, gc, cg)
               end if
            end if
         end if
      end do
      ld = kd

   end subroutine compute_con_and_dcon

!***********************************************************************
!> date: 97/12/01
!
! computation of the value and the gradient of the objective function.
!
!@note This routine was formerly called `pf1f01`.

   subroutine compute_obj_and_dobj(me, nf, x, gf, g, ff, f, kd, ld, iext)

      class(psqp_class), intent(inout) :: me
      integer,intent(in) :: nf       !! number of variables.
      real(wp),intent(in) :: x(nf)   !! x(nf)   vector of variables.
      real(wp),intent(out) :: gf(nf) !! gf(nf)  gradient of the model function.
      real(wp),intent(out) :: g(nf)  !! g(nf)   gradient of the objective function.
      real(wp),intent(out) :: ff     !! value of the model function.
      real(wp),intent(out) :: f      !! value of the objective function.
      integer,intent(in) :: kd       !! degree of required derivatives.
      integer,intent(inout) :: ld    !! degree of previously computed derivatives.
      integer,intent(in) :: iext     !! type of extremum.
                                     !!
                                     !! * `iext=0` -- minimum.
                                     !! * `iext=1` -- maximum.

      if (kd <= ld) return
      if (ld < 0) then
         me%nfv = me%nfv + 1
         call me%obj(nf, x, ff)
         if (iext <= 0) then
            f = ff
         else
            f = -ff
         end if
      end if
      if (kd >= 1) then
         if (ld < 1) then
            me%nfg = me%nfg + 1
            call me%dobj(nf, x, gf)
            if (iext > 0) call mxvneg(nf, gf, g)
         end if
      end if
      ld = kd

   end subroutine compute_obj_and_dobj

!***********************************************************************
!> date: 97/12/01
!
! dual range space quadratic programming method.
!
!@note This routine was formerly called `plqdb1`.

   subroutine dual_range_space_quad_prog(me, nf, nc, x, ix, xl, xu, cf, cfd, &
                                         ic, ica, cl, cu, cg, cr, cz, g, go, h, s, &
                                         mfp, kbf, kbc, idecf, &
                                         eta2, eta9, eps7, eps9, umax, gmax, n, iterq)

      class(psqp_class), intent(inout) :: me
      integer :: nf     !! number of variables.
      integer :: nc     !! number of linear constraints.
      integer :: ix(*)  !! ix(nf)  vector containing types of bounds.
      integer :: ic(*)  !! ic(nc)  vector containing types of constraints.
      integer :: ica(*) !! ica(nf)  vector containing indices of active constraints.
      integer :: mfp    !! type of feasible point.
                        !!
                        !! * mfp=1-arbitrary feasible point.
                        !! * mfp=2-optimum feasible point.
                        !! * mfp=3-repeated solution.
      integer :: kbf    !! specification of simple bounds.
                        !!
                        !! * kbf=0-no simple bounds.
                        !! * kbf=1-one sided simple bounds.
                        !! * kbf=2=two sided simple bounds.
      integer :: kbc    !! specification of linear constraints.
                        !!
                        !! * kbc=0 - no linear constraints.
                        !! * kbc=1 - one sided linear constraints.
                        !! * kbc=2 - two sided linear constraints.
      integer :: idecf  !! decomposition indicator.
                        !!
                        !! * idecf=0  - no decomposition.
                        !! * idecf=1  - gill-murray decomposition.
                        !! * idecf=9  - inversion.
                        !! * idecf=10 - diagonal matrix.
      integer :: n      !! dimension of the manifold defined by active constraints.
      integer :: iterq  !! type of feasible point.
                        !!
                        !! * iterq=1  - arbitrary feasible point.
                        !! * iterq=2  - optimum feasible point.
                        !! * iterq=-1 - feasible point does not exists.
                        !! * iterq=-2 - optimum feasible point does not exists.
      real(wp) :: x(*)  !! x(nf)   vector of variables.
      real(wp) :: xl(*) !! xl(nf)  vector containing lower bounds for variables.
      real(wp) :: xu(*) !! xu(nf)  vector containing upper bounds for variables.
      real(wp) :: cf(*) !! cf(nf)  vector containing values of the constraint functions.
      real(wp) :: cfd(*)!! cfd(nc)  vector containing increments of the constraint functions.
      real(wp) :: cl(*) !! cl(nc)  vector containing lower bounds for constraint functions.
      real(wp) :: cu(*) !! cu(nc)  vector containing upper bounds for constraint functions.
      real(wp) :: cg(*) !! cg(nf*nc)  matrix whose columns are normals of the linear constraints.
      real(wp) :: cr(*) !! cr(nf*(nf+1)/2)  triangular decomposition of kernel of the orthogonal projection.
      real(wp) :: cz(*) !! cz(nf)  vector of lagrange multipliers.
      real(wp) :: g(*)  !! g(nf)  gradient of the lagrangian function.
      real(wp) :: go(*) !! go(nf)  saved gradient of the objective function.
      real(wp) :: h(*)  !! h(nf*(nf+1)/2)  triangular decomposition or inversion
                        !! of the hessian matrix approximation.
      real(wp) :: s(*)  !! s(nf)  direction vector.
      real(wp) :: eta2  !! tolerance for positive definiteness of the hessian matrix.
      real(wp) :: eta9  !! maximum for real numbers.
      real(wp) :: eps7  !! tolerance for linear independence of constraints.
      real(wp) :: eps9  !! tolerance for activity of constraints.
      real(wp) :: umax  !! maximum absolute value of a negative lagrange multiplier.
      real(wp) :: gmax  !! maximum absolute value of a partial derivative.

      real(wp) :: con, temp, step, step1, step2, dmax, par, snorm
      integer :: nca, ncr, i, j, k, iold, jold, inew, jnew, knew, &
                 inf, ier, krem, kc, nred, nadd

      con = eta9
      if (idecf < 0) idecf = 1
      if (idecf == 0) then
         ! gill-murray decomposition
         temp = eta2
         call mxdpgf(nf, h, inf, temp, step)
         me%ndec = me%ndec + 1
         idecf = 1
      end if
      if (idecf >= 2 .and. idecf <= 8) then
         iterq = -10
         return
      end if

      ! initiation

      nred = 0
      jold = 0
      jnew = 0
      iterq = 0
      dmax = 0.0_wp
      if (mfp /= 3) then
         n = nf
         nca = 0
         ncr = 0
         if (kbf > 0) call mxvina(nf, ix)
         if (kbc > 0) call mxvina(nc, ic)
      end if

      outer: do

         ! direction determination

         call mxvneg(nf, go, s)
         do j = 1, nca
            kc = ica(j)
            if (kc > 0) then
               ! s = s + cz(j) * constraint_gradient(kc)
               if (me%use_sparse) then
                  ! Sparse: use sparse row operations
                  call me%jac_sparse%sparse_axpy_row(kc, cz(j), nf, s)
               else
                  call mxvdir(nf, cz(j), cg((kc - 1)*nf + 1), s, s)
               end if
            else
               k = -kc
               s(k) = s(k) + cz(j)
            end if
         end do
         call mxvcop(nf, s, g)
         if (idecf == 1) then
            call mxdpgb(nf, h, s, 0)
         else
            call mxdsmm(nf, h, g, s)
         end if
         if (iterq /= 3) then
            ! check of feasibility
            inew = 0
            par = 0.0_wp
            call determine_new_active_linear_constr(me, nf, nc, cf, cfd, ic, cl, cu, &
                                                    cg, s, eps9, par, kbc, inew, knew)
            call determine_new_active_simple_bound(nf, ix, x, xl, xu, s, kbf, inew, &
                                                   knew, eps9, par)
            if (inew == 0) then
               ! solution achieved
               call mxvneg(nf, g, g)
               iterq = 2
               return
            else
               snorm = 0.0_wp
            end if

            inner: do

               ier = 0

               ! stepsize determination

               call update_tri_decomp_general(me, nf, nc, n, ica, cg, cr, h, s, g, eps7, gmax, umax, &
                                              idecf, inew, nadd, ier, 1)
               call mxdprb(nca, cr, g, -1)
               if (knew < 0) call mxvneg(nca, g, g)

               ! primal stepsize

               if (ier /= 0) then
                  step1 = con
               else
                  step1 = -par/umax
               end if

               ! dual stepsize

               iold = 0
               step2 = con
               do j = 1, nca
                  kc = ica(j)
                  if (kc >= 0) then
                     k = ic(kc)
                  else
                     i = -kc
                     k = ix(i)
                  end if
                  if (k <= -5) then
                  elseif ((k == -1 .or. k == -3.) .and. g(j) <= 0.0_wp) then
                  elseif (.not. ((k == -2 .or. k == -4.) .and. g(j) >= 0.0_wp)) then
                     temp = cz(j)/g(j)
                     if (step2 > temp) then
                        iold = j
                        step2 = temp
                     end if
                  end if
               end do

               ! final stepsize

               step = min(step1, step2)
               if (step >= con) then
                  ! feasible solution does not exist
                  iterq = -1
                  return
               end if

               ! new lagrange multipliers

               dmax = step
               call mxvdir(nca, -step, g, cz, cz)
               snorm = snorm + sign(1, knew)*step
               par = par - (step/step1)*par
               if (step == step1) then
                  if (n <= 0) then
                     ! impossible situation
                     iterq = -5
                     return
                  end if

                  ! constraint addition

                  if (ier == 0) then
                     n = n - 1
                     nca = nca + 1
                     ncr = ncr + nca
                     cz(nca) = snorm
                  end if
                  if (inew > 0) then
                     kc = inew
                     call mxvinv(ic, kc, knew)
                  elseif (abs(knew) == 1) then
                     i = -inew
                     call mxvinv(ix, i, knew)
                  else
                     i = -inew
                     if (knew > 0) ix(i) = -3
                     if (knew < 0) ix(i) = -4
                  end if
                  nred = nred + 1
                  me%nadd = me%nadd + 1
                  jnew = inew
                  jold = 0
                  cycle outer
               end if

               ! constraint deletion

               do j = iold, nca - 1
                  cz(j) = cz(j + 1)
               end do
               call me%ops_after_constr_deletion(nf, nc, ix, ic, ica, cr, ic, g, n, iold, krem, ier)
               ncr = ncr - nca
               nca = nca - 1
               jold = iold
               jnew = 0
               if (kbc > 0) call mxvina(nc, ic)
               if (kbf > 0) call mxvina(nf, ix)
               do j = 1, nca
                  kc = ica(j)
                  if (kc > 0) then
                     ic(kc) = -ic(kc)
                  else
                     kc = -kc
                     ix(kc) = -ix(kc)
                  end if
               end do

            end do inner

         end if

         exit outer
      end do outer

   end subroutine dual_range_space_quad_prog

!***********************************************************************
!> date: 97/12/01
!
! triangular decomposition of kernel of the general projection
! is updated after constraint addition.
!
!@note This routine was formerly called `pladr1`.

   subroutine update_tri_decomp_general(me, nf, nc, n, ica, cg, cr, h, s, g, eps7, &
                                        gmax, umax, idecf, inew, nadd, ier, job)

      class(psqp_class), intent(in) :: me !! psqp class for sparse support
      integer :: nf       !! declared number of variables.
      integer :: nc       !! number of constraints.
      integer :: n        !! actual number of variables.
      integer :: ica(*)   !! ica(nf)  vector containing indices of active constraints.
      integer :: idecf    !! decomposition indicator.
                          !!
                          !! * idecf=0-no decomposition.
                          !! * idecf=1-gill-murray decomposition.
                          !! * idecf=9-inversion.
                          !! * idecf=10-diagonal matrix.
      integer :: inew     !! index of the new active constraint.
      integer :: nadd     !! number of constraint additions.
      integer :: ier      !! error indicator.
      integer :: job      !! specification of computation.
                          !! output vector g is not or is
                          !! computed in case when n<=0 if
                          !! job=0 or job=1 respectively.
      real(wp) :: cg(*)   !! cg(nf*nc)  matrix whose columns are normals of
                          !! the linear constraints.
      real(wp) :: cr(*)   !! cr(nf*(nf+1)/2)  triangular decomposition of
                          !! kernel of the orthogonal projection.
      real(wp) :: h(*)    !! h(nf*(nf+1)/2)  triangular decomposition or
                          !! inversion of the hessian matrix approximation.
      real(wp) :: s(*)    !! s(nf)  auxiliary vector.
      real(wp) :: g(*)    !! g(nf)  vector used in the dual range space
                          !! quadratic programming method.
      real(wp) :: eps7    !! tolerance for linear independence of constraints.
      real(wp) :: gmax    !! maximum absolute value of a partial derivative.
      real(wp) :: umax    !! maximum absolute value of a negative
                          !! lagrange multiplier.

      integer :: nca, ncr, jcg, j, k, l

      ier = 0
      if (job == 0 .and. n <= 0) ier = 2
      if (inew == 0) ier = 3
      if (idecf /= 1 .and. idecf /= 9) ier = -2
      if (ier /= 0) return
      nca = nf - n
      ncr = nca*(nca + 1)/2
      if (inew > 0) then
         jcg = (inew - 1)*nf + 1
         if (idecf == 1) then
            ! Get constraint gradient and apply H^-1
            if (me%use_sparse) then
               call me%jac_sparse%get_row(inew, s(1:nf))
            else
               call mxvcop(nf, cg(jcg), s)
            end if
            call mxdpgb(nf, h, s, 0)
         else
            ! Get constraint gradient and apply H
            if (me%use_sparse) then
               call me%jac_sparse%get_row(inew, g(1:nf))
               call mxdsmm(nf, h, g, s)
            else
               call mxdsmm(nf, h, cg(jcg), s)
            end if
         end if
         ! Compute dot product with constraint gradient
         gmax = me%constraint_gradient_dot(nf, nc, inew, cg, s)
      else
         k = -inew
         if (idecf == 1) then
            call mxvset(nf, 0.0_wp, s)
            s(k) = 1.0_wp
            call mxdpgb(nf, h, s, 0)
         else
            call mxdsmv(nf, h, s, k)
         end if
         gmax = s(k)
      end if
      do j = 1, nca
         l = ica(j)
         if (l > 0) then
            ! Dot product of constraint gradient with vector s
            g(j) = me%constraint_gradient_dot(nf, nc, l, cg, s)
         else
            l = -l
            g(j) = s(l)
         end if
      end do
      if (n == 0) then
         call mxdprb(nca, cr, g, 1)
         umax = 0.0_wp
         ier = 2
         return
      elseif (nca == 0) then
         umax = gmax
      else
         call mxdprb(nca, cr, g, 1)
         umax = gmax - mxvdot(nca, g, g)
         call mxvcop(nca, g, cr(ncr + 1))
      end if
      if (umax <= eps7*gmax) then
         ier = 1
         return
      else
         nca = nca + 1
         ncr = ncr + nca
         ica(nca) = inew
         cr(ncr) = sqrt(umax)
         if (job == 0) then
            n = n - 1
            nadd = nadd + 1
         end if
      end if

   end subroutine update_tri_decomp_general

!***********************************************************************
!> date: 97/12/01
!
! determination of the new active linear constraint.
!
!@note This routine was formerly called `plminn`.

   subroutine determine_new_active_linear_constr(me, nf, nc, cf, cfd, ic, cl, cu, &
                                                 cg, s, eps9, par, kbc, inew, knew)

      class(psqp_class), intent(in) :: me !! psqp class for sparse support
      integer :: nf      !! number of variables.
      integer :: nc      !! number of constraints.
      integer :: ic(*)   !! ic(nc)  vector containing types of constraints.
      integer :: kbc     !! specification of linear constraints.
                         !!
                         !! * kbc=0 - no linear constraints.
                         !! * kbc=1 - one sided linear constraints.
                         !! * kbc=2 - two sided linear constraints.
      integer :: inew    !! index of the new active constraint.
      integer :: knew    !! signum of the new active normal.
      real(wp) :: cf(*)  !! cf(nc)  vector containing values of the
                         !! constraint functions.
      real(wp) :: cfd(*) !! cfd(nc)  vector containing increments of
                         !! the constraint functions.
      real(wp) :: cl(*)  !! cl(nc)  vector containing lower bounds for
                         !! constraint functions.
      real(wp) :: cu(*)  !! cu(nc)  vector containing upper bounds for
                         !! constraint functions.
      real(wp) :: cg(*)  !! cg(nf*nc)  matrix whose columns are normals
                         !! of the linear constraints.
      real(wp) :: s(*)   !! s(nf)  direction vector.
      real(wp) :: eps9   !! tolerance for active constraints.
      real(wp) :: par    !! auxiliary variable.

      real(wp) :: temp, pom
      integer :: jcg, kc

      if (kbc > 0) then
         jcg = 1
         do kc = 1, nc
            if (ic(kc) > 0) then
               ! Compute dot product with constraint gradient
               temp = me%constraint_gradient_dot(nf, nc, kc, cg, s)
               cfd(kc) = temp
               temp = cf(kc) + temp
               if (ic(kc) == 1 .or. ic(kc) >= 3) then
                  pom = temp - cl(kc)
                  if (pom < min(par, -eps9*max(abs(cl(kc)), 1.0_wp))) then
                     inew = kc
                     knew = 1
                     par = pom
                  end if
               end if
               if (ic(kc) == 2 .or. ic(kc) >= 3) then
                  pom = cu(kc) - temp
                  if (pom < min(par, -eps9*max(abs(cu(kc)), 1.0_wp))) then
                     inew = kc
                     knew = -1
                     par = pom
                  end if
               end if
            end if
            jcg = jcg + nf
         end do
      end if

   end subroutine determine_new_active_linear_constr

!***********************************************************************
!> date: 91/12/01
!
! determination of the new active simple bound.
!
!@note This routine was formerly called `plmins`.

   subroutine determine_new_active_simple_bound(nf, ix, xo, xl, xu, s, kbf, &
                                                inew, knew, eps9, par)

      real(wp) :: eps9   !! tolerance for active constraints.
      real(wp) :: par    !! auxiliary variable.
      integer :: inew    !! index of the new active constraint.
      integer :: kbf     !! specification of simple bounds.
                         !!
                         !! * kbf=0-no simple bounds.
                         !! * kbf=1-one sided simple bounds.
                         !! * kbf=2=two sided simple bounds.
      integer :: knew    !! signum of the new normal.
      integer :: nf      !! declared number of variables.
      real(wp) :: s(*)   !! s(nf)  direction vector.
      real(wp) :: xl(*)  !! xl(nf)  vector containing lower bounds
                         !! for variables.
      real(wp) :: xo(*)  !! xo(nf)  saved vector of variables.
      real(wp) :: xu(*)  !! xu(nf)  vector containing upper bounds
                         !! for variables.
      integer :: ix(*)   !! ix(nf)  vector containing types of bounds.

      real(wp) :: pom, temp
      integer :: i

      if (kbf > 0) then
         do i = 1, nf
            if (ix(i) > 0) then
               temp = 1.0_wp
               if (ix(i) == 1 .or. ix(i) >= 3) then
                  pom = xo(i) + s(i)*temp - xl(i)
                  if (pom < min(par, -eps9*max(abs(xl(i)), temp))) then
                     inew = -i
                     knew = 1
                     par = pom
                  end if
               end if
               if (ix(i) == 2 .or. ix(i) >= 3) then
                  pom = xu(i) - s(i)*temp - xo(i)
                  if (pom < min(par, -eps9*max(abs(xu(i)), temp))) then
                     inew = -i
                     knew = -1
                     par = pom
                  end if
               end if
            end if
         end do
      end if

   end subroutine determine_new_active_simple_bound

!***********************************************************************
!> date: 97/12/01
!
! test on activity of a given simple bound.
!
!@note This routine was formerly called `plnews`.

   subroutine test_simple_bound(x, ix, xl, xu, eps9, i)

      real(wp),intent(in) :: x(*)    !! x(nf)  vector of variables.
      integer,intent(inout) :: ix(*) !! ix(nf)  vector containing types of bounds.
      real(wp),intent(in) :: xl(*)   !! xl(nf)  vector containing lower bounds for variables.
      real(wp),intent(in) :: xu(*)   !! xu(nf)  vector containing upper bounds for variables.
      real(wp),intent(in) :: eps9    !! tolerance for active constraints.
      integer,intent(in) :: i        !! index of tested simple bound.

      real(wp) :: temp

      temp = 1.0_wp
      if (ix(i) <= 0) then
      elseif (ix(i) == 1) then
         if (x(i) <= xl(i) + eps9*max(abs(xl(i)), temp)) then
            ix(i) = 11
         end if
      elseif (ix(i) == 2) then
         if (x(i) >= xu(i) - eps9*max(abs(xu(i)), temp)) then
            ix(i) = 12
         end if
      elseif (ix(i) == 3 .or. ix(i) == 4) then
         if (x(i) <= xl(i) + eps9*max(abs(xl(i)), temp)) then
            ix(i) = 13
         end if
         if (x(i) >= xu(i) - eps9*max(abs(xu(i)), temp)) then
            ix(i) = 14
         end if
      end if

   end subroutine test_simple_bound

!***********************************************************************
!> date: 98/12/01
!
! transformation of the incompatible quadratic programming subproblem.
!
!@note This routine was formerly called `plredl`.

   subroutine transform_incompatible_qp_subproblem(nc, cf, ic, cl, cu, kbc)

      integer :: nc     !! number of current linear constraints.
      integer :: ic(nc) !! ic(nc)  vector containing types of constraints.
      integer :: kbc    !! specification of linear constraints.
                        !!
                        !! * kbc=0-no linear constraints.
                        !! * kbc=1-one sided linear constraints.
                        !! * kbc=2=two sided linear constraints.
      real(wp) :: cf(*) !! cf(nf)  vector containing values of the constraint functions.
      real(wp) :: cl(*) !! cl(nc)  vector containing lower bounds for constraint functions.
      real(wp) :: cu(*) !! cu(nc)  vector containing upper bounds for constraint functions.

      real(wp) :: temp
      integer :: k, kc

      if (kbc > 0) then
         do kc = 1, nc
            k = ic(kc)
            if (abs(k) == 1 .or. abs(k) == 3 .or. abs(k) == 4) then
               temp = (cf(kc) - cl(kc))
               if (temp < 0) cf(kc) = cl(kc) + 0.1_wp*temp
            end if
            if (abs(k) == 2 .or. abs(k) == 3 .or. abs(k) == 4) then
               temp = (cf(kc) - cu(kc))
               if (temp > 0) cf(kc) = cu(kc) + 0.1_wp*temp
            end if
            if (abs(k) == 5 .or. abs(k) == 6) then
               temp = (cf(kc) - cl(kc))
               cf(kc) = cl(kc) + 0.1_wp*temp
            end if
         end do
      end if

   end subroutine transform_incompatible_qp_subproblem

!***********************************************************************
!> date: 91/12/01
!
! operations after constraint deletion.
!
!@note This routine was formerly called `plrmf0`.

   subroutine ops_after_constr_deletion(me, nf, nc, ix, ia, iaa, ar, &
                                        ic, s, n, iold, krem, ier)

      class(psqp_class), intent(inout) :: me
      integer :: ier    !! error indicator.
      integer :: iold   !! index of the old active constraint.
      integer :: krem   !! auxiliary variable.
      integer :: n      !! actual number of variables.
      integer :: nc     !! number of constraints.
      integer :: nf     !! declared number of variables.
      real(wp) :: ar(*) !! ar((nf+1)*(nf+2)/2)  triangular decomposition
                        !! of kernel of the orthogonal projection.
      real(wp) :: s(*)  !! s(nf+1)  auxiliary vector.
      integer :: ia(*)  !! ia(na)  vector containing types of deviations.
      integer :: iaa(*) !! iaa(nf+1)  vector containing indices of active
                        !! functions.
      integer :: ic(*)  !! ic(nc)  vector containing types of constraints.
      integer :: ix(*)  !! ix(nf)  vector containing types of bounds.

      integer :: l

      call update_tri_decomp_orthogonal(nf, iaa, ar, s, n, iold, krem, ier)
      n = n + 1
      me%nrem = me%nrem + 1
      l = iaa(nf - n + 1)
      if (l > nc) then
         l = l - nc
         ia(l) = -ia(l)
      elseif (l > 0) then
         ic(l) = -ic(l)
      else
         l = -l
         ix(l) = -ix(l)
      end if

   end subroutine ops_after_constr_deletion

!***********************************************************************
!> date: 91/12/01
!
! triangular decomposition of kernel of the orthogonal projection is
! updated after constraint deletion.
!
!@note This routine was formerly called `plrmr0`.

   subroutine update_tri_decomp_orthogonal(nf, ica, cr, g, n, iold, krem, ier)

      integer :: ier    !! error indicator.
      integer :: iold   !! index of the old active constraint.
      integer :: krem   !! auxiliary variable.
      integer :: n      !! actual number of variables.
      integer :: nf     !! declared number of variables.
      real(wp) :: cr(*) !! cr(nf*(nf+1)/2)  triangular decomposition
                        !! of kernel of the orthogonal projection.
      real(wp) :: g(*)  !! g(nf)  auxiliary vector.
      integer :: ica(*) !! ica(nf)  vector containing indices of active constraints.

      real(wp) :: ck, cl
      integer :: i, j, k, kc, l, nca

      nca = nf - n
      if (iold < nca) then
         k = iold*(iold - 1)/2
         kc = ica(iold)
         call mxvcop(iold, cr(k + 1), g)
         call mxvset(nca - iold, 0.0_wp, g(iold + 1))
         k = k + iold
         do i = iold + 1, nca
            k = k + i
            call mxvort(cr(k - 1), cr(k), ck, cl, ier)
            call mxvrot(g(i - 1), g(i), ck, cl, ier)
            l = k
            do j = i, nca - 1
               l = l + j
               call mxvrot(cr(l - 1), cr(l), ck, cl, ier)
            end do
         end do
         k = iold*(iold - 1)/2
         do i = iold, nca - 1
            l = k + i
            ica(i) = ica(i + 1)
            call mxvcop(i, cr(l + 1), cr(k + 1))
            k = l
         end do
         ica(nca) = kc
         call mxvcop(nca, g, cr(k + 1))
      end if
      krem = 1

   end subroutine update_tri_decomp_orthogonal

!***********************************************************************
!> date: 91/12/01
!
! extrapolation or interpolation for line search without directional
! derivatives.
!
!### Method
! extrapolation or interpolation with standard model functions.
!
!@note This routine was formerly called `pnint3`.

   subroutine line_search_interpolation(ro, rl, ru, ri, fo, fl, fu, fi, &
                                        po, r, mode, mtyp, merr)

      real(wp) :: fo   !! value of the objective function for r=ro.
      real(wp) :: fl   !! value of the objective function for r=rl.
      real(wp) :: fu   !! value of the objective function for r=ru.
      real(wp) :: fi   !! value of the objective function for r=ri.
      real(wp) :: po   !! initial value of the directional derivative.
      real(wp) :: r    !! value of the stepsize parameter obtained.
      real(wp) :: rl   !! lower value of the stepsize parameter.
      real(wp) :: ru   !! upper value of the stepsize parameter.
      real(wp) :: ri   !! inner value of the stepsize parameter.
      real(wp) :: ro   !! initial value of the stepsize parameter.
      integer :: merr  !! error indicator. merr=0 for normal return.
      integer :: mode  !! mode of line search.
      integer :: mtyp  !! method selection
                       !!
                       !! * mtyp=1 - bisection.
                       !! * mtyp=2 - two point quadratic interpolation.
                       !! * mtyp=2 - three point quadratic interpolation.

      real(wp) :: ai, al, au, den, dis
      integer :: ntyp
      logical :: l1, l2

      real(wp), parameter :: zero = 0.0_wp
      real(wp), parameter :: half = 0.5_wp
      real(wp), parameter :: one = 1.0_wp
      real(wp), parameter :: two = 2.0_wp
      real(wp), parameter :: three = 3.0_wp
      real(wp), parameter :: c1l = 1.1_wp
      real(wp), parameter :: c1u = 1000.0_wp
      real(wp), parameter :: c2l = 1.0e-2_wp
      real(wp), parameter :: c2u = 0.9_wp
      real(wp), parameter :: c3l = 1.0e-1_wp

      merr = 0
      if (mode <= 0) return
      if (po >= zero) then
         merr = 2
         return
      elseif (ru <= rl) then
         merr = 3
         return
      end if
      l1 = rl <= ro
      l2 = ri <= rl
      main: do ntyp = mtyp, 1, -1
         if (ntyp == 1) then
            ! bisection
            if (mode == 1) then
               r = two*ru
               return
            elseif (ri - rl <= ru - ri) then
               r = half*(ri + ru)
               return
            else
               r = half*(rl + ri)
               return
            end if
         elseif (ntyp == mtyp .and. l1) then
            if (.not. l2) ai = (fi - fo)/(ri*po)
            au = (fu - fo)/(ru*po)
         end if
         if (l1 .and. (ntyp == 2 .or. l2)) then
            ! two point quadratic extrapolation or interpolation
            if (au >= one) cycle main
            r = half*ru/(one - au)
         elseif (.not. l1 .or. .not. l2 .and. ntyp == 3) then
            ! three point quadratic extrapolation or interpolation
            al = (fi - fl)/(ri - rl)
            au = (fu - fi)/(ru - ri)
            den = au - al
            if (den <= zero) cycle main
            r = ri - half*(au*(ri - rl) + al*(ru - ri))/den
         elseif (l1 .and. .not. l2 .and. ntyp == 4) then
            ! three point cubic extrapolation or interpolation
            dis = (ai - one)*(ru/ri)
            den = (au - one)*(ri/ru) - dis
            dis = au + ai - den - two*(one + dis)
            dis = den*den - three*dis
            if (dis < zero) cycle main
            den = den + sqrt(dis)
            if (den == zero) cycle main
            r = (ru - ri)/den
         else
            cycle main
         end if
         if (mode == 1 .and. r > ru) then
            ! extrapolation accepted
            r = max(r, c1l*ru)
            r = min(r, c1u*ru)
            return
         elseif (mode == 2 .and. r > rl .and. r < ru) then
            ! interpolation accepted
            if (ri == zero .and. ntyp /= 4) then
               r = max(r, rl + c2l*(ru - rl))
            else
               r = max(r, rl + c3l*(ru - rl))
            end if
            r = min(r, rl + c2u*(ru - rl))
            if (r /= ri) return
         end if
      end do main

   end subroutine line_search_interpolation

!***********************************************************************
!> date: 97/12/01
!
! computation of value of the augmented lagrangian function.
!
!@note This routine was formerly called `pp0af8`.

   pure subroutine compute_augmented_lagrangian(nf, n, nc, cf, ic, ica, cl, cu, cz, rpf, fc, f)

      integer,intent(in) :: nf      !! number of variables.
      integer,intent(in) :: n       !! dimension of the constraint null space.
      integer,intent(in) :: nc      !! number of constraints.
      integer,intent(in) :: ic(*)   !! ic(nc)  vector containing types of constraints.
      integer,intent(in) :: ica(*)  !! ica(nf)  vector containing indices of active constraints.
      real(wp),intent(in) :: cf(*)  !! cf(nc+1)  vector containing values of the constraints.
      real(wp),intent(in) :: cl(*)  !! cl(nc)  vector containing lower bounds for constraint functions.
      real(wp),intent(in) :: cu(*)  !! cu(nc)  vector containing upper bounds for constraint functions.
      real(wp),intent(in) :: cz(*)  !! cz(nc)  vector of lagrange multipliers.
      real(wp),intent(in) :: rpf    !! penalty coefficient.
      real(wp),intent(out) :: fc    !! value of the penalty term.
      real(wp),intent(out) :: f     !! value of the penalty function.

      real(wp) :: pom, temp
      integer :: j, kc

      fc = 0.0_wp
      do kc = 1, nc
         if (ic(kc) > 0) then
            pom = 0.0_wp
            temp = cf(kc)
            if (ic(kc) == 1 .or. ic(kc) >= 3) pom = min(pom, temp - cl(kc))
            if (ic(kc) == 2 .or. ic(kc) >= 3) pom = min(pom, cu(kc) - temp)
            fc = fc + rpf*abs(pom)
         end if
      end do
      do j = 1, nf - n
         kc = ica(j)
         if (kc > 0) then
            pom = 0.0_wp
            temp = cf(kc)
            if (ic(kc) == 1 .or. ic(kc) == 3 .or. ic(kc) == 5) &
               pom = min(pom, temp - cl(kc))
            if (ic(kc) == 2 .or. ic(kc) == 4 .or. ic(kc) == 6) &
               pom = max(pom, temp - cu(kc))
            fc = fc - cz(j)*pom
         end if
      end do
      f = cf(nc + 1) + fc

   end subroutine compute_augmented_lagrangian

!***********************************************************************
!> date: 97/12/01
!
! computation of the new penalty parameters.
!
!@note This routine was formerly called `ppset2`.

   pure subroutine compute_new_penalty_parameters(nf, n, nc, ica, cz, cp)

      integer,intent(in) :: nf         !! declared number of variables.
      integer,intent(in) :: n          !! actual number of variables.
      integer,intent(in) :: nc         !! number of constraints.
      integer,intent(in) :: ica(*)     !! vector containing indices of active constraints.
      real(wp),intent(in) :: cz(*)     !! vector of lagrange multipliers.
      real(wp),intent(inout) :: cp(*)  !! vector containing penalty parameters.

      real(wp) :: temp
      integer :: j, l, kc

      do kc = 1, nc
         cp(kc) = 0.5_wp*cp(kc)
      end do
      do j = 1, nf - n
         l = ica(j)
         if (l > 0) then
            temp = abs(cz(j))
            cp(l) = max(temp, cp(l) + 0.5_wp*temp)
         end if
      end do

   end subroutine compute_new_penalty_parameters

!***********************************************************************
!> date: 97/12/01
!
!  extended line search without directional derivatives.
!
!### Method
! safeguarded extrapolation and interpolation with extended termination
! criteria.
!
!@note This routine was formerly called `ps0l02`.

   subroutine extended_line_search(me, r, ro, rp, f, fo, fp, po, pp, fmin, fmax, &
                                   rmin, rmax, tols, kd, ld, nit, kit, nred, mred, maxst, iest, &
                                   inits, iters, kters, mes, isys)

      class(psqp_class), intent(inout) :: me
      integer :: kd     !! degree of required dervatives.
      integer :: ld     !! degree of previously computed derivatives.
      integer :: nit    !! actual number of iterations.
      integer :: kit    !! number of the iteration after last restart.
      integer :: nred   !! actual number of extrapolations or interpolations.
      integer :: mred   !! maximum number of extrapolations or interpolations.
      integer :: maxst  !! maximum stepsize indicator. maxst=0 or maxst=1
                        !! if maximum stepsize was not or was reached.
      integer :: iest   !! lower bound specification. iest=0 or iest=1
                        !! if lower bound is not or is given.
      integer :: inits  !! choice of the initial stepsize.
                        !!
                        !! * inits=0 - initial stepsize is specified in the calling program.
                        !! * inits=1 - unit initial stepsize.
                        !! * inits=2 - combined unit and quadratically estimated initial stepsize.
                        !! * inits=3 - quadratically estimated initial stepsize.
      integer :: iters  !! termination indicator.
                        !!
                        !! * iters=0 - zero step.
                        !! * iters=1 - perfect line search.
                        !! * iters=2   goldstein stepsize.
                        !! * iters=3 - curry stepsize.
                        !! * iters=4 - extended curry stepsize.
                        !! * iters=5 - armijo stepsize.
                        !! * iters=6 - first stepsize.
                        !! * iters=7 - maximum stepsize.
                        !! * iters=8 - unbounded function.
                        !! * iters=-1 - mred reached.
                        !! * iters=-2 - positive directional derivative.
                        !! * iters=-3 - error in interpolation.
      integer :: kters  !! termination selection.
                        !!
                        !! * kters=1 - perfect line search.
                        !! * kters=2 - goldstein stepsize.
                        !! * kters=3 - curry stepsize.
                        !! * kters=4 - extended curry stepsize.
                        !! * kters=5 - armijo stepsize.
                        !! * kters=6 - first stepsize.
      integer :: mes    !! method selection.
                        !!
                        !! * mes=1 - bisection.
                        !! * mes=2 - quadratic interpolation (with one directional derivative).
                        !! * mes=3 - quadratic interpolation (with two directional derivatives).
                        !! * mes=4 - cubic interpolation.
                        !! * mes=5 - conic interpolation.
      integer :: isys   !! control parameter.
      real(wp) :: r     !! value of the stepsize parameter.
      real(wp) :: ro    !! initial value of the stepsize parameter.
      real(wp) :: rp    !! previous value of the stepsize parameter.
      real(wp) :: f     !! value of the objective function.
      real(wp) :: fo    !! initial value of the objective function.
      real(wp) :: fp    !! previous value of the objective function.
      real(wp) :: po    !! initial value of the directional derivative.
      real(wp) :: pp    !! previous value of the directional derivative.
      real(wp) :: fmin  !! lower bound for value of the objective function.
      real(wp) :: fmax  !! upper bound for value of the objective function.
      real(wp) :: rmin  !! minimum value of the stepsize parameter
      real(wp) :: rmax  !! maximum value of the stepsize parameter
      real(wp) :: tols  !! termination tolerance for line search
                        !! (in test on the change of the function value).

      real(wp) :: rtemp
      integer :: merr, init1
      logical :: l1, l2, l3, l4, l6, l7

      real(wp), parameter :: tol = 1.0d-4

      if (isys /= 1) then
         ! go to (1,3) isys+1
         me%mes1 = 2
         me%mes2 = 2
         iters = 0
         if (po >= 0.0_wp) then
            r = 0.0_wp
            iters = -2
            isys = 0
            return
         end if
         if (rmax <= 0.0_wp) then
            iters = 0
            isys = 0
            return
         end if
         ! initial stepsize selection
         if (inits > 0) then
            rtemp = fmin - f
         elseif (iest == 0) then
            rtemp = f - fp
         else
            rtemp = max(f - fp, 10.0_wp*(fmin - f))
         end if
         init1 = abs(inits)
         rp = 0.0_wp
         fp = fo
         pp = po
         if (init1 == 0) then
         elseif (init1 == 1 .or. inits >= 1 .and. iest == 0) then
            r = 1.0_wp
         elseif (init1 == 2) then
            r = min(1.0_wp, 4.0_wp*rtemp/po)
         elseif (init1 == 3) then
            r = min(1.0_wp, 2.0_wp*rtemp/po)
         elseif (init1 == 4) then
            r = 2.0_wp*rtemp/po
         end if
         rtemp = r
         r = max(r, rmin)
         r = min(r, rmax)
         me%mode = 0
         me%rl = 0.0_wp
         me%fl = fo
         me%ru = 0.0_wp
         me%fu = fo
         me%ri = 0.0_wp
         me%fi = fo
      elseif (iters /= 0) then
         isys = 0
         return
      else
         if (f <= fmin) then
            iters = 7
            isys = 0
            return
         else
            l1 = r <= rmin .and. nit /= kit
            l2 = r >= rmax
            l3 = f - fo <= tols*r*po .or. f - fmin <= (fo - fmin)/10.0_wp
            l4 = f - fo >= (1.0_wp - tols)*r*po .or. me%mes2 == 2 .and. me%mode == 2
            l6 = me%ru - me%rl <= tol*me%ru .and. me%mode == 2
            l7 = me%mes2 <= 2 .or. me%mode /= 0
            maxst = 0
            if (l2) maxst = 1
         end if
         ! test on termination
         if (l1 .and. .not. l3) then
            iters = 0
            isys = 0
            return
         elseif (l2 .and. .not. f >= me%fu) then
            iters = 7
            isys = 0
            return
         elseif (l6) then
            iters = 1
            isys = 0
            return
         elseif (l3 .and. l7 .and. kters == 5) then
            iters = 5
            isys = 0
            return
         elseif (l3 .and. l4 .and. l7 .and. &
                 (kters == 2 .or. kters == 3 .or. kters == 4)) then
            iters = 2
            isys = 0
            return
         elseif (kters < 0 .or. kters == 6 .and. l7) then
            iters = 6
            isys = 0
            return
         elseif (abs(nred) >= mred) then
            iters = -1
            isys = 0
            return
         else
            rp = r
            fp = f
            me%mode = max(me%mode, 1)
            me%mtyp = abs(mes)
            if (f >= fmax) me%mtyp = 1
         end if
         if (me%mode == 1) then
            ! interval change after extrapolation
            me%rl = me%ri
            me%fl = me%fi
            me%ri = me%ru
            me%fi = me%fu
            me%ru = r
            me%fu = f
            if (f >= me%fi) then
               nred = 0
               me%mode = 2
            elseif (me%mes1 == 1) then
               me%mtyp = 1
            end if
            ! interval change after interpolation
         elseif (r <= me%ri) then
            if (f <= me%fi) then
               me%ru = me%ri
               me%fu = me%fi
               me%ri = r
               me%fi = f
            else
               me%rl = r
               me%fl = f
            end if
         elseif (f <= me%fi) then
            me%rl = me%ri
            me%fl = me%fi
            me%ri = r
            me%fi = f
         else
            me%ru = r
            me%fu = f
         end if
      end if
      ! new stepsize selection (extrapolation or interpolation)
      call line_search_interpolation(ro, me%rl, me%ru, me%ri, fo, me%fl, me%fu, &
                                     me%fi, po, r, me%mode, me%mtyp, merr)
      if (merr > 0) then
         iters = -merr
         isys = 0
         return
      elseif (me%mode == 1) then
         nred = nred - 1
         r = min(r, rmax)
      elseif (me%mode == 2) then
         nred = nred + 1
      end if
      ! computation of the new function value
      kd = 0
      ld = -1
      isys = 1
   end subroutine extended_line_search

!***********************************************************************
!> date: 92/12/01
!
! variable metric update of a dense symmetric positive definite matrix
! using the factorization b=l*d*trans(l).
!
!### Method
! bfgs variable metric method.
!
!@note This routine was formerly called `pudbg1`.

   subroutine bfgs_variable_metric_update(n, h, g, s, xo, go, r, po, nit, &
                                          kit, iterh, met, met1, mec)

      real(wp) :: po  !! old value of the directional derivative.
      real(wp) :: r  !! value of the stepsize parameter.
      integer :: iterh   !! termination indicator.
                         !!
                         !! * iterh<0-bad decomposition.
                         !! * iterh=0-successful update.
                         !! * iterh>0-nonpositive parameters.
      integer :: kit   !! number of the iteration after last restart.
      integer :: met   !!
      integer :: met1   !! selection of self scaling.
                        !!
                        !! * met1=1-self scaling suppressed.
                        !! * met1=2 initial self scaling.
      integer :: mec   !! correction if the negative curvature occurs.
                       !!
                       !! * mec=1-correction suppressed.
                       !! * mec=2-powell's correction.
      integer :: n   !! actual number of variables.
      integer :: nit  !! actual number of iterations.
      real(wp) :: g(*)   !! g(nf)  gradient of the objective function.
      real(wp) :: go(*)   !! go(nf)  gradients difference.
      real(wp) :: h(*)   !! h(m)  factorization b=l*d*trans(l) of a positive
                         !! definite approximation of the hessian matrix.
      real(wp) :: s(*)   !! s(nf)  auxiliary vector.
      real(wp) :: xo(*)  !! xo(nf)  vectors of variables difference.

      real(wp) :: a, b, c, gam, par, den, dis
      logical :: l1, l3

      l1 = met1 >= 3 .or. met1 == 2 .and. nit == kit
      l3 = .not. l1

      ! determination of the parameters b, c
      b = mxvdot(n, xo, go)
      a = 0.0_wp
      if (l1) then
         call mxvcop(n, go, s)
         call mxdpgb(n, h, s, 1)
         a = mxdpgp(n, h, s, s)
         if (a <= 0.0_wp) then
            iterh = 1
            return
         end if
      end if
      call mxvdif(n, go, g, s)
      call mxvscl(n, r, s, s)
      c = -r*po
      if (c <= 0.0_wp) then
         iterh = 3
         return
      end if
      if (mec > 1) then
         if (b <= 1.0e-4_wp*c) then
            ! powell's correction
            dis = (1.0_wp - 0.1_wp)*c/(c - b)
            call mxvdif(n, go, s, go)
            call mxvdir(n, dis, go, s, go)
            b = c + dis*(b - c)
            if (l1) a = c + 2.0_wp*(1.0_wp - dis)*(b - c) + dis*dis*(a - c)
         end if
      elseif (b <= 1.0e-4_wp*c) then
         iterh = 2
         return
      end if
      if (l1) then
         ! determination of the parameter gam (self scaling)
         if (met == 1) then
            par = c/b
         elseif (a <= 0.0_wp) then
            par = c/b
         else
            par = sqrt(c/a)
         end if
         gam = par
         if (met1 > 1) then
            if (nit /= kit) l3 = gam < 0.5_wp .or. gam > 4.0_wp
         end if
      end if
      if (l3) then
         gam = 1.0_wp
         par = gam
      end if
      if (met == 1) then
         ! bfgs update
         call mxdpgu(n, h, par/b, go, xo)
         call mxdpgu(n, h, -1.0_wp/c, s, xo)
      else
         ! hoshino update
         den = par*b + c
         dis = 0.5_wp*b
         call mxvdir(n, par, go, s, s)
         call mxdpgu(n, h, par/dis, go, xo)
         call mxdpgu(n, h, -1.0_wp/den, s, xo)
      end if
      iterh = 0
      if (gam == 1.0_wp) return
      call mxdpgs(n, h, 1.0_wp/gam)

   end subroutine bfgs_variable_metric_update

!***********************************************************************
!> date: 91/12/01
!
! dual range space quadratic programming method for minimax
! approximation.
!
!@note This routine was formerly called `pytrnd`.

   subroutine dual_range_space_qp(me, nf, n, x, xo, ica, cg, cz, g, go, r, f, fo, &
                                  p, po, cmax, cmaxo, dmax, kd, ld, iters)

      class(psqp_class), intent(in) :: me  !! psqp class for sparse support
      integer,intent(in) :: nf  !! declared number of variables.
      integer,intent(inout) :: n  !! actual number of variables.
      integer,intent(in) :: ica(*)  !! ica(nf)  vector containing indices of active constraints.
      real(wp) :: x(*)  !! x(nf)  vector of variables.
      real(wp) :: xo(*)  !! xo(nf)  saved vector of variables.
      real(wp) :: cg(*)  !! cg(nf*nc)  matrix whose columns are normals of the linear constraints.
      real(wp) :: cz(*)  !! cz(nf)  vector of lagrange multipliers.
      real(wp) :: g(*)  !! g(nf)  gradient of the lagrangian function.
      real(wp) :: go(*)  !! go(nf)  saved gradient of the lagrangian function.
      real(wp) :: r  !! value of the stepsize parameter.
      real(wp) :: f  !! new value of the objective function.
      real(wp) :: fo  !! old value of the objective function.
      real(wp) :: p  !! new value of the directional derivative.
      real(wp) :: po  !! old value of the directional derivative.
      real(wp) :: cmax  !! value of the constraint violation.
      real(wp) :: cmaxo  !! saved value of the constraint violation.
      real(wp),intent(out) :: dmax  !! maximum relative difference of variables.
      integer :: kd  !!
      integer :: ld  !!
      integer :: iters  !! termination indicator for steplength determination.
                        !! iters=0 for zero step.

      integer :: i, j, l

      do j = 1, nf - n
         l = ica(j)
         if (l > 0) then
            ! g = g - cz(j) * constraint_gradient(l)
            if (me%use_sparse) then
               call me%jac_sparse%sparse_axpy_row(l, -cz(j), nf, g)
            else
               call mxvdir(nf, -cz(j), cg((l - 1)*nf + 1), g, g)
            end if
         else
            l = -l
            g(l) = g(l) - cz(j)
         end if
      end do
      if (iters > 0) then
         call mxvdif(nf, x, xo, xo)
         call mxvdif(nf, g, go, go)
         po = r*po
         p = r*p
      else
         f = fo
         p = po
         cmax = cmaxo
         call mxvsav(nf, x, xo)
         call mxvsav(nf, g, go)
         ld = kd
      end if
      dmax = 0.0_wp
      do i = 1, nf
         dmax = max(dmax, abs(xo(i))/max(abs(x(i)), 1.0_wp))
      end do
      n = nf
   end subroutine dual_range_space_qp

!***********************************************************************
!>
!  Set the sparsity pattern for the constraint Jacobian.
!  This must be called before using sparse mode.
!
!### Example
!```fortran
!  ! Jacobian has 3 constraints, 5 variables
!  ! Constraint 1: depends on vars 1, 3
!  ! Constraint 2: depends on vars 2, 4, 5
!  ! Constraint 3: depends on var 1
!  integer :: rows(6) = [1, 1, 2, 2, 2, 3]
!  integer :: cols(6) = [1, 3, 2, 4, 5, 1]
!  call solver%set_jacobian_pattern(3, 5, rows, cols, 6)
!```

   subroutine set_jacobian_pattern(me, nc, nf, rows, cols, nnz)

      class(psqp_class), intent(inout) :: me
      integer, intent(in) :: nc   !! number of constraints
      integer, intent(in) :: nf   !! number of variables
      integer, intent(in) :: rows(nnz) !! row indices (constraint indices, 1-based)
      integer, intent(in) :: cols(nnz) !! column indices (variable indices, 1-based)
      integer, intent(in) :: nnz  !! number of non-zeros

      ! Create sparse matrix structure from COO format
      if (allocated(me%jac_sparse)) deallocate(me%jac_sparse)
      allocate(me%jac_sparse)
      call me%jac_sparse%from_coo(nc, nf, rows, cols, nnz)
      me%use_sparse = .true.

   end subroutine set_jacobian_pattern

!***********************************************************************
!>
!  Set the sparse Jacobian callback function.
!
!  This should be called after set_jacobian_pattern to enable sparse mode.

   subroutine set_sparse_jacobian_callback(me, sparse_jac)

      class(psqp_class), intent(inout) :: me
      procedure(sparse_jac_func) :: sparse_jac  !! sparse Jacobian evaluation callback

      me%sparse_jac => sparse_jac

   end subroutine set_sparse_jacobian_callback

!***********************************************************************
!>
!  Get the gradient of a specific constraint (abstraction layer).
!  Works in both dense and sparse modes.

   subroutine get_constraint_gradient(me, nf, nc, kc, cg, gc)

      class(psqp_class), intent(in) :: me
      integer, intent(in) :: nf   !! number of variables
      integer, intent(in) :: nc   !! number of constraints
      integer, intent(in) :: kc   !! constraint index (1-based)
      real(wp), intent(in) :: cg(*) !! dense Jacobian (if dense mode)
      real(wp), intent(out) :: gc(nf) !! constraint gradient output

      if (me%use_sparse) then
         ! Extract row kc from sparse Jacobian
         call me%jac_sparse%get_row(kc, gc)
      else
         ! Dense: copy from cg array
         call mxvcop(nf, cg((kc - 1)*nf + 1), gc)
      end if

   end subroutine get_constraint_gradient

!***********************************************************************
!>
!  Compute dot product of constraint gradient with a vector (abstraction layer).
!  Works in both dense and sparse modes.

   function constraint_gradient_dot(me, nf, nc, kc, cg, v) result(dotprod)

      class(psqp_class), intent(in) :: me
      integer, intent(in) :: nf   !! number of variables
      integer, intent(in) :: nc   !! number of constraints
      integer, intent(in) :: kc   !! constraint index (1-based)
      real(wp), intent(in) :: cg(*) !! dense Jacobian (if dense mode)
      real(wp), intent(in) :: v(nf) !! vector to dot with
      real(wp) :: dotprod

      if (me%use_sparse) then
         ! Sparse row-vector dot product
         dotprod = me%jac_sparse%row_dot(kc, v)
      else
         ! Dense
         dotprod = mxvdot(nf, cg((kc - 1)*nf + 1), v)
      end if

   end function constraint_gradient_dot

!***********************************************************************
!>
!  Set the gradient of a specific constraint (abstraction layer).
!  Works in both dense and sparse modes.

   subroutine set_constraint_gradient(me, nf, nc, kc, gc, cg)

      class(psqp_class), intent(inout) :: me
      integer, intent(in) :: nf   !! number of variables
      integer, intent(in) :: nc   !! number of constraints
      integer, intent(in) :: kc   !! constraint index (1-based)
      real(wp), intent(in) :: gc(nf) !! constraint gradient to set
      real(wp), intent(inout) :: cg(*) !! dense Jacobian (if dense mode)

      if (me%use_sparse) then
         ! Set row kc in sparse Jacobian
         call me%jac_sparse%set_row(kc, gc)
      else
         ! Dense: copy to cg array
         call mxvcop(nf, gc, cg((kc - 1)*nf + 1))
      end if

   end subroutine set_constraint_gradient

!***********************************************************************
end module psqp_module
!***********************************************************************
