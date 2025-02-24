!*******************************************************************************
!> author: Jacob Williams
!
!  Test for the [[psqp_module]].
!
!  Same test case from SLSQP (see `slsqp_test_2`).

    program psqp_test_2

    use psqp_module, only: psqp_class, wp => psqp_wp

    implicit none

    integer,parameter :: nf = 3          !! number of optimization variables
    integer,parameter :: nc = 2          !! total number of constraints
    integer,parameter :: max_iter = 1000 !! maximum number of allowed iterations
    integer,parameter :: nb = 1          !! use simple bounds
    real(wp),dimension(nf),parameter :: xl = [-10.0_wp, -10.0_wp, -10.0_wp]  !! lower bounds
    real(wp),dimension(nf),parameter :: xu = [ 10.0_wp,  10.0_wp,  10.0_wp]  !! upper bounds
    integer,dimension(nc) :: ic = [5, 1]  !! ic(nc) vector containing types of constraints:
                                          !! c(1) = x(1)*x(2) - x(3) !equality constraint (==0)    -- cf(kc) == cl(kc).
                                          !! c(2) = x(3) - 1.0_wp    !inequality constraint (>=0)  -- lower constraint cl(kc) <= cf(kc).
    integer,dimension(6),parameter :: ipar = [  max_iter, &  ! maximum number of iterations.
                                                max_iter, &  ! maximum number of function evaluations.
                                                0, &  ! this parameter is not used in the subroutine psqp.
                                                0, &  ! this parameter is not used in the subroutine psqp.
                                                1, &  ! variable metric update used.
                                                1]    ! correction of the variable metric update if a negative curvature occurs.
    real(wp),dimension(5) :: rpar = [0.0_wp, 2*epsilon(1.0_wp), 2*epsilon(1.0_wp), 0.0_wp, 0.0_wp]    !! real parameters:
                                                                                !!
                                                                                !! * rpar(1)  maximum stepsize.
                                                                                !! * rpar(2)  tolerance for change of variables.
                                                                                !! * rpar(3)  tolerance for constraint violations.
                                                                                !! * rpar(4)  tolerance for the gradient of the lagrangian function.
                                                                                !! * rpar(5)  penalty coefficient.
    integer,parameter :: iprnt = 2  !! print specification:
                                    !!
                                    !! * iprnt=0      - no print.
                                    !! * abs(iprnt)=1 - print of final results.
                                    !! * abs(iprnt)=2 - print of final results and iterations.
                                    !! * iprnt>0      - basic final results.
                                    !! * iprnt<0      - extended final results.

    type(psqp_class)         :: solver      !! instantiate an slsqp solver
    real(wp),dimension(nf)   :: x           !! optimization variable vector
    integer                  :: istat       !! for solver status check
    logical                  :: status_ok   !! for initialization status check
    integer                  :: iterations  !! number of iterations by the solver
    integer,dimension(nf)    :: ix          !! bounds type
    real(wp)                 :: f           !! objective function value
    integer                  :: iterm
    real(wp),dimension(nc+1) :: cf          !! cf(nc+1) vector containing values of the constraint functions.
    real(wp),dimension(nc)   :: cl          !! cl(nc) vector containing lower bounds for constraint functions.
    real(wp),dimension(nc)   :: cu          !! cu(nc) vector containing upper bounds for constraint functions.
    real(wp)                 :: cmax        !! maximum constraint violation.
    real(wp)                 :: gmax        !! maximum partial derivative of the lagrangian function.

    x   = [1.0_wp, 2.0_wp, 3.0_wp]   ! initial guess
    ix  = [3, 3, 3]                  ! all have upper and lower bounds
    cl  = [0.0_wp, 0.0_wp]           ! lower bounds for constraints
    cu  = [10000.0_wp, 10000.0_wp]   ! upper bounds for constraints [NOT USED]

    call solver%psqpn(nf,nb,nc,x,ix,xl,xu,cf,ic,cl,cu,&
                        ipar,rpar,f,cmax,gmax,iprnt,iterm,&
                        obj,dobj,con,dcon)

    write(*,*) ''
    write(*,*) 'iterm = ', iterm
    write(*,*) 'x    = ', x
    write(*,*) 'f    = ', f
    write(*,*) 'cmax = ', cmax
    write(*,*) 'gmax = ', gmax
    write(*,*) ''

    ! Solution is: x = [1,1,1], f = 3

    contains

    subroutine obj(me,nf,x,ff)

        !! user supplied subroutine (calculation of ff)

        class(psqp_class),intent(inout) :: me
        real(wp) :: ff
        integer  :: nf
        real(wp) :: x(nf)

        ff = x(1)**2 + x(2)**2 + x(3)  !objective function

    end subroutine obj

    subroutine dobj(me,nf,x,gf)

        !! user supplied subroutine (calculation of gf)

        class(psqp_class),intent(inout) :: me
        integer  :: nf
        real(wp) :: gf(nf),x(nf)

        gf(1) = 2.0_wp*x(1)
        gf(2) = 2.0_wp*x(2)
        gf(3) = 1.0_wp

    end subroutine dobj

    subroutine con(me,nf,kc,x,fc)

        !! user supplied subroutine (calculation of fc)

        class(psqp_class),intent(inout) :: me
        real(wp) :: fc
        integer  :: kc,nf
        real(wp) :: x(nf)

        select case (kc)
        case(1)
            fc = x(1)*x(2) - x(3)       !equality constraint (==0)
        case(2)
            fc = x(3) - 1.0_wp          !inequality constraint (>=0)
        case default
            error stop 'invalid kc'
        end select

    end subroutine con

    subroutine dcon(me,nf,kc,x,gc)

        !! user supplied subroutine (calculation of gc)

        class(psqp_class),intent(inout) :: me
        integer  :: kc,nf
        real(wp) :: gc(nf),x(nf)

        select case (kc)
        case(1)
            gc(1) = x(2)
            gc(2) = x(1)
            gc(3) = -1.0_wp
        case(2)
            gc(1) = 0.0_wp
            gc(2) = 0.0_wp
            gc(3) = 1.0_wp
        case default
            error stop 'invalid kc'
        end select

    end subroutine dcon

    end program psqp_test_2
!*******************************************************************************