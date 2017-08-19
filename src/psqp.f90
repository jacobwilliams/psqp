!***********************************************************************
!>
!  module for psqp

    module psqp_module

    use matrix_routines

    implicit none

    private

    type,public :: psqp_class

        private

        ! these were formerly in the `stat` common block:
        integer,public :: nres = 0 !! number of restarts.
        integer,public :: ndec = 0 !! number of matrix decomposition.
        integer,public :: nrem = 0 !! number of constraint deletions.
        integer,public :: nadd = 0 !! number of constraint additions.
        integer,public :: nit  = 0 !! number of iterations.
        integer,public :: nfv  = 0 !! number of function evaluations.
        integer,public :: nfg  = 0 !! number of gradient evaluations.
        integer,public :: nfh  = 0 !! number of hessian evaluations.

        ! formerly saved variables in ps0l02:
        integer :: mtyp = 0
        integer :: mode = 0
        integer :: mes1 = 0
        integer :: mes2 = 0
        double precision :: rl = 0.0d0
        double precision :: fl = 0.0d0
        double precision :: ru = 0.0d0
        double precision :: fu = 0.0d0
        double precision :: ri = 0.0d0
        double precision :: fi = 0.0d0

        procedure(obj_func) ,pointer :: obj  => null() !! objective function
        procedure(dobj_func),pointer :: dobj => null() !! gradient of the objective function
        procedure(con_func) ,pointer :: con  => null() !! constraint function
        procedure(dcon_func),pointer :: dcon => null() !! gradient of the constraint function

    contains

        private

        procedure,public :: psqpn
        procedure,public :: psqp

        procedure :: pf1f01
        procedure :: plqdb1
        procedure :: plrmf0
        procedure :: pc1f01
        procedure :: ps0l02

    end type psqp_class

    abstract interface

        subroutine obj_func(me,nf,x,ff)
            !! objective function interface
            import
            implicit none
            class(psqp_class),intent(inout) :: me
            integer          :: nf    !! the number of variables
            double precision :: x(nf) !! a vector of variables
            double precision :: ff    !! the value of the objective function
        end subroutine obj_func

        subroutine dobj_func(me,nf,x,gf)
            !! gradient of the objective function interface
            import
            implicit none
            class(psqp_class),intent(inout) :: me
            integer          :: nf     !! the number of variables
            double precision :: x(nf)  !! a vector of variables
            double precision :: gf(nf) !! the gradient of the objective function
        end subroutine dobj_func

        subroutine con_func(me,nf,kc,x,fc)
            !! constraint function interface
            import
            implicit none
            class(psqp_class),intent(inout) :: me
            integer          :: nf      !! the number of variables
            integer          :: kc      !! the index of the constraint function
            double precision :: x(nf)   !! a vector of variables
            double precision :: fc      !! the value of the constraint function
        end subroutine con_func

        subroutine dcon_func(me,nf,kc,x,gc)
            !! gradient of the constraint function interface
            import
            implicit none
            class(psqp_class),intent(inout) :: me
            integer          :: nf     !! the number of variables
            integer          :: kc     !! the index of the constraint function
            double precision :: x(nf)  !! a vector of variables and
            double precision :: gc(nf) !! the gradient of the constraint function
        end subroutine dcon_func

    end interface

    contains
!***********************************************************************

!***********************************************************************
!> date: 97/01/22
!
! easy to use subroutine for general nonlinear programming problems.

      subroutine psqpn(me,nf,nb,nc,x,ix,xl,xu,cf,ic,cl,cu,ipar,rpar,f,gmax,&
                       cmax,iprnt,iterm,obj,dobj,con,dcon)

      implicit none

      class(psqp_class),intent(inout) :: me

      double precision :: f     !! value of the objective function.
      double precision :: cmax  !! maximum constraint violation.
      double precision :: gmax  !! maximum partial derivative of the lagrangian function.
      integer :: iprnt  !! print specification:
                        !!
                        !! * iprnt=0      - no print.
                        !! * abs(iprnt)=1 - print of final results.
                        !! * abs(iprnt)=2 - print of final results and iterations.
                        !! * iprnt>0      - basic final results.
                        !! * iprnt<0      - extended final results.
      integer,intent(out) :: iterm  !! variable that indicates the cause of termination.
                                    !!
                                    !! * iterm=1  - if abs(x-xo) was less than or equal to tolx in mtesx (usually two) subsequent iterations.
                                    !! * iterm=2  - if abs(f-fo) was less than or equal to tolf in mtesf (usually two) subsequent iterations.
                                    !! * iterm=3  - if f is less than or equal to tolb.
                                    !! * iterm=4  - if gmax is less than or equal to tolg.
                                    !! * iterm=11 - if nit exceeded mit. iterm=12-if nfv exceeded mfv.
                                    !! * iterm=13 - if nfg exceeded mfg. iterm<0-if the method failed.
                                    !! * iterm=-6, then the termination criterion has not been satisfied, but the point obtained if usually acceptable.
      integer,intent(in) :: nb     !! choice of simple bounds. nb=0-simple bounds suppressed.
                                   !! nb>0-simple bounds accepted.
      integer,intent(in) :: nc     !! number of linear constraints.
      integer,intent(in) :: nf     !! number of variables
      double precision :: cf(*)    !! vector containing values of the constraint functions.
      double precision :: cl(*)    !! vector containing lower bounds for constraint functions.
      double precision :: cu(*)    !! vector containing upper bounds for constraint functions.
      double precision :: rpar(5)  !! real parameters:
                                   !!
                                   !! * rpar(1)  maximum stepsize.
                                   !! * rpar(2)  tolerance for change of variables.
                                   !! * rpar(3)  tolerance for constraint violations.
                                   !! * rpar(4)  tolerance for the gradient of the lagrangian function.
                                   !! * rpar(5)  penalty coefficient.
      double precision :: x(*)    !! vector of variables.
      double precision :: xl(*)   !! vector containing lower bounds for variables.
      double precision :: xu(*)   !! vector containing upper bounds for variables.
      integer :: ic(*)  !! vector containing types of constraints:
                        !!
                        !! * ic(kc) = 0 -- constraint cf(kc) is not used.
                        !! * ic(kc) = 1 -- lower constraint cl(kc) <= cf(kc).
                        !! * ic(kc) = 2 -- upper constraint cf(kc) <= cu(kc).
                        !! * ic(kc) = 3 -- two side constraint cl(kc) <= cf(kc) <= cu(kc).
                        !! * ic(kc) = 5 -- equality constraint cf(kc) == cl(kc).
      integer :: ipar(6)    !! integer paremeters:
                            !!
                            !! * ipar(1)  maximum number of iterations.
                            !! * ipar(2)  maximum number of function evaluations.
                            !! * ipar(3)  this parameter is not used in the subroutine psqp.
                            !! * ipar(4)  this parameter is not used in the subroutine psqp.
                            !! * ipar(5)  variable metric update used.
                            !!   ipar(5)=1 - the bfgs update.
                            !!   ipar(5)=2 - the hoshino update.
                            !! * ipar(6)  correction of the variable metric update if a negative
                            !!   curvature occurs.
                            !!   ipar(6)=1 - no correction.
                            !!   ipar(6)=2 - powell's correction.
      integer,intent(in) :: ix(*)  !! vector containing types of bounds.
                                   !!
                                   !! * ix(i) = 0 -- variable x(i) is unbounded.
                                   !! * ix(i) = 1 -- lower bound xl(i) <= x(i).
                                   !! * ix(i) = 2 -- upper bound x(i) <= xu(i).
                                   !! * ix(i) = 3 -- two side bound xl(i) <= x(i) <= xu(i).
                                   !! * ix(i) = 5 -- variable x(i) is fixed.
      procedure(obj_func)  :: obj  !! computation of the value of the objective function
      procedure(dobj_func) :: dobj !! computation of the gradient of the objective function
      procedure(con_func)  :: con  !! computation of the value of the constraint function
      procedure(dcon_func) :: dcon !! computation of the gradient of the constraint function

      integer :: lcfd,lcfo,lcg,lcp,lcr,lcz,lg,lgc,lgf,lgo,lh,lia,ls,lxo
      integer,dimension(:),allocatable :: ia
      double precision,dimension(:),allocatable :: ra

      ! set the functions:
      me%obj  => obj
      me%dobj => dobj
      me%con  => con
      me%dcon => dcon

      allocate (ia(nf),ra((nf+nc+8)*nf+3*nc+1))

      lcg = 1
      lcfo = lcg + nf*nc
      lcfd = lcfo + nc + 1
      lgc = lcfd + nc
      lcr = lgc + nf
      lcz = lcr + nf*(nf+1)/2
      lcp = lcz + nf
      lgf = lcp + nc
      lg = lgf + nf
      lh = lg + nf
      ls = lh + nf*(nf+1)/2
      lxo = ls + nf
      lgo = lxo + nf
      lia = 1
      call me%psqp(nf,nb,nc,x,ix,xl,xu,cf,ic,cl,cu,ra,ra(lcfo),ra(lcfd), &
                ra(lgc),ia,ra(lcr),ra(lcz),ra(lcp),ra(lgf),ra(lg),ra(lh),&
                ra(ls),ra(lxo),ra(lgo),rpar(1),rpar(2),rpar(3),rpar(4),  &
                rpar(5),cmax,gmax,f,ipar(1),ipar(2),ipar(5),ipar(6),     &
                iprnt,iterm)
      deallocate (ia,ra)

    end subroutine psqpn
!***********************************************************************

!***********************************************************************
!> date: 97/01/22
!
! recursive quadratic programming method with the bfgs variable metric
! update for general nonlinear programming problems.
!
! parameters :
!  ii  nf  number of variables.
!  ii  nb  choice of simple bounds. nb=0-simple bounds suppressed.
!         nb>0-simple bounds accepted.
!  ii  nc  number of linear constraints.
!  ri  x(nf)  vector of variables.
!  ii  ix(nf)  vector containing types of bounds. ix(i)=0-variable
!         x(i) is unbounded. ix(i)=1-lover bound xl(i)<=x(i).
!         ix(i)=2-upper bound x(i)<=xu(i). ix(i)=3-two side bound
!         xl(i)<=x(i)<=xu(i). ix(i)=5-variable x(i) is fixed.
!  ri  xl(nf)  vector containing lower bounds for variables.
!  ri  xu(nf)  vector containing upper bounds for variables.
!  ro  cf(nc+1)  vector containing values of the constraint functions.
!  ii  ic(nc)  vector containing types of constraints.
!         ic(kc)=0-constraint cf(kc) is not used. ic(kc)=1-lover
!         constraint cl(kc)<=cf(kc). ic(kc)=2-upper constraint
!         cf(kc)<=cu(kc). ic(kc)=3-two side constraint
!         cl(kc)<=cf(kc)<=cu(kc). ic(kc)=5-equality constraint
!         cf(kc).eq.cl(kc).
!  ri  cl(nc)  vector containing lower bounds for constraint functions.
!  ri  cu(nc)  vector containing upper bounds for constraint functions.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ra  cfo(nc)  vector containing saved values of the constraint
!         functions.
!  ra  cfd(nc)  vector containing increments of the constraint
!         functions.
!  ra  gc(nf)  gradient of the selected constraint function.
!  io  ica(nf)  vector containing indices of active constraints.
!  ro  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ro  cz(nf)  vector of lagrange multipliers.
!  ro  gf(nf)  gradient of the model function.
!  ro  g(nf)  gradient of the objective function.
!  ru  h(nf*(nf+1)/2)  triangular decomposition or inversion of the
!         hessian matrix approximation.
!  ro  s(nf)  direction vector.
!  ru  xo(nf)  vectors of variables difference.
!  ri  go(nf)  gradients difference.
!  ri  xmax  maximum stepsize.
!  ri  tolx  tolerance for change of variables.
!  ri  tolc  tolerance for constraint violations.
!  ri  tolg  tolerance for the gradient of the lagrangian function.
!  ri  rpf  penalty coefficient.
!  ro  cmax  maximum constraint violation.
!  ro  gmax  maximum partial derivative of the lagrangian function.
!  ro  f  value of the objective function.
!         functions.
!  ii  mit  maximum number of iterations.
!  ii  mfv  maximum number of function evaluations.
!  ii  met  variable metric update used. met=1-the bfgs update.
!         met=2-the hoshino update.
!  ii  mec  correction if the negative curvature occurs.
!         mec=1-correction suppressed. mec=2-powell's correction.
!  ii  iprnt  print specification. iprnt=0-no print.
!         abs(iprnt)=1-print of final results.
!         abs(iprnt)=2-print of final results and iterations.
!         iprnt>0-basic final results. iprnt<0-extended final
!         results.
!  io  iterm  variable that indicates the cause of termination.
!         iterm=1-if abs(x-xo) was less than or equal to tolx in
!                   mtesx (usually two) subsequebt iterations.
!         iterm=2-if abs(f-fo) was less than or equal to tolf in
!                   mtesf (usually two) subsequebt iterations.
!         iterm=3-if f is less than or equal to tolb.
!         iterm=4-if gmax is less than or equal to tolg.
!         iterm=11-if nit exceeded mit. iterm=12-if nfv exceeded mfv.
!         iterm=13-if nfg exceeded mfg. iterm<0-if the method failed.
!         if iterm=-6, then the termination criterion has not been
!         satisfied, but the point obtained if usually acceptable.
!
!### Method
! recursive quadratic programming method with the bfgs variable metric
! update.

      subroutine psqp(me,nf,nb,nc,x,ix,xl,xu,cf,ic,cl,cu,cg,cfo,cfd,gc,ica,&
                      cr,cz,cp,gf,g,h,s,xo,go,xmax,tolx,tolc,tolg,rpf,  &
                      cmax,gmax,f,mit,mfv,met,mec,iprnt,iterm)
      implicit none

      class(psqp_class),intent(inout) :: me

      double precision f , cmax , gmax , rpf , tolc , told , tolg ,     &
                       tols , tolx , xmax
      integer iprnt , iterm , met , met1 , mec , mes , mfv , mit , nb , &
              nc , nf
      double precision cf(*) , cfd(*) , cfo(*) , cg(*) , cl(*) , cp(*) ,&
                       cr(*) , cz(*) , cu(*) , g(*) , gc(*) , gf(*) ,   &
                       go(*) , h(*) , s(*) , x(*) , xl(*) , xo(*) ,     &
                       xu(*)
      integer ic(*) , ica(*) , ix(*)
      double precision alf1 , alf2 , cmaxo , dmax , eps7 , eps9 , eta0 ,&
                       eta2 , eta9 , fmax , fmin , fo , gnorm , p , po ,&
                       r , rmax , rmin , ro , snorm , tolb , tolf ,     &
                       umax , rp , fp , pp , ff , fc
      integer i , idecf , iext , irest , iterd , iterl , iterh , iterq ,&
              iters , kbc , kbf , kc , kd , kit , ld , mred , mtesf ,   &
              mtesx , n , k , ntesx , iest , inits , kters , maxst ,    &
              isys , mfp , nred , ipom , lds

      if ( abs(iprnt)>1 ) write (6,'(1x,''entry to psqp :'')')
!
!     initiation
!
      kbf = 0
      kbc = 0
      if ( nb>0 ) kbf = 2
      if ( nc>0 ) kbc = 2
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
      eta0 = 1.0d-15
      eta2 = 1.0d-15
      eta9 = 1.0d60
      eps7 = 1.0d-15
      eps9 = 1.0d-8
      alf1 = 1.0d-10
      alf2 = 1.0d10
      fmax = 1.0d60
      fmin = -fmax
      tolb = -fmax
      dmax = eta9
      tolf = 1.0d-16
      if ( xmax<=0.0d0 ) xmax = 1.0d+16
      if ( tolx<=0.0d0 ) tolx = 1.0d-16
      if ( tolg<=0.0d0 ) tolg = 1.0d-6
      if ( tolc<=0.0d0 ) tolc = 1.0d-6
      told = 1.0d-8
      tols = 1.0d-4
      if ( rpf<=0.0d0 ) rpf = 1.0d-4
      if ( met<=0 ) met = 1
      met1 = 2
      if ( mec<=0 ) mec = 2
      mes = 1
      if ( mit<=0 ) mit = 1000
      if ( mfv<=0 ) mfv = 2000
      kd = 1
      ld = -1
      kit = 0
      call mxvset(nc,0.0d0,cp)
!
!     initial operations with simple bounds
!
      if ( kbf>0 ) then
         do i = 1 , nf
            if ( (ix(i)==3 .or. ix(i)==4) .and. xu(i)<=xl(i) ) then
               xu(i) = xl(i)
               ix(i) = 5
            elseif ( ix(i)==5 .or. ix(i)==6 ) then
               xl(i) = x(i)
               xu(i) = x(i)
               ix(i) = 5
            endif
            if ( ix(i)==1 .or. ix(i)==3 ) x(i) = max(x(i),xl(i))
            if ( ix(i)==2 .or. ix(i)==3 ) x(i) = min(x(i),xu(i))
         enddo
      endif
!     initial operations with general constraints
!
      if ( kbc>0 ) then
         k = 0
         do kc = 1 , nc
            if ( (ic(kc)==3 .or. ic(kc)==4) .and. cu(kc)<=cl(kc) ) then
               cu(kc) = cl(kc)
               ic(kc) = 5
            elseif ( ic(kc)==5 .or. ic(kc)==6 ) then
               cu(kc) = cl(kc)
               ic(kc) = 5
            endif
            k = k + nf
         enddo
      endif
      if ( kbf>0 ) then
         do i = 1 , nf
            if ( ix(i)>=5 ) ix(i) = -ix(i)
            if ( ix(i)<=0 ) then
            elseif ( (ix(i)==1 .or. ix(i)==3) .and. x(i)<=xl(i) ) then
               x(i) = xl(i)
            elseif ( (ix(i)==2 .or. ix(i)==3) .and. x(i)>=xu(i) ) then
               x(i) = xu(i)
            endif
            call plnews(x,ix,xl,xu,eps9,i,iterl)
            if ( ix(i)>10 ) ix(i) = 10 - ix(i)
         enddo
      endif
      fo = fmin
      gmax = eta9
      dmax = eta9
 100  lds = ld
      call me%pf1f01(nf,x,gf,gf,ff,f,kd,ld,iext)
      ld = lds
      call me%pc1f01(nf,nc,x,fc,cf,cl,cu,ic,gc,cg,cmax,kd,ld)
      cf(nc+1) = f
      if ( abs(iprnt)>1 ) &
        write (6,'(1x,''nit='',i9,2x,''nfv='',i9,2x,''nfg='',i9,2x,&
                ''f='',g13.6,2x,''c='',e8.1,2x,''g='',e8.1)') &
                me%nit , me%nfv , me%nfg , f , cmax , gmax
!
!     start of the iteration with tests for termination.
!
      if ( iterm<0 ) goto 500
      if ( iters/=0 ) then
         if ( f<=tolb ) then
            iterm = 3
            goto 500
         endif
         if ( dmax<=tolx ) then
            iterm = 1
            ntesx = ntesx + 1
            if ( ntesx>=mtesx ) goto 500
         else
            ntesx = 0
         endif
      endif
      if ( me%nit>=mit ) then
         iterm = 11
         goto 500
      endif
      if ( me%nfv>=mfv ) then
         iterm = 12
         goto 500
      endif
      iterm = 0
      me%nit = me%nit + 1
!
!     restart
!
 200  n = nf
      if ( irest>0 ) then
         call mxdsmi(n,h)
         ld = min(ld,1)
         idecf = 1
         if ( kit<me%nit ) then
            me%nres = me%nres + 1
            kit = me%nit
         else
            iterm = -10
            if ( iters<0 ) iterm = iters - 5
            goto 500
         endif
      endif
!
!     direction determination using a quadratic programming procedure
!
      call mxvcop(nc+1,cf,cfo)
      mfp = 2
      ipom = 0
 300  call me%plqdb1(nf,nc,x,ix,xl,xu,cf,cfd,ic,ica,cl,cu,cg,cr,cz,g,gf,h, &
                  s,mfp,kbf,kbc,idecf,eta2,eta9,eps7,eps9,umax,gmax,n,  &
                  iterq)
      if ( iterq<0 ) then
         if ( ipom<10 ) then
            ipom = ipom + 1
            call plredl(nc,cf,ic,cl,cu,kbc)
            goto 300
         endif
         iterd = iterq - 10
         goto 400
      endif
      ipom = 0
      iterd = 1
      gmax = mxvmax(nf,g)
      gnorm = sqrt(mxvdot(nf,g,g))
      snorm = sqrt(mxvdot(nf,s,s))
 400  if ( iterd<0 ) iterm = iterd
      if ( iterm==0 ) then
         call mxvcop(nc+1,cfo,cf)
!
!     test for sufficient descent
!
         p = mxvdot(nf,g,s)
         irest = 1
         if ( snorm<=0.0d0 ) then
         elseif ( p+told*gnorm*snorm<=0.0d0 ) then
            irest = 0
         endif
         if ( irest/=0 ) goto 200
         nred = 0
         rmin = alf1*gnorm/snorm
         rmax = min(alf2*gnorm/snorm,xmax/snorm)
         if ( gmax<=tolg .and. cmax<=tolc ) then
            iterm = 4
            goto 500
         endif
         call ppset2(nf,n,nc,ica,cz,cp)
         call mxvina(nc,ic)
         call pp0af8(nf,n,nc,cf,ic,ica,cl,cu,cz,rpf,fc,f)
!
!     preparation of line search
!
         ro = 0.0d0
         fo = f
         po = p
         cmaxo = cmax
         call mxvcop(nf,x,xo)
         call mxvcop(nf,g,go)
         call mxvcop(nf,gf,cr)
         call mxvcop(nc+1,cf,cfo)
!
!     line search without directional derivatives
!
 450     call me%ps0l02(r,ro,rp,f,fo,fp,po,pp,fmin,fmax,rmin,rmax,tols,kd, &
                     ld,me%nit,kit,nred,mred,maxst,iest,inits,iters,kters, &
                     mes,isys)
         if ( isys==0 ) then
            kd = 1
!
!     decision after unsuccessful line search
!
            if ( iters<=0 ) then
               r = 0.0d0
               f = fo
               p = po
               call mxvcop(nf,xo,x)
               call mxvcop(nf,cr,gf)
               call mxvcop(nc+1,cfo,cf)
               irest = 1
               ld = kd
               goto 200
            endif
!
!     computation of the value and the gradient of the objective
!     function together with the values and the gradients of the
!     approximated functions
!
            if ( kd>ld ) then
               lds = ld
               call me%pf1f01(nf,x,gf,gf,ff,f,kd,ld,iext)
               ld = lds
               call me%pc1f01(nf,nc,x,fc,cf,cl,cu,ic,gc,cg,cmax,kd,ld)
            endif
!
!     preparation of variable metric update
!
            call mxvcop(nf,gf,g)
            call pytrnd(nf,n,x,xo,ica,cg,cz,g,go,r,f,fo,p,po,cmax,cmaxo,&
                        dmax,kd,ld,iters)
!
!     variable metric update
!
            call pudbg1(n,h,g,s,xo,go,r,po,me%nit,kit,iterh,met,met1,mec)
!      if (mer.gt.0.and.iterh.gt.0) irest=1
!
!     end of the iteration
!
            goto 100
         else
!      go to (11174,11172) isys+1
            call mxvdir(nf,r,s,xo,x)
            lds = ld
            call me%pf1f01(nf,x,gf,g,ff,f,kd,ld,iext)
            ld = lds
            call me%pc1f01(nf,nc,x,fc,cf,cl,cu,ic,gc,cg,cmax,kd,ld)
            cf(nc+1) = f
            call pp0af8(nf,n,nc,cf,ic,ica,cl,cu,cz,rpf,fc,f)
            goto 450
         endif
      endif

 500  if ( iprnt>1 .or. iprnt<0 ) write (6,'(1x,''exit from psqp :'')')
      if ( iprnt/=0 ) &
         write (6,'(1x,''nit='',i4,2x,''nfv='',i4,2x,''nfg='',i4,2x,''f='',&
                  g13.6,2x,''c='',e8.1,2x,''g='',e8.1,2x,''iterm='',i3)') &
                  me%nit , me%nfv , me%nfg , f , cmax , gmax , iterm
      if ( iprnt<0 ) write (6,'(1x,''x='',5(g14.7,1x):/(3x,5(g14.7,1x)))') (x(i),i=1,nf)
      end subroutine psqp

!***********************************************************************
!> date: 97/12/01
!
! computation of the value and the gradient of the constraint function.
!
! parameters :
!  ii  nf  number of variables.
!  ii  nc  number of constraints.
!  ri  x(nf)  vector of variables.
!  ri  fc  value of the selected constraint function.
!  ri  cf(nc)  vector containing values of constraint functions.
!  ri  cl(nc)  vector containing lower bounds for constraint functions.
!  ri  cu(nc)  vector containing upper bounds for constraint functions.
!  ii  ic(nc)  vector containing types of constraints.
!  ri  gc(nf)  gradient of the selected constraint function.
!  ri  cg(nf*nc)  matrix whose columns are gradients of constraint
!         functions.
!  ro  cmax  maximum constraint violation.
!  ii  kd  degree of required dervatives.
!  ii  ld  degree of previously computed derivatives.

      subroutine pc1f01(me,nf,nc,x,fc,cf,cl,cu,ic,gc,cg,cmax,kd,ld)

      implicit none

      class(psqp_class),intent(inout) :: me
      double precision fc , cmax
      integer kd , ld , nc , nf
      double precision cf(*) , cg(*) , cl(*) , cu(*) , gc(*) , x(*)
      integer ic(*)
      double precision pom , temp
      integer kc
      if ( kd<=ld ) return
      if ( ld<0 ) cmax = 0.0d0
      do kc = 1 , nc
         if ( kd>=0 ) then
            if ( ld>=0 ) then
               fc = cf(kc)
               goto 20
            else
               call me%con(nf,kc,x,fc)
               cf(kc) = fc
            endif
            if ( ic(kc)>0 ) then
               pom = 0.0d0
               temp = cf(kc)
               if ( ic(kc)==1 .or. ic(kc)>=3 ) &
                    pom = min(pom,temp-cl(kc))
               if ( ic(kc)==2 .or. ic(kc)>=3 ) &
                    pom = min(pom,cu(kc)-temp)
               if ( pom<0.0d0 ) cmax = max(cmax,-pom)
            endif
 20         if ( kd>=1 ) then
               if ( ld>=1 ) then
                  call mxvcop(nf,cg((kc-1)*nf+1),gc)
               else
                  call me%dcon(nf,kc,x,gc)
                  call mxvcop(nf,gc,cg((kc-1)*nf+1))
               endif
            endif
         endif
      enddo
      ld = kd
      end subroutine pc1f01

!***********************************************************************
!> date: 97/12/01
!
! computation of the value and the gradient of the objective function.
!
! parameters:
!  ii  nf  number of variables.
!  ri  x(nf)  vector of variables.
!  ri  gf(nf)  gradient of the model function.
!  ro  g(nf)  gradient of the objective function.
!  ri  ff  value of the model function.
!  ro  f  value of the objective function.
!  ii  kd  degree of required derivatives.
!  ii  ld  degree of previously computed derivatives.
!  ii  iext  type of extremum. iext=0-minimum. iext=1-maximum.

      subroutine pf1f01(me,nf,x,gf,g,ff,f,kd,ld,iext)

      implicit none

      class(psqp_class),intent(inout) :: me
      double precision f , ff
      integer iext , kd , ld , nf
      double precision gf(*) , g(*) , x(*)

      if ( kd<=ld ) return
      if ( ld<0 ) then
         me%nfv = me%nfv + 1
         call me%obj(nf,x,ff)
         if ( iext<=0 ) then
            f = ff
         else
            f = -ff
         endif
      endif
      if ( kd>=1 ) then
         if ( ld<1 ) then
            me%nfg = me%nfg + 1
            call me%dobj(nf,x,gf)
            if ( iext>0 ) call mxvneg(nf,gf,g)
         endif
      endif
      ld = kd
      end subroutine pf1f01

!***********************************************************************
!> date: 97/12/01
!
! new linear constraint or a new simple bound is added to the
! active set.
!
! parameters :
!  ii  nf  declared number of variables.
!  iu  n  actual number of variables.
!  iu  ica(nf)  vector containing indices of active constraints.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ru  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ru  cz(nf*nf)  matrix whose columns are basic vectors from the
!         current reduced subspace.
!  ra  s(nf)  auxiliary vector.
!  ri  eps7  tolerance for linear independence of constraints.
!  ro  gmax  maximum absolute value of a partial derivative.
!  ro  umax  maximum absolute value of a negative lagrange multiplier.
!  ii  inew  index of the new active constraint.
!  iu  nadd  number of constraint additions.
!  io  ier  error indicator.

      subroutine pladb0(nf,n,ica,cg,cr,cz,s,eps7,gmax,umax,inew,nadd,ier)
      implicit none
      integer nf , n , ica(*) , inew , nadd , ier
      double precision cg(*) , cr(*) , cz(*) , s(*) , eps7 , gmax , umax
      double precision ck , cl
      integer k , l , n1
      call pladr0(nf,n,ica,cg,cr,s,eps7,gmax,umax,inew,nadd,ier)
      if ( ier/=0 ) return
      if ( n>0 ) then
         n1 = n + 1
         if ( inew>0 ) then
            call mxdrmm(nf,n1,cz,cg((inew-1)*nf+1),s)
         else
            call mxdrmv(nf,n1,cz,s,-inew)
         endif
         do l = 1 , n
            k = l + 1
            call mxvort(s(k),s(l),ck,cl,ier)
            call mxdrgr(nf,cz,k,l,ck,cl,ier)
            if ( ier<0 ) return
         enddo
      endif
      ier = 0
      end subroutine pladb0

!***********************************************************************
!> date: 97/12/01
!
! new linear constraint or a new simple bound is added to the active
! set. transformed hessian matrix approximation or its inversion
! is updated.
!
! parameters :
!  ii  nf  declared number of variables.
!  iu  n  actual number of variables.
!  iu  ica(nf)  vector containing indices of active constraints.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ru  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ru  cz(nf*nf)  matrix whose columns are basic vectors from the
!         current reduced subspace.
!  ru  h(nf*(nf+1)/2)  transformed hessian matrix approximation or
!         its inversion.
!  ra  s(nf)  auxiliary vector.
!  ri  eps7  tolerance for linear independence of constraints.
!  ro  gmax  maximum absolute value of a partial derivative.
!  ro  umax  maximum absolute value of a negative lagrange multiplier.
!  iu  idecf  decomposition indicator. idecf=0-no decomposition.
!         idecf=1-gill-murray decomposition. idecf=9-inversion.
!         idecf=10-diagonal matrix.
!  ii  inew  index of the new active constraint.
!  iu  nadd  number of constraint additions.
!  io  ier  error indicator.

      subroutine pladb4(nf,n,ica,cg,cr,cz,h,s,eps7,gmax,umax,idecf,inew,&
                        nadd,ier)
      implicit none
      integer nf , n , ica(*) , idecf , inew , nadd , ier
      double precision cg(*) , cr(*) , cz(*) , h(*) , s(*) , eps7 ,     &
                       gmax , umax
      double precision ck , cl
      integer i , j , k , l , n1
      if ( idecf/=0 .and. idecf/=9 ) then
         ier = -2
         return
      endif
      call pladr0(nf,n,ica,cg,cr,s,eps7,gmax,umax,inew,nadd,ier)
      if ( ier/=0 ) return
      if ( n>0 ) then
         n1 = n + 1
         if ( inew>0 ) then
            call mxdrmm(nf,n1,cz,cg((inew-1)*nf+1),s)
         else
            call mxdrmv(nf,n1,cz,s,-inew)
         endif
         do l = 1 , n
            k = l + 1
            call mxvort(s(k),s(l),ck,cl,ier)
            call mxdrgr(nf,cz,k,l,ck,cl,ier)
            call mxdsmr(n1,h,k,l,ck,cl,ier)
            if ( ier<0 ) return
         enddo
         if ( idecf==9 ) then
            l = n*(n+1)/2
            if ( h(l+n1)/=0.0d0 ) then
               cl = 1.0d0/h(l+n1)
               k = 0
               do i = 1 , n
                  ck = cl*h(l+i)
                  do j = 1 , i
                     k = k + 1
                     h(k) = h(k) - ck*h(l+j)
                  enddo
               enddo
            endif
         endif
      endif
      ier = 0
      end subroutine pladb4

!***********************************************************************
!> date: 97/12/01
!
! triangular decomposition of kernel of the orthogonal projection
! is updated after constraint addition.
!
! parameters :
!  ii  nf  declared number of variables.
!  iu  n  actual number of variables.
!  iu  ica(nf)  vector containing indices of active constraints.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ru  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ra  s(nf)  auxiliary vector.
!  ri  eps7  tolerance for linear independence of constraints.
!  ro  gmax  maximum absolute value of a partial derivative.
!  ro  umax  maximum absolute value of a negative lagrange multiplier.
!  ii  inew  index of the new active constraint.
!  iu  nadd  number of constraint additions.
!  io  ier  error indicator.

      subroutine pladr0(nf,n,ica,cg,cr,s,eps7,gmax,umax,inew,nadd,ier)
      implicit none
      integer nf , n , ica(*) , inew , nadd , ier
      double precision cg(*) , cr(*) , s(*) , eps7 , gmax , umax
      integer nca , ncr , i , j , k , l
      ier = 0
      if ( n<=0 ) ier = 2
      if ( inew==0 ) ier = 3
      if ( ier/=0 ) return
      nca = nf - n
      ncr = nca*(nca+1)/2
      if ( inew>0 ) then
         call mxvcop(nf,cg((inew-1)*nf+1),s)
         gmax = mxvdot(nf,cg((inew-1)*nf+1),s)
         do j = 1 , nca
            l = ica(j)
            if ( l>0 ) then
               cr(ncr+j) = mxvdot(nf,cg((l-1)*nf+1),s)
            else
               i = -l
               cr(ncr+j) = s(i)
            endif
         enddo
      else
         k = -inew
         gmax = 1.0d0
         do j = 1 , nca
            l = ica(j)
            if ( l>0 ) then
               cr(ncr+j) = cg((l-1)*nf+k)*gmax
            else
               cr(ncr+j) = 0.0d0
            endif
         enddo
      endif
      if ( nca==0 ) then
         umax = gmax
      else
         call mxdprb(nca,cr,cr(ncr+1),1)
         umax = gmax - mxvdot(nca,cr(ncr+1),cr(ncr+1))
      endif
      if ( umax<=eps7*gmax ) then
         ier = 1
         return
      else
         n = n - 1
         nca = nca + 1
         ncr = ncr + nca
         ica(nca) = inew
         cr(ncr) = sqrt(umax)
         nadd = nadd + 1
      endif
      end subroutine pladr0

!***********************************************************************
!> date: 97/12/01
!
! dual range space quadratic programming method.
!
! parameters :
!  ii  nf  number of variables.
!  io  n  dimension of the manifold defined by active constraints.
!  ii  nc  number of linear constraints.
!  ru  x(nf)   vector of variables.
!  ii  ix(nf)  vector containing types of bounds.
!  ii  ixa(nf)  vector containing information on trust region activity.
!  ri  xn(nf)  vector of scaling factors.
!  ri  xl(nf)  vector containing lower bounds for variables.
!  ri  xu(nf)  vector containing upper bounds for variables.
!  ri  cf(nf)  vector containing values of the constraint functions.
!  ro  cfd(nc)  vector containing increments of the constraint
!            functions.
!  ii  ic(nc)  vector containing types of constraints.
!  ii  ica(nf)  vector containing indices of active constraints.
!  ri  cl(nc)  vector containing lower bounds for constraint functions.
!  ri  cu(nc)  vector containing upper bounds for constraint functions.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ri  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ro  cz(nf)  vector of lagrange multipliers.
!  ro  g(nf)  gradient of the lagrangian function.
!  ro  go(nf)  saved gradient of the objective function.
!  ru  h(nf*(nf+1)/2)  triangular decomposition or inversion of the
!         hessian matrix approximation.
!  ri  s(nf)  direction vector.
!  ri  eta2  tolerance for positive definiteness of the hessian matrix.
!  ri  eta9  maximum for real numbers.
!  ri  eps7  tolerance for linear independence of constraints.
!  ri  eps9  tolerance for activity of constraints.
!  ru  xdel  trust region bound.
!  ro  umax  maximum absolute value of a negative lagrange multiplier.
!  ro  gmax  maximum absolute value of a partial derivative.
!  ii  mfp  type of feasible point. mfp=1-arbitrary feasible point.
!         mfp=2-optimum feasible point. mfp=3-repeated solution.
!
! common data :
!  ii  kbf  specification of simple bounds. kbf=0-no simple bounds.
!         kbf=1-one sided simple bounds. kbf=2=two sided simple bounds.
!  ii  kbc  specification of linear constraints. kbc=0-no linear
!         constraints. kbc=1-one sided linear constraints. kbc=2=two
!         sided linear constraints.
!  ii  normf  scaling specification. normf=0-no scaling performed.
!         normf=1-scaling factors are determined automatically.
!         normf=2-scaling factors are supplied by user.
!  iu  idecf  decomposition indicator. idecf=0-no decomposition.
!         idecf=1-gill-murray decomposition. idecf=9-inversion.
!         idecf=10-diagonal matrix.
!  iu  ndecf  number of decompositions.
!  io  iterq  type of feasible point. iterq=1-arbitrary feasible point.
!         iterq=2-optimum feasible point. iterq=-1 feasible point does
!         not exists. iterq=-2 optimum feasible point does not exists.

      subroutine plqdb1(me,nf,nc,x,ix,xl,xu,cf,cfd,ic,ica,cl,cu,cg,cr,cz,g,&
                        go,h,s,mfp,kbf,kbc,idecf,eta2,eta9,eps7,eps9,   &
                        umax,gmax,n,iterq)
      implicit none

      class(psqp_class),intent(inout) :: me

      integer nf , nc , ix(*) , ic(*) , ica(*) , mfp , kbf , kbc ,      &
              idecf , n , iterq
      double precision x(*) , xl(*) , xu(*) , cf(*) , cfd(*) , cl(*) ,  &
                       cu(*) , cg(*) , cr(*) , cz(*) , g(*) , go(*) ,   &
                       h(*) , s(*) , eta2 , eta9 , eps7 , eps9 , umax , &
                       gmax
      double precision con , temp , step , step1 , step2 , dmax , par , &
                       snorm
      integer nca , ncr , i , j , k , iold , jold , inew , jnew , knew ,&
              inf , ier , krem , kc , nred

      con = eta9
      if ( idecf<0 ) idecf = 1
      if ( idecf==0 ) then
!
!     gill-murray decomposition
!
         temp = eta2
         call mxdpgf(nf,h,inf,temp,step)
         me%ndec = me%ndec + 1
         idecf = 1
      endif
      if ( idecf>=2 .and. idecf<=8 ) then
         iterq = -10
         return
      endif
!
!     initiation
!
      nred = 0
      jold = 0
      jnew = 0
      iterq = 0
      dmax = 0.0d0
      if ( mfp/=3 ) then
         n = nf
         nca = 0
         ncr = 0
         if ( kbf>0 ) call mxvina(nf,ix)
         if ( kbc>0 ) call mxvina(nc,ic)
      endif
!
!     direction determination
!
 100  call mxvneg(nf,go,s)
      do j = 1 , nca
         kc = ica(j)
         if ( kc>0 ) then
            call mxvdir(nf,cz(j),cg((kc-1)*nf+1),s,s)
         else
            k = -kc
            s(k) = s(k) + cz(j)
         endif
      enddo
      call mxvcop(nf,s,g)
      if ( idecf==1 ) then
         call mxdpgb(nf,h,s,0)
      else
         call mxdsmm(nf,h,g,s)
      endif
      if ( iterq/=3 ) then
!
!     check of feasibility
!
         inew = 0
         par = 0.0d0
         call plminn(nf,nc,cf,cfd,ic,cl,cu,cg,s,eps9,par,kbc,inew,knew)
         call plmins(nf,ix,x,xl,xu,s,kbf,inew,knew,eps9,par)
         if ( inew==0 ) then
!
!     solution achieved
!
            call mxvneg(nf,g,g)
            iterq = 2
            return
         else
            snorm = 0.0d0
         endif
 150     ier = 0
!
!     stepsize determination
!
         call pladr1(nf,n,ica,cg,cr,h,s,g,eps7,gmax,umax,idecf,inew,me%nadd,ier,1)
         call mxdprb(nca,cr,g,-1)
         if ( knew<0 ) call mxvneg(nca,g,g)
!
!     primal stepsize
!
         if ( ier/=0 ) then
            step1 = con
         else
            step1 = -par/umax
         endif
!
!     dual stepsize
!
         iold = 0
         step2 = con
         do j = 1 , nca
            kc = ica(j)
            if ( kc>=0 ) then
               k = ic(kc)
            else
               i = -kc
               k = ix(i)
            endif
            if ( k<=-5 ) then
            elseif ( (k==-1 .or. k==-3.) .and. g(j)<=0.0d0 ) then
            elseif ( .not.((k==-2 .or. k==-4.) .and. g(j)>=0.0d0) ) then
               temp = cz(j)/g(j)
               if ( step2>temp ) then
                  iold = j
                  step2 = temp
               endif
            endif
         enddo
!
!     final stepsize
!
         step = min(step1,step2)
         if ( step>=con ) then
!
!     feasible solution does not exist
!
            iterq = -1
            return
         endif
!
!     new lagrange multipliers
!
         dmax = step
         call mxvdir(nca,-step,g,cz,cz)
         snorm = snorm + sign(1,knew)*step
         par = par - (step/step1)*par
         if ( step==step1 ) then
            if ( n<=0 ) then
!
!     impossible situation
!
               iterq = -5
               return
            endif
!
!     constraint addition
!
            if ( ier==0 ) then
               n = n - 1
               nca = nca + 1
               ncr = ncr + nca
               cz(nca) = snorm
            endif
            if ( inew>0 ) then
               kc = inew
               call mxvinv(ic,kc,knew)
            elseif ( abs(knew)==1 ) then
               i = -inew
               call mxvinv(ix,i,knew)
            else
               i = -inew
               if ( knew>0 ) ix(i) = -3
               if ( knew<0 ) ix(i) = -4
            endif
            nred = nred + 1
            me%nadd = me%nadd + 1
            jnew = inew
            jold = 0
            goto 100
         else
!
!     constraint deletion
!
            do j = iold , nca - 1
               cz(j) = cz(j+1)
            enddo
            call me%plrmf0(nf,nc,ix,ic,ica,cr,ic,g,n,iold,krem,ier)
            ncr = ncr - nca
            nca = nca - 1
            jold = iold
            jnew = 0
            if ( kbc>0 ) call mxvina(nc,ic)
            if ( kbf>0 ) call mxvina(nf,ix)
            do j = 1 , nca
               kc = ica(j)
               if ( kc>0 ) then
                  ic(kc) = -ic(kc)
               else
                  kc = -kc
                  ix(kc) = -ix(kc)
               endif
            enddo
            goto 150
         endif
      endif
      end subroutine plqdb1

!***********************************************************************
!> date: 97/12/01
!
! triangular decomposition of kernel of the general projection
! is updated after constraint addition.
!
! parameters :
!  ii  nf  declared number of variables.
!  iu  n  actual number of variables.
!  iu  ica(nf)  vector containing indices of active constraints.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ri  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ru  h(nf*(nf+1)/2)  triangular decomposition or inversion of the
!         hessian matrix approximation.
!  ra  s(nf)  auxiliary vector.
!  ro  g(nf)  vector used in the dual range space quadratic programming
!         method.
!  ri  eps7  tolerance for linear independence of constraints.
!  ro  gmax  maximum absolute value of a partial derivative.
!  ro  umax  maximum absolute value of a negative lagrange multiplier.
!  ro  e  auxiliary variable.
!  ri  t  auxiliary variable.
!  iu  idecf  decomposition indicator. idecf=0-no decomposition.
!         idecf=1-gill-murray decomposition. idecf=9-inversion.
!         idecf=10-diagonal matrix.
!  ii  inew  index of the new active constraint.
!  iu  nadd  number of constraint additions.
!  io  ier  error indicator.
!  ii  job  specification of computation. output vector g is not or is
!         computed in case when n<=0 if job=0 or job=1 respectively.

      subroutine pladr1(nf,n,ica,cg,cr,h,s,g,eps7,gmax,umax,idecf,inew, &
                        nadd,ier,job)
      implicit none
      integer nf , n , ica(*) , idecf , inew , nadd , ier , job
      double precision cg(*) , cr(*) , h(*) , s(*) , g(*) , eps7 ,      &
                       gmax , umax
      integer nca , ncr , jcg , j , k , l

      ier = 0
      if ( job==0 .and. n<=0 ) ier = 2
      if ( inew==0 ) ier = 3
      if ( idecf/=1 .and. idecf/=9 ) ier = -2
      if ( ier/=0 ) return
      nca = nf - n
      ncr = nca*(nca+1)/2
      if ( inew>0 ) then
         jcg = (inew-1)*nf + 1
         if ( idecf==1 ) then
            call mxvcop(nf,cg(jcg),s)
            call mxdpgb(nf,h,s,0)
         else
            call mxdsmm(nf,h,cg(jcg),s)
         endif
         gmax = mxvdot(nf,cg(jcg),s)
      else
         k = -inew
         if ( idecf==1 ) then
            call mxvset(nf,0.0d0,s)
            s(k) = 1.0d0
            call mxdpgb(nf,h,s,0)
         else
            call mxdsmv(nf,h,s,k)
         endif
         gmax = s(k)
      endif
      do j = 1 , nca
         l = ica(j)
         if ( l>0 ) then
            g(j) = mxvdot(nf,cg((l-1)*nf+1),s)
         else
            l = -l
            g(j) = s(l)
         endif
      enddo
      if ( n==0 ) then
         call mxdprb(nca,cr,g,1)
         umax = 0.0d0
         ier = 2
         return
      elseif ( nca==0 ) then
         umax = gmax
      else
         call mxdprb(nca,cr,g,1)
         umax = gmax - mxvdot(nca,g,g)
         call mxvcop(nca,g,cr(ncr+1))
      endif
      if ( umax<=eps7*gmax ) then
         ier = 1
         return
      else
         nca = nca + 1
         ncr = ncr + nca
         ica(nca) = inew
         cr(ncr) = sqrt(umax)
         if ( job==0 ) then
            n = n - 1
            nadd = nadd + 1
         endif
      endif
      end subroutine pladr1

!***********************************************************************
!> date: 97/12/01
!
! determination of the new active linear constraint.
!
! parameters :
!  ii  nf  number of variables.
!  ii  nc  number of constraints.
!  ri  cf(nc)  vector containing values of the constraint functions.
!  ro  cfd(nc)  vector containing increments of the constraint
!            functions.
!  ii  ic(nc)  vector containing types of constraints.
!  ri  cl(nc)  vector containing lower bounds for constraint functions.
!  ri  cu(nc)  vector containing upper bounds for constraint functions.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ri  s(nf)  direction vector.
!  ri  eps9  tolerance for active constraints.
!  ra  par  auxiliary variable.
!  ii  kbc  specification of linear constraints. kbc=0-no linear
!         constraints. kbc=1-one sided linear constraints. kbc=2=two
!         sided linear constraints.
!  io  inew  index of the new active constraint.
!  io  knew  signum of the new active normal.

      subroutine plminn(nf,nc,cf,cfd,ic,cl,cu,cg,s,eps9,par,kbc,inew,   &
                        knew)
      implicit none
      integer nf , nc , ic(*) , kbc , inew , knew
      double precision cf(*) , cfd(*) , cl(*) , cu(*) , cg(*) , s(*) ,  &
                       eps9 , par
      double precision temp , pom !, mxvdot
      integer jcg , kc
      if ( kbc>0 ) then
         jcg = 1
         do kc = 1 , nc
            if ( ic(kc)>0 ) then
               temp = mxvdot(nf,cg(jcg),s)
               cfd(kc) = temp
               temp = cf(kc) + temp
               if ( ic(kc)==1 .or. ic(kc)>=3 ) then
                  pom = temp - cl(kc)
                  if ( pom<min(par,-eps9*max(abs(cl(kc)),1.0d0)) ) then
                     inew = kc
                     knew = 1
                     par = pom
                  endif
               endif
               if ( ic(kc)==2 .or. ic(kc)>=3 ) then
                  pom = cu(kc) - temp
                  if ( pom<min(par,-eps9*max(abs(cu(kc)),1.0d0)) ) then
                     inew = kc
                     knew = -1
                     par = pom
                  endif
               endif
            endif
            jcg = jcg + nf
         enddo
      endif
      end subroutine plminn

!***********************************************************************
!> date: 91/12/01
!
! determination of the new active simple bound.
!
! parameters :
!  ii  nf declared number of variables.
!  ii  ix(nf)  vector containing types of bounds.
!  ri  xo(nf)  saved vector of variables.
!  ri  xl(nf)  vector containing lower bounds for variables.
!  ri  xu(nf)  vector containing upper bounds for variables.
!  ri  s(nf)  direction vector.
!  ii  kbf  specification of simple bounds. kbf=0-no simple bounds.
!         kbf=1-one sided simple bounds. kbf=2=two sided simple bounds.
!  io  inew  index of the new active constraint.
!  io  knew  signum of the new normal.
!  ri  eps9  tolerance for active constraints.
!  ra  par  auxiliary variable.
!
      subroutine plmins(nf,ix,xo,xl,xu,s,kbf,inew,knew,eps9,par)
      implicit none
      double precision eps9 , par
      integer inew , kbf , knew , nf
      double precision s(*) , xl(*) , xo(*) , xu(*)
      integer ix(*)
      double precision pom , temp
      integer i
      if ( kbf>0 ) then
         do i = 1 , nf
            if ( ix(i)>0 ) then
               temp = 1.0d0
               if ( ix(i)==1 .or. ix(i)>=3 ) then
                  pom = xo(i) + s(i)*temp - xl(i)
                  if ( pom<min(par,-eps9*max(abs(xl(i)),temp)) ) then
                     inew = -i
                     knew = 1
                     par = pom
                  endif
               endif
               if ( ix(i)==2 .or. ix(i)>=3 ) then
                  pom = xu(i) - s(i)*temp - xo(i)
                  if ( pom<min(par,-eps9*max(abs(xu(i)),temp)) ) then
                     inew = -i
                     knew = -1
                     par = pom
                  endif
               endif
            endif
         enddo
      endif
      end subroutine plmins

!***********************************************************************
!> date: 97/12/01
!
! test on activity of a given simple bound.
!
! parameters :
!  ri  x(nf)  vector of variables.
!  iu  ix(nf)  vector containing types of bounds.
!  ri  xl(nf)  vector containing lower bounds for variables.
!  ri  xu(nf)  vector containing upper bounds for variables.
!  ri  eps9  tolerance for active constraints.
!  ii  i  index of tested simple bound.
!  io  inew  index of the new active constraint.
!
      subroutine plnews(x,ix,xl,xu,eps9,i,inew)
      implicit none
      integer ix(*) , i , inew
      double precision x(*) , xl(*) , xu(*) , eps9
      double precision temp
      temp = 1.0d0
      if ( ix(i)<=0 ) then
      elseif ( ix(i)==1 ) then
         if ( x(i)<=xl(i)+eps9*max(abs(xl(i)),temp) ) then
            ix(i) = 11
            inew = -i
         endif
      elseif ( ix(i)==2 ) then
         if ( x(i)>=xu(i)-eps9*max(abs(xu(i)),temp) ) then
            ix(i) = 12
            inew = -i
         endif
      elseif ( ix(i)==3 .or. ix(i)==4 ) then
         if ( x(i)<=xl(i)+eps9*max(abs(xl(i)),temp) ) then
            ix(i) = 13
            inew = -i
         endif
         if ( x(i)>=xu(i)-eps9*max(abs(xu(i)),temp) ) then
            ix(i) = 14
            inew = -i
         endif
      endif
      end subroutine plnews

!***********************************************************************
!> date: 98/12/01
!
! transformation of the incompatible quadratic programming subproblem.
!
! parameters :
!  ii  nc  number of current linear constraints.
!  ri  cf(nf)  vector containing values of the constraint functions.
!  ii  ic(nc)  vector containing types of constraints.
!  ri  cl(nc)  vector containing lower bounds for constraint functions.
!  ri  cu(nc)  vector containing upper bounds for constraint functions.
!  ii  kbc  specification of linear constraints. kbc=0-no linear
!         constraints. kbc=1-one sided linear constraints. kbc=2=two
!         sided linear constraints.
!
      subroutine plredl(nc,cf,ic,cl,cu,kbc)
      implicit none
      integer nc , ic(nc) , kbc , k
      double precision cf(*) , cl(*) , cu(*)
      double precision temp
      integer kc
      if ( kbc>0 ) then
         do kc = 1 , nc
            k = ic(kc)
            if ( abs(k)==1 .or. abs(k)==3 .or. abs(k)==4 ) then
               temp = (cf(kc)-cl(kc))
               if ( temp<0 ) cf(kc) = cl(kc) + 0.1d0*temp
            endif
            if ( abs(k)==2 .or. abs(k)==3 .or. abs(k)==4 ) then
               temp = (cf(kc)-cu(kc))
               if ( temp>0 ) cf(kc) = cu(kc) + 0.1d0*temp
            endif
            if ( abs(k)==5 .or. abs(k)==6 ) then
               temp = (cf(kc)-cl(kc))
               cf(kc) = cl(kc) + 0.1d0*temp
            endif
         enddo
      endif
      end subroutine plredl

!***********************************************************************
!> date: 91/12/01
!
! operations after constraint deletion.
!
! parameters :
!  ii  nf  declared number of variables.
!  ii  nc  number of constraints.
!  ii  ix(nf)  vector containing types of bounds.
!  ii  ia(na)  vector containing types of deviations.
!  iu  iaa(nf+1)  vector containing indices of active functions.
!  ru  ar((nf+1)*(nf+2)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ii  ic(nc)  vector containing types of constraints.
!  ra  s(nf+1)  auxiliary vector.
!  ii  n  actual number of variables.
!  ii  iold  index of the old active constraint.
!  io  krem  auxiliary variable.
!  io  ier  error indicator.

      subroutine plrmf0(me,nf,nc,ix,ia,iaa,ar,ic,s,n,iold,krem,ier)
      implicit none

      class(psqp_class),intent(inout) :: me

      integer ier , iold , krem , n , nc , nf
      double precision ar(*) , s(*)
      integer ia(*) , iaa(*) , ic(*) , ix(*)
      integer l

      call plrmr0(nf,iaa,ar,s,n,iold,krem,ier)
      n = n + 1
      me%nrem = me%nrem + 1
      l = iaa(nf-n+1)
      if ( l>nc ) then
         l = l - nc
         ia(l) = -ia(l)
      elseif ( l>0 ) then
         ic(l) = -ic(l)
      else
         l = -l
         ix(l) = -ix(l)
      endif
      end subroutine plrmf0

!***********************************************************************
!> date: 91/12/01
!
! triangular decomposition of kernel of the orthogonal projection is
! updated after constraint deletion.
!
! parameters :
!  ii  nf  declared number of variables.
!  iu  ica(nf)  vector containing indices of active constraints.
!  ru  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ra  g(nf)  auxiliary vector.
!  ii  n  actual number of variables.
!  ii  iold  index of the old active constraint.
!  io  krem  auxiliary variable.
!  io  ier  error indicator.

      subroutine plrmr0(nf,ica,cr,g,n,iold,krem,ier)
      implicit none
      integer ier , iold , krem , n , nf
      double precision cr(*) , g(*)
      integer ica(*)
      double precision ck , cl
      integer i , j , k , kc , l , nca
      nca = nf - n
      if ( iold<nca ) then
         k = iold*(iold-1)/2
         kc = ica(iold)
         call mxvcop(iold,cr(k+1),g)
         call mxvset(nca-iold,0.0d0,g(iold+1))
         k = k + iold
         do i = iold + 1 , nca
            k = k + i
            call mxvort(cr(k-1),cr(k),ck,cl,ier)
            call mxvrot(g(i-1),g(i),ck,cl,ier)
            l = k
            do j = i , nca - 1
               l = l + j
               call mxvrot(cr(l-1),cr(l),ck,cl,ier)
            enddo
         enddo
         k = iold*(iold-1)/2
         do i = iold , nca - 1
            l = k + i
            ica(i) = ica(i+1)
            call mxvcop(i,cr(l+1),cr(k+1))
            k = l
         enddo
         ica(nca) = kc
         call mxvcop(nca,g,cr(k+1))
      endif
      krem = 1
      end subroutine plrmr0

!***********************************************************************
!> date: 97/12/01
!
! determination of initial values of the constraint functions.
!
! parameters :
!  ii  nf  number of variables.
!  ii  nc  number of current linear constraints.
!  ri  x(nf)  vector of variables.
!  ri  xo(nf)  saved vector of variables.
!  ru  cf(nf)  vector containing values of the constraint funcyions.
!  ii  ic(nc)  vector containing types of constraints.
!  ri  cg(nf*mcl)  matrix whose columns are normals of the linear
!         constraints.
!  ra  s(nf)  auxiliary vector.

      subroutine plsetc(nf,nc,x,xo,cf,ic,cg,s)
      implicit none
      integer nf , nc , ic(*)
      double precision x(*) , xo(*) , cf(*) , cg(*) , s(*)
      integer jcg , kc

      call mxvdif(nf,x,xo,s)
      jcg = 0
      do kc = 1 , nc
         if ( ic(kc)/=0 ) cf(kc) = cf(kc) + mxvdot(nf,s,cg(jcg+1))
         jcg = jcg + nf
      enddo
      end subroutine plsetc

!***********************************************************************
!> date: 97/12/01
!
! gradient determination in the first phase of lp subroutine.
!
! parameters :
!  ii  nf  declared number of variables.
!  ii  nc  number of constraints.
!  ii  ic(nc)  vector containing types of constraints.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ro  g(nf)  gradient of the objective function.
!  io  inew  index of the new active constraint.

      subroutine plsetg(nf,nc,ic,cg,g,inew)
      implicit none
      integer nf , nc , ic(*) , inew
      double precision cg(*) , g(*)
      integer kc
      call mxvset(nf,0.0d0,g)
      inew = 0
      do kc = 1 , nc
         if ( ic(kc)>=-10 ) then
         elseif ( ic(kc)==-11 .or. ic(kc)==-13 .or. ic(kc)==-15 ) then
            call mxvdir(nf,-1.0d0,cg((kc-1)*nf+1),g,g)
            inew = 1
         elseif ( ic(kc)==-12 .or. ic(kc)==-14 .or. ic(kc)==-16 ) then
            call mxvdir(nf,1.0d0,cg((kc-1)*nf+1),g,g)
            inew = 1
         endif
      enddo
      end subroutine plsetg

!***********************************************************************
!> date: 97/12/01
!
! maximum absolute value of the negative lagrange multiplier is
! computed.
!
! parameters :
!  ii  nf  declared number of variables.
!  ii  n  actual number of variables.
!  ii  nc  number of linearized constraints.
!  ii  ix(nf)  vector containing types of bounds.
!  ii  ia(na)  vector containing types of deviations.
!  ii  iaa(nf+1)  vector containing indices of active functions.
!  ri  az(nf+1)  vector of lagrange multipliers.
!  ii  ic(nc)  vector containing types of constraints.
!  ri  eps7  tolerance for linear and quadratic programming.
!  ro  umax  maximum absolute value of the negative lagrange multiplier.
!  io  iold  index of the removed constraint.
!
      subroutine pltlag(nf,n,nc,ix,ia,iaa,az,ic,eps7,umax,iold)
      implicit none
      integer nf , n , nc , ix(*) , ia(*) , iaa(*) , ic(*) , iold
      double precision az(*) , eps7 , umax
      double precision temp
      integer naa , j , k , l
      iold = 0
      umax = 0.0d0
      naa = nf - n
      do j = 1 , naa
         temp = az(j)
         l = iaa(j)
         if ( l>nc ) then
            l = l - nc
            k = ia(l)
         elseif ( l>0 ) then
            k = ic(l)
         else
            l = -l
            k = ix(l)
         endif
         if ( k<=-5 ) then
         elseif ( (k==-1 .or. k==-3) .and. umax+temp>=0.0d0 ) then
         elseif ( .not.((k==-2 .or. k==-4) .and. umax-temp>=0.0d0) )    &
                  then
            iold = j
            umax = abs(temp)
         endif
      enddo
      if ( umax<=eps7 ) iold = 0
      end subroutine pltlag

!***********************************************************************
!> date: 97/12/01
!
! gradient of the objective function is scaled and reduced. lagrange
! multipliers are determined. test values gmax and umax are computed.
!
! parameters :
!  ii  nf  declared number of variables.
!  ii  n  actual number of variables.
!  ii  nc  number of current linear constraints.
!  ii  ix(nf)  vector containing types of bounds.
!  ii  ic(nc)  vector containing types of constraints.
!  ii  ica(nf)  vector containing indices of active constraints.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ri  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ru  cz(nf*nf)  matrix whose columns are basic vectors from the
!         current reduced subspace. vector cz(1,nf) contains lagrange
!         multipliers being determined.
!  ri  g(nf)  gradient of the objective function.
!  ro  gn(nf)  transformed gradient of the objective function if it is
!         nonzero.
!  ri  eps7  tolerance for linear and quadratic programming.
!  ro  gmax  norm of the transformed gradient.
!  ro  umax  maximum absolute value of the negative lagrange multiplier.
!  io  iold  index of the removed constraint.

      subroutine pltrbg(nf,n,nc,ix,ic,ica,cg,cr,cz,g,gn,eps7,gmax,umax, &
                        iold)
      implicit none
      integer nf , n , nc , ix(*) , ic(*) , ica(*) , iold
      double precision cg(*) , cr(*) , cz(*) , g(*) , gn(*) , eps7 ,    &
                       gmax , umax
      integer nca , ncz

      gmax = 0.0d0
      if ( n>0 ) then
         call mxdrmm(nf,n,cz,g,gn)
         gmax = mxvmax(n,gn)
      endif
      if ( nf>n .and. gmax<=eps7 ) then
         nca = nf - n
         ncz = n*nf
         call plvlag(nf,n,nc,ica,cg,cg,g,cz(ncz+1))
         call mxdprb(nca,cr,cz(ncz+1),0)
         call pltlag(nf,n,nc,ix,ic,ica,cz(ncz+1),ic,eps7,umax,iold)
         if ( umax<=eps7 ) iold = 0
         call mxvset(n,0.0d0,gn)
         gmax = 0.0d0
      else
         iold = 0
         umax = 0.0d0
      endif
      end subroutine pltrbg

!***********************************************************************
!> date: 97/12/01
!
! gradient of the objective function is premultiplied by transpose
! of the matrix whose columns are normals of current active constraints
! and gradients of current active functions.
!
! parameters :
!  ii  nf  declared number of variables.
!  ii  n  actual number of variables.
!  ii  nc  number of linearized constraints.
!  ii  iaa(nf+1)  vector containing indices of active functions.
!  ri  ag(nf*na)  vector containing scaling parameters.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ri  g(nf)  gradient of the objective function.
!  ro  gn(nf+1)  output vector.

      subroutine plvlag(nf,n,nc,iaa,ag,cg,g,gn)
      implicit none
      integer nf , n , nc , iaa(*)
      double precision ag(*) , cg(*) , g(*) , gn(*)
      integer naa , j , l

      naa = nf - n
      do j = 1 , naa
         l = iaa(j)
         if ( l>nc ) then
            l = l - nc
            gn(j) = mxvdot(nf,ag((l-1)*nf+1),g)
         elseif ( l>0 ) then
            gn(j) = mxvdot(nf,cg((l-1)*nf+1),g)
         else
            l = -l
            gn(j) = g(l)
         endif
      enddo
      end subroutine plvlag

!***********************************************************************
!> date: 91/12/01
!
! extrapolation or interpolation for line search with directional
! derivatives.
!
! parameters :
!  ri  rl  lower value of the stepsize parameter.
!  ri  ru  upper value of the stepsize parameter.
!  ri  fl  value of the objective function for r=rl.
!  ri  fu  value of the objective function for r=ru.
!  ri  pl  directional derivative for r=rl.
!  ri  pu  directional derivative for r=ru.
!  ro  r  value of the stepsize parameter obtained.
!  ii  mode  mode of line search.
!  ii  mtyp  method selection. mtyp=1-bisection. mtyp=2-quadratic
!         interpolation (with one directional derivative).
!         mtyp=3-quadratic interpolation (with two directional
!         derivatives). mtyp=4-cubic interpolation. mtyp=5-conic
!         interpolation.
!  io  merr  error indicator. merr=0 for normal return.
!
!### Method
! extrapolation or interpolation with standard model functions.
!
      subroutine pnint1(rl,ru,fl,fu,pl,pu,r,mode,mtyp,merr)
      implicit none
      double precision rl , ru , fl , fu , pl , pu , r
      integer mode , mtyp , merr , ntyp
      double precision a , b , c , d , dis , den
      double precision c1l , c1u , c2l , c2u , c3l
      parameter (c1l=1.1d0,c1u=1.0d3,c2l=1.0d-2,c2u=0.9d0,c3l=0.1d0)
      merr = 0
      if ( mode<=0 ) return
      if ( pl>=0.0d0 ) then
         merr = 2
         return
      elseif ( ru<=rl ) then
         merr = 3
         return
      endif
      do ntyp = mtyp , 1 , -1
         if ( ntyp==1 ) then
!
!     bisection
!
            if ( mode==1 ) then
               r = 4.0d0*ru
               return
            else
               r = 0.5d0*(rl+ru)
               return
            endif
         elseif ( ntyp==mtyp ) then
            a = (fu-fl)/(pl*(ru-rl))
            b = pu/pl
         endif
         if ( ntyp==2 ) then
!
!     quadratic extrapolation or interpolation with one directional
!     derivative
!
            den = 2.0d0*(1.0d0-a)
         elseif ( ntyp==3 ) then
!
!     quadratic extrapolation or interpolation with two directional
!     derivatives
!
            den = 1.0d0 - b
         elseif ( ntyp==4 ) then
!
!     cubic extrapolation or interpolation
!
            c = b - 2.0d0*a + 1.0d0
            d = b - 3.0d0*a + 2.0d0
            dis = d*d - 3.0d0*c
            if ( dis<0.0d0 ) goto 100
            den = d + sqrt(dis)
         elseif ( ntyp==5 ) then
!
!     conic extrapolation or interpolation
!
            dis = a*a - b
            if ( dis<0.0d0 ) goto 100
            den = a + sqrt(dis)
            if ( den<=0.0d0 ) goto 100
            den = 1.0d0 - b*(1.0d0/den)**3
         endif
         if ( mode==1 .and. den>0.0d0 .and. den<1.0d0 ) then
!
!     extrapolation accepted
!
            r = rl + (ru-rl)/den
            r = max(r,c1l*ru)
            r = min(r,c1u*ru)
            return
         elseif ( mode==2 .and. den>1.0d0 ) then
!
!     interpolation accepted
!
            r = rl + (ru-rl)/den
            if ( rl==0.0d0 ) then
               r = max(r,rl+c2l*(ru-rl))
            else
               r = max(r,rl+c3l*(ru-rl))
            endif
            r = min(r,rl+c2u*(ru-rl))
            return
         endif
 100  enddo
      end subroutine pnint1

!***********************************************************************
!> date: 91/12/01
!
! extrapolation or interpolation for line search without directional
! derivatives.
!
! parameters :
!  ri  ro  initial value of the stepsize parameter.
!  ri  rl  lower value of the stepsize parameter.
!  ri  ru  upper value of the stepsize parameter.
!  ri  ri  inner value of the stepsize parameter.
!  ri  fo  value of the objective function for r=ro.
!  ri  fl  value of the objective function for r=rl.
!  ri  fu  value of the objective function for r=ru.
!  ri  fi  value of the objective function for r=ri.
!  ro  po  initial value of the directional derivative.
!  ro  r  value of the stepsize parameter obtained.
!  ii  mode  mode of line search.
!  ii  mtyp  method selection. mtyp=1-bisection. mtyp=2-two point
!         quadratic interpolation. mtyp=2-three point quadratic
!         interpolation.
!  io  merr  error indicator. merr=0 for normal return.
!
!### Method
! extrapolation or interpolation with standard model functions.
!
      subroutine pnint3(ro,rl,ru,ri,fo,fl,fu,fi,po,r,mode,mtyp,merr)
      implicit none
      double precision zero , half , one , two , three , c1l , c1u ,    &
                       c2l , c2u , c3l
      parameter (zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0,three=3.0d0, &
                 c1l=1.1d0,c1u=1.0d3,c2l=1.0d-2,c2u=0.9d0,c3l=1.0d-1)
      double precision fi , fl , fo , fu , po , r , ri , rl , ro , ru
      integer merr , mode , mtyp
      double precision ai , al , au , den , dis
      integer ntyp
      logical l1 , l2
      merr = 0
      if ( mode<=0 ) return
      if ( po>=zero ) then
         merr = 2
         return

      elseif ( ru<=rl ) then
         merr = 3
         return
      endif
      l1 = rl<=ro
      l2 = ri<=rl
      do ntyp = mtyp , 1 , -1
         if ( ntyp==1 ) then
!
!     bisection
!
            if ( mode==1 ) then
               r = two*ru
               return
            elseif ( ri-rl<=ru-ri ) then
               r = half*(ri+ru)
               return
            else
               r = half*(rl+ri)
               return
            endif
         elseif ( ntyp==mtyp .and. l1 ) then
            if ( .not.l2 ) ai = (fi-fo)/(ri*po)
            au = (fu-fo)/(ru*po)
         endif
         if ( l1 .and. (ntyp==2 .or. l2) ) then
!
!     two point quadratic extrapolation or interpolation
!
            if ( au>=one ) goto 100
            r = half*ru/(one-au)
         elseif ( .not.l1 .or. .not.l2 .and. ntyp==3 ) then
!
!     three point quadratic extrapolation or interpolation
!
            al = (fi-fl)/(ri-rl)
            au = (fu-fi)/(ru-ri)
            den = au - al
            if ( den<=zero ) goto 100
            r = ri - half*(au*(ri-rl)+al*(ru-ri))/den
         elseif ( l1 .and. .not.l2 .and. ntyp==4 ) then
!
!     three point cubic extrapolation or interpolation
!
            dis = (ai-one)*(ru/ri)
            den = (au-one)*(ri/ru) - dis
            dis = au + ai - den - two*(one+dis)
            dis = den*den - three*dis
            if ( dis<zero ) goto 100
            den = den + sqrt(dis)
            if ( den==zero ) goto 100
            r = (ru-ri)/den
         else
            goto 100
         endif
         if ( mode==1 .and. r>ru ) then
!
!     extrapolation accepted
!
            r = max(r,c1l*ru)
            r = min(r,c1u*ru)
            return
         elseif ( mode==2 .and. r>rl .and. r<ru ) then
!
!     interpolation accepted
!
            if ( ri==zero .and. ntyp/=4 ) then
               r = max(r,rl+c2l*(ru-rl))
            else
               r = max(r,rl+c3l*(ru-rl))
            endif
            r = min(r,rl+c2u*(ru-rl))
            if ( r/=ri ) return
         endif
 100  enddo
      end subroutine pnint3

!***********************************************************************
!> date: 97/12/01
!
! computation of value of the augmented lagrangian function.
!
! parameters :
!  ii  nf  number of variables.
!  ii  n  dimension of the constraint null space.
!  ii  nc  number of constraints.
!  ri  cf(nc+1)  vector containing values of the constraints.
!  ii  ic(nc)  vector containing types of constraints.
!  ii  ica(nf)  vector containing indices of active constraints.
!  ri  cl(nc)  vector containing lower bounds for constraint functions.
!  ri  cu(nc)  vector containing upper bounds for constraint functions.
!  ri  cz(nc)  vector of lagrange multipliers.
!  ri  rpf  penalty coefficient.
!  ro  fc  value of the penalty term.
!  ro  f  value of the penalty function.
!
      subroutine pp0af8(nf,n,nc,cf,ic,ica,cl,cu,cz,rpf,fc,f)
      implicit none
      integer nf , n , nc , ic(*) , ica(*)
      double precision cf(*) , cl(*) , cu(*) , cz(*) , rpf , fc , f
      double precision pom , temp
      integer j , kc
      fc = 0.0d0
      do kc = 1 , nc
         if ( ic(kc)>0 ) then
            pom = 0.0d0
            temp = cf(kc)
            if ( ic(kc)==1 .or. ic(kc)>=3 ) pom = min(pom,temp-cl(kc))
            if ( ic(kc)==2 .or. ic(kc)>=3 ) pom = min(pom,cu(kc)-temp)
            fc = fc + rpf*abs(pom)
         endif
      enddo
      do j = 1 , nf - n
         kc = ica(j)
         if ( kc>0 ) then
            pom = 0.0d0
            temp = cf(kc)
            if ( ic(kc)==1 .or. ic(kc)==3 .or. ic(kc)==5 )              &
                 pom = min(pom,temp-cl(kc))
            if ( ic(kc)==2 .or. ic(kc)==4 .or. ic(kc)==6 )              &
                 pom = max(pom,temp-cu(kc))
            fc = fc - cz(j)*pom
         endif
      enddo
      f = cf(nc+1) + fc
      end subroutine pp0af8

!***********************************************************************
!> date: 97/12/01
!
! computation of the new penalty parameters.
!
! parameters :
!  ii  nf  declared number of variables.
!  ii  n  actual number of variables.
!  ii  nc  number of constraints.
!  ii  ica(nf)  vector containing indices of active constraints.
!  ri  cz(nf)  vector of lagrange multipliers.
!  ri  cp(nc)  vector containing penalty parameters.
!
      subroutine ppset2(nf,n,nc,ica,cz,cp)
      implicit none
      integer nf , n , nc , ica(*)
      double precision cz(*) , cp(*)
      double precision temp
      integer j , l , kc
      do kc = 1 , nc
         cp(kc) = 0.5d0*cp(kc)
      enddo
      do j = 1 , nf - n
         l = ica(j)
         if ( l>0 ) then
            temp = abs(cz(j))
            cp(l) = max(temp,cp(l)+0.5d0*temp)
         endif
      enddo
      end subroutine ppset2

!***********************************************************************
!> date: 97/12/01
!
!  extended line search without directional derivatives.
!
! parameters :
!  ro  r  value of the stepsize parameter.
!  ro  ro  initial value of the stepsize parameter.
!  ro  rp  previous value of the stepsize parameter.
!  ro  f  value of the objective function.
!  ri  fo  initial value of the objective function.
!  ro  fp  previous value of the objective function.
!  ri  po  initial value of the directional derivative.
!  ro  pp  previous value of the directional derivative.
!  ri  fmin  lower bound for value of the objective function.
!  ri  fmax  upper bound for value of the objective function.
!  ri  rmin  minimum value of the stepsize parameter
!  ri  rmax  maximum value of the stepsize parameter
!  ri  tols  termination tolerance for line search (in test on the
!         change of the function value).
!  ii  kd  degree of required dervatives.
!  io  ld  degree of previously computed derivatives.
!  ii  nit  actual number of iterations.
!  ii  kit  number of the iteration after last restart.
!  io  nred  actual number of extrapolations or interpolations.
!  ii  mred  maximum number of extrapolations or interpolations.
!  io  maxst  maximum stepsize indicator. maxst=0 or maxst=1 if maximum
!         stepsize was not or was reached.
!  ii  iest  lower bound specification. iest=0 or iest=1 if lower bound
!         is not or is given.
!  ii  inits  choice of the initial stepsize. inits=0-initial stepsize
!         is specified in the calling program. inits=1-unit initial
!         stepsize. inits=2-combined unit and quadratically estimated
!         initial stepsize. inits=3-quadratically estimated initial
!         stepsize.
!  io  iters  termination indicator. iters=0-zero step. iters=1-perfect
!         line search. iters=2 goldstein stepsize. iters=3-curry
!         stepsize. iters=4-extended curry stepsize.
!         iters=5-armijo stepsize. iters=6-first stepsize.
!         iters=7-maximum stepsize. iters=8-unbounded function.
!         iters=-1-mred reached. iters=-2-positive directional
!         derivative. iters=-3-error in interpolation.
!  ii  kters  termination selection. kters=1-perfect line search.
!         kters=2-goldstein stepsize. kters=3-curry stepsize.
!         kters=4-extended curry stepsize. kters=5-armijo stepsize.
!         kters=6-first stepsize.
!  ii  mes  method selection. mes=1-bisection. mes=2-quadratic
!         interpolation (with one directional derivative).
!         mes=3-quadratic interpolation (with two directional
!         derivatives). mes=4-cubic interpolation. mes=5-conic
!         interpolation.
!  iu  isys  control parameter.
!
!### Method
! safeguarded extrapolation and interpolation with extended termination
! criteria.

      subroutine ps0l02(me,r,ro,rp,f,fo,fp,po,pp,fmin,fmax,rmin,rmax,tols, &
                        kd,ld,nit,kit,nred,mred,maxst,iest,inits,iters, &
                        kters,mes,isys)
      implicit none

      class(psqp_class),intent(inout) :: me

      integer kd , ld , nit , kit , nred , mred , maxst , iest , inits ,&
              iters , kters , mes , isys
      double precision r , ro , rp , f , fo , fp , po , pp , fmin , &
                       fmax , rmin , rmax , tols
      double precision rtemp
      integer merr , init1
      logical l1 , l2 , l3 , l4 , l6 , l7

      double precision,parameter :: tol = 1.0d-4

      if ( isys/=1 ) then
!      go to (1,3) isys+1
         me%mes1 = 2
         me%mes2 = 2
         iters = 0
         if ( po>=0.0d0 ) then
            r = 0.0d0
            iters = -2
            isys = 0
            return
         endif
         if ( rmax<=0.0d0 ) then
            iters = 0
            isys = 0
            return
         endif
!
!     initial stepsize selection
!
         if ( inits>0 ) then
            rtemp = fmin - f
         elseif ( iest==0 ) then
            rtemp = f - fp
         else
            rtemp = max(f-fp,1.0d1*(fmin-f))
         endif
         init1 = abs(inits)
         rp = 0.0d0
         fp = fo
         pp = po
         if ( init1==0 ) then
         elseif ( init1==1 .or. inits>=1 .and. iest==0 ) then
            r = 1.0d0
         elseif ( init1==2 ) then
            r = min(1.0d0,4.0d0*rtemp/po)
         elseif ( init1==3 ) then
            r = min(1.0d0,2.0d0*rtemp/po)
         elseif ( init1==4 ) then
            r = 2.0d0*rtemp/po
         endif
         rtemp = r
         r = max(r,rmin)
         r = min(r,rmax)
         me%mode = 0
         me%rl = 0.0d0
         me%fl = fo
         me%ru = 0.0d0
         me%fu = fo
         me%ri = 0.0d0
         me%fi = fo
      elseif ( iters/=0 ) then
         isys = 0
         return
      else
         if ( f<=fmin ) then
            iters = 7
            isys = 0
            return
         else
            l1 = r<=rmin .and. nit/=kit
            l2 = r>=rmax
            l3 = f - fo<=tols*r*po .or. f - fmin<=(fo-fmin)/1.0d1
            l4 = f - fo>=(1.0d0-tols)*r*po .or. me%mes2==2 .and. me%mode==2
            l6 = me%ru - me%rl<=tol*me%ru .and. me%mode==2
            l7 = me%mes2<=2 .or. me%mode/=0
            maxst = 0
            if ( l2 ) maxst = 1
         endif
!
!     test on termination
!
         if ( l1 .and. .not.l3 ) then
            iters = 0
            isys = 0
            return
         elseif ( l2 .and. .not.f>=me%fu ) then
            iters = 7
            isys = 0
            return
         elseif ( l6 ) then
            iters = 1
            isys = 0
            return
         elseif ( l3 .and. l7 .and. kters==5 ) then
            iters = 5
            isys = 0
            return
         elseif ( l3 .and. l4 .and. l7 .and. &
                  (kters==2 .or. kters==3 .or. kters==4) ) then
            iters = 2
            isys = 0
            return
         elseif ( kters<0 .or. kters==6 .and. l7 ) then
            iters = 6
            isys = 0
            return
         elseif ( abs(nred)>=mred ) then
            iters = -1
            isys = 0
            return
         else
            rp = r
            fp = f
            me%mode = max(me%mode,1)
            me%mtyp = abs(mes)
            if ( f>=fmax ) me%mtyp = 1
         endif
         if ( me%mode==1 ) then
!
!     interval change after extrapolation
!
            me%rl = me%ri
            me%fl = me%fi
            me%ri = me%ru
            me%fi = me%fu
            me%ru = r
            me%fu = f
            if ( f>=me%fi ) then
               nred = 0
               me%mode = 2
            elseif ( me%mes1==1 ) then
               me%mtyp = 1
            endif
!
!     interval change after interpolation
!
         elseif ( r<=me%ri ) then
            if ( f<=me%fi ) then
               me%ru = me%ri
               me%fu = me%fi
               me%ri = r
               me%fi = f
            else
               me%rl = r
               me%fl = f
            endif
         elseif ( f<=me%fi ) then
            me%rl = me%ri
            me%fl = me%fi
            me%ri = r
            me%fi = f
         else
            me%ru = r
            me%fu = f
         endif
      endif
!
!     new stepsize selection (extrapolation or interpolation)
!
      call pnint3(ro,me%rl,me%ru,me%ri,fo,me%fl,me%fu,me%fi,po,r,me%mode,me%mtyp,merr)
      if ( merr>0 ) then
         iters = -merr
         isys = 0
         return
      elseif ( me%mode==1 ) then
         nred = nred - 1
         r = min(r,rmax)
      elseif ( me%mode==2 ) then
         nred = nred + 1
      endif
!
!     computation of the new function value
!
      kd = 0
      ld = -1
      isys = 1
      end subroutine ps0l02

!***********************************************************************
!> date: 92/12/01
!
! variable metric update of a dense symmetric positive definite matrix
! using the factorization b=l*d*trans(l).
!
! parameters :
!  ii  n  actual number of variables.
!  ru  h(m)  factorization b=l*d*trans(l) of a positive
!         definite approximation of the hessian matrix.
!  ri  g(nf)  gradient of the objective function.
!  ra  s(nf)  auxiliary vector.
!  ru  xo(nf)  vectors of variables difference.
!  ri  go(nf)  gradients difference.
!  ri  r  value of the stepsize parameter.
!  ri  po  old value of the directional derivative.
!  ii  nit  actual number of iterations.
!  ii  kit  number of the iteration after last restart.
!  io  iterh  termination indicator. iterh<0-bad decomposition.
!         iterh=0-successful update. iterh>0-nonpositive parameters.
!  ii  met1  selection of self scaling. met1=1-self scaling suppressed.
!         met1=2 initial self scaling.
!  ii  mec  correction if the negative curvature occurs.
!         mec=1-correction suppressed. mec=2-powell's correction.
!
!### Method
! bfgs variable metric method.
!
      subroutine pudbg1(n,h,g,s,xo,go,r,po,nit,kit,iterh,met,met1,mec)
      implicit none
      double precision po , r
      integer iterh , kit , met , met1 , mec , n , nit
      double precision g(*) , go(*) , h(*) , s(*) , xo(*)
      double precision a , b , c , gam , par , den , dis
      logical l1 , l3

      l1 = met1>=3 .or. met1==2 .and. nit==kit
      l3 = .not.l1
!
!     determination of the parameters b, c
!
      b = mxvdot(n,xo,go)
      a = 0.0d0
      if ( l1 ) then
         call mxvcop(n,go,s)
         call mxdpgb(n,h,s,1)
         a = mxdpgp(n,h,s,s)
         if ( a<=0.0d0 ) then
            iterh = 1
            return
         endif
      endif
      call mxvdif(n,go,g,s)
      call mxvscl(n,r,s,s)
      c = -r*po
      if ( c<=0.0d0 ) then
         iterh = 3
         return
      endif
      if ( mec>1 ) then
         if ( b<=1.0d-4*c ) then
!
!     powell's correction
!
            dis = (1.0d0-0.1d0)*c/(c-b)
            call mxvdif(n,go,s,go)
            call mxvdir(n,dis,go,s,go)
            b = c + dis*(b-c)
            if ( l1 ) a = c + 2.0d0*(1.0d0-dis)*(b-c) + dis*dis*(a-c)
         endif
      elseif ( b<=1.0d-4*c ) then
         iterh = 2
         return
      endif
      if ( l1 ) then
!
!     determination of the parameter gam (self scaling)
!
         if ( met==1 ) then
            par = c/b
         elseif ( a<=0.0d0 ) then
            par = c/b
         else
            par = sqrt(c/a)
         endif
         gam = par
         if ( met1>1 ) then
            if ( nit/=kit ) l3 = gam<0.5d0 .or. gam>4.0d0
         endif
      endif
      if ( l3 ) then
         gam = 1.0d0
         par = gam
      endif
      if ( met==1 ) then
!
!     bfgs update
!
         call mxdpgu(n,h,par/b,go,xo)
         call mxdpgu(n,h,-1.0d0/c,s,xo)
      else
!
!     hoshino update
!
         den = par*b + c
         dis = 0.5d0*b
         call mxvdir(n,par,go,s,s)
         call mxdpgu(n,h,par/dis,go,xo)
         call mxdpgu(n,h,-1.0d0/den,s,xo)
      endif
      iterh = 0
      if ( gam==1.0d0 ) return
      call mxdpgs(n,h,1.0d0/gam)
      end subroutine pudbg1

!***********************************************************************
!> date: 91/12/01
!
! dual range space quadratic programming method for minimax
! approximation.
!
! parameters :
!  ii  nf declared number of variables.
!  ii  n  actual number of variables.
!  ii  nc number of constraints.
!  ri  x(nf)  vector of variables.
!  ri  xn(nf)  vector of scaling factors.
!  ro  xo(nf)  saved vector of variables.
!  ii  ica(nf)  vector containing indices of active constraints.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ro  cz(nf)  vector of lagrange multipliers.
!  ro  czs(nf)  saved vector of lagrange multipliers.
!  ri  g(nf)  gradient of the lagrangian function.
!  ri  go(nf)  saved gradient of the lagrangian function.
!  ro  r  value of the stepsize parameter.
!  ro  f  new value of the objective function.
!  ri  fo  old value of the objective function.
!  ro  p  new value of the directional derivative.
!  ri  po  old value of the directional derivative.
!  ri  cmax  value of the constraint violation.
!  ro  cmaxo  saved value of the constraint violation.
!  ro  dmax  maximum relative difference of variables.
!
! common data :
!  ii  normf  scaling specification. normf=0-no scaling performed.
!         normf=1-scaling factors are determined automatically.
!         normf=2-scaling factors are supplied by user.
!  ii  iters  termination indicator for steplength determination.
!         iters=0 for zero step.

      subroutine pytrnd(nf,n,x,xo,ica,cg,cz,g,go,r,f,fo,p,po,cmax,cmaxo,&
                        dmax,kd,ld,iters)
      implicit none
      integer nf , n , kd , ld , iters
      integer ica(*)
      double precision x(*) , xo(*) , cg(*) , cz(*) , g(*) , go(*) , r ,&
                       f , fo , p , po , cmax , cmaxo , dmax
      integer i , j , l
      do j = 1 , nf - n
         l = ica(j)
         if ( l>0 ) then
            call mxvdir(nf,-cz(j),cg((l-1)*nf+1),g,g)
         else
            l = -l
            g(l) = g(l) - cz(j)
         endif
      enddo
      if ( iters>0 ) then
         call mxvdif(nf,x,xo,xo)
         call mxvdif(nf,g,go,go)
         po = r*po
         p = r*p
      else
         f = fo
         p = po
         cmax = cmaxo
         call mxvsav(nf,x,xo)
         call mxvsav(nf,g,go)
         ld = kd
      endif
      dmax = 0.0d0
      do i = 1 , nf
         dmax = max(dmax,abs(xo(i))/max(abs(x(i)),1.0d0))
      enddo
      n = nf
      end subroutine pytrnd

!***********************************************************************
    end module psqp_module
!***********************************************************************
