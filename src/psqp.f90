!***********************************************************************
!>
!  module for psqp

  module psqp_module

  use matrix_routines

  implicit none

  private

  public :: psqpn

  contains
!***********************************************************************

!***********************************************************************
! subroutine psqpn              all systems                   97/01/22
! purpose :
! easy to use subroutine for general nonlinear programming problems.
!
! parameters :
!  ii  nf  number of variables.
!  ii  nb  choice of simple bounds. nb=0-simple bounds suppressed.
!         nb>0-simple bounds accepted.
!  ii  nc  number of linear constraints.
!  ri  x(nf)  vector of variables.
!  ii  ix(nf)  vector containing types of bounds. ix(i)=0-variable
!         x(i) is unbounded. ix(i)=1-lover bound xl(i).le.x(i).
!         ix(i)=2-upper bound x(i).le.xu(i). ix(i)=3-two side bound
!         xl(i).le.x(i).le.xu(i). ix(i)=5-variable x(i) is fixed.
!  ri  xl(nf)  vector containing lower bounds for variables.
!  ri  xu(nf)  vector containing upper bounds for variables.
!  ri  cf(nc+1)  vector containing values of the constraint functions.
!  ii  ic(nc)  vector containing types of constraints.
!         ic(kc)=0-constraint cf(kc) is not used. ic(kc)=1-lover
!         constraint cl(kc).le.cf(kc). ic(kc)=2-upper constraint
!         cf(kc).le.cu(kc). ic(kc)=3-two side constraint
!         cl(kc).le.cf(kc).le.cu(kc). ic(kc)=5-equality constraint
!         cf(kc).eq.cl(kc).
!  ri  cl(nc)  vector containing lower bounds for constraint functions.
!  ri  cu(nc)  vector containing upper bounds for constraint functions.
!  ii  ipar(6)  integer paremeters:
!      ipar(1)  maximum number of iterations.
!      ipar(2)  maximum number of function evaluations.
!      ipar(3)  this parameter is not used in the subroutine psqp.
!      ipar(4)  this parameter is not used in the subroutine psqp.
!      ipar(5)  variable metric update used. ipar(5)=1-the bfgs update.
!         ipar(5)-the hoshino update.
!      ipar(6)  correction of the variable metric update if a negative
!         curvature occurs. ipar(6)=1-no correction. ipar(6)=2-powell's
!         correction.
!  ri  rpar(5)  real parameters:
!      rpar(1)  maximum stepsize.
!      rpar(2)  tolerance for change of variables.
!      rpar(3)  tolerance for constraint violations.
!      rpar(4)  tolerance for the gradient of the lagrangian function.
!      rpar(5)  penalty coefficient.
!  ro  f  value of the objective function.
!  ro  gmax  maximum partial derivative of the lagrangian function.
!  ro  cmax  maximum constraint violation.
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
! variables in common /stat/ (statistics) :
!  io  nres  number of restarts.
!  io  ndec  number of matrix decomposition.
!  io  nrem  number of constraint deletions.
!  io  nadd  number of constraint additions.
!  io  nit  number of iterations.
!  io  nfv  number of function evaluations.
!  io  nfg  number of gradient evaluations.
!  io  nfh  number of hessian evaluations.
!
! subprograms used :
!  s   psqp  recursive quadratic programming method with the bfgs
!         variable metric update.
!
! external subroutines :
!  se  obj  computation of the value of the objective function.
!         calling sequence: call obj(nf,x,ff) where nf is the number
!         of variables, x(nf) is a vector of variables and ff is the
!         value of the objective function.
!  se  dobj  computation of the gradient of the objective function.
!         calling sequence: call dobj(nf,x,gf) where nf is the number
!         of variables, x(nf) is a vector of variables and gc(nf) is
!         the gradient of the objective function.
!  se  con  computation of the value of the constraint function.
!         calling sequence: call con(nf,kc,x,fc) where nf is the
!         number of variables, kc is the index of the constraint
!         function, x(nf) is a vector of variables and fc is the
!         value of the constraint function.
!  se  dcon  computation of the gradient of the constraint function.
!         calling sequence: call dcon(nf,kc,x,gc) where nf is the
!         number of variables, kc is the index of the constraint
!         function, x(nf) is a vector of variables and gc(nf) is the
!         gradient of the constraint function.
!
      subroutine psqpn(nf,nb,nc,x,ix,xl,xu,cf,ic,cl,cu,ipar,rpar,f,gmax,&
                       cmax,iprnt,iterm)
      implicit none
!
!     pointers for auxiliary arrays
!
      double precision f , cmax , gmax
      integer iprnt , iterm , nb , nc , nf
      double precision cf(*) , cl(*) , cu(*) , rpar(5) , x(*) , xl(*) , &
                       xu(*)
      integer ic(*) , ipar(6) , ix(*)
      integer nadd , ndec , nfg , nfh , nfv , nit , nrem , nres
      integer lcfd , lcfo , lcg , lcp , lcr , lcz , lg , lgc , lgf ,    &
              lgo , lh , lia , ls , lxo
      common /stat  / nres , ndec , nrem , nadd , nit , nfv , nfg , nfh
      integer ia(:)
      double precision ra(:)
      allocatable ia , ra
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
      call psqp(nf,nb,nc,x,ix,xl,xu,cf,ic,cl,cu,ra,ra(lcfo),ra(lcfd),   &
                ra(lgc),ia,ra(lcr),ra(lcz),ra(lcp),ra(lgf),ra(lg),ra(lh)&
                ,ra(ls),ra(lxo),ra(lgo),rpar(1),rpar(2),rpar(3),rpar(4),&
                rpar(5),cmax,gmax,f,ipar(1),ipar(2),ipar(5),ipar(6),    &
                iprnt,iterm)
      deallocate (ia,ra)
      end subroutine psqpn

!***********************************************************************
! subroutine psqp               all systems                   97/01/22
! purpose :
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
!         x(i) is unbounded. ix(i)=1-lover bound xl(i).le.x(i).
!         ix(i)=2-upper bound x(i).le.xu(i). ix(i)=3-two side bound
!         xl(i).le.x(i).le.xu(i). ix(i)=5-variable x(i) is fixed.
!  ri  xl(nf)  vector containing lower bounds for variables.
!  ri  xu(nf)  vector containing upper bounds for variables.
!  ro  cf(nc+1)  vector containing values of the constraint functions.
!  ii  ic(nc)  vector containing types of constraints.
!         ic(kc)=0-constraint cf(kc) is not used. ic(kc)=1-lover
!         constraint cl(kc).le.cf(kc). ic(kc)=2-upper constraint
!         cf(kc).le.cu(kc). ic(kc)=3-two side constraint
!         cl(kc).le.cf(kc).le.cu(kc). ic(kc)=5-equality constraint
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
!  ii  mit  maximun number of iterations.
!  ii  mfv  maximun number of function evaluations.
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
! variables in common /stat/ (statistics) :
!  io  nres  number of restarts.
!  io  ndec  number of matrix decomposition.
!  io  nrem  number of constraint deletions.
!  io  nadd  number of constraint additions.
!  io  nit  number of iterations.
!  io  nfv  number of function evaluations.
!  io  nfg  number of gradient evaluations.
!  io  nfh  number of hessian evaluations.
!
! subprograms used :
!  s   pc1f01  computation of the value and the gradient of the
!         constraint function.
!  s   pf1f01  computation of the value and the gradient of the
!         objective function.
!  s   plqdb1  general quadratic programming subroutine based on the
!         goldfarb-idnani dual method.
!  s   plnews  identification of active simple bounds.
!  s   plredl  transformation of the incompatible quadratic programming
!         subproblem.
!  s   pp0af8  computation of value of the augmented lagrangian
!         function.
!  s   ppset2  computation of the new penalty parameters.
!  s   ps0l02  line search using only function values.
!  s   pytrnd  determination of differences for variable metric
!         updates.
!  s   pudbg1  variable metric update after gill-murray decomposition.
!  s   mxdsmi  symmetric matrix is replaced by the unit matrix.
!  s   mxvdir  vector augmented by the scaled vector.
!  rf  mxvdot  dot product of two vectors.
!  s   mxvcop  copying of a vector.
!  s   mxvina  absolute values of elements of an integer vector.
!  rf  mxvmax  l-infinity norm of a vector.
!  s   mxvset  initiation of a vector.
!
! external subroutines :
!  se  obj  computation of the value of the objective function.
!         calling sequence: call obj(nf,x,ff) where nf is the number
!         of variables, x(nf) is a vector of variables and ff is the
!         value of the objective function.
!  se  dobj  computation of the gradient of the objective function.
!         calling sequence: call dobj(nf,x,gf) where nf is the number
!         of variables, x(nf) is a vector of variables and gc(nf) is
!         the gradient of the objective function.
!  se  con  computation of the value of the constraint function.
!         calling sequence: call con(nf,kc,x,fc) where nf is the
!         number of variables, kc is the index of the constraint
!         function, x(nf) is a vector of variables and fc is the
!         value of the constraint function.
!  se  dcon  computation of the gradient of the constraint function.
!         calling sequence: call dcon(nf,kc,x,gc) where nf is the
!         number of variables, kc is the index of the constraint
!         function, x(nf) is a vector of variables and gc(nf) is the
!         gradient of the constraint function.
!
! method :
! recursive quadratic programming method with the bfgs variable metric
! update.
!
      subroutine psqp(nf,nb,nc,x,ix,xl,xu,cf,ic,cl,cu,cg,cfo,cfd,gc,ica,&
                      cr,cz,cp,gf,g,h,s,xo,go,xmax,tolx,tolc,tolg,rpf,  &
                      cmax,gmax,f,mit,mfv,met,mec,iprnt,iterm)
      implicit none
      double precision f , cmax , gmax , rpf , tolc , told , tolg ,     &
                       tols , tolx , xmax
      integer iprnt , iterm , met , met1 , mec , mes , mfv , mit , nb , &
              nc , nf
      double precision cf(*) , cfd(*) , cfo(*) , cg(*) , cl(*) , cp(*) ,&
                       cr(*) , cz(*) , cu(*) , g(*) , gc(*) , gf(*) ,   &
                       go(*) , h(*) , s(*) , x(*) , xl(*) , xo(*) ,     &
                       xu(*)
      integer ic(*) , ica(*) , ix(*)
      integer nadd , ndec , nfg , nfh , nfv , nit , nrem , nres
      double precision alf1 , alf2 , cmaxo , dmax , eps7 , eps9 , eta0 ,&
                       eta2 , eta9 , fmax , fmin , fo , gnorm , p , po ,&
                       r , rmax , rmin , ro , snorm , tolb , tolf ,     &
                       umax , rp , fp , pp , ff , fc
      integer i , idecf , iext , irest , iterd , iterl , iterh , iterq ,&
              iters , kbc , kbf , kc , kd , kit , ld , mred , mtesf ,   &
              mtesx , n , k , ntesx , iest , inits , kters , maxst ,    &
              isys , mfp , nred , ipom , lds
      !double precision mxvdot , mxvmax
      common /stat  / nres , ndec , nrem , nadd , nit , nfv , nfg , nfh
      if ( abs(iprnt)>1 ) write (6,'(1x,''entry to psqp :'')')
!
!     initiation
!
      kbf = 0
      kbc = 0
      if ( nb>0 ) kbf = 2
      if ( nc>0 ) kbc = 2
      nres = 0
      ndec = 0
      nrem = 0
      nadd = 0
      nit = 0
      nfv = 0
      nfg = 0
      nfh = 0
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
      call pf1f01(nf,x,gf,gf,ff,f,kd,ld,iext)
      ld = lds
      call pc1f01(nf,nc,x,fc,cf,cl,cu,ic,gc,cg,cmax,kd,ld)
      cf(nc+1) = f
      if ( abs(iprnt)>1 ) write (6,                                     &
      '(1x,''nit='',i4,2x,''nfv='',i4,2x,''nfg='',i4,2x,       ''f='',g1&
      2.6,2x,''c='',e7.1,2x,''g='',e7.1)') nit , nfv , nfg , f , cmax , &
      gmax
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
      if ( nit>=mit ) then
         iterm = 11
         goto 500
      endif
      if ( nfv>=mfv ) then
         iterm = 12
         goto 500
      endif
      iterm = 0
      nit = nit + 1
!
!     restart
!
 200  n = nf
      if ( irest>0 ) then
         call mxdsmi(n,h)
         ld = min(ld,1)
         idecf = 1
         if ( kit<nit ) then
            nres = nres + 1
            kit = nit
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
 300  call plqdb1(nf,nc,x,ix,xl,xu,cf,cfd,ic,ica,cl,cu,cg,cr,cz,g,gf,h, &
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
 450     call ps0l02(r,ro,rp,f,fo,fp,po,pp,fmin,fmax,rmin,rmax,tols,kd, &
                     ld,nit,kit,nred,mred,maxst,iest,inits,iters,kters, &
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
               call pf1f01(nf,x,gf,gf,ff,f,kd,ld,iext)
               ld = lds
               call pc1f01(nf,nc,x,fc,cf,cl,cu,ic,gc,cg,cmax,kd,ld)
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
            call pudbg1(n,h,g,s,xo,go,r,po,nit,kit,iterh,met,met1,mec)
!      if (mer.gt.0.and.iterh.gt.0) irest=1
!
!     end of the iteration
!
            goto 100
         else
!      go to (11174,11172) isys+1
            call mxvdir(nf,r,s,xo,x)
            lds = ld
            call pf1f01(nf,x,gf,g,ff,f,kd,ld,iext)
            ld = lds
            call pc1f01(nf,nc,x,fc,cf,cl,cu,ic,gc,cg,cmax,kd,ld)
            cf(nc+1) = f
            call pp0af8(nf,n,nc,cf,ic,ica,cl,cu,cz,rpf,fc,f)
            goto 450
         endif
      endif

 500  if ( iprnt>1 .or. iprnt<0 ) write (6,'(1x,''exit from psqp :'')')
      if ( iprnt/=0 ) &
         write (6,'(1x,''nit='',i4,2x,''nfv='',i4,2x,''nfg='',i4,2x,''f='',&
                  g12.6,2x,''c='',e7.1,2x,''g='',e7.1,2x,''iterm='',i3)') &
                  nit , nfv , nfg , f , cmax , gmax , iterm
      if ( iprnt<0 ) write (6,'(1x,''x='',5(g14.7,1x):/(3x,5(g14.7,1x)))') (x(i),i=1,nf)
      end subroutine psqp

! subroutine pa0gs1             all systems                 97/12/01
! 97/12/01 lu : original version
!
! purpose:
! numerical computation of the gradient of the model function.
!
! parameters :
!  ii  n  number of variables.
!  ii  ka  indef of the approximated function.
!  ri  x(n)  vector of variables.
!  ro  ga(n)  gradient of the approximated function.
!  ri  fa  value of the approximated function.
!  ri  eta1  precision of the computed values.
!  iu  nav  number of approximated function evaluations.
!
      subroutine pa0gs1(n,ka,x,ga,fa,eta1,nav)
      implicit none
      integer n , ka , nav
      double precision x(*) , ga(*) , fa , eta1
      double precision xstep , xtemp , ftemp , eta
      integer ivar
      eta = sqrt(eta1)
      ftemp = fa
      do ivar = 1 , n
!
!     step selection
!
         xstep = 1.0d0
         xstep = eta*max(abs(x(ivar)),xstep)*sign(1.0d0,x(ivar))
         xtemp = x(ivar)
         x(ivar) = x(ivar) + xstep
         xstep = x(ivar) - xtemp
         nav = nav + 1
         call fun(n,ka,x,fa)
!
!     numerical differentiation
!
         ga(ivar) = (fa-ftemp)/xstep
         x(ivar) = xtemp
      enddo
      fa = ftemp
      end subroutine pa0gs1

! subroutine pa1sq1             all systems                 97/12/01
! 97/12/01 lu : original version
!
! purpose :
! computation of the value and the gradient of the objective function
! which is defined as a sum of squares.
!
! parameters:
!  ii  n  number of variables.
!  ri  x(n)  vector of variables.
!  ro  f  value of the objective function.
!  ro  af(n)  values of the approximated functions.
!  ri  ga(nf)  gradient of the approximated function.
!  ri  ag(n*n)  rectangular matrix which is used for the direction
!         vector determination.
!  ro  g(nf)  gradient of the objective function.
!  ri  eta1  precision of the computes function values.
!  ii  kda  degree of computed derivatives.
!  ii  kd  degree of required dervatives.
!  io  ld  degree of previously computed derivatives.
!  iu  nfv  number of objective function values computed.
!  iu  nfg  number of objective function gradients computed.
!
! subprograms used :
!  s   mxvcop  copying of a vector.
!  s   mxvdir  vector augmented by the scaled vector.
!  s   mxvset  initiation of a vector.
!
      subroutine pa1sq1(n,x,f,af,ga,ag,g,eta1,kda,kd,ld,nfv,nfg)
      implicit none
      integer n , kda , kd , ld , nfv , nfg
      double precision x(*) , f , af(*) , ga(*) , ag(*) , g(*) , eta1
      integer ka , nav
      double precision fa
      if ( kd<=ld ) return
      if ( kd>=0 .and. ld<0 ) then
         f = 0.0d0
         nfv = nfv + 1
      endif
      if ( kd>=1 .and. ld<1 ) then
         call mxvset(n,0.0d0,g)
         if ( kda>0 ) nfg = nfg + 1
      endif
      nav = 0
      do ka = 1 , n
         if ( kd>=0 ) then
            if ( ld>=0 ) then
               fa = af(ka)
               goto 20
            else
               call fun(n,ka,x,fa)
               af(ka) = fa
            endif
            if ( ld<0 ) f = f + fa*fa
 20         if ( kd>=1 ) then
               if ( kda>0 ) then
                  call dfun(n,ka,x,ga)
               else
                  call pa0gs1(n,ka,x,ga,fa,eta1,nav)
               endif
               call mxvdir(n,fa,ga,g,g)
               call mxvcop(n,ga,ag((ka-1)*n+1))
            endif
         endif
      enddo
      nfv = nfv + nav/n
      if ( kd>=0 .and. ld<0 ) f = 0.5d0*f
      ld = kd
      end subroutine pa1sq1

! subroutine pa2sq1             all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
!  computation of the value and the gradient and the hessian matrix
!  of the objective function which is defined as a sum of squares.
!
! parameters:
!  ii  nf  number of variables.
!  ii  na  number of approximated functions.
!  ri  ga(nf)  gradient of the approximated function.
!  ro  g(nf)  gradient of the objective function.
!  ri  ha(nf*(nf+1)/2)  hessian matrix of the approximated function.
!  ro  h(nf*(nf+1)/2)  hessian matrix of the objective function.
!  ri  fa  value of the selected function.
!  ro  f  value of the objective function.
!
! common data :
!  ii  kda  degree of analytically computed derivatives.
!  ii  kd  degree of required dervatives.
!  io  ld  degree of previously computed derivatives.
!  iu  nav  number of approximated function values computed.
!  iu  nag  number of approximated function gradients computed.
!  iu  nah  number of approximated function hessian matrices computed.
!  iu  nfv  number of objective function values computed.
!  iu  nfg  number of objective function gradients computed.
!  iu  nfh  number of objective function hessian matrices computed.
!  iu  idecf  decomposition indicator. idecf=0-no decomposition.
!
! status variables :
!  ns,isp,tss
!
! subprograms used :
!  s  uypro1  prologue.
!  s  uyepi1  epilogue.
!  s  uyset0  status definition.
!  s  mxvset  initiation of a vector.
!  s  mxvdir  vector augmented by the scaled vector.
!  s  mxdsmo  initiation of a dense symmetric matrix.
!  s  mxdsma  dense symmetric matrix augmented by the scaled dense
!         symmetric matrix.
!  s  mxdsmu  correction of a dense symmetric matrix.
!
      subroutine pa2sq1(nf,na,x,f,af,ga,g,h,eta1,kda,kd,ld,nfv,nfg)
      implicit none
      integer nf , na , kda , kd , ld , nfv , nfg
      double precision x(*) , f , af(*) , ga(*) , g(*) , h(*) , eta1
      integer ka , nav
      double precision fa
      if ( kd<=ld ) return
      if ( kd>=0 .and. ld<0 ) then
         f = 0.0d0
         nfv = nfv + 1
      endif
      if ( kd>=1 .and. ld<1 ) then
         call mxvset(nf,0.0d0,g)
         if ( kda>0 ) nfg = nfg + 1
      endif
      if ( kd>=2 .and. ld<2 ) call mxvset(nf*(nf+1)/2,0.0d0,h)
      nav = 0
      do ka = 1 , na
         if ( kd>=0 ) then
            if ( ld>=0 ) then
               fa = af(ka)
               goto 20
            else
               call fun(nf,ka,x,fa)
               af(ka) = fa
            endif
            if ( ld<0 ) f = f + fa*fa
 20         if ( kd>=1 ) then
               if ( kda>0 ) then
                  call dfun(nf,ka,x,ga)
               else
                  call pa0gs1(nf,ka,x,ga,fa,eta1,nav)
               endif
               call mxvdir(nf,fa,ga,g,g)
               if ( kd>=2 ) call mxdsmu(nf,h,1.0d0,ga)
            endif
         endif
      enddo
      nfv = nfv + nav/na
      if ( kd>=0 .and. ld<0 ) f = 0.5d0*f
      ld = kd
      end subroutine pa2sq1

! subroutine pc1f01             all systems                 97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  s   mxvcop  copying of a vector.
!
      subroutine pc1f01(nf,nc,x,fc,cf,cl,cu,ic,gc,cg,cmax,kd,ld)
      implicit none
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
               call con(nf,kc,x,fc)
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
                  call dcon(nf,kc,x,gc)
                  call mxvcop(nf,gc,cg((kc-1)*nf+1))
               endif
            endif
         endif
      enddo
      ld = kd
      end subroutine pc1f01

! subroutine pf1f01                all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  s   mxvcop  copying of a vector.
!  s   mxvneg  copying of a vector with change of the sign.
!
      subroutine pf1f01(nf,x,gf,g,ff,f,kd,ld,iext)
      implicit none
      double precision f , ff
      integer iext , kd , ld , nf
      double precision gf(*) , g(*) , x(*)
      integer nadd , ndec , nfg , nfh , nfv , nit , nrem , nres
      common /stat  / nres , ndec , nrem , nadd , nit , nfv , nfg , nfh
      if ( kd<=ld ) return
      if ( ld<0 ) then
         nfv = nfv + 1
         call obj(nf,x,ff)
         if ( iext<=0 ) then
            f = ff
         else
            f = -ff
         endif
      endif
      if ( kd>=1 ) then
         if ( ld<1 ) then
            nfg = nfg + 1
            call dobj(nf,x,gf)
            if ( iext>0 ) call mxvneg(nf,gf,g)
         endif
      endif
      ld = kd
      end subroutine pf1f01

! subroutine pladb0               all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  s   pladr0  correction of kernel of the orthogonal projection
!         after constraint addition.
!  s   mxdrmm  premultiplication of a vector by a rowwise stored dense
!         rectangular matrix.
!  s   mxdrmv  copy of the selected column of a rowwise stored dense
!         rectangular matrix.
!  s   mxdrgr  plane rotation of a transposed dense rectangular matrix.
!  s   mxvort  determination of an elementary orthogonal matrix for
!         plane rotation.
!
      subroutine pladb0(nf,n,ica,cg,cr,cz,s,eps7,gmax,umax,inew,nadd,   &
                        ier)
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

! subroutine pladb4               all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  s   pladr0  correction of kernel of the orthogonal projection
!         after constraint addition.
!  s   mxdrmm  premultiplication of a vector by a rowwise stored dense
!         rectangular matrix.
!  s   mxdrmv  copy of the selected column of a rowwise stored dense
!         rectangular matrix.
!  s   mxdrgr  plane rotation of a transposed dense rectangular matrix.
!         rectangular matrix.
!  s   mxdsmr  plane rotation of a dense symmetric matrix.
!  s   mxvort  determination of an elementary orthogonal matrix for
!         plane rotation.
!
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

! subroutine pladr0               all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  s   mxsprb  sparse back substitution.
!  s   mxvcop  copying of a vector.
!  rf  mxvdot  dot product of two vectors.
!
      subroutine pladr0(nf,n,ica,cg,cr,s,eps7,gmax,umax,inew,nadd,ier)
      implicit none
      integer nf , n , ica(*) , inew , nadd , ier
      double precision cg(*) , cr(*) , s(*) , eps7 , gmax , umax
      !double precision mxvdot
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

! subroutine pllpb1             all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
! determination of the initial feasible point and the linear programming
! subroutine.
!
! parameters :
!  ii  nf  number of variables.
!  ii  nc  number of linear constraints.
!  ru  x(nf)  vector of variables.
!  ii  ix(nf)  vector containing types of bounds.
!  ro  xo(nf)  saved vector of variables.
!  ri  xl(nf)  vector containing lower bounds for variables.
!  ri  xu(nf)  vector containing upper bounds for variables.
!  ru  cf(nf)  vector containing values of the constraint functions.
!  ra  cfd(nf)  vector containing increments of the constraint
!         functions.
!  ii  ic(nc)  vector containing types of constraints.
!  io  ica(nf)  vector containing indices of active constraints.
!  ri  cl(nc)  vector containing lower bounds for constraint functions.
!  ri  cu(nc)  vector containing upper bounds for constraint functions.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ro  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ro  cz(nf*nf)  matrix whose columns are basic vectors from the
!         current reduced subspace.
!  ri  g(nf)  gradient of the objective function.
!  ro  go(nf)  saved gradient of the objective function.
!  ra  s(nf)  direction vector.
!  ii  mfp  type of feasible point. mfp=1-arbitrary feasible point.
!         mfp=2-optimum feasible point. mfp=3-repeated solution.
!  ii  kbf  specification of simple bounds. kbf=0-no simple bounds.
!         kbf=1-one sided simple bounds. kbf=2=two sided simple bounds.
!  ii  kbc  specification of linear constraints. kbc=0-no linear
!         constraints. kbc=1-one sided linear constraints. kbc=2=two
!         sided linear constraints.
!  ri  eta9  maximum for real numbers.
!  ri  eps7  tolerance for linear independence of constraints.
!  ri  eps9  tolerance for activity of constraints.
!  ro  umax  maximum absolute value of a negative lagrange multiplier.
!  ro  gmax  maximum absolute value of a partial derivative.
!  io  n  dimension of the manifold defined by active constraints.
!  io  iterl  type of feasible point. iterl=1-arbitrary feasible point.
!         iterl=2-optimum feasible point. iterl=-1 feasible point does
!         not exists. iterl=-2 optimum feasible point does not exists.
!
! subprograms used :
!  s   plinit  determination of initial point satisfying simple bounds.
!  s   plmaxl  maximum stepsize using linear constraints.
!  s   plmaxs  maximum stepsize using simple bounds.
!  s   plmaxt  maximum stepsize using trust region bounds.
!  s   plnewl  identification of active linear constraints.
!  s   plnews  identification of active simple bounds.
!  s   plnewt  identification of active trust region bounds.
!  s   pldirl  new values of constraint functions.
!  s   pldirs  new values of variables.
!  s   plsetc  initial values of constraint functions.
!  s   plsetg  determination of the first phase gradient vector.
!  s   pltrbg  determination of lagrange multipliers and computation
!  s   pladb0  constraint addition.
!  s   plrmb0  constraint deletion.
!  s   mxdcmm  premultiplication of a vector by a dense rectangular
!         matrix stored by columns.
!  s   mxdrmm  premultiplication of a vector by a dense rectangular
!         matrix stored by rows.
!  s   mxdsmi  determination of the initial unit dense symmetric
!         matrix.
!  s   mxvcop  copying of a vector.
!  s   mxvina  absolute values of elements of an integer vector.
!  s   mxvinc  update of an integer vector.
!  s   mxvind  change of the integer vector for constraint addition.
!  s   mxvint  change of the integer vector for trust region bound
!         addition.
!  s   mxvmul  diagonal premultiplication of a vector.
!  s   mxvneg  copying of a vector with change of the sign.
!  s   mxvset  initiation of a vector.
!
      subroutine pllpb1(nf,nc,x,ix,xo,xl,xu,cf,cfd,ic,ica,cl,cu,cg,cr,  &
                        cz,g,go,s,mfp,kbf,kbc,eta9,eps7,eps9,umax,gmax, &
                        n,iterl)
      implicit none
      integer nf , nc , ix(*) , ic(*) , ica(*) , mfp , kbf , kbc , n ,  &
              iterl
      double precision x(*) , xo(*) , xl(*) , xu(*) , cf(*) , cfd(*) ,  &
                       cl(*) , cu(*) , cg(*) , cr(*) , cz(*) , g(*) ,   &
                       go(*) , s(*) , eta9 , eps7 , eps9 , umax , gmax
      double precision pom , con , dmax
      integer nca , ncr , ncz , ipom , i , k , iold , inew , ier ,      &
              krem , kc , nred
      integer nres , ndec , nrem , nadd , nit , nfv , nfg , nfh
      common /stat  / nres , ndec , nrem , nadd , nit , nfv , nfg , nfh
      con = eta9
!
!     initiation
!
      call mxvcop(nf,x,xo)
      ipom = 0
      nred = 0
      krem = 0
      iterl = 1
      dmax = 0.0d0
      if ( mfp==3 ) goto 200
      if ( kbf>0 ) call mxvina(nf,ix)
!
!     shift of variables for satisfying simple bounds
!
      call plinit(nf,x,ix,xl,xu,eps9,kbf,inew,iterl)
      if ( iterl<0 ) return
      n = 0
      nca = 0
      ncz = 0
      do i = 1 , nf
         if ( kbf>0 .and. ix(i)<0 ) then
            nca = nca + 1
            ica(nca) = -i
         else
            n = n + 1
            call mxvset(nf,0.0d0,cz(ncz+1))
            cz(ncz+i) = 1.0d0
            ncz = ncz + nf
         endif
      enddo
      call mxdsmi(nca,cr)
      if ( nc>0 ) then
         call mxdrmm(nf,nc,cg,x,cf)
!
!     addition of active constraints and initial check of feasibility
!
         call mxvina(nc,ic)
!      if (nf.gt.n) call plsetc(nf,nc,x,xo,cf,ic,cg,s)
         do kc = 1 , nc
            if ( ic(kc)/=0 ) then
               inew = 0
               call plnewl(kc,cf,ic,cl,cu,eps9,inew)
               call pladb0(nf,n,ica,cg,cr,cz,s,eps7,gmax,umax,inew,nadd,&
                           ier)
               call mxvind(ic,kc,ier)
               if ( ic(kc)<-10 ) ipom = 1
            endif
         enddo
      endif
 100  if ( ipom==1 ) then
!
!     check of feasibility and update of the first phase objective
!     function
!
         call plsetg(nf,nc,ic,cg,go,inew)
         if ( inew==0 ) ipom = 0
      endif
      if ( ipom==0 .and. iterl==0 ) then
!
!     feasibility achieved
!
         iterl = 1
         call mxvcop(nf,g,go)
         if ( mfp==1 ) return
      endif
!
!     lagrange multipliers and reduced gradient determination
!
 200  call pltrbg(nf,n,nc,ix,ic,ica,cg,cr,cz,go,s,eps7,gmax,umax,iold)
      inew = 0
      if ( gmax/=0.0d0 ) then
!
!     direction determination
!
         nca = nf - n
         ncr = nca*(nca+1)/2
         call mxdcmm(nf,n,cz,s,cr(ncr+1))
         call mxvneg(nf,cr(ncr+1),s)
!
!     stepsize selection
!
         pom = con
         call plmaxl(nf,nc,cf,cfd,ic,cl,cu,cg,s,pom,kbc,krem,inew)
         call plmaxs(nf,x,ix,xl,xu,s,pom,kbf,krem,inew)
         if ( inew/=0 ) then
!
!     step realization
!
            call pldirs(nf,x,ix,s,pom,kbf)
            call pldirl(nc,cf,cfd,ic,pom,kbc)
!
!     constraint addition
!
            if ( inew>0 ) then
               kc = inew
               inew = 0
               call plnewl(kc,cf,ic,cl,cu,eps9,inew)
               call pladb0(nf,n,ica,cg,cr,cz,s,eps7,gmax,umax,inew,nadd,&
                           ier)
               call mxvind(ic,kc,ier)
            elseif ( inew+nf>=0 ) then
               i = -inew
               inew = 0
               call plnews(x,ix,xl,xu,eps9,i,inew)
               call pladb0(nf,n,ica,cg,cr,cz,s,eps7,gmax,umax,inew,nadd,&
                           ier)
               call mxvind(ix,i,ier)
            endif
            dmax = pom
            if ( dmax>0.0d0 ) nred = nred + 1
            goto 100
         elseif ( ipom==0 ) then
!
!     bounded solution does not exist
!
            iterl = -2
         else
!
!     feasible solution does not exist
!
            iterl = -3
         endif
!
!     optimum on a linear manifold obtained
!
      elseif ( iold/=0 ) then
!
!     constraint deletion
!
         call plrmb0(nf,n,ica,cg,cr,cz,go,s,iold,krem,nrem,ier)
         kc = ica(nf-n+1)
         if ( kc>0 ) then
            ic(kc) = -ic(kc)
         else
            k = -kc
            ix(k) = -ix(k)
         endif
         dmax = 0.0d0
         goto 200
      elseif ( ipom==0 ) then
!
!     optimal solution achieved
!
         iterl = 2
      else
         ipom = 0
         do kc = 1 , nc
            if ( ic(kc)<-10 ) then
               inew = 0
               call plnewl(kc,cf,ic,cl,cu,eps9,inew)
               if ( ic(kc)<-10 ) ipom = 1
            endif
         enddo
         if ( ipom==0 ) then
!
!     optimal solution achieved
!
            call mxvcop(nf,go,g)
            iterl = 2
         else
!
!     feasible solution does not exist
!
            call mxvcop(nf,go,g)
            iterl = -1
         endif
      endif
      end subroutine pllpb1

! subroutine plrmb0               all systems                92/12/01
! 92/12/01 lu : original version
!
! purpose :
! old linear constraint or an old simple bound is removed from the
! active set.
!
! parameters :
!  ii  nf  declared number of variables.
!  ii  n  actual number of variables.
!  ii  ica(nf)  vector containing indices of active constraints.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ru  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ru  cz(nf*nf)  matrix whose columns are basic vectors from the
!         current reduced subspace.
!  ri  g(nf)  gradient of the objective function.
!  ru  gn(nf)  transformed gradient of the objective function.
!  ii  iold  index of the old active constraint.
!  io  krem  auxiliary variable.
!  iu  nrem number of constraint deletion.
!  io  ier  error indicator.
!
! subprograms used :
!  s   plrmr0  correction of kernel of the orthogonal projection
!         after constraint deletion.
!  s   mxdprb  back substitution.
!  s   mxvcop  copying of a vector.
!  s   mxvdir  vector augmented by the scaled vector.
!  rf  mxvdot  dot product of two vectors.
!  s   mxvmul  diagonal premultiplication of a vector.
!  s   mxvset  initiation of a vector.
!
      subroutine plrmb0(nf,n,ica,cg,cr,cz,g,gn,iold,krem,nrem,ier)
      implicit none
      integer nf , n , ica(*) , iold , krem , nrem , ier
      double precision cg(*) , cr(*) , cz(*) , g(*) , gn(*)
      !double precision mxvdot
      integer nca , ncr , ncz , i , j , kc
      ier = 0
      if ( n==nf ) ier = 2
      if ( iold==0 ) ier = 3
      if ( ier/=0 ) return
      nca = nf - n
      ncr = nca*(nca-1)/2
      ncz = n*nf
      call plrmr0(nf,ica,cr,cz(ncz+1),n,iold,krem,ier)
      call mxvset(nca,0.0d0,cz(ncz+1))
      cz(ncz+nca) = 1.0d0
      call mxdprb(nca,cr,cz(ncz+1),-1)
      call mxvcop(nca,cz(ncz+1),cr(ncr+1))
      n = n + 1
      call mxvset(nf,0.0d0,cz(ncz+1))
      do j = 1 , nca
         kc = ica(j)
         if ( kc>0 ) then
            call mxvdir(nf,cr(ncr+j),cg((kc-1)*nf+1),cz(ncz+1),cz(ncz+1)&
                        )
         else
            i = -kc
            cz(ncz+i) = cz(ncz+i) + cr(ncr+j)
         endif
      enddo
      gn(n) = mxvdot(nf,cz(ncz+1),g)
      nrem = nrem + 1
      ier = 0
      end subroutine plrmb0

! subroutine plqdb1             all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  s   plmins  determination of the new active simple bound.
!  s   plminl  determination of the new active linear constraint.
!  s   plmint  determination of the new active trust region bound.
!  s   pladr1  addition of a new active constraint.
!  s   plrmr0  constrain deletion.
!  s   plsob1  transformation of the local solution to the solution
!         of the original qp problem.
!  s   mxdpgf  gill-murray decomposition of a dense symmetric matrix.
!  s   mxdpgb  back substitution after gill-murray decomposition.
!  s   mxdprb  back substitution.
!  s   mxdsmm  matrix vector product.
!  s   mxvcop  copying of a vector.
!  s   mxvdir  vector augmented by the scaled vector.
!  s   mxvina  absolute values of elements of an integer vector.
!  s   mxvinv  change of an integer vector after constraint addition.
!  s   mxvneg  copying of a vector with change of the sign.
!
      subroutine plqdb1(nf,nc,x,ix,xl,xu,cf,cfd,ic,ica,cl,cu,cg,cr,cz,g,&
                        go,h,s,mfp,kbf,kbc,idecf,eta2,eta9,eps7,eps9,   &
                        umax,gmax,n,iterq)
      implicit none
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
      integer nres , ndec , nrem , nadd , nit , nfv , nfg , nfh
      common /stat  / nres , ndec , nrem , nadd , nit , nfv , nfg , nfh
      con = eta9
      if ( idecf<0 ) idecf = 1
      if ( idecf==0 ) then
!
!     gill-murray decomposition
!
         temp = eta2
         call mxdpgf(nf,h,inf,temp,step)
         ndec = ndec + 1
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
         call pladr1(nf,n,ica,cg,cr,h,s,g,eps7,gmax,umax,idecf,inew,    &
                     nadd,ier,1)
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
            nadd = nadd + 1
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
            call plrmf0(nf,nc,ix,ic,ica,cr,ic,g,n,iold,krem,ier)
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

! subroutine pladr1               all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
!         computed in case when n.le.0 if job=0 or job=1 respectively.
!
! subprograms used :
!  s   mxdpgb  back substitution.
!  s   mxdprb  back substitution.
!  s   mxdsmm  matrix-vector product.
!  s   mxdsmv  copying of a row of dense symmetric matrix.
!  s   mxvcop  copying of a vector.
!  rf  mxvdot  dot product of two vectors.
!
      subroutine pladr1(nf,n,ica,cg,cr,h,s,g,eps7,gmax,umax,idecf,inew, &
                        nadd,ier,job)
      implicit none
      integer nf , n , ica(*) , idecf , inew , nadd , ier , job
      double precision cg(*) , cr(*) , h(*) , s(*) , g(*) , eps7 ,      &
                       gmax , umax
      !double precision mxvdot
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

! subroutine pldirl               all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
! determination of the new values of the constraint functions.
!
! parameters :
!  ii  nc  number of constraints.
!  ru  cf(nf)  vector containing values of the constraint functions.
!  ri  cfd(nf)  vector containing increments of the constraint
!         functions.
!  ii  ic(nc)  vector containing types of constraints.
!  ri  step  current stepsize.
!  ii  kbc  specification of linear constraints. kbc=0-no linear
!         constraints. kbc=1-one sided linear constraints. kbc=2=two
!         sided linear constraints.
!
      subroutine pldirl(nc,cf,cfd,ic,step,kbc)
      implicit none
      integer nc , ic(*) , kbc
      double precision cf(*) , cfd(*) , step
      integer kc
      if ( kbc>0 ) then
         do kc = 1 , nc
            if ( ic(kc)>=0 .and. ic(kc)<=10 ) then
               cf(kc) = cf(kc) + step*cfd(kc)
            elseif ( ic(kc)<-10 ) then
               cf(kc) = cf(kc) + step*cfd(kc)
            endif
         enddo
      endif
      end subroutine pldirl

! subroutine pldirs               all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
! determination of the new vector of variables.
!
! parameters :
!  ii  nf  number of variables.
!  ru  x(nf)  vector of variables.
!  ii  ix(nf)  vector containing types of bounds.
!  ri  s(nf)  direction vector.
!  ri  step  current stepsize.
!  ii  kbf  specification of simple bounds. kbf=0-no simple bounds.
!         kbf=1-one sided simple bounds. kbf=2=two sided simple bounds.
!
      subroutine pldirs(nf,x,ix,s,step,kbf)
      implicit none
      integer nf , ix(*) , kbf
      double precision x(*) , s(*) , step
      integer i
      do i = 1 , nf
         if ( kbf<=0 ) then
            x(i) = x(i) + step*s(i)
         elseif ( ix(i)>=0 .and. ix(i)<=10 ) then
            x(i) = x(i) + step*s(i)
         elseif ( ix(i)<-10 ) then
            x(i) = x(i) + step*s(i)
         endif
      enddo
      end subroutine pldirs

! subroutine plinit             all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
! determination of the initial point which satisfies simple bounds.
!
! parameters :
!  ii  nf  number of variables.
!  ru  x(nf)  vector of variables.
!  ii  ix(nf)  vector containing types of bounds.
!  ri  xl(nf)  vector containing lower bounds for variables.
!  ri  xu(nf)  vector containing upper bounds for variables.
!  ri  eps9  tolerance for active constraints.
!  io  inew  index of the new active constraint.
!  io  ind  indicator. if ind.ne.0 then trust region bounds cannot
!         be satisfied.
!
! subprograms used :
!  s   plnews  test on activity of a given simple bound.
!
      subroutine plinit(nf,x,ix,xl,xu,eps9,kbf,inew,ind)
      implicit none
      integer nf , ix(*) , kbf , inew , ind
      double precision x(*) , xl(*) , xu(*) , eps9
      integer i
      ind = 0
      if ( kbf>0 ) then
         do i = 1 , nf
            call plnews(x,ix,xl,xu,eps9,i,inew)
            if ( ix(i)<5 ) then
            elseif ( ix(i)==5 ) then
               ix(i) = -5
            elseif ( ix(i)==11 .or. ix(i)==13 ) then
               x(i) = xl(i)
               ix(i) = 10 - ix(i)
            elseif ( ix(i)==12 .or. ix(i)==14 ) then
               x(i) = xu(i)
               ix(i) = 10 - ix(i)
            endif
         enddo
      endif
      end subroutine plinit

! subroutine plmaxl               all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
! determination of the maximum stepsize using linear constraints.
!
! parameters :
!  ii  nf  declared number of variables.
!  ii  nc  number of current linear constraints.
!  ri  cf(nf)  vector containing values of the constraint funcyions.
!  ro  cfd(nf)  vector containing increments of the constraint
!         functions.
!  ii  ic(nc)  vector containing types of constraints.
!  ri  cl(nc)  vector containing lower bounds for constraint functions.
!  ri  cu(nc)  vector containing upper bounds for constraint functions.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ri  s(nf)  direction vector.
!  ro  step  maximum stepsize.
!  ii  kbc  specification of linear constraints. kbc=0-no linear
!         constraints. kbc=1-one sided linear constraints. kbc=2=two
!         sided linear constraints.
!  ii  krem  indication of linearly dependent gradients.
!  io  inew  index of the new active function.
!
! subprograms used :
!  rf  mxvdot  dot product of two vectors.
!
      subroutine plmaxl(nf,nc,cf,cfd,ic,cl,cu,cg,s,step,kbc,krem,inew)
      implicit none
      integer nf , nc , ic(*) , kbc , krem , inew
      double precision cf(*) , cfd(*) , cl(*) , cu(*) , cg(*) , s(*) ,  &
                       step
      double precision temp !, mxvdot
      integer jcg , kc
      if ( kbc>0 ) then
         jcg = 1
         do kc = 1 , nc
            if ( krem>0 .and. ic(kc)>10 ) ic(kc) = ic(kc) - 10
            if ( ic(kc)>0 .and. ic(kc)<=10 ) then
               temp = mxvdot(nf,cg(jcg),s)
               cfd(kc) = temp
               if ( temp<0.0d0 ) then
                  if ( ic(kc)==1 .or. ic(kc)>=3 ) then
                     temp = (cl(kc)-cf(kc))/temp
                     if ( temp<=step ) then
                        inew = kc
                        step = temp
                     endif
                  endif
               elseif ( temp>0.0d0 ) then
                  if ( ic(kc)==2 .or. ic(kc)>=3 ) then
                     temp = (cu(kc)-cf(kc))/temp
                     if ( temp<=step ) then
                        inew = kc
                        step = temp
                     endif
                  endif
               endif
            elseif ( ic(kc)<-10 ) then
               temp = mxvdot(nf,cg(jcg),s)
               cfd(kc) = temp
               if ( temp>0.0d0 ) then
                  if ( ic(kc)==-11 .or. ic(kc)==-13 .or. ic(kc)==-15 )  &
                       then
                     temp = (cl(kc)-cf(kc))/temp
                     if ( temp<=step ) then
                        inew = kc
                        step = temp
                     endif
                  endif
               elseif ( temp<0.0d0 ) then
                  if ( ic(kc)==-12 .or. ic(kc)==-14 .or. ic(kc)==-16 )  &
                       then
                     temp = (cu(kc)-cf(kc))/temp
                     if ( temp<=step ) then
                        inew = kc
                        step = temp
                     endif
                  endif
               endif
            endif
            jcg = jcg + nf
         enddo
      endif
      end subroutine plmaxl

! subroutine plmaxs               all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
! determination of the maximum stepsize using the simple bounds
! for variables.
!
! parameters :
!  ii  nf  number of variables.
!  ri  x(nf)  vector of variables.
!  ii  ix(nf)  vector containing types of bounds.
!  ri  xl(nf)  vector containing lower bounds for variables.
!  ri  xu(nf)  vector containing upper bounds for variables.
!  ri  s(nf)  direction vector.
!  ro  step  maximum stepsize.
!  ii  kbf  specification of simple bounds. kbf=0-no simple bounds.
!         kbf=1-one sided simple bounds. kbf=2=two sided simple bounds.
!  io  krem  indication of linearly dependent gradients.
!  io  inew  index of the new active constraint.
!
      subroutine plmaxs(nf,x,ix,xl,xu,s,step,kbf,krem,inew)
      implicit none
      integer nf , ix(*) , kbf , krem , inew
      double precision x(*) , xl(*) , xu(*) , s(*) , step
      double precision temp
      integer i
      if ( kbf>0 ) then
         do i = 1 , nf
            if ( krem>0 .and. ix(i)>10 ) ix(i) = ix(i) - 10
            if ( ix(i)>0 .and. ix(i)<=10 ) then
               if ( s(i)<0.0d0 ) then
                  if ( ix(i)==1 .or. ix(i)>=3 ) then
                     temp = (xl(i)-x(i))/s(i)
                     if ( temp<=step ) then
                        inew = -i
                        step = temp
                     endif
                  endif
               elseif ( s(i)>0.0d0 ) then
                  if ( ix(i)==2 .or. ix(i)>=3 ) then
                     temp = (xu(i)-x(i))/s(i)
                     if ( temp<=step ) then
                        inew = -i
                        step = temp
                     endif
                  endif
               endif
            endif
         enddo
      endif
      krem = 0
      end subroutine plmaxs

! subroutine plnewl             all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
! test on activity of a given linear constraint.
!
! parameters :
!  ii  kc  index of a given constraint.
!  ri  cf(nc)  vector containing values of the constraint functions.
!  iu  ic(nc)  vector containing types of constraints.
!  ri  cl(nc)  vector containing lower bounds for constraint functions.
!  ri  cu(nc)  vector containing upper bounds for constraint functions.
!  ri  eps9  tolerance for active constraints.
!  io  inew  index of the new active constraint.
!
      subroutine plnewl(kc,cf,ic,cl,cu,eps9,inew)
      implicit none
      integer kc , ic(*) , inew
      double precision cf(*) , cl(*) , cu(*) , eps9
      double precision temp
      if ( ic(kc)<-10 ) ic(kc) = -ic(kc) - 10
      if ( ic(kc)<=0 ) then
      elseif ( ic(kc)==1 ) then
         temp = eps9*max(abs(cl(kc)),1.0d0)
         if ( cf(kc)>cl(kc)+temp ) then
         elseif ( cf(kc)>=cl(kc)-temp ) then
            ic(kc) = 11
            inew = kc
         else
            ic(kc) = -11
         endif
      elseif ( ic(kc)==2 ) then
         temp = eps9*max(abs(cu(kc)),1.0d0)
         if ( cf(kc)<cu(kc)-temp ) then
         elseif ( cf(kc)<=cu(kc)+temp ) then
            ic(kc) = 12
            inew = kc
         else
            ic(kc) = -12
         endif
      elseif ( ic(kc)==3 .or. ic(kc)==4 ) then
         temp = eps9*max(abs(cl(kc)),1.0d0)
         if ( cf(kc)>cl(kc)+temp ) then
            temp = eps9*max(abs(cu(kc)),1.0d0)
            if ( cf(kc)<cu(kc)-temp ) then
            elseif ( cf(kc)<=cu(kc)+temp ) then
               ic(kc) = 14
               inew = kc
            else
               ic(kc) = -14
            endif
         elseif ( cf(kc)>=cl(kc)-temp ) then
            ic(kc) = 13
            inew = kc
         else
            ic(kc) = -13
         endif
      elseif ( ic(kc)==5 .or. ic(kc)==6 ) then
         temp = eps9*max(abs(cl(kc)),1.0d0)
         if ( cf(kc)>cl(kc)+temp ) then
            temp = eps9*max(abs(cu(kc)),1.0d0)
            if ( cf(kc)<cu(kc)-temp ) then
            elseif ( cf(kc)<=cu(kc)+temp ) then
               ic(kc) = 16
               inew = kc
            else
               ic(kc) = -16
            endif
         elseif ( cf(kc)>=cl(kc)-temp ) then
            ic(kc) = 15
            inew = kc
         else
            ic(kc) = -15
         endif
      endif
      end subroutine plnewl

! subroutine plminn             all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  rf  mxvdot  dot product of two vectors.
!
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

! subroutine plmins             all systems                   91/12/01
! 91/12/01 lu : original version
!
! purpose :
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

! subroutine plnews             all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
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

! subroutine plredl               all systems                   98/12/01
! 98/12/01 lu : original version
!
! purpose :
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

! subroutine plrmf0             all systems                   91/12/01
! 91/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  s   plrmr0  correction of kernel of the orthogonal projection
!         after constraint deletion.
!
      subroutine plrmf0(nf,nc,ix,ia,iaa,ar,ic,s,n,iold,krem,ier)
      implicit none
      integer ier , iold , krem , n , nc , nf
      double precision ar(*) , s(*)
      integer ia(*) , iaa(*) , ic(*) , ix(*)
      integer nadd , ndec , nfg , nfh , nfv , nit , nrem , nres
      integer l
      common /stat  / nres , ndec , nrem , nadd , nit , nfv , nfg , nfh
      call plrmr0(nf,iaa,ar,s,n,iold,krem,ier)
      n = n + 1
      nrem = nrem + 1
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

! subroutine plrmr0               all systems                91/12/01
! 91/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  s   mxvcop  copying of a vector.
!  s   mxvort  determination of an elementary orthogonal matrix for
!         plane rotation.
!  s   mxvrot  plane rotation of a vector.
!  s   mxvset  initiation of a vector.
!
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

! subroutine plsetc             all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  s   mxvdif  difference of two vectors.
!  rf  mxvdot  dot product of two vectors.
!  s   mxvmul  diagonal premultiplication of a vector.
!
      subroutine plsetc(nf,nc,x,xo,cf,ic,cg,s)
      implicit none
      integer nf , nc , ic(*)
      double precision x(*) , xo(*) , cf(*) , cg(*) , s(*)
      !double precision mxvdot
      integer jcg , kc
      call mxvdif(nf,x,xo,s)
      jcg = 0
      do kc = 1 , nc
         if ( ic(kc)/=0 ) cf(kc) = cf(kc) + mxvdot(nf,s,cg(jcg+1))
         jcg = jcg + nf
      enddo
      end subroutine plsetc

! subroutine plsetg             all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  s   mxvdir  vector augmented by the scaled vector.
!  s   mxvset  initiation of a vector.
!
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

! subroutine pltlag               all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
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

! subroutine pltrbg               all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  s   plvlag  gradient is premultiplied by the matrix whose columns
!         are normals of the active constraints.
!  s   pltlag  computation of the maximum absolute value of the negative
!         lagrange multiplier.
!  s   mxdrmm  premultiplication of a vector by a rowwise stored dense
!         rectangular matrix.
!  s   mxdprb  back substitution after a choleski decomposition.
!  rf  mxvmax  l-infinity norm of a vector.
!  s   mxvset  initiation of a vector.
!
      subroutine pltrbg(nf,n,nc,ix,ic,ica,cg,cr,cz,g,gn,eps7,gmax,umax, &
                        iold)
      implicit none
      integer nf , n , nc , ix(*) , ic(*) , ica(*) , iold
      double precision cg(*) , cr(*) , cz(*) , g(*) , gn(*) , eps7 ,    &
                       gmax , umax
      !double precision mxvmax
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

! subroutine plvlag               all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  rf  mxvdot  dot product of two vectors.
!
      subroutine plvlag(nf,n,nc,iaa,ag,cg,g,gn)
      implicit none
      integer nf , n , nc , iaa(*)
      double precision ag(*) , cg(*) , g(*) , gn(*)
      !double precision mxvdot
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

! subroutine pnint1                all systems                91/12/01
! 91/12/01 lu : original version
!
! purpose :
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
! method :
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

! subroutine pnint3                all systems                91/12/01
! 91/12/01 lu : original version
!
! purpose :
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
! method :
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

! subroutine pnstep                all systems                89/12/01
! 89/01/01 lu : original version
!
! purpose :
! determination of a scaling factor for the boundary step.
!
! parameters :
!  ri  del  maximum stepsize.
!  ri  a  input parameter.
!  ri  b  input parameter.
!  ri  c  input parameter.
!  ro  alf  scaling factor for the boundary step such that
!         a**2+2*b*alf+c*alf**2=del**2.
!
      subroutine pnstep(del,a,b,c,alf)
      implicit none
      double precision del , a , b , c , alf
      double precision den , dis
      double precision zero
      parameter (zero=0.0d0)
      alf = zero
      den = (del+a)*(del-a)
      if ( den<=zero ) return
      dis = b*b + c*den
      if ( b>=zero ) then
         alf = den/(sqrt(dis)+b)
      else
         alf = (sqrt(dis)-b)/c
      endif
      end subroutine pnstep

! subroutine pp0af8             all systems                 97/12/01
! 97/12/01 lu : original version
!
! purpose :
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

! subroutine ppset2             all systems                   97/12/01
! 97/12/01 lu : original version
!
! purpose :
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

! subroutine ps0g01                all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
! simple search with trust region update.
!
! parameters :
!  ro  r  value of the stepsize parameter.
!  ro  f  value of the objective function.
!  ri  fo  initial value of the objective function.
!  ri  po  initial value of the directional derivative.
!  ri  pp  quadratic part of the predicted function value.
!  ru  xdel  trust region bound.
!  ro  xdelo  previous trust region bound.
!  ri  xmax maximum stepsize.
!  ri  rmax  maximum value of the stepsize parameter.
!  ri  snorm  euclidean norm of the direction vector.
!  ri  bet1  lower bound for stepsize reduction.
!  ri  bet2  upper bound for stepsize reduction.
!  ri  gam1  lower bound for stepsize expansion.
!  ri  gam2  upper bound for stepsize expansion.
!  ri  eps4  first tolerance for ratio df/dfpred. step bound is
!         decreased if df/dfpred<eps4.
!  ri  eps5  second tolerance for ratio df/dfpred. step bound is
!         increased if it is active and df/dfpred>eps5.
!  ii  kd  degree of required dervatives.
!  io  ld  degree of previously computed derivatives.
!  iu  idir indicator for direction determination.
!         idir=0-basic determination. idir=1-determination
!         after stepsize reduction. idir=2-determination after
!         stepsize expansion.
!  io  iters  termination indicator. iters=0-zero step. iters=1-step
!         bound was decreased. iters=2-step bound was unchanged.
!         iters=3-step bound was increased. iters=6-first stepsize.
!  ii  iterd termination indicator. iterd<0-bad decomposition.
!         iterd=0-descent direction. iterd=1-newton like step.
!         iterd=2-inexact newton like step. iterd=3-boundary step.
!         iterd=4-direction with the negative curvature.
!         iterd=5-marquardt step.
!  io  maxst  maximum stepsize indicator. maxst=0 or maxst=1 if maximum
!         stepsize was not or was reached.
!  io  nred  actual number of extrapolations or interpolations.
!  ii  mred  maximum number of extrapolations or interpolations.
!  ii  kters  termination selection. kters=1-normal termination.
!         kters=6-first stepsize.
!  ii  mes1  switch for extrapolation. mes1=1-constant increasing of
!         the interval. mes1=2-extrapolation specified by the parameter
!         mes. mes1=3 suppressed extrapolation.
!  ii  mes2  switch for termination. mes2=1-normal termination.
!         mes2=2-termination after at least two steps (asymptotically
!         perfect line search).
!  ii  mes3  safeguard against rounding errors. mes3=0-safeguard
!         suppressed. mes3=1-first level of safeguard. mes3=2-second
!         level of safeguard.
!  iu  isys  control parameter.
!
! common data :
!
! method :
! g.a.schultz, r.b.schnabel, r.h.byrd: a family of trust-region-based
! algorithms for unconstrained minimization with strong global
! convergence properties, siam j. numer.anal. 22 (1985) pp. 47-67.
!
      subroutine ps0g01(r,f,fo,po,pp,xdel,xdelo,xmax,rmax,snorm,bet1,   &
                        bet2,gam1,gam2,eps4,eps5,kd,ld,idir,iters,iterd,&
                        maxst,nred,mred,kters,mes1,mes2,mes3,isys)
      implicit none
      integer kd , ld , idir , iters , iterd , maxst , nred , mred ,    &
              kters , mes1 , mes2 , mes3 , isys
      double precision r , f , fo , po , pp , xdel , xdelo , xmax ,     &
                       rmax , snorm , bet1 , bet2 , gam1 , gam2 , eps4 ,&
                       eps5
      double precision df , dfpred
      integer nred1 , nred2
      save nred1 , nred2
      if ( isys==1 ) then
         if ( kters<0 .or. kters>5 ) then
            iters = 6
         else
            df = fo - f
            dfpred = -r*(po+r*pp)
            if ( df<eps4*dfpred ) then
!
!     step is too large, it has to be reduced
!
               if ( mes1==1 ) then
                  xdel = bet2*snorm
               elseif ( mes1==2 ) then
                  xdel = bet2*min(0.5d0*xdel,snorm)
               else
                  xdel = 0.5d0*po*snorm/(po+df)
                  xdel = max(xdel,bet1*snorm)
                  xdel = min(xdel,bet2*snorm)
               endif
               iters = 1
               if ( mes3<=1 ) then
                  nred2 = nred2 + 1
               else
                  if ( iterd>2 ) nred2 = nred2 + 1
               endif
            elseif ( df<=eps5*dfpred ) then
!
!     step is suitable
!
               iters = 2
            else
!
!     step is too small, it has to be enlarged
!
               if ( mes2==2 ) then
                  xdel = max(xdel,gam1*snorm)
               elseif ( iterd>2 ) then
                  xdel = gam1*xdel
               endif
               iters = 3
            endif
            xdel = min(xdel,xmax,gam2*snorm)
            if ( fo<=f ) then
               if ( nred1>=mred ) then
                  iters = -1
               else
                  idir = 1
                  iters = 0
                  nred1 = nred1 + 1
               endif
            endif
         endif
         maxst = 0
         if ( xdel>=xmax ) maxst = 1
         if ( mes3==0 ) then
            nred = nred1
         else
            nred = nred2
         endif
         isys = 0
         return
      endif
!      go to (1,2) isys+1
      if ( idir==0 ) then
         nred1 = 0
         nred2 = 0
      endif
      idir = 0
      xdelo = xdel
!
!     computation of the new function value
!
      r = min(1.0d0,rmax)
      kd = 0
      ld = -1
      isys = 1
      end subroutine ps0g01

! subroutine ps0l02                all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
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
! subprogram used :
!  s   pnint3  extrapolation or interpolation without directional
!         derivatives.
!
! method :
! safeguarded extrapolation and interpolation with extended termination
! criteria.
!
      subroutine ps0l02(r,ro,rp,f,fo,fp,po,pp,fmin,fmax,rmin,rmax,tols, &
                        kd,ld,nit,kit,nred,mred,maxst,iest,inits,iters, &
                        kters,mes,isys)
      implicit none
      integer kd , ld , nit , kit , nred , mred , maxst , iest , inits ,&
              iters , kters , mes , isys
      double precision r , ro , rp , f , fo , fp , po , pp , fmin ,     &
                       fmax , rmin , rmax , tols
      double precision rl , fl , ru , fu , ri , fi , rtemp , tol
      integer mtyp , merr , mode , init1 , mes1 , mes2
      logical l1 , l2 , l3 , l4 , l6 , l7
      parameter (tol=1.0d-4)
      save mtyp , mode , mes1 , mes2
      save rl , fl , ru , fu , ri , fi
      if ( isys/=1 ) then
!      go to (1,3) isys+1
         mes1 = 2
         mes2 = 2
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
         mode = 0
         rl = 0.0d0
         fl = fo
         ru = 0.0d0
         fu = fo
         ri = 0.0d0
         fi = fo
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
            l4 = f - fo>=(1.0d0-tols)*r*po .or. mes2==2 .and. mode==2
            l6 = ru - rl<=tol*ru .and. mode==2
            l7 = mes2<=2 .or. mode/=0
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
         elseif ( l2 .and. .not.f>=fu ) then
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
         elseif ( l3 .and. l4 .and. l7 .and.                            &
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
            mode = max(mode,1)
            mtyp = abs(mes)
            if ( f>=fmax ) mtyp = 1
         endif
         if ( mode==1 ) then
!
!     interval change after extrapolation
!
            rl = ri
            fl = fi
            ri = ru
            fi = fu
            ru = r
            fu = f
            if ( f>=fi ) then
               nred = 0
               mode = 2
            elseif ( mes1==1 ) then
               mtyp = 1
            endif
!
!     interval change after interpolation
!
         elseif ( r<=ri ) then
            if ( f<=fi ) then
               ru = ri
               fu = fi
               ri = r
               fi = f
            else
               rl = r
               fl = f
            endif
         elseif ( f<=fi ) then
            rl = ri
            fl = fi
            ri = r
            fi = f
         else
            ru = r
            fu = f
         endif
      endif
!
!     new stepsize selection (extrapolation or interpolation)
!
      call pnint3(ro,rl,ru,ri,fo,fl,fu,fi,po,r,mode,mtyp,merr)
      if ( merr>0 ) then
         iters = -merr
         isys = 0
         return
      elseif ( mode==1 ) then
         nred = nred - 1
         r = min(r,rmax)
      elseif ( mode==2 ) then
         nred = nred + 1
      endif
!
!     computation of the new function value
!
      kd = 0
      ld = -1
      isys = 1
      end subroutine ps0l02

! subroutine ps1l01                all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
!  standard line search with directional derivatives.
!
! parameters :
!  ro  r  value of the stepsize parameter.
!  ro  rp  previous value of the stepsize parameter.
!  ro  f  value of the objective function.
!  ri  fo  initial value of the objective function.
!  ro  fp  previous value of the objective function.
!  ro  p  value of the directional derivative.
!  ri  po  initial value of the directional derivative.
!  ro  pp  previous value of the directional derivative.
!  ri  fmin  lower bound for value of the objective function.
!  ri  fmax  upper bound for value of the objective function.
!  ri  rmin  minimum value of the stepsize parameter
!  ri  rmax  maximum value of the stepsize parameter
!  ri  tols  termination tolerance for line search (in test on the
!         change of the function value).
!  ri  tolp  termination tolerance for line search (in test on the
!         change of the directional derivative).
!  ro  par1  parameter for controlled scaling of variable metric
!         updates.
!  ro  par2  parameter for controlled scaling of variable metric
!         updates.
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
! subprogram used :
!  s   pnint1  extrapolation or interpolation with directional
!         derivatives.
!
! method :
! safeguarded extrapolation and interpolation with standard termination
! criteria.
!
      subroutine ps1l01(r,rp,f,fo,fp,p,po,pp,fmin,fmax,rmin,rmax,tols,  &
                        tolp,par1,par2,kd,ld,nit,kit,nred,mred,maxst,   &
                        iest,inits,iters,kters,mes,isys)
      implicit none
      integer kd , ld , nit , kit , nred , mred , maxst , iest , inits ,&
              iters , kters , mes , isys
      double precision r , rp , f , fo , fp , p , po , pp , fmin ,      &
                       fmax , rmin , rmax , tols , tolp , par1 , par2
      double precision rl , fl , pl , ru , fu , pu , rtemp
      integer mtyp , merr , mode , init1 , mes1 , mes2 , mes3
      logical l1 , l2 , l3 , l5 , l7 , m1 , m2 , m3
      double precision con , con1
      parameter (con=1.0d-2,con1=1.0d-13)
      save mtyp , mode , mes1 , mes2 , mes3
      save rl , fl , pl , ru , fu , pu
      if ( isys==1 ) then
         if ( mode==0 ) then
            par1 = p/po
            par2 = f - fo
         endif
         if ( iters/=0 ) then
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
               l3 = f - fo<=tols*r*po
               l5 = p>=tolp*po .or. mes2==2 .and. mode==2
               l7 = mes2<=2 .or. mode/=0
               m1 = .false.
               m2 = .false.
               m3 = l3
               if ( mes3>=1 ) then
                  m1 = abs(p)<=con*abs(po) .and. fo - f>=(con1/con)     &
                       *abs(fo)
                  l3 = l3 .or. m1
               endif
               if ( mes3>=2 ) then
                  m2 = abs(p)<=0.5d0*abs(po) .and. abs(fo-f)            &
                       <=2.0d0*con1*abs(fo)
                  l3 = l3 .or. m2
               endif
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
            elseif ( l2 .and. l3 .and. .not.l5 ) then
               iters = 7
               isys = 0
               return
            elseif ( m3 .and. mes1==3 ) then
               iters = 5
               isys = 0
               return
            elseif ( l3 .and. l5 .and. l7 ) then
               iters = 4
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
               pp = p
               mode = max(mode,1)
               mtyp = abs(mes)
               if ( f>=fmax ) mtyp = 1
            endif
            if ( mode==1 ) then
!
!     interval change after extrapolation
!
               rl = ru
               fl = fu
               pl = pu
               ru = r
               fu = f
               pu = p
               if ( .not.l3 ) then
                  nred = 0
                  mode = 2
               elseif ( mes1==1 ) then
                  mtyp = 1
               endif
!
!     interval change after interpolation
!
            elseif ( .not.l3 ) then
               ru = r
               fu = f
               pu = p
            else
               rl = r
               fl = f
               pl = p
            endif
         endif
      else
!      go to (1,3) isys+1
         mes1 = 2
         mes2 = 2
         mes3 = 2
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
         r = max(r,rmin)
         r = min(r,rmax)
         mode = 0
         ru = 0.0d0
         fu = fo
         pu = po
      endif
!
!     new stepsize selection (extrapolation or interpolation)
!
      call pnint1(rl,ru,fl,fu,pl,pu,r,mode,mtyp,merr)
      if ( merr>0 ) then
         iters = -merr
         isys = 0
         return
      elseif ( mode==1 ) then
         nred = nred - 1
         r = min(r,rmax)
      elseif ( mode==2 ) then
         nred = nred + 1
      endif
!
!     computation of the new function value and the new directional
!     derivative
!
      kd = 1
      ld = -1
      isys = 1
      end subroutine ps1l01

! subroutine pudbg1                all systems                92/12/01
! 92/12/01 lu : original version
!
! purpose :
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
! subprograms used :
!  s   mxdpgu  correction of a dense symmetric positive definite
!         matrix in the factored form b=l*d*trans(l).
!  s   mxdpgs  scaling of a dense symmetric positive definite matrix
!         in the factored form b=l*d*trans(l).
!  s   mxvdif  difference of two vectors.
!  rf  mxvdot  dot product of vectors.
!  s   mxvscl  scaling of a vector.
!
! method :
! bfgs variable metric method.
!
      subroutine pudbg1(n,h,g,s,xo,go,r,po,nit,kit,iterh,met,met1,mec)
      implicit none
      double precision po , r
      integer iterh , kit , met , met1 , mec , n , nit
      double precision g(*) , go(*) , h(*) , s(*) , xo(*)
      double precision a , b , c , gam , par , den , dis
      logical l1 , l3
      !double precision mxvdot , mxdpgp
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

! subroutine pudbi1                all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
! variable metric updates.
!
! parameters :
!  ii  n  number of variables.
!  ru  h(n*(n+1)/2)  updated approximation of the inverse hessian
!         matrix.
!  ra  s(n)  auxiliary vector.
!  ri  xo(n)  vector of variables difference.
!  ri  go(n)  gradient difference.
!  ro  r  value of the stepsize parameter.
!  ri  po  initial value of the directional derivative.
!  ro  par1  parameter for control scaling.
!  ro  par2  parameter for control scaling.
!  ro  f  value of the objective function.
!  ri  fo  initial value of the objective function.
!  ri  p  current value of the directional derivative.
!  ii  nit  number of iterations.
!  ii  kit  index of the iteration with the last restart.
!  ii  met  variable metric update.
!  ii  met1  scaling strategy.
!  ii  met2  correction rule.
!  iu  idecf  decomposition indicator.
!  ii  iterd  termination indicator. iterd<0-bad decomposition.
!         iterd=0-descent direction. iterd=1-newton like step.
!         iterd=2-inexact newton like step. iterd=3-boundary step.
!         iterd=4-direction with the negative curvature.
!         iterd=5-marquardt step.
!  io  iterh  update indicator. iterh=0-successful update.
!         iterh>0-unsuccessful update.
!
! method :
! various variable metric updates including bfgs update.
!
      subroutine pudbi1(n,h,s,xo,go,r,po,par1,par2,f,fo,p,nit,kit,met,  &
                        met1,met2,idecf,iterd,iterh)
      implicit none
      integer n , nit , kit , met , met1 , met2 , idecf , iterd , iterh
      double precision h(*) , s(*) , xo(*) , go(*) , r , po
      double precision par1 , par2
      double precision f , fo , p
      double precision aa , cc
      !double precision mxvdot
      double precision dis , pom , pom3 , pom4 , a , b , c , gam , rho ,&
                       par
      double precision den
      logical l1 , l2 , l3
      if ( met>0 ) then
         if ( idecf/=9 ) then
            iterh = -1
            goto 100
         endif
         l1 = abs(met1)>=3 .or. abs(met1)==2 .and. nit==kit
         l3 = .not.l1
!
!     determination of the parameters a, b, c
!
         b = mxvdot(n,xo,go)
         if ( b<=0.0d0 ) then
            iterh = 2
            goto 100
         endif
         call mxdsmm(n,h,go,s)
         a = mxvdot(n,go,s)
         if ( a<=0.0d0 ) then
            iterh = 1
            goto 100
         endif
         if ( .not.(met/=1 .or. l1) ) then
            c = 0.0d0
         elseif ( iterd/=1 ) then
            c = 0.0d0
         else
            c = -r*po
            if ( c<=0.0d0 ) then
               iterh = 3
               goto 100
            endif
         endif
!
!     determination of the parameter rho (nonquadratic properties)
!
         if ( met2==1 ) then
            rho = 1.0d0
         elseif ( fo-f+p==0 ) then
            rho = 1.0d0
         else
            rho = 0.5d0*b/(fo-f+p)
         endif
         if ( rho<=1.0d-2 ) rho = 1.0d0
         if ( rho>=1.0d2 ) rho = 1.0d0
         aa = a/b
         cc = c/b
         if ( l1 ) then
!
!     determination of the parameter gam (self scaling)
!
            par = a/b
            pom3 = 0.7d0
            pom4 = 6.0d0
            gam = rho/par
            if ( nit/=kit ) then
               if ( met1==3 ) then
                  l2 = par2<=0.0d0
                  l3 = l2 .and. abs(par1)<=0.2d0
                  l3 = l3 .or. (.not.l2 .and. gam>1.0d0)
                  l3 = l3 .or. (l2 .and. par1<0.0d0 .and. gam>1.0d0)
                  l3 = l3 .or. (l2 .and. par1>0.0d0 .and. gam<1.0d0)
                  l3 = l3 .or. gam<pom3
                  l3 = l3 .or. gam>pom4
               elseif ( met1==4 ) then
                  l3 = gam<pom3 .or. gam>pom4
               endif
            endif
         endif
         if ( l3 ) then
            gam = 1.0d0
            par = rho/gam
         endif
         if ( met/=1 ) then
!
!     new update
!
            pom = 1.0d0/(aa*cc)
            if ( pom<1.0d0 ) then
               pom = max(1.0d-15,(sqrt(c/a)-pom)/(1.0d0-pom))
!
!     general update
!
               den = par + pom*aa
               dis = pom/den
               call mxdsmu(n,h,(par*dis-1.0d0)/a,s)
               call mxvdir(n,-dis,s,xo,s)
               call mxdsmu(n,h,den/b,s)
               goto 50
            endif
         endif
!
!     bfgs update
!
         pom = 1.0d0
         dis = par + aa
         call mxvdir(n,-dis,xo,s,xo)
         dis = 1.0d0/(b*dis)
         call mxdsmu(n,h,dis,xo)
         call mxdsmu(n,h,-dis,s)
!
!     scaling
!
 50      if ( gam/=1.0d0 ) call mxdsms(n,h,gam)
      endif
 100  iterh = 0
      end subroutine pudbi1

! subroutine pudbm2                all systems                92/12/01
! 92/12/01 lu : original version
!
! purpose :
! variable metric update of a dense symmetric positive definite matrix.
!
! parameters :
!  ru  h(m)  positive definite approximation of the hessian
!         matrix.
!  ri  g(nf)  gradient of the objective function.
!  ra  s(nf)  auxiliary vector.
!  ru  xo(nf)  vectors of variables difference.
!  ri  go(nf)  gradients difference.
!
! common data :
!  ii  nf declared number of variables.
!  ii  n  actual number of variables.
!  ii  m  number of nonzero elements of the matrix.
!  ii  met  method selection. met=1-bfgs update. met=2-dfp update.
!         met=3-hoshino update.
!  ii  met1  selection of self scaling.  met1=1-self scaling suppressed.
!         met1=2-initial self scaling.
!  ii  met2  selection of the line search model. met2=1-quadratic model.
!         met2=2 use of taylor expansion.
!  ii  met3  method correction. met3=1-no correction.
!         met3=2-powell's correction.
!  ii  iterd  termination indicator. iterd<0-bad decomposition.
!         iterd=0-descent direction. iterd=1-newton like step.
!         iterd=2-inexact newton like step. iterd=3-boundary step.
!         iterd=4-direction with the negative curvature.
!         iterd=5-marquardt step.
!  io  iterh  termination indicator. iterh<0-bad decomposition.
!         iterh=0-successful update. iterh>0-nonpositive parameters.
!  ii  idecf  decomposition indicator. idecf=0-no decomposition.
!         idecf=1-gill-murray decomposition. idecf=2-bunch-parlett
!         decomposition. idecf=3-inversion.
!  ii  itran  transformation indicator. itran=0 or itran=1 if
!         transformation is not or is used.
!  ii  nit  actual number of iterations.
!  ii  kit  number of the iteration after last restart.
!  ri  r  value of the stepsize parameter.
!  ri  f  new value of the objective function.
!  ri  fo  old value of the objective function.
!  ri  p  new value of the directional derivative.
!  ri  po  old value of the directional derivative.
!  to  tuxx  text information on the correction used.
!
! subprograms used :
!  s   mxdsmm  matrix-vector product.
!  s   mxdsmu  correction of a dense symmetric matrix.
!  s   mxdsms  scaling of a dense symmetric matrix.
!  s   mxvdif  difference of two vectors.
!  s   mxvdir  vector augmented by the scaled vector.
!  rf  mxvdot  dot product of vectors.
!  s   mxvneg  copying of a vector with the change of the sign.
!  s   mxvscl  scaling of a vector.
!  s   uou1d1  print of entry to variable metric update.
!  s   uou1d2  print of exit from variable metric update.
!
! method :
! basic variable metric methods.
!
      subroutine pudbm2(nf,n,h,hh,s,xo,go,so,fo,par,met1,met3,idecf,    &
                        iterh)
      implicit none
      integer nf , n , met1 , met3 , idecf , iterh
      double precision h(nf*(nf+1)/2) , hh(nf*(nf+1)/2) , s(nf) , xo(nf)&
                       , go(nf) , so(nf) , fo , par
      double precision den , a , b , c , gam , pom !, mxvdot
      logical l1
      double precision con
      parameter (con=1.0d-8)
      if ( idecf/=0 ) then
         iterh = -1
         return
      endif
      l1 = met1>=2
!
!     determination of the parameters b, c
!
      call mxdsmm(n,h,xo,s)
      call mxvdif(n,go,s,so)
      if ( met3==2 ) call mxvscl(n,1.0d0/sqrt(fo),so,so)
      b = mxvdot(n,xo,so)
      if ( b<=0.0d0 ) l1 = .false.
      a = 0.0d0
      call mxdsmm(n,hh,xo,s)
      c = mxvdot(n,xo,s)
      if ( c<=0.0d0 ) l1 = .false.
      if ( l1 ) then
!
!     determination of the parameter gam (self scaling)
!
         gam = c/b
      else
         gam = 1.0d0
      endif
      par = gam
!
!     rank one update
!
      den = par*b - c
      if ( abs(den)<=con*max(con,abs(par*b),abs(c)) ) then
         if ( b>0.0d0 .and. c>0.0d0 ) then
!
!     bfgs update
!
            pom = 0.0d0
            call mxdsmu(n,hh,par/b,so)
            if ( c>0.0d0 ) call mxdsmu(n,hh,-1.0d0/c,s)
            goto 100
         else
            iterh = 4
            return
         endif
      endif
      pom = par*b/den
      call mxvdir(n,-par,so,s,s)
      call mxdsmu(n,hh,1.0d0/den,s)
 100  iterh = 0
      if ( gam/=1.0d0 ) call mxdsms(n,hh,1.0d0/gam)
      end subroutine pudbm2

! subroutine pudbq1                all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
! broyden good update of a rectangular matrix after the qr
! decomposition.
!
! parameters :
!  ii  n  number of variables.
!  ii  na  number of equations.
!  ru  h(n*(n+1)/2)  updated upper triangular matrix.
!  ri  eta2  parameter which controls a nonsingularity
!  ru  ag(n*na)  updated rectangular matrix.
!  ra  s(n)  auxiliary vector.
!  ri  xo(n)  vector of variables difference.
!  ri  afo(na)  right hand sides difference.
!  ii  met  variable metric update.
!  io  iterh  update indicator. iterh=0-successful update.
!         iterh>0-unsuccessful update.
!  iu  ideca  decomposition indicator.
!  ii  ndeca  number of decompositions.
!
! method :
! various variable metric updates including bfgs update.
!
      subroutine pudbq1(n,na,h,eta2,ag,s,xo,afo,met,iterh,ideca,ndeca)
      implicit none
      integer n , na , met , inf , iterh , ideca , ndeca
      double precision h(*) , eta2 , ag(*) , s(*) , xo(*) , afo(*)
      double precision den !, mxvdot
      if ( met<=0 ) return
      if ( ideca==0 ) then
!
!     qr decomposition
!
         den = eta2
         call mxdrqf(n,na,ag,h)
         call mxdprc(n,h,inf,den)
         ndeca = ndeca + 1
         ideca = 1
      elseif ( ideca/=1 ) then
         iterh = -1
         return
      endif
!
!     the good broyden update
!
      den = mxvdot(n,xo,xo)
      if ( den<=0.0d0 ) then
         iterh = 2
         return
      endif
      call mxvcop(n,xo,s)
      call mxvneg(n,xo,xo)
      call mxdprm(n,h,xo,1)
      call mxdrmd(n,na,ag,xo,1.0d0,afo,afo)
      call mxdrqu(n,na,ag,h,1.0d0/den,afo,s,xo,inf)
      iterh = 0
      end subroutine pudbq1

! subroutine pudfm1                all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
! variable metric update of a dense symmetric positive definite matrix.
!
! parameters :
!  ii  n  actual number of variables.
!  ru  b(m)  positive definite approximation of the hessian matrix.
!  ri  g(nf)  gradient of the objective function.
!  ra  s(nf)  auxiliary vector.
!  ru  xo(nf)  vectors of variables difference.
!  ri  go(nf)  gradients difference.
!  ri  f  current value of the objective function.
!  ri  fo  previous value of the objective function.
!  ri  eta5  tolerance for a hybrid method.
!  ri  eta9  maximum for real numbers.
!  ii  ipom1  method indicator.
!  io  ipom2  indicator for scaling.
!  ii  met  method selection. met=0-no update. met=1-bfgs update.
!  ii  met1  selection of self scaling.  met1=1-self scaling suppressed.
!         met1=2 self scaling in the first iteration after restart.
!         met1=3-self scaling in each iteration.
!
! common data :
!  io  iterh  termination indicator. iterh<0-bad decomposition.
!         iterh=0-successful update. iterh>0-nonpositive parameters.
!  ii  idecf  decomposition indicator. idecf=0-no decomposition.
!  to  tuxx  text information on the correction used.
!
! subprograms used :
!  s   mxdsmm  matrix-vector product.
!  s   mxdsmu  correction of a dense symmetric matrix.
!  s   mxdsms  scaling of a dense symmetric matrix.
!  rf  mxvdot  dot product of vectors.
!  s   uoerr1  error mesages.
!  s   uou1d1  print of entry to variable metric update.
!  s   uou1d2  print of exit from variable metric update.
!
! method :
!  fletcher's combination of the gauss-newton and the bfgs methods.
!
      subroutine pudfm1(n,b,s,xo,go,f,fo,eta5,ipom1,ipom2,met1,idecf,   &
                        iterh)
      implicit none
      integer n , ipom1 , ipom2 , met1 , idecf , iterh
      double precision b(n*(n+1)/2) , s(n) , xo(n) , go(n) , f , fo ,   &
                       eta5
      !double precision mxvdot
      double precision ab , bb , cb , gam , par
      logical l1
      if ( idecf/=0 ) then
         iterh = -1
         return
      endif
      par = 0.0d0
!
!     determination of the parameters a,b,c
!
      bb = mxvdot(n,xo,go)
      if ( bb<=0.0d0 ) then
         iterh = 2
         ipom1 = 0
         return
      endif
      ab = 0.0d0
      call mxdsmm(n,b,xo,s)
      cb = mxvdot(n,xo,s)
      if ( cb<=0.0d0 ) then
         iterh = 3
         return
      endif
      l1 = met1==4 .or. met1==3 .and. ipom2>=1 .or. met1==2 .and.       &
           ipom2==1
      if ( fo-f>=eta5*fo ) then
         ipom1 = 0
      else
         ipom1 = 1
      endif
      if ( l1 ) then
!
!     determination of the parameter gam (self scaling)
!
         gam = cb/bb
      else
         gam = 1.0d0
      endif
!
!     bfgs update
!
      call mxdsmu(n,b,gam/bb,go)
      call mxdsmu(n,b,-1.0d0/cb,s)
      iterh = 0
      ipom2 = 0
      if ( l1 ) call mxdsms(n,b,1.0d0/gam)
      end subroutine pudfm1

! subroutine pudrv1                all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
! driver for hybrid quasi-newton updates.
!
! parameters:
!  ri  r  value of the stepsize parameter.
!  ri  fo  previous value of the objective function.
!  ri  f  current value of the objective function.
!  ri  po  previous value of the directional derivative.
!  ii  ipom1  update selection.
!  ii  ipom2  method selection.
!  io  nred  actual number of extrapolations or interpolations.
!  ii  irest  restart specification. if irest=0 does not hold then a
!         restart is performed.
!
      subroutine pudrv1(r,fo,f,po,ipom1,ipom2,nred,irest)
      implicit none
      integer ipom1 , ipom2 , nred , irest
      double precision r , fo , f , po
      double precision pom
      double precision con2
      parameter (con2=1.0d-2)
      pom = (fo-f)/fo
      select case (ipom2)
      case (2)
         irest = 1
         if ( pom>=con2 ) then
            ipom1 = 0
         elseif ( f-fo<=r*po ) then
            ipom1 = 0
         else
            ipom1 = 1
            irest = 0
         endif
      case (3)
         irest = 1
         if ( nred<=0 ) then
            if ( ipom1/=1 ) then
               ipom1 = 2
               irest = 0
            else
               ipom1 = 0
            endif
         elseif ( pom>=con2 ) then
            ipom1 = 0
         elseif ( ipom1/=2 ) then
            ipom1 = 1
            irest = 0
         else
            ipom1 = 0
         endif
      case (4)
         irest = 1
         ipom1 = 0
      case default
         irest = 1
         if ( nred<=0 ) then
            ipom1 = 2
            irest = 0
         else
            ipom1 = 0
         endif
      end select
      end subroutine pudrv1

! subroutine pudsd2                all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
! initiation of a dense symmetric positive definite matrix
!
! parameters :
!  ii  n  actual number of variables.
!  ru  h(n*(n+1)/2)  positive definite approximation of the hessian
!         matrix
!  ru  b(n*(n+1)/2)  positive definite approximation of the hessian
!         matrix
!  ri  f  current value of the objective function.
!  ri  fo  previous value of the objective function.
!  ri  eta5  tolerance for a hybrid method.
!  ii  met3  type of structured update. met3=1-standard structured
!         update. met3=2-totally structured update.
!
! common data :
!  ru  ran  random number.
!  ii  idecf  decomposition indicator. idecf=0-no decomposition.
!  ii  irest  restart specification. if irest=0 does not hold then a
!         restart is performed.
!  ii  idir indicator of direction determination. idir=0-basic
!         determination. idir=1-determination after stepsize
!         reduction. idir=2-determination after stepsize expansion.
!  iu  ld  degree of previously computed derivatives.
!  to  tuxx   text information on the restart used.
!
! subprograms used :
!  s   mxdsmi  generation of the unit matrix.
!  s   mxdsdo  initiation of a diagonal matrix.
!  s   mxdsma  dense symmetric matrix augmented
!              by the scaled dense symmetric matrix.
!  s   mxdsmo  a scalar is set to all elements
!              of a dense symmetric matrix.
!  s   mxdsms  scaling of a dense symmetric matrix.
!  s   uyset3  definition of the restart variables.
!
      subroutine pudsd2(n,h,b,f,fo,eta5,met3,ld,idir,idecf,irest,ind)
      implicit none
      integer n , met3 , ld , idir , idecf , irest , ind
      double precision h(n*(n+1)/2) , b(n*(n+1)/2) , f , fo , eta5
      integer iudsd
      save iudsd
      ind = 0
      if ( irest<0 ) then
         call mxdsmi(n,b)
         if ( f<1.0d0 ) call mxdsms(n,b,sqrt(f))
         iudsd = 1
      elseif ( irest==0 ) then
         if ( idir<=0 ) then
            if ( fo-f<=eta5*fo ) then
               if ( met3==2 ) then
                  call mxdsma(n,sqrt(f),b,h,h)
               else
                  call mxdsma(n,1.0d0,b,h,h)
               endif
               ld = min(ld,1)
               iudsd = 0
            endif
         endif
      elseif ( iudsd==0 ) then
         if ( met3==2 ) then
            call mxdsma(n,-sqrt(f),b,h,h)
         else
            call mxdsma(n,-1.0d0,b,h,h)
         endif
         call mxdsmi(n,b)
         if ( f<1.0d0 ) call mxdsms(n,b,sqrt(f))
         iudsd = 1
      else
         call mxdsmi(n,h)
         ld = min(ld,1)
         iudsd = 1
         ind = 1
      endif
      idecf = 0
      end subroutine pudsd2

! subroutine pudsd3                all systems                97/12/01
! 97/12/01 lu : original version
!
! purpose :
! initiation of a dense symmetric positive definite matrix
!
! parameters :
!  ii  n  actual number of variables.
!  ru  h(n*(n+1)/2)  factorization h=l*d*trans(l) of a positive
!         semidefinite hessian matrix.
!  ru  b(n*(n+1)/2)  factorization b=l*d*trans(l) of a positive
!         definite approximation of the hessian matrix.
!  iu  ipom1  method indicator.
!  iu  ipom2  indicator for scaling.
!
! common data :
!  ru  ran  random number.
!  ii  idecf  decomposition indicator. idecf=0-no decomposition.
!  ii  irest  restart specification. if irest=0 does not hold then a
!         restart is performed.
!  ii  iters  termination indicator. iters=0-zero step.
!  iu  ld  degree of previously computed derivatives.
!
! subprograms used :
!  s   mxdsmi  generation of the unit matrix.
!  s   mxdsdo  initiation of a diagonal matrix.
!  s   uudsmc  copying of a dense symmetric matrix.
!  rf  unran1  random number generator.
!  s   uyset3  definition of the restart variables.
!
      subroutine pudsd3(n,h,b,ipom1,ipom2,ld,idecf,iters,irest,ind)
      implicit none
      integer n , ipom1 , ipom2 , ld , idecf , iters , irest , ind
      double precision h(n*(n+1)/2) , b(n*(n+1)/2)
      integer kdecf
      save kdecf
      ind = 0
      if ( .not.(irest==0 .or. irest>0 .and. ipom1==0) ) then
         call mxdsmi(n,b)
         ipom2 = 1
         kdecf = -1
      endif
      if ( ipom1==1 ) then
         if ( iters>0 .or. irest>0 ) then
            call mxdsmc(n,b,h)
            ld = min(ld,1)
            if ( ipom2==1 ) idecf = kdecf
            if ( irest>0 ) ind = 1
         endif
      elseif ( irest>0 ) then
         ipom1 = 1
         call mxdsmc(n,b,h)
         ld = min(ld,1)
         if ( ipom2==1 ) idecf = kdecf
      endif
      end subroutine pudsd3

! subroutine pyadb4             all systems                   98/12/01
! 98/12/01 lu : original version
!
! purpose :
! new linear constraints or new simple bounds are added to the active
! set. gill-murray factorization of the transformed hessian matrix
! approximation is updated.
!
! parameters :
!  ii  nf  declared number of variables.
!  iu  n  actual number of variables.
!  ii  nc  number of linearized constraints.
!  ri  x(nf)  vector of variables.
!  iu  ix(nf)  vector containing types of bounds.
!  ri  xl(nf)  vector containing lower bounds for variables.
!  ri  xu(nf)  vector containing upper bounds for variables.
!  ri  cf(nc)  vector containing values of the constraint functions.
!  ri  cfd(nc) vector containing increments of the constraint functions.
!  iu  ic(nc)  vector containing types of constraints.
!  iu  ica(nf)  vector containing indices of active constraints.
!  ri  cl(nc)  vector containing lower bounds for constraint functions.
!  ri  cu(nc)  vector containing upper bounds for constraint functions.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ru  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ru  cz(nf*nf)  matrix whose columns are basic vectors from the
!         current reduced subspace.
!  ru  h(nf*(nf+1)/2)  gill-murray factorization of the transformed
!         hessian matrix approximation.
!  ra  s(nf)  auxiliary vector.
!  ri  r  value of the stepsize parameter.
!  ri  eps7  tolerance for linear independence of constraints.
!  ri  eps9  tolerance for active constraints.
!  ro  gmax  maximum absolute value of a partial derivative.
!  ro  umax  maximum absolute value of a negative lagrange multiplier.
!  ii  kbf  type of simple bounds. kbf=0-no simple bounds. kbf=1-one
!         sided simple bounds. kbf=2-two sided simple bounds.
!  ii  kbc  type of constraints. kbc=0-no constraints. kbc=1-constraints
!         with one sided bounds. kbc=2-constraints with two sided
!         bounds.
!  iu  inew  index of the new active constraint.
!  io  ier  error indicator.
!  io  iterm  termination indicator.
!
! common data :
!  iu  nadd  number of constraint additions.
!
! subprograms used :
!  s   pladb4  addition of a new active constraint.
!  s   plnews  identification of active upper bounds.
!  s   plnewl  identification of active linear constrainrs.
!  s   pldirl  new values of constraint functions.
!  s   mxvind  change of the integer vector for constraint addition.
!
      subroutine pyadb4(nf,n,nc,x,ix,xl,xu,cf,cfd,ic,ica,cl,cu,cg,cr,cz,&
                        h,s,r,eps7,eps9,gmax,umax,kbf,kbc,inew,ier,     &
                        iterm)
      implicit none
      integer nf , n , nc , ix(*) , ic(*) , ica(*) , kbf , kbc , inew , &
              ier , iterm
      double precision x(*) , xl(*) , xu(*) , cf(*) , cfd(*) , cl(*) ,  &
                       cu(*) , cg(*) , cr(*) , cz(*) , h(*) , s(*) , r ,&
                       eps7 , eps9 , gmax , umax
      integer i , j , k , l , ij , ik , kc , kj , kk , ll
      double precision den , temp
      integer nres , ndec , nrem , nadd , nit , nfv , nfg , nfh
      common /stat  / nres , ndec , nrem , nadd , nit , nfv , nfg , nfh
      if ( kbc>0 ) then
         if ( r/=0.0d0 ) call pldirl(nc,cf,cfd,ic,r,kbc)
         if ( inew/=0 ) then
            if ( kbf>0 ) then
               do i = 1 , nf
                  inew = 0
                  call plnews(x,ix,xl,xu,eps9,i,inew)
                  call pladb4(nf,n,ica,cg,cr,cz,h,s,eps7,gmax,umax,9,   &
                              inew,nadd,ier)
                  call mxvind(ix,i,ier)
                  if ( ier<0 ) then
                     iterm = -15
                     return
                  endif
               enddo
            endif
            do kc = 1 , nc
               inew = 0
               call plnewl(kc,cf,ic,cl,cu,eps9,inew)
               call pladb4(nf,n,ica,cg,cr,cz,h,s,eps7,gmax,umax,9,inew, &
                           nadd,ier)
               call mxvind(ic,kc,ier)
               if ( ier<0 ) then
                  iterm = -15
                  return
               endif
            enddo
         endif
      elseif ( kbf>0 ) then
         k = 0
         do l = 1 , nf
            if ( ix(l)>=0 ) k = k + 1
            inew = 0
            call plnews(x,ix,xl,xu,eps9,l,inew)
            if ( inew/=0 ) then
               ix(l) = 10 - ix(l)
               kk = k*(k-1)/2
               den = h(kk+k)
               if ( den/=0.0d0 ) then
                  ij = 0
                  kj = kk
                  do j = 1 , n
                     if ( j<=k ) then
                        kj = kj + 1
                     else
                        kj = kj + j - 1
                     endif
                     if ( j/=k ) temp = h(kj)/den
                     ik = kk
                     do i = 1 , j
                        if ( i<=k ) then
                           ik = ik + 1
                        else
                           ik = ik + i - 1
                        endif
                        ij = ij + 1
                        if ( i/=k .and. j/=k ) h(ij) = h(ij)            &
                             + temp*h(ik)
                     enddo
                  enddo
               endif
               ll = kk + k
               do i = k + 1 , n
                  do j = 1 , i
                     ll = ll + 1
                     if ( j/=k ) then
                        kk = kk + 1
                        h(kk) = h(ll)
                     endif
                  enddo
               enddo
               n = n - 1
            endif
         enddo
      endif
      end subroutine pyadb4

! subroutine pyfut1                all systems                98/12/01
! 98/12/01 lu : original version
!
! purpose :
! termination criteria and test on restart.
!
! parameters :
!  ii  n  actual number of variables.
!  ri  f  new value of the objective function.
!  ri  fo  old value of the objective function.
!  ri  umax  maximun absolute value of the negative lagrange multiplier.
!  ro  gmax  norm of the transformed gradient.
!  ri  dmax  maximum relative difference of variables.
!  ri  tolx  lower bound for steplength.
!  ri  tolf  lower bound for function decrease.
!  ri  tolb  lower bound for function value.
!  ri  tolg  lower bound for gradient.
!  ii  kd  degree of required derivatives.
!  iu  nit  actual number of iterations.
!  ii  kit  number of the iteration after restart.
!  ii  mit  maximum number of iterations.
!  iu  nfv  actual number of computed function values.
!  ii  mfv  maximum number of computed function values.
!  iu  nfg  actual number of computed gradient values.
!  ii  mfg  maximum number of computed gradient values.
!  iu  ntesx  actual number of tests on steplength.
!  ii  mtesx  maximum number of tests on steplength.
!  iu  ntesf  actual number of tests on function decrease.
!  ii  mtesf  maximum number of tests on function decrease.
!  ii  ites  system varible which specifies termination. if ites=0
!         then termination is suppressed.
!  ii  ires1  restart specification. restart is performed after
!         ires1*n+ires2 iterations.
!  ii  ires2  restart specification. restart is performed after
!         ires1*n+ires2 iterations.
!  iu  irest  restart indicator. restart is performed if irest>0.
!  ii  iters  termination indicator for steplength determination.
!         iters=0 for zero step.
!  io  iterm  termination indicator. iterm=1-termination after mtesx
!         unsufficient steplengths. iterm=2-termination after mtesf
!         unsufficient function decreases. iterm=3-termination on lower
!         bound for function value. iterm=4-termination on lower bound
!         for gradient. iterm=11-termination after maximum number of
!         iterations. iterm=12-termination after maximum number of
!         computed function values.
!
      subroutine pyfut1(n,f,fo,umax,gmax,dmax,tolx,tolf,tolb,tolg,kd,   &
                        nit,kit,mit,nfv,mfv,nfg,mfg,ntesx,mtesx,ntesf,  &
                        mtesf,ites,ires1,ires2,irest,iters,iterm)
      implicit none
      integer n , kd , nit , kit , mit , nfv , mfv , nfg , mfg , ntesx ,&
              mtesx , ntesf , mtesf , ites , ires1 , ires2 , irest ,    &
              iters , iterm
      double precision f , fo , umax , gmax , dmax , tolx , tolf ,      &
                       tolg , tolb
      double precision temp
      if ( iterm<0 ) return
      if ( ites>0 ) then
         if ( iters/=0 ) then
            if ( nit<=0 ) fo = f + min(sqrt(abs(f)),abs(f)/1.0d1)
            if ( f<=tolb ) then
               iterm = 3
               return
            endif
            if ( kd>0 ) then
               if ( gmax<=tolg .and. umax<=tolg ) then
                  iterm = 4
                  return
               endif
            endif
            if ( nit<=0 ) then
               ntesx = 0
               ntesf = 0
            endif
            if ( dmax<=tolx ) then
               iterm = 1
               ntesx = ntesx + 1
               if ( ntesx>=mtesx ) return
            else
               ntesx = 0
            endif
            temp = abs(fo-f)/max(abs(f),1.0d0)
            if ( temp<=tolf ) then
               iterm = 2
               ntesf = ntesf + 1
               if ( ntesf>=mtesf ) return
            else
               ntesf = 0
            endif
         endif
         if ( nit>=mit ) then
            iterm = 11
            return
         endif
         if ( nfv>=mfv ) then
            iterm = 12
            return
         endif
         if ( nfg>=mfg ) then
            iterm = 13
            return
         endif
      endif
      iterm = 0
      if ( n>0 .and. nit-kit>=ires1*n+ires2 ) irest = max(irest,1)
      nit = nit + 1
      end subroutine pyfut1

! subroutine pyrmb1               all systems                98/12/01
! 98/12/01 lu : original version
!
! purpose :
! old linear constraint or an old simple bound is removed from the
! active set. transformed gradient of the objective function and
! transformed hessian matrix approximation are updated.
!
! parameters :
!  ii  nf  declared number of variables.
!  ii  n  actual number of variables.
!  iu  ix(nf)  vector containing types of bounds.
!  iu  ic(nc)  vector containing types of constraints.
!  iu  ica(nf)  vector containing indices of active constraints.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ru  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ru  cz(nf*nf)  matrix whose columns are basic vectors from the
!         current reduced subspace.
!  ri  g(nf)  gradient of the objective function.
!  ru  gn(nf)  transformed gradient of the objective function.
!  ru  h(nf*(nf+1)/2)  transformed hessian matrix approximation.
!  ri  eps8  tolerance for constraint to be removed.
!  ri  umax  maximun absolute value of the negative lagrange multiplier.
!  ri  gmax  norm of the transformed gradient.
!  ii  kbf  specification of simple bounds. kbf=0-no simple bounds.
!         kbf=1-one sided simple bounds. kbf=2=two sided simple bounds.
!  ii  kbc  specification of linear constraints. kbc=0-no linear
!         constraints. kbc=1-one sided linear constraints. kbc=2=two
!         sided linear constraints.
!  ii  iold  index of the removed constraint.
!  ia  kold  auxiliary variable.
!  ia  krem  auxiliary variable.
!  io  ier  error indicator.
!  io  iterm  termination indicator.
!
! common data :
!  iu  nrem  number of constraint deletions.
!
! subprograms used :
!  s   plrmb0  constraint deletion.
!  s   mxvset  initiation of a vector.
!
      subroutine pyrmb1(nf,n,ix,ic,ica,cg,cr,cz,g,gn,h,eps8,umax,gmax,  &
                        kbf,kbc,iold,kold,krem,ier,iterm)
      implicit none
      integer nf , n , ix(*) , ic(*) , ica(*) , kbf , kbc , iold ,      &
              kold , krem , ier , iterm
      double precision cg(*) , cr(*) , cz(*) , g(*) , gn(*) , h(*) ,    &
                       eps8 , umax , gmax
      integer i , j , k , kc , l
      integer nres , ndec , nrem , nadd , nit , nfv , nfg , nfh
      common /stat  / nres , ndec , nrem , nadd , nit , nfv , nfg , nfh
      if ( kbc>0 ) then
         if ( umax>eps8*gmax ) then
            call plrmb0(nf,n,ica,cg,cr,cz,g,gn,iold,krem,nrem,ier)
            if ( ier<0 ) then
               iterm = -16
            elseif ( ier>0 ) then
               iold = 0
            else
               k = n*(n-1)/2
               call mxvset(n,0.0d0,h(k+1))
               h(k+n) = 1.0d0
               kc = ica(nf-n+1)
               if ( kc>0 ) then
                  ic(kc) = -ic(kc)
               else
                  k = -kc
                  ix(k) = -ix(k)
               endif
            endif
         else
            iold = 0
         endif
      elseif ( kbf>0 ) then
         if ( umax>eps8*gmax ) then
            ix(iold) = min(abs(ix(iold)),3)
            do i = n , kold , -1
               gn(i+1) = gn(i)
            enddo
            gn(kold) = g(iold)
            n = n + 1
            k = n*(n-1)/2
            l = k + n
            do i = n , kold , -1
               do j = i , 1 , -1
                  if ( i/=kold .and. j/=kold ) then
                     h(l) = h(k)
                     k = k - 1
                     l = l - 1
                  elseif ( i==kold .and. j==kold ) then
                     h(l) = 1.0d0
                     l = l - 1
                  else
                     h(l) = 0.0d0
                     l = l - 1
                  endif
               enddo
            enddo
         else
            iold = 0
            kold = 0
         endif
      endif
      end subroutine pyrmb1

! subroutine pytrbd             all systems                   98/12/01
! 98/12/01 lu : original version
!
! purpose :
! vectors of variables difference and gradients difference are computed
! and transformed. test value dmax is determined.
!
! parameters :
!  ii  nf declared number of variables.
!  ii  n  actual number of variables.
!  ri  x(nf)  vector of variables.
!  ru  xo(nf)  vectors of variables difference.
!  ri  g(nf)  gradient of the objective function.
!  ru  go(nf)  gradients difference.
!  ri  cz(nf*nf)  matrix whose columns are basic vectors from current
!         reduced subspace.
!  ru  sn(nf)  transformed direction vector.
!  ri  r  value of the stepsize parameter.
!  ru  f  new value of the objective function.
!  ri  fo  old value of the objective function.
!  ru  p  new value of the directional derivative.
!  ru  po  old value of the directional derivative.
!  ro  dmax  maximum relative difference of variables.
!  ii  iters  termination indicator for steplength determination.
!         iters=0 for zero step.
!  ii  kbf  specification of simple bounds. kbf=0-no simple bounds.
!         kbf=1-one sided simple bounds. kbf=2=two sided simple bounds.
!  ii  kbc  specification of linear constraints. kbc=0-no linear
!         constraints. kbc=1-one sided linear constraints. kbc=2=two
!         sided linear constraints.
!
! subprograms used :
!  s   mxdrmm  premultiplication of a vector by transpose of a dense
!         rectangular matrix.
!  s   mxvcop  copying of a vector.
!  s   mxvdif  difference of two vectors.
!  s   mxvmul  diagonal premultiplication of a vector.
!  s   mxvsav  difference of two vectors with copying and saving the
!         substracted one.
!  s   mxvscl  scaling of a vector.
!
      subroutine pytrbd(nf,n,x,ix,xo,g,go,cz,sn,r,f,fo,p,po,dmax,iters, &
                        kbf,kbc)
      implicit none
      integer nf , n , ix(*) , iters , kbf , kbc
      double precision x(*) , xo(*) , g(*) , go(*) , cz(*) , sn(*) , r ,&
                       f , fo , p , po , dmax
      integer i , k
      if ( iters>0 ) then
         call mxvdif(nf,x,xo,xo)
         call mxvdif(nf,g,go,go)
         po = r*po
         p = r*p
      else
         f = fo
         p = po
         call mxvsav(nf,x,xo)
         call mxvsav(nf,g,go)
      endif
      dmax = 0.0d0
      if ( kbc>0 ) then
         do i = 1 , nf
            dmax = max(dmax,abs(xo(i))/max(abs(x(i)),1.0d0))
         enddo
         if ( n>0 ) then
            call mxvscl(n,r,sn,xo)
            call mxvcop(nf,go,sn)
            call mxdrmm(nf,n,cz,sn,go)
         endif
      elseif ( kbf>0 ) then
         k = 0
         do i = 1 , nf
            if ( ix(i)>=0 ) then
               k = k + 1
               dmax = max(dmax,abs(xo(i))/max(abs(x(i)),1.0d0))
               xo(k) = xo(i)
               go(k) = go(i)
            endif
         enddo
      else
         do i = 1 , nf
            dmax = max(dmax,abs(xo(i))/max(abs(x(i)),1.0d0))
         enddo
      endif
      end subroutine pytrbd

! subroutine pytrbg               all systems                98/12/01
! 98/12/01 lu : original version
!
! purpose :
! gradient of the objective function is scaled and reduced.
! test values gmax and umax are computed.
!
! parameters :
!  ii  nf  declared number of variables.
!  ii  n  actual number of variables.
!  ii  ix(nf)  vector containing types of bounds.
!  ii  ic(nc)  vector containing types of constraints.
!  ii  ica(nf)  vector containing indices of active constraints.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ri  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ru  cz(nf*nf)  matrix whose columns are basic vectors from the
!         current reduced subspace.
!  ri  g(nf)  gradient of the objective function.
!  ro  gn(nf)  transformed gradient of the objective function.
!  ri  eps7  tolerance for linear independence of constraints.
!  ro  umax  maximum absolute value of the negative lagrange multiplier.
!  ro  gmax  norm of the transformed gradient.
!  ii  kbf  specification of simple bounds. kbf=0-no simple bounds.
!         kbf=1-one sided simple bounds. kbf=2=two sided simple bounds.
!  ii  kbc  specification of linear constraints. kbc=0-no linear
!         constraints. kbc=1-one sided linear constraints. kbc=2=two
!         sided linear constraints.
!  ii  iold  index of the removed constraint.
!  ia  kold  auxiliary variable.
!
! subprograms used :
!  s   mxdrmm  premultiplication of a vector by a rowwise stored dense
!         rectangular matrix.
!  s   mxdprb  back substitution.
!  s   mxvcop  copying of a vector.
!  rf  mxvdot  dot product of two vectors.
!  rf  mxvmax  l-infinity norm of a vector.
!  s   mxvmul  diagonal premultiplication of a vector.
!
      subroutine pytrbg(nf,n,ix,ic,ica,cg,cr,cz,g,gn,umax,gmax,kbf,kbc, &
                        iold,kold)
      implicit none
      integer nf , n , ix(*) , ic(*) , ica(*) , kbf , kbc , iold , kold
      double precision cg(*) , cr(*) , cz(*) , g(*) , gn(*) , umax ,    &
                       gmax
      double precision temp !, mxvmax , mxvdot
      integer nca , ncz , i , j , k , kc
      iold = 0
      kold = 0
      umax = 0.0d0
      gmax = 0.0d0
      if ( kbc>0 ) then
         if ( nf>n ) then
            nca = nf - n
            ncz = n*nf
            call mxvcop(nf,g,gn)
            do j = 1 , nca
               k = ica(j)
               if ( k>0 ) then
                  cz(ncz+j) = mxvdot(nf,cg((k-1)*nf+1),gn)
               else
                  i = -k
                  cz(ncz+j) = gn(i)
               endif
            enddo
            call mxdprb(nca,cr,cz(ncz+1),0)
            do j = 1 , nca
               temp = cz(ncz+j)
               kc = ica(j)
               if ( kc>0 ) then
                  k = ic(kc)
               else
                  i = -kc
                  k = ix(i)
               endif
               if ( k<=-5 ) then
               elseif ( (k==-1 .or. k==-3) .and. umax+temp>=0.0d0 ) then
               elseif ( .not.((k==-2 .or. k==-4) .and. umax-temp>=0.0d0)&
                        ) then
                  iold = j
                  umax = abs(temp)
               endif
            enddo
         endif
         if ( n>0 ) then
            call mxdrmm(nf,n,cz,g,gn)
            gmax = mxvmax(n,gn)
         endif
      elseif ( kbf>0 ) then
         j = 0
         iold = 0
         kold = 0
         do i = 1 , nf
            temp = g(i)
            k = ix(i)
            if ( k>=0 ) then
               j = j + 1
               gn(j) = temp
               gmax = max(gmax,abs(temp))
            elseif ( k<=-5 ) then
            elseif ( (k==-1 .or. k==-3) .and. umax+temp>=0.0d0 ) then
            elseif ( .not.((k==-2 .or. k==-4) .and. umax-temp>=0.0d0) ) &
                     then
               iold = i
               kold = j + 1
               umax = abs(temp)
            endif
         enddo
         n = j
      else
         do i = 1 , nf
            temp = g(i)
            gmax = max(gmax,abs(temp))
         enddo
         n = nf
      endif
      end subroutine pytrbg

! subroutine pytrbh               all systems                98/12/01
! 98/12/01 lu : original version
!
! purpose :
! hessian matrix of the objective function or its approximation is
! scaled and reduced.
!
! parameters :
!  ii  nf  declared number of variables.
!  ii  n  actual number of variables.
!  ri  cr(nf*(nf+1)/2)  triangular decomposition of kernel of the
!         orthogonal projection.
!  ri  cz(nf*nf)  matrix whose columns are basic vectors from the
!         current reduced subspace.
!  ri  h(nf*(nf+1)/2)  hessian matrix or its approximation.
!  ra  s(nf)  auxiliary vector.
!  ii  kbf  specification of simple bounds. kbf=0-no simple bounds.
!         kbf=1-one sided simple bounds. kbf=2=two sided simple bounds.
!  ii  kbc  specification of linear constraints. kbc=0-no linear
!         constraints. kbc=1-one sided linear constraints. kbc=2=two
!         sided linear constraints.
!  ii  ld  degree of previously computed derivatives.
!  ii  iters  termination indicator for steplength determination.
!
! subprograms used :
!  s   mxdsmm  matrix vector product.
!  s   mxvcop  copying of a vector.
!  rf  mxvdot  dot product of two vectors.
!
      subroutine pytrbh(nf,n,ix,cr,cz,h,s,kbf,kbc,ld,iters)
      implicit none
      integer nf , n , ix(*) , kbf , kbc , ld , iters
      double precision cr(*) , cz(*) , h(*) , s(*)
      !double precision mxvdot
      integer nca , ncr , icz , jcz , i , j , k , l
      if ( ld/=2 .or. iters==0 ) return
      if ( kbc>0 ) then
         if ( n<=0 ) return
         nca = nf - n
         ncr = nca*(nca+1)/2
         k = ncr
         jcz = 1
         do j = 1 , n
            call mxdsmm(nf,h,cz(jcz),s)
            icz = 1
            do i = 1 , j
               k = k + 1
               cr(k) = mxvdot(nf,cz(icz),s)
               icz = icz + nf
            enddo
            jcz = jcz + nf
         enddo
         call mxvcop(n*(n+1)/2,cr(ncr+1),h)
      elseif ( kbf>0 ) then
         k = 0
         l = 0
         do i = 1 , nf
            do j = 1 , i
               k = k + 1
               if ( ix(i)>=0 .and. ix(j)>=0 ) then
                  l = l + 1
                  h(l) = h(k)
               endif
            enddo
         enddo
      endif
      end subroutine pytrbh

! subroutine pytrbs               all systems                98/12/01
! 98/12/01 lu : original version
!
! purpose :
! scaled and reduced direction vector is back transformed.
! vectors x,g and values f,p are saved.
!
! parameters :
!  ii  nf  declared number of variables.
!  iu  n  actual number of variables.
!  ii  nc  number of linearized constraints.
!  ri  x(nf)  vector of variables.
!  ii  ix(nf)  vector containing types of bounds.
!  ro  xo(nf)  saved vector of variables.
!  ri  xl(nf)  vector containing lower bounds for variables.
!  ri  xu(nf)  vector containing upper bounds for variables.
!  ri  g(nf)  gradient of the objective function.
!  ro  go(nf)  saved gradient of the objective function.
!  ri  cf(nf)  vector containing values of the constraint funcyions.
!  ro  cfd(nf)  vector containing increments of the constraint
!         functions.
!  ii  ic(nc)  vector containing types of constraints.
!  ri  cl(nc)  vector containing lower bounds for constraint functions.
!  ri  cu(nc)  vector containing upper bounds for constraint functions.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ri  cz(nf*nf)  matrix whose columns are basic vectors from the
!         current reduced subspace.
!  ri  sn(nf)  transformed direction vector.
!  ro  s(nf)  direction vector.
!  ro  ro  saved value of the stepsize parameter.
!  ro  fp  previous value of the objective function.
!  ru  fo  saved value of the objective function.
!  ri  f  value of the objective function.
!  ro  po  saved value of the directional derivative.
!  ri  p  value of the directional derivative.
!  ru  rmax  maximum value of the stepsize parameter.
!  ii  kbf  specification of simple bounds. kbf=0-no simple bounds.
!         kbf=1-one sided simple bounds. kbf=2=two sided simple bounds.
!  ii  kbc  specification of linear constraints. kbc=0-no linear
!         constraints. kbc=1-one sided linear constraints. kbc=2=two
!         sided linear constraints.
!  io  krem  indication of linearly dependent gradients.
!  io  inew  index of the new active function.
!
! subprograms used :
!  s   plmaxs  determination of the maximum stepsize using simple
!         bounds.
!  s   plmaxl  determination of the maximum stepsize using linear
!         constraints.
!  s   mxdcmm  matrix vector product.
!  s   mxvcop  copying of a vector.
!  s   mxvset  initiation of a vector.
!
      subroutine pytrbs(nf,n,nc,x,ix,xo,xl,xu,g,go,cf,cfd,ic,cl,cu,cg,  &
                        cz,sn,s,ro,fp,fo,f,po,p,rmax,kbf,kbc,krem,inew)
      implicit none
      integer nf , n , nc , ix(*) , ic(*) , kbf , kbc , krem , inew
      double precision x(*) , xo(*) , xl(*) , xu(*) , g(*) , go(*) ,    &
                       cf(*) , cfd(*) , cl(*) , cu(*) , cg(*) , cz(*) , &
                       sn(*) , s(*) , ro , fp , fo , f , po , p , rmax
      integer i , k
      fp = fo
      ro = 0.0d0
      fo = f
      po = p
      call mxvcop(nf,x,xo)
      call mxvcop(nf,g,go)
      if ( kbc>0 ) then
         if ( n>0 ) then
            call mxdcmm(nf,n,cz,sn,s)
            inew = 0
            call plmaxl(nf,nc,cf,cfd,ic,cl,cu,cg,s,rmax,kbc,krem,inew)
            call plmaxs(nf,x,ix,xl,xu,s,rmax,kbf,krem,inew)
         else
            call mxvset(nf,0.0d0,s)
         endif
      elseif ( kbf>0 ) then
         k = n + 1
         do i = nf , 1 , -1
            if ( ix(i)<0 ) then
               s(i) = 0.0d0
            else
               k = k - 1
               s(i) = sn(k)
            endif
         enddo
         inew = 0
         call plmaxs(nf,x,ix,xl,xu,s,rmax,kbf,krem,inew)
      endif
      end subroutine pytrbs

! subroutine pytrfd             all systems                   90/12/01
! 90/12/01 lu : original version
!
! purpose :
! preparation of variable metric update.
!
! parameters :
!  ii  nf  declared number of variables.
!  ii  nc  number of constraints.
!  ri  x(nf)  vector of variables.
!  ru  xo(nf)  saved vector of variables.
!  ii  iaa(nf+1)  vector containing indices of active functions.
!  ri  ag(nf*na)  matrix whose columns are gradients of the linear
!          approximated functions.
!  ri  az(nf+1)  vector of lagrange multipliers.
!  ri  cg(nf*nc)  matrix whose columns are normals of the linear
!         constraints.
!  ri  g(nf)  gradient of the lagrangian function.
!  ru  go(nf)  saved gradient of the lagrangian function.
!  ii  n  actual number of variables.
!  ii  kd  degree of required dervatives.
!  iu  ld  degree of previously computed derivatives.
!  ru  r  value of the stepsize parameter.
!  ru  f  value of the objective function.
!  ri  fo  saved value of the objective function.
!  ru  p  value of the directional derivative.
!  ru  po  saved value of the directional derivative.
!  ro  dmax  relative stepsize.
!  io  iters  termination indicator. iters=0-zero step. iters=1-perfect
!         line search. iters=2 goldstein stepsize. iters=3-curry
!         stepsize. iters=4-extended curry stepsize.
!         iters=5-armijo stepsize. iters=6-first stepsize.
!         iters=7-maximum stepsize. iters=8-unbounded function.
!         iters=-1-mred reached. iters=-2-positive directional
!         derivative. iters=-3-error in interpolation.
!
! subprograms used :
!  s   mxvcop  copying of a vector.
!  s   mxvdif  difference of two vectors.
!  s   mxvdir  vector augmented by the scaled vector.
!  s   mxvset  initiation of a vector.
!  s   mxvsav  difference of two vectors with copying and saving the
!         substracted one.
!
      subroutine pytrfd(nf,nc,x,xo,iaa,ag,az,cg,g,go,n,kd,ld,r,f,fo,p,  &
                        po,dmax,iters)
      implicit none
      double precision dmax , f , fo , p , po , r
      integer iters , kd , ld , n , nc , nf
      double precision ag(*) , az(*) , cg(*) , g(*) , go(*) , x(*) ,    &
                       xo(*)
      integer iaa(*)
      integer i , j , l
      call mxvset(nf,0.0d0,g)
      do j = 1 , nf - n
         l = iaa(j)
         if ( l>nc ) then
            l = l - nc
            call mxvdir(nf,-az(j),ag((l-1)*nf+1),g,g)
         elseif ( l>0 ) then
            call mxvdir(nf,-az(j),cg((l-1)*nf+1),g,g)
         else
            l = -l
            g(l) = g(l) - az(j)
         endif
      enddo
      if ( iters>0 ) then
         call mxvdif(nf,x,xo,xo)
         call mxvdif(nf,g,go,go)
         po = r*po
         p = r*p
      else
         r = 0.0d0
         f = fo
         p = po
         call mxvsav(nf,x,xo)
         call mxvsav(nf,g,go)
         ld = kd
      endif
      dmax = 0.0d0
      do i = 1 , nf
         dmax = max(dmax,abs(xo(i))/max(abs(x(i)),1.0d0))
      enddo
      n = nf
      end subroutine pytrfd

! subroutine pytrnd             all systems                   91/12/01
! 91/12/01 lu : original version
!
! purpose :
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
!
! subprograms used :
!  s   mxvcop  copying of a vector.
!  s   mxvdir  vector augmented by the scaled vector.
!  s   mxvset  initiation of a vector.
!  s   mxvsav  difference of two vectors with copying and saving the
!         substracted one.
!
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

! subroutine pytrud             all systems                   98/12/01
! 98/12/01 lu : original version
!
! purpose :
! vectors of variables difference and gradients difference are computed
! and scaled. test value dmax is determined.
!
! parameters :
!  ii  nf declared number of variables.
!  ri  x(nf)  vector of variables.
!  ru  xo(nf)  vectors of variables difference.
!  ri  g(nf)  gradient of the objective function.
!  ru  go(nf)  gradients difference.
!  ro  r  value of the stepsize parameter.
!  ro  f  new value of the objective function.
!  ri  fo  old value of the objective function.
!  ro  p  new value of the directional derivative.
!  ri  po  old value of the directional derivative.
!  ro  dmax  maximum relative difference of variables.
!  ii  kd  degree of required dervatives.
!  io  ld  degree of previously computed derivatives.
!  ii  iters  termination indicator for steplength determination.
!         iters=0 for zero step.
!
! subprograms used :
!  s   pyset1  degree definition of the computed derivatives.
!  s   mxvdif  difference of two vectors.
!  s   mxvsav  difference of two vectors with copying and saving the
!         substracted one.
!
      subroutine pytrud(nf,x,xo,g,go,r,f,fo,p,po,dmax,kd,ld,iters)
      implicit none
      integer nf , kd , ld , iters
      double precision x(*) , xo(*) , g(*) , go(*) , r , f , fo , p ,   &
                       po , dmax
      integer i
      if ( iters>0 ) then
         call mxvdif(nf,x,xo,xo)
         call mxvdif(nf,g,go,go)
         po = r*po
         p = r*p
      else
         f = fo
         p = po
         call mxvsav(nf,x,xo)
         call mxvsav(nf,g,go)
         ld = kd
      endif
      dmax = 0.0d0
      do i = 1 , nf
         dmax = max(dmax,abs(xo(i))/max(abs(x(i)),1.0d0))
      enddo
      end subroutine pytrud

! subroutine pytruf             all systems                   98/12/01
! 98/12/01 lu : original version
!
! purpose :
! vectors of variables difference and right hand sides difference are
! computed and scaled. test value dmax is determined.
!
! parameters :
!  ii  nf declared number of variables.
!  ii  na number of approximated functions.
!  ri  x(nf)  vector of variables.
!  ru  xo(nf)  vectors of variables difference.
!  ri  af(na)  vector of right hand sides.
!  ri  afo(na)  vector of right hand sides difference.
!  ro  r  value of the stepsize parameter.
!  ro  f  new value of the objective function.
!  ri  fo  old value of the objective function.
!  ro  p  new value of the directional derivative.
!  ri  po  old value of the directional derivative.
!  ro  dmax  maximum relative difference of variables.
!  ii  kd  degree of required dervatives.
!  io  ld  degree of previously computed derivatives.
!  ii  iters  termination indicator for steplength determination.
!         iters=0 for zero step.
!
! subprograms used :
!  s   pyset1  degree definition of the computed derivatives.
!  s   mxvdif  difference of two vectors.
!  s   mxvsav  difference of two vectors with copying and saving the
!         substracted one.
!
      subroutine pytruf(nf,na,x,xo,af,afo,r,f,fo,p,po,dmax,kd,ld,iters)
      implicit none
      integer nf , na , kd , ld , iters
      double precision x(*) , xo(*) , af(*) , afo(*) , r , f , fo , p , &
                       po , dmax
      integer i
      if ( iters>0 ) then
         call mxvdif(nf,x,xo,xo)
         call mxvdif(na,af,afo,afo)
         po = r*po
         p = r*p
      else
         r = 0.0d0
         f = fo
         p = po
         call mxvsav(nf,x,xo)
         call mxvsav(na,af,afo)
         ld = kd
      endif
      dmax = 0.0d0
      do i = 1 , nf
         dmax = max(dmax,abs(xo(i))/max(abs(x(i)),1.0d0))
      enddo
      end subroutine pytruf





!***********************************************************************
    end module psqp_module
!***********************************************************************
