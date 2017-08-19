!***********************************************************************
!     program tsqpn
!
!     test program for the subroutine pminu

    program main

      use psqp_module, only: psqp_class

      type(psqp_class) :: solver

      integer          :: next,nfv,nit
      double precision :: f,fmin,cmax,gmax
      integer          :: i,ierr,iprnt,iterm,itime,nb,nc,nf
      double precision :: cf(200),cl(200),cu(200),rpar(5),x(40),xl(40),xu(40)
      integer          :: ic(200),ix(40),ipar(6)
      integer          :: niter,nfval,nsucc

      niter=0
      nfval=0
      nsucc=0
      call tytim1(itime)

      ! loop for 34 test problems

      do next = 1,34

          ! choice of integer and real parameters

          do i = 1,6
              ipar(i) = 0
          end do
          do i = 1,5
              rpar(i) = 0.0d0
          end do
          ipar(5)=2
          iprnt = 1

          ! problem dimension

          nf = 20
          nb = 20
          nc = 30

          ! initiation of x and choice of rpar(1)

          call tind07(nf,nc,x,ix,xl,xu,ic,cl,cu,fmin,rpar(1),next,ierr)
          if (ierr==0) then
            rpar(1)=0.0d0

            ! solution

            call solver%psqpn(nf,nb,nc,x,ix,xl,xu,cf,ic,cl,cu,&
                              ipar,rpar,f,cmax,gmax,iprnt,iterm,&
                              obj,dobj,con,dcon)
            niter=niter+solver%nit
            nfval=nfval+solver%nfv
            if (iterm>0 .and. iterm<9) nsucc=nsucc+1
          end if
      end do

      write(*,'(A,I5,3X,A,I5,3X,A,I5)') &
            ' niter =',niter,' nfval =',nfval,' nsucc =',nsucc
      call tytim2(itime)

    contains

      subroutine obj(me,nf,x,ff)

      !! user supplied subroutine (calculation of ff)

      class(psqp_class),intent(inout) :: me

      double precision ff
      integer          nf
      double precision x(nf)
      external         tffu07
      call tffu07(nf,x,ff,next)
      end subroutine obj

      subroutine dobj(me,nf,x,gf)

      !! user supplied subroutine (calculation of gf)

      class(psqp_class),intent(inout) :: me

      integer          nf
      double precision gf(nf),x(nf)
      external         tfgu07
      call tfgu07(nf,x,gf,next)
      end subroutine dobj

      subroutine con(me,nf,kc,x,fc)

      !! user supplied subroutine (calculation of fc)

      class(psqp_class),intent(inout) :: me

      double precision fc
      integer          kc,nf
      double precision x(nf)
      external         tcfu07
      call tcfu07(nf,kc,x,fc,next)
      end subroutine con

      subroutine dcon(me,nf,kc,x,gc)

      !! user supplied subroutine (calculation of gc)

      class(psqp_class),intent(inout) :: me

      integer          kc,nf
      double precision gc(nf),x(nf)
      external         tcgu07
      call tcgu07(nf,kc,x,gc,next)
      end subroutine dcon

      ! ! empty subroutines

      ! subroutine fun(nf,ka,x,fa)
      ! integer ka,nf
      ! double precision fa,x(*)
      ! ka=nf
      ! fa=x(1)
      ! end subroutine fun

      ! subroutine dfun(nf,ka,x,ga)
      ! integer ka,nf
      ! double precision ga(*),x(*)
      ! ka=nf
      ! ga(1)=x(1)
      ! end subroutine dfun

    end program main


