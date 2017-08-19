!***********************************************************************
!>
!  Matrix routines.
!
!### History
! * Original version: LU, 1991

   module matrix_routines

      use kind_module, only: wp

      implicit none

      public

   contains
!***********************************************************************

!***********************************************************************
!> date: 91/12/01
!
! multiplication of a columnwise stored dense rectangular matrix a
! by a vector x.

      subroutine mxdcmm(n,m,a,x,y)

      implicit none

      integer :: n !! number of rows of the matrix a.
      integer :: m !! number of columns of the matrix a.
      real(wp) :: a(*) !! a(n*m) rectangular matrix stored columnwise
                       !! in the one-dimensional array.
      real(wp) :: x(*) !! x(m)  input vector.
      real(wp) :: y(*) !! y(n)  output vector equal to a*x.

      integer :: j , k

      call mxvset(n,0.0_wp,y)
      k = 0
      do j = 1 , m
         call mxvdir(n,x(j),a(k+1),y,y)
         k = k + n
      end do
      end subroutine mxdcmm

!***********************************************************************
!> date: 91/12/01
!
! solution of a system of linear equations with a dense symmetric
! positive definite matrix a+e using the factorization a+e=l*d*trans(l)
! obtained by the subroutine mxdpgf.
!
!### Method
! back substitution

      subroutine mxdpgb(n,a,x,job)

      implicit none

      integer :: job  !! option:
                      !!
                      !! * if job=0 then x:=(a+e)**(-1)*x.
                      !! * if job>0 then x:=l**(-1)*x.
                      !! * if job<0 then x:=trans(l)**(-1)*x.
      integer :: n  !! order of the matrix a.
      real(wp) :: a(*)  !! a(n*(n+1)/2) factorization a+e=l*d*trans(l)
                                !! obtained by the subroutine mxdpgf.
      real(wp) :: x(*)  !! x(n)  on input the right hand side of a
                                !! system of linear equations. on output the
                                !! solution of a system of linear equations.

      integer :: i , ii , ij , j
      if ( job>=0 ) then
!
!     phase 1 : x:=l**(-1)*x
!
         ij = 0
         do i = 1 , n
            do j = 1 , i - 1
               ij = ij + 1
               x(i) = x(i) - a(ij)*x(j)
            end do
            ij = ij + 1
         end do
      endif
      if ( job==0 ) then
!
!     phase 2 : x:=d**(-1)*x
!
         ii = 0
         do i = 1 , n
            ii = ii + i
            x(i) = x(i)/a(ii)
         end do
      endif
      if ( job<=0 ) then
!
!     phase 3 : x:=trans(l)**(-1)*x
!
         ii = n*(n-1)/2
         do i = n - 1 , 1 , -1
            ij = ii
            do j = i + 1 , n
               ij = ij + j - 1
               x(i) = x(i) - a(ij)*x(j)
            end do
            ii = ii - i
         end do
      endif
      end subroutine mxdpgb

!***********************************************************************
!> date: 91/12/01
!
! computation of a direction of negative curvature with respect to a
! dense symmetric matrix a using the factorization a+e=l*d*trans(l)
!         obtained by the subroutine mxdpgf.
!
!### Method
! p.e.gill, w.murray : newton type methods for unconstrained and
! linearly constrained optimization, math. programming 28 (1974)
! pp. 311-350.

      subroutine mxdpgd(n,a,x,inf)

      implicit none

      integer :: n  !! order of the matrix a.
      integer :: inf  !! information obtained in the factorization process. the
                      !! direction of negative curvature exists only if inf>0.
      real(wp) :: a(n*(n+1)/2)  !! factorization a+e=l*d*trans(l) obtained by the
                                !! subroutine mxdpgf.
      real(wp) :: x(n)  !! computed direction of negative curvature (i.e.
                        !! trans(x)*a*x<0) if it exists.

      integer :: i , j , ii , ij

      real(wp),parameter :: zero = 0.0_wp
      real(wp),parameter :: one = 1.0_wp

!     right hand side formation

      do i = 1 , n
         x(i) = zero
      end do
      if ( inf<=0 ) return
      x(inf) = one

!     back substitution

      ii = inf*(inf-1)/2
      do i = inf - 1 , 1 , -1
         ij = ii
         do j = i + 1 , inf
            ij = ij + j - 1
            x(i) = x(i) - a(ij)*x(j)
         end do
         ii = ii - i
      end do

      end subroutine mxdpgd

!***********************************************************************
!> date: 89/12/01
!
! factorization a+e=l*d*trans(l) of a dense symmetric positive definite
! matrix a+e where d and e are diagonal positive definite matrices and
! l is a lower triangular matrix. if a is sufficiently positive
! definite then e=0.
!
!### Method
!  * p.e.gill, w.murray : newton type methods for unconstrained and
!    linearly constrained optimization, math. programming 28 (1974)
!    pp. 311-350.

      subroutine mxdpgf(n,a,inf,alf,tau)

      implicit none

      real(wp) :: alf  !! on input a desired tolerance for positive definiteness. on
                       !! output the most negative diagonal element used in the
                       !! factorization process (if inf>0).
      real(wp) :: tau  !! maximum diagonal element of the matrix e.
      integer :: inf  !! an information obtained in the factorization process. if
                      !! inf=0 then a is sufficiently positive definite and e=0. if
                      !! inf<0 then a is not sufficiently positive definite and e>0. if
                      !! inf>0 then a is indefinite and inf is an index of the
                      !! most negative diagonal element used in the factorization
                      !! process.
      integer :: n  !! order of the matrix a.
      real(wp) :: a(*)  !! a(n*(n+1)/2)  on input a given dense symmetric (usually positive
                        !! definite) matrix a stored in the packed form. on output the
                        !! computed factorization a+e=l*d*trans(l).

      real(wp) :: bet , del , gam , rho , sig , tol
      integer :: i , ij , ik , j , k , kj , kk , l

      l = 0
      inf = 0
      tol = alf
!
!     estimation of the matrix norm
!
      alf = 0.0_wp
      bet = 0.0_wp
      gam = 0.0_wp
      tau = 0.0_wp
      kk = 0
      do k = 1 , n
         kk = kk + k
         bet = max(bet,abs(a(kk)))
         kj = kk
         do j = k + 1 , n
            kj = kj + j - 1
            gam = max(gam,abs(a(kj)))
         end do
      end do
      bet = max(tol,bet,gam/n)
!      del = tol*bet
      del = tol*max(bet,1.0_wp)
      kk = 0
      do k = 1 , n
         kk = kk + k
!
!     determination of a diagonal correction
!
         sig = a(kk)
         if ( alf>sig ) then
            alf = sig
            l = k
         endif
         gam = 0.0_wp
         kj = kk
         do j = k + 1 , n
            kj = kj + j - 1
            gam = max(gam,abs(a(kj)))
         end do
         gam = gam*gam
         rho = max(abs(sig),gam/bet,del)
         if ( tau<rho-sig ) then
            tau = rho - sig
            inf = -1
         endif
!
!     gaussian elimination
!
         a(kk) = rho
         kj = kk
         do j = k + 1 , n
            kj = kj + j - 1
            gam = a(kj)
            a(kj) = gam/rho
            ik = kk
            ij = kj
            do i = k + 1 , j
               ik = ik + i - 1
               ij = ij + 1
               a(ij) = a(ij) - a(ik)*gam
            end do
         end do
      end do
      if ( l>0 .and. abs(alf)>del ) inf = l
      end subroutine mxdpgf

!***********************************************************************
!> date: 91/12/01
!
! estimation of the minimum eigenvalue and the corresponding eigenvector
! of a dense symmetric positive definite matrix a+e using the
! factorization a+e=l*d*trans(l) obtained by the subroutine mxdpgf.
!
!### Method
!  * a.k.cline, c.b.moler, g.w.stewart, j.h.wilkinson : an estimate
!    for the condition number of a matrix. siam j. numer. anal. 16
!    (1979) 368-373.

      subroutine mxdpgn(n,a,x,alf,job)

      implicit none

      integer :: n  !! order of the matrix a.
      integer :: job  !! option:
                      !!
                      !! * if job=0 then only estimated eigenvalue is computed.
                      !! * if job>0 then both estimated eigenvalue and
                      !!   estimated eigenvector are computed by job iterations.
      real(wp) :: a(n*(n+1)/2)  !! a(n*(n+1)/2) factorization a+e=l*d*trans(l)
                                !! obtained by the subroutine mxdpgf.
      real(wp) :: x(n)  !! estimated eigenvector.
      real(wp) :: alf  !! estimated eigenvalue.

      real(wp) :: xp , xm , fp , fm
      integer :: i , k , ik , kk

      real(wp),parameter :: zero = 0.0_wp
      real(wp),parameter :: one = 1.0_wp
!
!     computation of the vector v with possible maximum norm such
!     that  l*d**(1/2)*v=u  where u has elements +1 or -1
!
      do i = 1 , n
         x(i) = zero
      end do
      kk = 0
      do k = 1 , n
         kk = kk + k
         xp = -x(k) + one
         xm = -x(k) - one
         fp = abs(xp)
         fm = abs(xm)
         ik = kk
         do i = k + 1 , n
            ik = ik + i - 1
            fp = fp + abs(x(i)+a(ik)*xp)
            fm = fm + abs(x(i)+a(ik)*xm)
         end do
         if ( fp>=fm ) then
            x(k) = xp
            ik = kk
            do i = k + 1 , n
               ik = ik + i - 1
               x(i) = x(i) + a(ik)*xp
            end do
         else
            x(k) = xm
            ik = kk
            do i = k + 1 , n
               ik = ik + i - 1
               x(i) = x(i) + a(ik)*xm
            end do
         endif
      end do
!
!     computation of the vector x such that
!     d**(1/2)*trans(l)*x=v
!
      fm = zero
      kk = 0
      do k = 1 , n
         kk = kk + k
         if ( job<=0 ) then
            fp = sqrt(a(kk))
            x(k) = x(k)/fp
            fm = fm + x(k)*x(k)
         else
            x(k) = x(k)/a(kk)
         endif
      end do
      fp = dble(n)
      if ( job<=0 ) then
!
!     first estimation of the minimum eigenvalue by the formula
!     alf=(trans(u)*u)/(trans(v)*v)
!
         alf = fp/fm
         return
      endif
      fm = zero
      kk = n*(n+1)/2
      do k = n , 1 , -1
         ik = kk
         do i = k + 1 , n
            ik = ik + i - 1
            x(k) = x(k) - a(ik)*x(i)
         end do
         fm = fm + x(k)*x(k)
         kk = kk - k
      end do
      fm = sqrt(fm)
      if ( job<=1 ) then
!
!     second estimation of the minimum eigenvalue by the formula
!     alf=sqrt(trans(u)*u)/sqrt(trans(x)*x)
!
         alf = sqrt(fp)/fm
      else
!
!     inverse iterations
!
         do k = 2 , job
!
!     scaling the vector x by its norm
!
            do i = 1 , n
               x(i) = x(i)/fm
            end do
            call mxdpgb(n,a,x,0)
            fm = sqrt(mxvdot(n,x,x))
         end do
         alf = one/fm
      endif
!
!     scaling the vector x by its norm
!
      do i = 1 , n
         x(i) = x(i)/fm
      end do
      end subroutine mxdpgn

!***********************************************************************
!> date: 91/12/01
!
! computation of the number mxdpgp=trans(x)*d**(-1)*y where d is a
! diagonal matrix in the factorization a+e=l*d*trans(l) obtained by the
! subroutine mxdpgf.

      function mxdpgp(n,a,x,y)

      implicit none

      integer :: n  !! order of the matrix a.
      real(wp) :: a(*)  !! a(n*(n+1)/2) factorization a+e=l*d*trans(l)
                        !! obtained by the subroutine mxdpgf.
      real(wp) :: x(*)  !! input vector.
      real(wp) :: y(*)  !! input vector.
      real(wp) :: mxdpgp  !! computed number mxdpgp=trans(x)*d**(-1)*y.

      real(wp) :: temp
      integer :: i , j

      j = 0
      temp = 0.0_wp
      do i = 1 , n
         j = j + i
         temp = temp + x(i)*y(i)/a(j)
      end do
      mxdpgp = temp
      end function mxdpgp

!***********************************************************************
!> date: 91/12/01
!
! scaling of a dense symmetric positive definite matrix a+e using the
! factorization a+e=l*d*trans(l) obtained by the subroutine mxdpgf.

      subroutine mxdpgs(n,a,alf)

      implicit none

      real(wp) :: alf  !! scaling factor.
      integer :: n  !! order of the matrix a.
      real(wp) :: a(*)  !! a(n*(n+1)/2) factorization a+e=l*d*trans(l)
                        !! obtained by the subroutine mxdpgf.

      integer :: i , j

      j = 0
      do i = 1 , n
         j = j + i
         a(j) = a(j)*alf
      end do

      end subroutine mxdpgs

!***********************************************************************
!> date: 89/12/01
!
! correction of a dense symmetric positive definite matrix a+e in the
! factored form a+e=l*d*trans(l) obtained by the subroutine mxdpgf.
! the correction is defined as a+e:=a+e+alf*x*trans(x) where alf is a
! given scaling factor and x is a given vector.
!
!### Method
!  * p.e.gill, w.murray, m.saunders: methods for computing and modifying
!    the ldv factors of a matrix, math. of comp. 29 (1974) pp. 1051-1077.

      subroutine mxdpgu(n,a,alf,x,y)

      implicit none

      real(wp) :: alf  !! scaling factor in the correction term.
      integer :: n  !! order of the matrix a.
      real(wp) :: a(*)  !! a(n*(n+1)/2) factorization a+e=l*d*trans(l)
                        !! obtained by the subroutine mxdpgf.
      real(wp) :: x(*)  !! vector in the correction term.
      real(wp) :: y(*)  !! auxiliary vector.

      real(wp) :: alfr
      real(wp) :: b , d , p , r , t , to
      integer :: i , ii , ij , j

      real(wp), parameter :: zero = 0.0_wp
      real(wp), parameter :: one = 1.0_wp
      real(wp), parameter :: four = 4.0_wp
      real(wp), parameter :: con = 1.0e-8_wp

      if ( alf>=zero ) then
!
!     forward correction in case when the scaling factor is nonnegative
!
         alfr = sqrt(alf)
         call mxvscl(n,alfr,x,y)
         to = one
         ii = 0
         do i = 1 , n
            ii = ii + i
            d = a(ii)
            p = y(i)
            t = to + p*p/d
            r = to/t
            a(ii) = d/r
            b = p/(d*t)
            if ( a(ii)<=four*d ) then
!
!     an easy formula for limited diagonal element
!
               ij = ii
               do j = i + 1 , n
                  ij = ij + j - 1
                  d = a(ij)
                  y(j) = y(j) - p*d
                  a(ij) = d + b*y(j)
               end do
            else
!
!     a more complicate but numerically stable formula for unlimited
!     diagonal element
!
               ij = ii
               do j = i + 1 , n
                  ij = ij + j - 1
                  d = a(ij)
                  a(ij) = r*d + b*y(j)
                  y(j) = y(j) - p*d
               end do
            endif
            to = t
         end do
      else
!
!     backward correction in case when the scaling factor is negative
!
         alfr = sqrt(-alf)
         call mxvscl(n,alfr,x,y)
         to = one
         ij = 0
         do i = 1 , n
            d = y(i)
            do j = 1 , i - 1
               ij = ij + 1
               d = d - a(ij)*y(j)
            end do
            y(i) = d
            ij = ij + 1
            to = to - d*d/a(ij)
         end do
         if ( to<=zero ) to = con
         ii = n*(n+1)/2
         do i = n , 1 , -1
            d = a(ii)
            p = y(i)
            t = to + p*p/d
            a(ii) = d*to/t
            b = -p/(d*to)
            to = t
            ij = ii
            do j = i + 1 , n
               ij = ij + j - 1
               d = a(ij)
               a(ij) = d + b*y(j)
               y(j) = y(j) + p*d
            end do
            ii = ii - i
         end do
      endif
      end subroutine mxdpgu

!***********************************************************************
!> date: 89/12/01
!
! solution of a system of linear equations with a dense symmetric
! positive definite matrix a using the factorization a=trans(r)*r.
!
!### Method
! back substitution

      subroutine mxdprb(n,a,x,job)

      implicit none

      integer :: job    !! option:
                        !!
                        !! * if job=0 then x:=a**(-1)*x.
                        !! * if job>0 then x:=trans(r)**(-1)*x.
                        !! * if job<0 then x:=r**(-1)*x.
      integer :: n  !! order of the matrix a.
      real(wp) :: a(*)  !! a(n*(n+1)/2) factorization a=trans(r)*r.
      real(wp) :: x(*)  !! x(n)  on input the right hand side of a system of linear
                        !! equations. on output the solution of a system of linear
                        !! equations.

      integer :: i , ii , ij , j

      if ( job>=0 ) then
!
!     phase 1 : x:=trans(r)**(-1)*x
!
         ij = 0
         do i = 1 , n
            do j = 1 , i - 1
               ij = ij + 1
               x(i) = x(i) - a(ij)*x(j)
            end do
            ij = ij + 1
            x(i) = x(i)/a(ij)
         end do
      endif
      if ( job<=0 ) then
!
!     phase 2 : x:=r**(-1)*x
!
         ii = n*(n+1)/2
         do i = n , 1 , -1
            ij = ii
            do j = i + 1 , n
               ij = ij + j - 1
               x(i) = x(i) - a(ij)*x(j)
            end do
            x(i) = x(i)/a(ii)
            ii = ii - i
         end do
      endif
      end subroutine mxdprb

!***********************************************************************
!> date: 92/12/01
!
! correction of a singular dense symmetric positive semidefinite matrix
! a decomposed as a=trans(r)*r.

      subroutine mxdprc(n,a,inf,tol)

      implicit none

      integer :: n  !! order of the matrix a.
      integer :: inf  !! an information obtained in the correction process:
                      !!
                      !! * if inf=0 then a is sufficiently positive definite.
                      !! * if inf<0 then a is not sufficiently positive definite.
      real(wp) :: a(n*(n+1)/2)  !! dense symmetric matrix stored in the packed form.
      real(wp) :: tol  !! desired tolerance for positive definiteness.

      real(wp) :: tol1 , temp
      integer :: l , i

      inf = 0
      tol1 = sqrt(tol)
      temp = tol1
      do i = 1 , n*(n+1)/2
         temp = max(temp,abs(a(i)))
      end do
      temp = temp*tol1
      l = 0
      do i = 1 , n
         l = l + i
         if ( abs(a(l))<=temp ) then
            a(l) = sign(temp,a(l))
            inf = -1
         endif
      end do

      end subroutine mxdprc

!***********************************************************************
!> date: 91/12/01
!
! multiplication of a given vector x by a dense symmetric positive
! definite matrix a using the factorization a=trans(r)*r.

      subroutine mxdprm(n,a,x,job)

      implicit none

      integer :: n  !! order of the matrix a.
      integer :: job  !! option:
                      !!
                      !! * if job=0 then x:=a*x.
                      !! * if job>0 then x:=r*x.
                      !! * if job<0 then x:=trans(r)*x.
      real(wp) :: a(n*(n+1)/2)  !! factorization a=trans(r)*r.
      real(wp) :: x(n)  !! on input the given vector.
                        !! on output the result of multiplication.

      integer :: i , j , ii , ij

      if ( job>=0 ) then
!
!     phase 1 : x:=r*x
!
         ii = 0
         do i = 1 , n
            ii = ii + i
            x(i) = a(ii)*x(i)
            ij = ii
            do j = i + 1 , n
               ij = ij + j - 1
               x(i) = x(i) + a(ij)*x(j)
            end do
         end do
      endif
      if ( job<=0 ) then
!
!     phase 2 : x:=trans(r)*x
!
         ij = n*(n+1)/2
         do i = n , 1 , -1
            x(i) = a(ij)*x(i)
            do j = i - 1 , 1 , -1
               ij = ij - 1
               x(i) = x(i) + a(ij)*x(j)
            end do
            ij = ij - 1
         end do
      endif
      end subroutine mxdprm

!***********************************************************************
!> date: 91/12/01
!
! plane rotation is applied to a rowwise stored dense rectangular
! matrix a.

      subroutine mxdrgr(n,a,k,l,ck,cl,ier)

      implicit none

      integer :: n  !! number of columns of the matrix a.
      integer :: k  !! first index of the plane rotation.
      integer :: l  !! second index of the plane rotation.
      integer :: ier  !! type of the plane rotation:
                      !!
                      !! * ier=0-general plane rotation.
                      !! * ier=1-permutation.
                      !! * ier=2-transformation suppressed.
      real(wp) :: a(*)  !! a(m*n)  rectangular matrix stored rowwise in the
                        !! one-dimensional array.
      real(wp) :: ck  !! diagonal element of the elementary orthogonal matrix.
      real(wp) :: cl  !! off-diagonal element of the elementary orthogonal matrix.

      integer :: i , ik , il

      if ( ier/=0 .and. ier/=1 ) return
      ik = (k-1)*n
      il = (l-1)*n
      do i = 1 , n
         ik = ik + 1
         il = il + 1
         call mxvrot(a(ik),a(il),ck,cl,ier)
      end do

      end subroutine mxdrgr

!***********************************************************************
!> date: 91/12/01
!
! multiplication of a rowwise stored dense rectangular matrix a by
! a vector x and addition of a scaled vector alf*y.

      subroutine mxdrmd(n,m,a,x,alf,y,z)

      implicit none

      integer :: n  !! number of columns of the matrix a.
      integer :: m  !! number of rows of the matrix a.
      real(wp) :: a(m*n)  !! rectangular matrix stored rowwise in the one-dimensional array.
      real(wp) :: x(n)  !! input vector.
      real(wp) :: alf  !! scaling factor.
      real(wp) :: y(m)  !! input vector.
      real(wp) :: z(m)  !! output vector equal to a*x+alf*y.

      real(wp) :: temp
      integer :: i , j , k

      k = 0
      do j = 1 , m
         temp = alf*y(j)
         do i = 1 , n
            temp = temp + a(k+i)*x(i)
         end do
         z(j) = temp
         k = k + n
      end do

      end subroutine mxdrmd

!***********************************************************************
!> date: 91/12/01
!
! rowwise stored dense rectangular matrix a is set to be a part of the
! unit matrix.

      subroutine mxdrmi(n,m,a)

      implicit none

      integer :: n  !! number of columns of the matrix a.
      integer :: m  !! number of rows of the matrix a.
      real(wp) :: a(m*n)  !! rectangular matrix stored rowwise in
                          !! the one-dimensional array.
                          !! this matrix is set to trans([i,0]).

      integer :: i , j , k
      real(wp) :: zero , one
      parameter (zero=0.0_wp,one=1.0_wp)
      k = 0
      do j = 1 , m
         do i = 1 , n
            a(i+k) = zero
            if ( i==j ) a(i+k) = one
         end do
         k = k + n
      end do

      end subroutine mxdrmi

!***********************************************************************
!> date: 91/12/01
!
! multiplication of a rowwise stored dense rectangular matrix a by
! a vector x.

      subroutine mxdrmm(n,m,a,x,y)

      implicit none

      integer :: n  !! number of columns of the matrix a.
      integer :: m  !! number of rows of the matrix a.
      real(wp) :: a(*)  !! a(m*n)  rectangular matrix stored rowwise
                        !! in the one-dimensional array.
      real(wp) :: x(*)  !! x(n)  input vector.
      real(wp) :: y(*)  !! y(m)  output vector equal to a*x.

      real(wp) :: temp
      integer :: i , j , k

      k = 0
      do j = 1 , m
         temp = 0.0_wp
         do i = 1 , n
            temp = temp + a(k+i)*x(i)
         end do
         y(j) = temp
         k = k + n
      end do
      end subroutine mxdrmm

!***********************************************************************
!> date: 91/12/01
!
! euclidean norm of a part of the i-th column of a rowwise stored dense
! rectangular matrix a is computed.

      function mxdrmn(n,m,a,i,j)

      implicit none

      integer :: n  !! number of columns of the matrix a.
      integer :: m  !! number of rows of the matrix a.
      real(wp) :: a(m*n)  !! rectangular matrix stored rowwise in the one-dimensional array.
      integer :: i  !! index of the column whose norm is computed.
      integer :: j  !! index of the first element from which the norm is computed.
      real(wp) :: mxdrmn

      real(wp) :: pom , den
      integer :: k , l

      real(wp),parameter :: zero = 0.0_wp

      den = zero
      l = (j-1)*n
      do k = j , m
         den = max(den,abs(a(l+i)))
         l = l + n
      end do
      pom = zero
      if ( den>zero ) then
         l = (j-1)*n
         do k = j , m
            pom = pom + (a(l+i)/den)**2
            l = l + n
         end do
      endif
      mxdrmn = den*sqrt(pom)

      end function mxdrmn

!***********************************************************************
!> date: 91/12/01
!
! k-th column of a rowwise stored dense rectangular matrix a is copied
! to the vector x.

      subroutine mxdrmv(n,m,a,x,k)

      implicit none

      integer :: n  !! number of columns of the matrix a.
      integer :: m  !! number of rows of the matrix a.
      real(wp) :: a(*)  !! a(m*n)  rectangular matrix stored rowwise
                        !! in the one-dimensional array.
      real(wp) :: x(*)  !! x(m)  output vector such that x(j)=a(j,k) for all j.
      integer :: k  !! index of the row being copied to the output vector.

      integer :: i , j

      if ( k<1 .or. k>n ) return
      i = k
      do j = 1 , m
         x(j) = a(i)
         i = i + n
      end do

      end subroutine mxdrmv

!***********************************************************************
!> date: 92/12/01
!
! qr decomposition of rowwise stored dense rectangular matrix q using
! householder transformations without pivoting.
!
!### Method
!  * p.a.bussinger, g.h.golub : linear least squares solution by
!    householder transformation. numer. math. 7 (1965) 269-276.

      subroutine mxdrqf(n,m,q,r)

      implicit none

      integer :: n  !! number of columns of the matrix q.
      integer :: m  !! number of rows of the matrix q.
      real(wp) :: q(m*n)   !! rectangular matrix stored rowwise in the one-dimensional array.
      real(wp) :: r(n*(n+1)/2)  !! upper triangular matrix stored in the packed form.

      real(wp) :: alf , pom
      integer :: i , j , k , l , jp , kp , nm

      real(wp),parameter :: zero=0.0_wp
      real(wp),parameter :: one=1.0_wp

      nm = min(n,m)
!
!     qr decomposition
!
      l = 0
      kp = 0
      do k = 1 , nm
         pom = mxdrmn(n,m,q,k,k)
         if ( q(kp+k)/=zero ) pom = sign(pom,q(kp+k))
         r(l+k) = -pom
         jp = 0
         do j = 1 , k - 1
            r(l+j) = q(jp+k)
            q(jp+k) = zero
            jp = jp + n
         end do
         if ( pom/=zero ) then
!
!     householder transformation
!
            do j = k , m
               q(jp+k) = q(jp+k)/pom
               jp = jp + n
            end do
            q(kp+k) = q(kp+k) + one
            do i = k + 1 , n
               alf = zero
               jp = kp
               do j = k , m
                  alf = alf + q(jp+k)*q(jp+i)
                  jp = jp + n
               end do
               alf = alf/q(kp+k)
               jp = kp
               do j = k , m
                  q(jp+i) = q(jp+i) - alf*q(jp+k)
                  jp = jp + n
               end do
            end do
         endif
         l = l + k
         kp = kp + n
      end do
!
!     explicit formulation of the orthogonal matrix
!
      kp = n*n
      do k = n , 1 , -1
         kp = kp - n
         if ( q(kp+k)/=zero ) then
            do i = k + 1 , n
               alf = zero
               jp = kp
               do j = k , m
                  alf = alf + q(jp+k)*q(jp+i)
                  jp = jp + n
               end do
               alf = alf/q(kp+k)
               jp = kp
               do j = k , m
                  q(jp+i) = q(jp+i) - alf*q(jp+k)
                  jp = jp + n
               end do
            end do
            jp = kp
            do j = k , m
               q(jp+k) = -q(jp+k)
               jp = jp + n
            end do
         endif
         q(kp+k) = q(kp+k) + one
      end do

      end subroutine mxdrqf

!***********************************************************************
!> date: 92/12/01
!
! update of a qr decomposition. this qr decomposition is updated
! by the rule q*r:=q*r+alf*x*trans(y).
!
!### Method
!  * j.w.daniel, w.b.gragg, l.kaufman, g.w.steward : reorthogonalization
!    and stable algorithms for updating the gram-schmidt qr factorization.
!    mathematics of computation 30 (1976) 772-795.

      subroutine mxdrqu(n,m,q,r,alf,x,y,z,inf)

      implicit none

      integer :: n             !! number of columns of the matrix q.
      integer :: m             !! number of rows of the matrix q.
      real(wp) :: q(m*n)       !! rectangular matrix stored rowwise in
                               !! the one-dimensional array (part of the
                               !! orthogonal matrix).
      real(wp) :: r(n*(n+1)/2) !! upper triangular matrix stored in a
                               !! packed form.
      real(wp) :: alf          !! scalar parameter.
      real(wp) :: x(m)         !! input vector.
      real(wp) :: y(n)         !! input vector.
      real(wp) :: z(n)         !! auxiliary vector.
      integer :: inf           !! information. if inf=0 then x lies in
                               !! the column space of q. if inf=1 then x
                               !! does not lie in the column space of q.

      real(wp) :: ck , cl , zk , zl
      integer :: j , k , l , kj , kk , ier

      real(wp),parameter :: one = 1.0_wp
      real(wp),parameter :: con = 1.0e-10_wp

      inf = 0
      kk = n*(n+1)/2
!
!     computation of the vector trans(q)*x
!
      call mxdcmm(n,m,q,x,z)
      if ( m>n ) then
!
!     if x does not lie in the column space of q we have to use
!     a subproblem whose dimension is by one greater (inf=1).
!
         zk = mxvnor(m,x)
         call mxdrmd(n,m,q,z,-one,x,x)
         zl = mxvnor(m,x)
         if ( zl>con*zk ) then
            inf = 1
            call mxvscl(m,-one/zl,x,x)
            call mxvort(z(n),zl,ck,cl,ier)
            if ( ier==0 .or. ier==1 ) then
               call mxvrot(r(kk),zl,ck,cl,ier)
               kj = n
               do j = 1 , m
                  call mxvrot(q(kj),x(j),ck,cl,ier)
                  kj = kj + n
               end do
            endif
         endif
      endif
!
!     application of plane rotations to the vector z so that
!     trans(q1)*z=e1 where q1 is an orthogonal matrix (accumulation of
!     the plane rotations) and e1 is the first column of the unit
!     matrix. at the same time both the upper hessenberg matrix
!     trans(q1)*r and the orthogonal matrix q*q1 are constructed so that
!     q*q1*r1=q*q1*(trans(q1)*r+alf*e1*trans(y)) where r1 is an upper
!     hessenberg matrix.
!
      do l = n , 2 , -1
         k = l - 1
         kk = kk - l
         call mxvort(z(k),z(l),ck,cl,ier)
         if ( ier==0 .or. ier==1 ) then
            call mxvrot(r(kk),z(l),ck,cl,ier)
            kj = kk
            do j = l , n
               kj = kj + j - 1
               call mxvrot(r(kj),r(kj+1),ck,cl,ier)
            end do
            kj = k
            do j = 1 , m
               call mxvrot(q(kj),q(kj+1),ck,cl,ier)
               kj = kj + n
            end do
         endif
      end do
      z(1) = alf*z(1)
      kj = 1
      do j = 1 , n
         r(kj) = r(kj) + z(1)*y(j)
         kj = kj + j
      end do
!
!     application of plane rotations to the upper hessenberg matrix r1
!     given above so that r2=trans(q2)*r1 where q2 is an orthogonal
!     matrix (accumulation of the plane rotations) and r2 is an upper
!     triangular matrix. we obtain the new qr decomposition q*q1*q2*r2.
!
      kk = 1
      do l = 2 , n
         k = l - 1
         call mxvort(r(kk),z(l),ck,cl,ier)
         if ( ier==0 .or. ier==1 ) then
            kj = kk
            do j = l , n
               kj = kj + j - 1
               call mxvrot(r(kj),r(kj+1),ck,cl,ier)
            end do
            kj = k
            do j = 1 , m
               call mxvrot(q(kj),q(kj+1),ck,cl,ier)
               kj = kj + n
            end do
         endif
         kk = kk + l
      end do
!
!     back transformation of the greater subproblem if inf=1.
!
      if ( inf==1 ) then
         call mxvort(r(kk),zl,ck,cl,ier)
         if ( ier==0 .or. ier==1 ) then
            kj = n
            do j = 1 , m
               call mxvrot(q(kj),x(j),ck,cl,ier)
               kj = kj + n
            end do
         endif
      endif
      end subroutine mxdrqu

!***********************************************************************
!> date: 91/12/01
!
! a dense symmetric matrix a is augmented by the scaled unit matrix
! such that a:=a+alf*i (i is the unit matrix of order n).

      subroutine mxdsda(n,a,alf)

      implicit none

      integer :: n      !! order of the matrix a.
      real(wp) :: a(*)  !! a(n*(n+1)/2)  dense symmetric matrix
                        !! stored in the packed form.
      real(wp) :: alf   !! scaling factor.

      integer :: i , j

      j = 0
      do i = 1 , n
         j = j + i
         a(j) = a(j) + alf
      end do
      end subroutine mxdsda

!***********************************************************************
!> date: 91/12/01
!
! determination of the minimum diagonal element of a dense symmetric
! matrix.

      function mxdsdl(n,a,inf)

      implicit none

      integer :: n  !! order of the matrix a.
      integer :: inf  !! index of the minimum diagonal element
                      !! of the matrix a.
      real(wp) :: a(n*(n+1)/2)   !! dense symmetric matrix stored
                                 !! in the packed form.
      real(wp) :: mxdsdl  !! minimum diagonal element of the matrix a.

      real(wp) :: temp
      integer :: i , j

      j = 1
      inf = 1
      temp = a(1)
      do i = 2 , n
         j = j + i
         if ( temp>a(j) ) then
            inf = j
            temp = a(j)
         endif
      end do
      mxdsdl = temp

      end function mxdsdl

!***********************************************************************
!> date: 91/12/01
!
! dense symmetric matrix augmented by the scaled dense symmetric matrix.

      subroutine mxdsma(n,alf,a,b,c)

      implicit none

      integer :: n  !! order of the matrices.
      real(wp) :: alf  !! scaling factor.
      real(wp) :: a(n*(n+1)/2)  !! input matrix.
      real(wp) :: b(n*(n+1)/2)  !! input matrix.
      real(wp) :: c(n*(n+1)/2)  !! output matrix where c:=b+alf*a.

      integer :: i

      do i = 1 , n*(n+1)/2
         c(i) = b(i) + alf*a(i)
      end do

      end subroutine mxdsma

!***********************************************************************
!> date: 91/12/01
!
! copying of a dense symmetric matrix.

      subroutine mxdsmc(n,a,b)

      implicit none

      integer :: n  !! order of the matrices a and b.
      real(wp) :: a(n*(n+1)/2)  !! dense symmetric matrix stored in the packed form.
      real(wp) :: b(n*(n+1)/2)  !! dense symmetric matrix stored in the packed form
                                !! where b:=a.

      integer :: m , i

      m = n*(n+1)/2
      do i = 1 , m
         b(i) = a(i)
      end do

      end subroutine mxdsmc

!***********************************************************************
!> date: 91/12/01
!
!  gershgorin bounds for eigenvalues of a dense symmetric matrix.
!  amin <= any eigenvalue of a <= amax.

      subroutine mxdsmg(n,a,amin,amax)

      implicit none

      integer :: n  !! dimension of the matrix a.
      real(wp) :: a(n*(n+1)/2)  !! dense symmetric matrix stored
                                !! in the packed form.
      real(wp) :: amin  !! lower bound for eigenvalues of a.
      real(wp) :: amax  !! upper bound for eigenvalues of a.

      real(wp) :: temp
      integer :: i , j , k , l

      real(wp),parameter :: zero = 0.0_wp

      amax = a(1)
      amin = a(1)
      k = 0
      do i = 1 , n
         temp = zero
         l = k
         do j = 1 , i - 1
            l = l + 1
            temp = temp + abs(a(l))
         end do
         l = l + 1
         do j = i + 1 , n
            l = l + j - 1
            temp = temp + abs(a(l))
         end do
         k = k + i
         amax = max(amax,a(k)+temp)
         amin = min(amin,a(k)-temp)
      end do

      end subroutine mxdsmg

!***********************************************************************
!> date: 88/12/01
!
! dense symmetric matrix a is set to the unit matrix with the same
! order.

      subroutine mxdsmi(n,a)

      implicit none

      integer :: n  !! order of the matrix a.
      real(wp) :: a(*)  !! a(n*(n+1)/2)  dense symmetric matrix
                        !! stored in the packed form which is set
                        !! to the unit matrix (i.e. a:=i).

      integer :: i , m

      m = n*(n+1)/2
      do i = 1 , m
         a(i) = 0.0_wp
      end do
      m = 0
      do i = 1 , n
         m = m + i
         a(m) = 1.0_wp
      end do

      end subroutine mxdsmi

!***********************************************************************
!> date: 89/12/01
!
! multiplication of a dense symmetric matrix a by a vector x.

      subroutine mxdsmm(n,a,x,y)

      implicit none

      integer :: n      !! order of the matrix a.
      real(wp) :: a(*)  !! a(n*(n+1)/2)  dense symmetric matrix stored in the packed form.
      real(wp) :: x(*)  !! x(n)  input vector.
      real(wp) :: y(*)  !! y(n)  output vector equal to  a*x.

      real(wp) :: temp
      integer :: i , j , k , l

      k = 0
      do i = 1 , n
         temp = 0.0_wp
         l = k
         do j = 1 , i
            l = l + 1
            temp = temp + a(l)*x(j)
         end do
         do j = i + 1 , n
            l = l + j - 1
            temp = temp + a(l)*x(j)
         end do
         y(i) = temp
         k = k + i
      end do

      end subroutine mxdsmm

!***********************************************************************
!> date: 91/12/01
!
! value of a quadratic form with a dense symmetric matrix a.

      function mxdsmq(n,a,x,y)

      implicit none

      integer :: n  !! order of the matrix a.
      real(wp) :: a(*)  !! a(n*(n+1)/2)  dense symmetric matrix stored in the packed form.
      real(wp) :: x(*)  !! x(n)  given vector.
      real(wp) :: y(*)  !! y(n)  given vector.
      real(wp) :: mxdsmq  !! value of the quadratic form mxdsmq=trans(x)*a*y.

      real(wp) :: temp , temp1 , temp2
      integer :: i , j , k

      real(wp), parameter :: zero = 0.0_wp

      temp = zero
      k = 0
      do i = 1 , n
         temp1 = zero
         temp2 = zero
         do j = 1 , i - 1
            k = k + 1
            temp1 = temp1 + a(k)*x(j)
            temp2 = temp2 + a(k)*y(j)
         end do
         k = k + 1
         temp = temp + x(i)*(temp2+a(k)*y(i)) + y(i)*temp1
      end do
      mxdsmq = temp

      end function mxdsmq

!***********************************************************************
!> date: 92/12/01
!
! plane rotation is applied to a dense symmetric matrix a. the case
! k=l+1 is required.

      subroutine mxdsmr(n,a,k,l,ck,cl,ier)

      implicit none

      integer :: n  !! order of the matrix a.
      integer :: k  !! first index of plane rotation.
      integer :: l  !! second index of plane rotation.
      integer :: ier  !! information on the transformation:
                      !!
                      !! * ier<0-k or l out of range.
                      !! * ier=0-plane rotation.
                      !! * ier=1-permutation.
                      !! * ier=2-transformation suppressed.
      real(wp) :: a(*)  !! a(n*(n+1)/2) dense symmetric matrix stored in the packed form.
      real(wp) :: ck  !! diagonal element of the elementary orthogonal matrix.
      real(wp) :: cl  !! off-diagonal element of the elementary orthogonal matrix.

      real(wp) :: akk , akl , all , ckk , ckl , cll
      integer :: j , kj , lj , kk , kl , ll

      if ( ier/=0 .and. ier/=1 ) return
      if ( k/=l+1 ) then
         ier = -1
         return
      endif
      lj = l*(l-1)/2
      do j = 1 , n
         if ( j<=l ) then
            lj = lj + 1
            kj = lj + l
         else
            lj = lj + j - 1
            kj = lj + 1
         endif
         if ( j/=k .and. j/=l ) call mxvrot(a(kj),a(lj),ck,cl,ier)
      end do
      if ( ier==0 ) then
         ckk = ck**2
         ckl = ck*cl
         cll = cl**2
         ll = l*(l+1)/2
         kl = ll + l
         kk = ll + k
         akl = (ckl+ckl)*a(kl)
         akk = ckk*a(kk) + cll*a(ll) + akl
         all = cll*a(kk) + ckk*a(ll) - akl
         a(kl) = (cll-ckk)*a(kl) + ckl*(a(kk)-a(ll))
         a(kk) = akk
         a(ll) = all
      else
         ll = l*(l+1)/2
         kk = ll + k
         akk = a(kk)
         a(kk) = a(ll)
         a(ll) = akk
      endif

      end subroutine mxdsmr

!***********************************************************************
!> date: 91/12/01
!
! scaling of a dense symmetric matrix.

      subroutine mxdsms(n,a,alf)

      implicit none

      integer :: n  !! order of the matrix a.
      real(wp) :: a(n*(n+1)/2)  !! dense symmetric matrix stored
                                !! in the packed form which is scaled
                                !! by the value alf (i.e. a:=alf*a).
      real(wp) :: alf  !! scaling factor.

      integer :: i , m

      m = n*(n+1)/2
      do i = 1 , m
         a(i) = a(i)*alf
      end do

      end subroutine mxdsms

!***********************************************************************
!> date: 89/12/01
!
! update of a dense symmetric matrix a. this update is defined as
! a:=a+alf*x*trans(x) where alf is a given scaling factor and x is
! a given vector.

      subroutine mxdsmu(n,a,alf,x)

      implicit none

      integer :: n  !! order of the matrix a.
      real(wp) :: a(n*(n+1)/2)  !! dense symmetric matrix stored in the packed form.
      real(wp) :: alf  !! scaling factor in the correction term.
      real(wp) :: x(n)  !! vector in the correction term.

      real(wp) :: temp
      integer :: i , j , k

      k = 0
      do i = 1 , n
         temp = alf*x(i)
         do j = 1 , i
            k = k + 1
            a(k) = a(k) + temp*x(j)
         end do
      end do

      end subroutine mxdsmu

!***********************************************************************
!> date: 91/12/01
!
! k-th row of a dense symmetric matrix a is copied to the vector x.

      subroutine mxdsmv(n,a,x,k)

      implicit none

      integer :: k      !! index of copied row.
      integer :: n      !! order of the matrix a.
      real(wp) :: a(*)  !! a(n*(n+1)/2)  dense symmetric matrix
                        !! stored in the packed form.
      real(wp) :: x(*)  !! x(n)  output vector.

      integer :: i , l

      l = k*(k-1)/2
      do i = 1 , n
         if ( i<=k ) then
            l = l + 1
         else
            l = l + i - 1
         endif
         x(i) = a(l)
      end do

      end subroutine mxdsmv

!***********************************************************************
!> date: 88/12/01
!
! copying of a vector.

      subroutine mxvcop(n,x,y)

      implicit none

      integer :: n  !! vector dimension.
      real(wp) :: x(*)  !! x(n)  input vector.
      real(wp) :: y(*)  !! y(n)  output vector where y:= x.

      integer :: i

      do i = 1 , n
         y(i) = x(i)
      end do

      end subroutine mxvcop

!***********************************************************************
!> date: 88/12/01
!
! vector difference.

      subroutine mxvdif(n,x,y,z)

      implicit none

      integer :: n  !! vector dimension.
      real(wp) :: x(*)  !! x(n)  input vector.
      real(wp) :: y(*)  !! y(n)  input vector.
      real(wp) :: z(*)  !! z(n)  output vector where z:= x - y.

      integer :: i

      do i = 1 , n
         z(i) = x(i) - y(i)
      end do

      end subroutine mxvdif

!***********************************************************************
!> date: 91/12/01
!
! vector augmented by the scaled vector.

      subroutine mxvdir(n,a,x,y,z)

      implicit none

      integer :: n !! vector dimension.
      real(wp) :: a !! scaling factor.
      real(wp) :: x(*) !! x(n)  input vector.
      real(wp) :: y(*) !! y(n)  input vector.
      real(wp) :: z(*) !! z(n)  output vector where z:= y + a*x.

      integer :: i

      do i = 1 , n
         z(i) = y(i) + a*x(i)
      end do

      end subroutine mxvdir

!***********************************************************************
!> date: 91/12/01
!
! dot product of two vectors.

      function mxvdot(n,x,y)

      implicit none

      integer :: n  !!vector dimension.
      real(wp) :: x(*)  !! x(n)  input vector.
      real(wp) :: y(*)  !! y(n)  input vector.
      real(wp) :: mxvdot  !! value of dot product mxvdot=trans(x)*y.

      real(wp) :: temp
      integer :: i

      temp = 0.0_wp
      do i = 1 , n
         temp = temp + x(i)*y(i)
      end do
      mxvdot = temp

      end function mxvdot

!***********************************************************************
!> date: 90/12/01
!
! elements of the integer vector are replaced by their absolute values.

      subroutine mxvina(n,ix)

      implicit none

      integer :: n  !! dimension of the integer vector.
      integer :: ix(*)  !! vector which is updated so that
                        !! ix(i):=abs(ix(i)) for all i.

      integer :: i

      do i = 1 , n
         ix(i) = abs(ix(i))
         if ( ix(i)>10 ) ix(i) = ix(i) - 10
      end do

      end subroutine mxvina

!***********************************************************************
!> date: 91/12/01
!
! change of the integer vector element for the constraint addition.

      subroutine mxvind(ix,i,job)

      implicit none

      integer :: ix(*)   !! ix(n)  integer vector.
      integer :: i  !! index of the changed element.
      integer :: job  !! change specification. is job.eq.0 then ix(i)=10-ix(i).

      if ( job==0 ) ix(i) = 10 - ix(i)

      end subroutine mxvind

!***********************************************************************
!> date: 90/12/01
!
! initiation of the integer vector.

      subroutine mxvins(n,ip,ix)

      implicit none

      integer :: n  !! dimension of the integer vector.
      integer :: ip  !! integer parameter.
      integer :: ix(n)  !! integer vector such that ix(i)=ip for all i.

      integer :: i

      do i = 1 , n
         ix(i) = ip
      end do

      end subroutine mxvins

!***********************************************************************
!> date: 91/12/01
!
! change of the integer vector element for the constraint addition.

      subroutine mxvinv(ix,i,job)

      implicit none

      integer :: i      !! index of the changed element.
      integer :: job    !! change specification
      integer :: ix(*)  !! ix(n)  integer vector.

      if ( (ix(i)==3 .or. ix(i)==5) .and. job<0 ) ix(i) = ix(i) + 1
      if ( (ix(i)==4 .or. ix(i)==6) .and. job>0 ) ix(i) = ix(i) - 1
      ix(i) = -ix(i)

      end subroutine mxvinv

!***********************************************************************
!> date: 91/12/01
!
! l-infinity norm of a vector.

      function mxvmax(n,x)

      implicit none

      integer :: n  !! vector dimension.
      real(wp) :: x(*)  !! x(n)  input vector.
      real(wp) :: mxvmax  !! l-infinity norm of the vector x.

      integer :: i

      mxvmax = 0.0_wp
      do i = 1 , n
         mxvmax = max(mxvmax,abs(x(i)))
      end do

      end function mxvmax

!***********************************************************************
!> date: 88/12/01
!
! change the signs of vector elements.

      subroutine mxvneg(n,x,y)

      implicit none

      integer :: n  !! vector dimension.
      real(wp) :: x(*)  !! x(n)  input vector.
      real(wp) :: y(*)  !! y(n)  output vector where y:= - x.

      integer :: i

      do i = 1 , n
         y(i) = -x(i)
      end do

      end subroutine mxvneg

!***********************************************************************
!> date: 91/12/01
!
! euclidean norm of a vector.

      function mxvnor(n,x)

      implicit none

      integer :: n        !! vector dimension.
      real(wp) :: x(*)    !! x(n)  input vector.
      real(wp) :: mxvnor  !! euclidean norm of x.

      real(wp) :: pom , den
      integer :: i

      real(wp),parameter :: zero = 0.0_wp

      den = zero
      do i = 1 , n
         den = max(den,abs(x(i)))
      end do
      pom = zero
      if ( den>zero ) then
         do i = 1 , n
            pom = pom + (x(i)/den)**2
         end do
      endif
      mxvnor = den*sqrt(pom)
      end function mxvnor

!***********************************************************************
!> date: 91/12/01
!
! determination of an elementary orthogonal matrix for plane rotation.
!
! parameters :
!  ru  xk  first value for plane rotation (xk is transformed to
!         sqrt(xk**2+xl**2))
!  ru  xl  second value for plane rotation (xl is transformed to
!         zero)
!  ro  ck  diagonal element of the elementary orthogonal matrix.
!  ro  cl  off-diagonal element of the elementary orthogonal matrix.
!  io  ier  information on the transformation. ier=0-general plane
!         rotation. ier=1-permutation. ier=2-transformation suppressed.

      subroutine mxvort(xk,xl,ck,cl,ier)

      implicit none

      real(wp) :: ck , cl , xk , xl
      integer :: ier
      real(wp) :: den , pom

      if ( xl==0.0_wp ) then
         ier = 2
      elseif ( xk==0.0_wp ) then
         xk = xl
         xl = 0.0_wp
         ier = 1
      else
         if ( abs(xk)>=abs(xl) ) then
            pom = xl/xk
            den = sqrt(1.0_wp+pom*pom)
            ck = 1.0_wp/den
            cl = pom/den
            xk = xk*den
         else
            pom = xk/xl
            den = sqrt(1.0_wp+pom*pom)
            cl = 1.0_wp/den
            ck = pom/den
            xk = xl*den
         endif
         xl = 0.0_wp
         ier = 0
      endif

      end subroutine mxvort

!***********************************************************************
!> date: 91/12/01
!
! plane rotation is applied to two values.

      subroutine mxvrot(xk,xl,ck,cl,ier)

      implicit none

      real(wp) :: ck  !! diagonal element of the elementary orthogonal matrix.
      real(wp) :: cl  !! off-diagonal element of the elementary orthogonal matrix.
      real(wp) :: xk  !! first value for plane rotation.
      real(wp) :: xl  !! second value for plane rotation.
      integer :: ier  !! information on the transformation:
                      !!
                      !! * ier=0-general plane rotation.
                      !! * ier=1-permutation.
                      !! * ier=2-transformation suppressed.

      real(wp) :: yk , yl

      if ( ier==0 ) then
         yk = xk
         yl = xl
         xk = ck*yk + cl*yl
         xl = cl*yk - ck*yl
      elseif ( ier==1 ) then
         yk = xk
         xk = xl
         xl = yk
      endif

      end subroutine mxvrot

!***********************************************************************
!> date: 91/12/01
!
! difference of two vectors returned in the subtracted one.

      subroutine mxvsav(n,x,y)

      implicit none

      integer :: n  !! vector dimension.
      real(wp) :: x(*)  !! x(n)  input vector.
      real(wp) :: y(*)  !! y(n)  update vector where y:= x - y.

      real(wp) :: temp
      integer :: i

      do i = 1 , n
         temp = y(i)
         y(i) = x(i) - y(i)
         x(i) = temp
      end do

      end subroutine mxvsav

!***********************************************************************
!> date: 88/12/01
!
! scaling of a vector.

      subroutine mxvscl(n,a,x,y)

      implicit none

      integer :: n  !! vector dimension.
      real(wp) :: a  !! scaling factor.
      real(wp) :: x(*)  !! x(n)  input vector.
      real(wp) :: y(*)  !! y(n)  output vector where y:= a*x.

      integer :: i

      do i = 1 , n
         y(i) = a*x(i)
      end do

      end subroutine mxvscl

!***********************************************************************
!> date: 88/12/01
!
! a scalar is set to all the elements of a vector.

      subroutine mxvset(n,a,x)

      implicit none

      real(wp) :: a  !! initial value.
      integer :: n  !! vector dimension.
      real(wp) :: x(*)  !! x(n)  output vector such that x(i)=a for all i.

      integer :: i
      do i = 1 , n
         x(i) = a
      end do

      end subroutine mxvset

!***********************************************************************
!> date: 88/12/01
!
! sum of two vectors.

      subroutine mxvsum(n,x,y,z)

      implicit none

      integer :: n  !! vector dimension.
      real(wp) :: x(*)  !! x(n)  input vector.
      real(wp) :: y(*)  !! y(n)  input vector.
      real(wp) :: z(*)  !! z(n)  output vector where z:= x + y.

      integer :: i

      do i = 1 , n
         z(i) = x(i) + y(i)
      end do

      end subroutine mxvsum

    end module matrix_routines
