!***********************************************************************
!>
!  Matrix routines.
!
!### History
! * Original version: LU, 1991

   module psqp_matrix_module

      use psqp_kind_module, only: wp => psqp_wp

      implicit none

      public

   contains
!***********************************************************************

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

      integer :: job  !! option
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

      integer :: job    !! option
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
! determination of an elementary orthogonal matrix for plane rotation.

      subroutine mxvort(xk,xl,ck,cl,ier)

      implicit none

      real(wp) :: ck !! diagonal element of the elementary orthogonal matrix.
      real(wp) :: cl !! off-diagonal element of the elementary orthogonal matrix.
      real(wp) :: xk !! first value for plane rotation
                     !! (xk is transformed to sqrt(xk**2+xl**2))
      real(wp) :: xl !! second value for plane rotation
                     !! (xl is transformed to zero)
      integer :: ier !! information on the transformation.
                     !!
                     !! * ier=0-general plane rotation.
                     !! * ier=1-permutation.
                     !! * ier=2-transformation suppressed.

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

    end module psqp_matrix_module
