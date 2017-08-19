!***********************************************************************
!>
!  Matrix routines.
!
!### History
! * Original version: LU, 1991

   module matrix_routines

      implicit none

      public

   contains

!***********************************************************************
!> date: 91/12/01
!
! multiplication of a columnwise stored dense rectangular matrix a
! by a vector x.

      subroutine mxdcmm(n,m,a,x,y)

      implicit none

      integer :: n !! number of rows of the matrix a.
      integer :: m !! number of columns of the matrix a.
      double precision :: a(*) !! a(n*m) rectangular matrix stored columnwise in the one-dimensional array.
      double precision :: x(*) !! x(m)  input vector.
      double precision :: y(*) !! y(n)  output vector equal to a*x.

      integer :: j , k

      call mxvset(n,0.0d0,y)
      k = 0
      do j = 1 , m
         call mxvdir(n,x(j),a(k+1),y,y)
         k = k + n
      enddo
      end subroutine mxdcmm

!***********************************************************************
!> date: 91/12/01
!
! solution of a system of linear equations with a dense symmetric
! positive definite matrix a+e using the factorization a+e=l*d*trans(l)
! obtained by the subroutine mxdpgf.
!
! parameters :
!  ii  n order of the matrix a.
!  ri  a(n*(n+1)/2) factorization a+e=l*d*trans(l) obtained by the
!         subroutine mxdpgf.
!  ru  x(n)  on input the right hand side of a system of linear
!         equations. on output the solution of a system of linear
!         equations.
!  ii  job  option. if job=0 then x:=(a+e)**(-1)*x. if job>0 then
!         x:=l**(-1)*x. if job<0 then x:=trans(l)**(-1)*x.
!
!### Method
! back substitution

      subroutine mxdpgb(n,a,x,job)
      implicit none
      integer :: job , n
      double precision :: a(*) , x(*)
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
            enddo
            ij = ij + 1
         enddo
      endif
      if ( job==0 ) then
!
!     phase 2 : x:=d**(-1)*x
!
         ii = 0
         do i = 1 , n
            ii = ii + i
            x(i) = x(i)/a(ii)
         enddo
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
            enddo
            ii = ii - i
         enddo
      endif
      end subroutine mxdpgb

!***********************************************************************
!> date: 91/12/01
!
! computation of a direction of negative curvature with respect to a
! dense symmetric matrix a using the factorization a+e=l*d*trans(l)
!         obtained by the subroutine mxdpgf.
!
! parameters :
!  ii  n order of the matrix a.
!  ri  a(n*(n+1)/2) factorization a+e=l*d*trans(l) obtained by the
!         subroutine mxdpgf.
!  ro  x(n)  computed direction of negative curvature (i.e.
!         trans(x)*a*x<0) if it exists.
!  ii  inf  information obtained in the factorization process. the
!         direction of negative curvature exists only if inf>0.
!
!### Method
! p.e.gill, w.murray : newton type methods for unconstrained and
! linearly constrained optimization, math. programming 28 (1974)
! pp. 311-350.

      subroutine mxdpgd(n,a,x,inf)
      implicit none
      integer :: n , inf
      double precision :: a(n*(n+1)/2) , x(n)
      integer :: i , j , ii , ij
      double precision :: zero , one
      parameter (zero=0.0d0,one=1.0d0)
!
!     right hand side formation
!
      do i = 1 , n
         x(i) = zero
      enddo
      if ( inf<=0 ) return
      x(inf) = one
!
!     back substitution
!
      ii = inf*(inf-1)/2
      do i = inf - 1 , 1 , -1
         ij = ii
         do j = i + 1 , inf
            ij = ij + j - 1
            x(i) = x(i) - a(ij)*x(j)
         enddo
         ii = ii - i
      enddo
      end subroutine mxdpgd

!***********************************************************************
!> date: 89/12/01
!
! factorization a+e=l*d*trans(l) of a dense symmetric positive definite
! matrix a+e where d and e are diagonal positive definite matrices and
! l is a lower triangular matrix. if a is sufficiently positive
! definite then e=0.
!
! parameters :
!  ii  n order of the matrix a.
!  ru  a(n*(n+1)/2)  on input a given dense symmetric (usually positive
!         definite) matrix a stored in the packed form. on output the
!         computed factorization a+e=l*d*trans(l).
!  io  inf  an information obtained in the factorization process. if
!         inf=0 then a is sufficiently positive definite and e=0. if
!         inf<0 then a is not sufficiently positive definite and e>0. if
!         inf>0 then a is indefinite and inf is an index of the
!         most negative diagonal element used in the factorization
!         process.
!  ru  alf  on input a desired tolerance for positive definiteness. on
!         output the most negative diagonal element used in the
!         factorization process (if inf>0).
!  ro  tau  maximum diagonal element of the matrix e.
!
!### Method
! p.e.gill, w.murray : newton type methods for unconstrained and
! linearly constrained optimization, math. programming 28 (1974)
! pp. 311-350.

      subroutine mxdpgf(n,a,inf,alf,tau)
      implicit none
      double precision :: alf , tau
      integer :: inf , n
      double precision :: a(*)
      double precision :: bet , del , gam , rho , sig , tol
      integer :: i , ij , ik , j , k , kj , kk , l
      l = 0
      inf = 0
      tol = alf
!
!     estimation of the matrix norm
!
      alf = 0.0d0
      bet = 0.0d0
      gam = 0.0d0
      tau = 0.0d0
      kk = 0
      do k = 1 , n
         kk = kk + k
         bet = max(bet,abs(a(kk)))
         kj = kk
         do j = k + 1 , n
            kj = kj + j - 1
            gam = max(gam,abs(a(kj)))
         enddo
      enddo
      bet = max(tol,bet,gam/n)
!      del = tol*bet
      del = tol*max(bet,1.0d0)
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
         gam = 0.0d0
         kj = kk
         do j = k + 1 , n
            kj = kj + j - 1
            gam = max(gam,abs(a(kj)))
         enddo
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
            enddo
         enddo
      enddo
      if ( l>0 .and. abs(alf)>del ) inf = l
      end subroutine mxdpgf

!***********************************************************************
!> date: 91/12/01
!
! estimation of the minimum eigenvalue and the corresponding eigenvector
! of a dense symmetric positive definite matrix a+e using the
! factorization a+e=l*d*trans(l) obtained by the subroutine mxdpgf.
!
! parameters :
!  ii  n order of the matrix a.
!  ri  a(n*(n+1)/2) factorization a+e=l*d*trans(l) obtained by the
!         subroutine mxdpgf.
!  ro  x(n)  estimated eigenvector.
!  ro  alf  estimated eigenvalue.
!  ii  job  option. if job=0 then only estimated eigenvalue is
!         computed. is job>0 then both estimated eigenvalue and
!         estimated eigenvector are computed by job iterations.
!
!### Method
! a.k.cline, c.b.moler, g.w.stewart, j.h.wilkinson : an estimate
! for the condition number of a matrix. siam j. numer. anal. 16
! (1979) 368-373.

      subroutine mxdpgn(n,a,x,alf,job)
      implicit none
      integer :: n , job
      double precision :: a(n*(n+1)/2) , x(n) , alf
      double precision :: xp , xm , fp , fm !, mxvdot
      integer :: i , k , ik , kk
      double precision :: zero , one
      parameter (zero=0.0d0,one=1.0d0)
!
!     computation of the vector v with possible maximum norm such
!     that  l*d**(1/2)*v=u  where u has elements +1 or -1
!
      do i = 1 , n
         x(i) = zero
      enddo
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
         enddo
         if ( fp>=fm ) then
            x(k) = xp
            ik = kk
            do i = k + 1 , n
               ik = ik + i - 1
               x(i) = x(i) + a(ik)*xp
            enddo
         else
            x(k) = xm
            ik = kk
            do i = k + 1 , n
               ik = ik + i - 1
               x(i) = x(i) + a(ik)*xm
            enddo
         endif
      enddo
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
      enddo
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
         enddo
         fm = fm + x(k)*x(k)
         kk = kk - k
      enddo
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
            enddo
            call mxdpgb(n,a,x,0)
            fm = sqrt(mxvdot(n,x,x))
         enddo
         alf = one/fm
      endif
!
!     scaling the vector x by its norm
!
      do i = 1 , n
         x(i) = x(i)/fm
      enddo
      end subroutine mxdpgn

! function mxdpgp                  all systems                91/12/01
!
! computation of the number mxdpgp=trans(x)*d**(-1)*y where d is a
! diagonal matrix in the factorization a+e=l*d*trans(l) obtained by the
! subroutine mxdpgf.
!
! parameters :
!  ii  n order of the matrix a.
!  ri  a(n*(n+1)/2) factorization a+e=l*d*trans(l) obtained by the
!         subroutine mxdpgf.
!  ri  x(n)  input vector.
!  ri  y(n)  input vector.
!  rr  mxdpgp  computed number mxdpgp=trans(x)*d**(-1)*y.
!
      function mxdpgp(n,a,x,y)
      implicit none
      integer :: n
      double precision :: a(*) , x(*) , y(*) , mxdpgp
      double precision :: temp
      integer :: i , j
      j = 0
      temp = 0.0d0
      do i = 1 , n
         j = j + i
         temp = temp + x(i)*y(i)/a(j)
      enddo
      mxdpgp = temp
      end function mxdpgp

!***********************************************************************
!> date: 91/12/01
!
! scaling of a dense symmetric positive definite matrix a+e using the
! factorization a+e=l*d*trans(l) obtained by the subroutine mxdpgf.
!
! parameters :
!  ii  n order of the matrix a.
!  ru  a(n*(n+1)/2) factorization a+e=l*d*trans(l) obtained by the
!         subroutine mxdpgf.
!  ri  alf  scaling factor.

      subroutine mxdpgs(n,a,alf)
      implicit none
      double precision :: alf
      integer :: n
      double precision :: a(*)
      integer :: i , j
      j = 0
      do i = 1 , n
         j = j + i
         a(j) = a(j)*alf
      enddo
      end subroutine mxdpgs

!***********************************************************************
!> date: 89/12/01
!
! correction of a dense symmetric positive definite matrix a+e in the
! factored form a+e=l*d*trans(l) obtained by the subroutine mxdpgf.
! the correction is defined as a+e:=a+e+alf*x*trans(x) where alf is a
! given scaling factor and x is a given vector.
!
! parameters :
!  ii  n order of the matrix a.
!  ru  a(n*(n+1)/2) factorization a+e=l*d*trans(l) obtained by the
!         subroutine mxdpgf.
!  ri  alf  scaling factor in the correction term.
!  ri  x(n)  vector in the correction term.
!  ra  y(n) auxiliary vector.
!
!### Method
! p.e.gill, w.murray, m.saunders: methods for computing and modifying
! the ldv factors of a matrix, math. of comp. 29 (1974) pp. 1051-1077.

      subroutine mxdpgu(n,a,alf,x,y)
      implicit none
      double precision :: zero , one , four , con
      parameter (zero=0.0d0,one=1.0d0,four=4.0d0,con=1.0d-8)
      double precision :: alf , alfr
      integer :: n
      double precision :: a(*) , x(*) , y(*)
      double precision :: b , d , p , r , t , to
      integer :: i , ii , ij , j
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
               enddo
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
               enddo
            endif
            to = t
         enddo
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
            enddo
            y(i) = d
            ij = ij + 1
            to = to - d*d/a(ij)
         enddo
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
            enddo
            ii = ii - i
         enddo
      endif
      end subroutine mxdpgu

!***********************************************************************
!> date: 89/12/01
!
! solution of a system of linear equations with a dense symmetric
! positive definite matrix a using the factorization a=trans(r)*r.
!
! parameters :
!  ii  n order of the matrix a.
!  ri  a(n*(n+1)/2) factorization a=trans(r)*r.
!  ru  x(n)  on input the right hand side of a system of linear
!         equations. on output the solution of a system of linear
!         equations.
!  ii  job  option. if job=0 then x:=a**(-1)*x. if job>0 then
!         x:=trans(r)**(-1)*x. if job<0 then x:=r**(-1)*x.
!
!### Method
! back substitution

      subroutine mxdprb(n,a,x,job)
      implicit none
      integer :: job , n
      double precision :: a(*) , x(*)
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
            enddo
            ij = ij + 1
            x(i) = x(i)/a(ij)
         enddo
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
            enddo
            x(i) = x(i)/a(ii)
            ii = ii - i
         enddo
      endif
      end subroutine mxdprb

!***********************************************************************
!> date: 92/12/01
!
! correction of a singular dense symmetric positive semidefinite matrix
! a decomposed as a=trans(r)*r.
!
! parameters :
!  ii  n  order of the matrix a.
!  ru  a(n*(n+1)/2)  dense symmetric matrix stored in the packed form.
!  io  inf  an information obtained in the correction process. if
!         inf=0 then a is sufficiently positive definite. if
!         inf<0 then a is not sufficiently positive definite.
!         process.
!  ri  tol  desired tolerance for positive definiteness.

      subroutine mxdprc(n,a,inf,tol)
      implicit none
      integer :: n , inf
      double precision :: a(n*(n+1)/2) , tol
      double precision :: tol1 , temp
      integer :: l , i
      inf = 0
      tol1 = sqrt(tol)
      temp = tol1
      do i = 1 , n*(n+1)/2
         temp = max(temp,abs(a(i)))
      enddo
      temp = temp*tol1
      l = 0
      do i = 1 , n
         l = l + i
         if ( abs(a(l))<=temp ) then
            a(l) = sign(temp,a(l))
            inf = -1
         endif
      enddo
      end subroutine mxdprc

!***********************************************************************
!> date: 91/12/01
!
! multiplication of a given vector x by a dense symmetric positive
! definite matrix a using the factorization a=trans(r)*r.
!
! parameters :
!  ii  n order of the matrix a.
!  ri  a(n*(n+1)/2) factorization a=trans(r)*r.
!  ru  x(n)  on input the given vector. on output the result of
!         multiplication.
!  ii  job  option. if job=0 then x:=a*x. if job>0 then x:=r*x.
!         if job<0 then x:=trans(r)*x.

      subroutine mxdprm(n,a,x,job)
      implicit none
      integer :: n , job
      double precision :: a(n*(n+1)/2) , x(n)
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
            enddo
         enddo
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
            enddo
            ij = ij - 1
         enddo
      endif
      end subroutine mxdprm

!***********************************************************************
!> date: 91/12/01
!
! plane rotation is applied to a rowwise stored dense rectangular
! matrix a.
!
! parameters :
!  ii  n  number of columns of the matrix a.
!  ru  a(m*n)  rectangular matrix stored rowwise in the
!         one-dimensional array.
!  ii  k  first index of the plane rotation.
!  ii  l  second index of the plane rotation.
!  ri  ck  diagonal element of the elementary orthogonal matrix.
!  ri  cl  off-diagonal element of the elementary orthogonal matrix.
!  ii  ier  type of the plane rotation. ier=0-general plane rotation.
!         ier=1-permutation. ier=2-transformation suppressed.

      subroutine mxdrgr(n,a,k,l,ck,cl,ier)
      implicit none
      integer :: n , k , l , ier
      double precision :: a(*) , ck , cl
      integer :: i , ik , il
      if ( ier/=0 .and. ier/=1 ) return
      ik = (k-1)*n
      il = (l-1)*n
      do i = 1 , n
         ik = ik + 1
         il = il + 1
         call mxvrot(a(ik),a(il),ck,cl,ier)
      enddo
      end subroutine mxdrgr

!***********************************************************************
!> date: 91/12/01
!
! multiplication of a rowwise stored dense rectangular matrix a by
! a vector x and addition of a scaled vector alf*y.
!
! parameters :
!  ii  n  number of columns of the matrix a.
!  ii  m  number of rows of the matrix a.
!  ri  a(m*n)  rectangular matrix stored rowwise in the
!         one-dimensional array.
!  ri  x(n)  input vector.
!  ri  alf  scaling factor.
!  ri  y(m)  input vector.
!  ro  z(m)  output vector equal to a*x+alf*y.

      subroutine mxdrmd(n,m,a,x,alf,y,z)
      implicit none
      integer :: n , m
      double precision :: a(m*n) , x(n) , alf , y(m) , z(m)
      double precision :: temp
      integer :: i , j , k
      k = 0
      do j = 1 , m
         temp = alf*y(j)
         do i = 1 , n
            temp = temp + a(k+i)*x(i)
         enddo
         z(j) = temp
         k = k + n
      enddo
      end subroutine mxdrmd

!***********************************************************************
!> date: 91/12/01
!
! rowwise stored dense rectangular matrix a is set to be a part of the
! unit matrix.
!
! parameters :
!  ii  n  number of columns of the matrix a.
!  ii  m  number of rows of the matrix a.
!  ro  a(m*n)  rectangular matrix stored rowwise in the one-dimensional
!          array. this matrix is set to trans([i,0]).

      subroutine mxdrmi(n,m,a)
      implicit none
      integer :: n , m
      double precision :: a(m*n)
      integer :: i , j , k
      double precision :: zero , one
      parameter (zero=0.0d0,one=1.0d0)
      k = 0
      do j = 1 , m
         do i = 1 , n
            a(i+k) = zero
            if ( i==j ) a(i+k) = one
         enddo
         k = k + n
      enddo
      end subroutine mxdrmi

!***********************************************************************
!> date: 91/12/01
!
! multiplication of a rowwise stored dense rectangular matrix a by
! a vector x.
!
! parameters :
!  ii  n  number of columns of the matrix a.
!  ii  m  number of rows of the matrix a.
!  ri  a(m*n)  rectangular matrix stored rowwise in the
!         one-dimensional array.
!  ri  x(n)  input vector.
!  ro  y(m)  output vector equal to a*x.

      subroutine mxdrmm(n,m,a,x,y)
      implicit none
      integer :: n , m
      double precision :: a(*) , x(*) , y(*)
      double precision :: temp
      integer :: i , j , k
      k = 0
      do j = 1 , m
         temp = 0.0d0
         do i = 1 , n
            temp = temp + a(k+i)*x(i)
         enddo
         y(j) = temp
         k = k + n
      enddo
      end subroutine mxdrmm

!***********************************************************************
!> date: 91/12/01
!
! euclidean norm of a part of the i-th column of a rowwise stored dense
! rectangular matrix a is computed.
!
! parameters :
!  ii  n  number of columns of the matrix a.
!  ii  m  number of rows of the matrix a.
!  ri  a(m*n)  rectangular matrix stored rowwise in the
!         one-dimensional array.
!  ii  i  index of the column whose norm is computed.
!  ii  j  index of the first element from which the norm is computed.

      function mxdrmn(n,m,a,i,j)
      implicit none
      integer :: n , m , i , j
      double precision :: a(m*n) , mxdrmn
      double precision :: pom , den
      integer :: k , l
      double precision :: zero
      parameter (zero=0.0d0)
      den = zero
      l = (j-1)*n
      do k = j , m
         den = max(den,abs(a(l+i)))
         l = l + n
      enddo
      pom = zero
      if ( den>zero ) then
         l = (j-1)*n
         do k = j , m
            pom = pom + (a(l+i)/den)**2
            l = l + n
         enddo
      endif
      mxdrmn = den*sqrt(pom)
      end function mxdrmn

!***********************************************************************
!> date: 91/12/01
!
! k-th column of a rowwise stored dense rectangular matrix a is copied
! to the vector x.
!
! parameters :
!  ii  n  number of columns of the matrix a.
!  ii  m  number of rows of the matrix a.
!  ri  a(m*n)  rectangular matrix stored rowwise in the
!         one-dimensional array.
!  ro  x(m)  output vector such that x(j)=a(j,k) for all j.
!  ii  k  index of the row being copied to the output vector.

      subroutine mxdrmv(n,m,a,x,k)
      implicit none
      integer :: n , m , k
      double precision :: a(*) , x(*)
      integer :: i , j
      if ( k<1 .or. k>n ) return
      i = k
      do j = 1 , m
         x(j) = a(i)
         i = i + n
      enddo
      end subroutine mxdrmv

!***********************************************************************
!> date: 92/12/01
!
! qr decomposition of rowwise stored dense rectangular matrix q using
! householder transformations without pivoting.
!
! parameters :
!  ii  n  number of columns of the matrix q.
!  ii  m  number of rows of the matrix q.
!  ru  q(m*n)  rectangular matrix stored rowwise in the
!         one-dimensional array.
!  ro  r(n*(n+1)/2)  upper triangular matrix stored in the packed form.
!
!### Method
! p.a.bussinger, g.h.golub : linear least squares solution by
! householder transformation. numer. math. 7 (1965) 269-276.

      subroutine mxdrqf(n,m,q,r)
      implicit none
      integer :: n , m
      double precision :: q(m*n) , r(n*(n+1)/2)
      double precision :: alf , pom
      integer :: i , j , k , l , jp , kp , nm
      double precision :: zero , one
      parameter (zero=0.0d0,one=1.0d0)
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
         enddo
         if ( pom/=zero ) then
!
!     householder transformation
!
            do j = k , m
               q(jp+k) = q(jp+k)/pom
               jp = jp + n
            enddo
            q(kp+k) = q(kp+k) + one
            do i = k + 1 , n
               alf = zero
               jp = kp
               do j = k , m
                  alf = alf + q(jp+k)*q(jp+i)
                  jp = jp + n
               enddo
               alf = alf/q(kp+k)
               jp = kp
               do j = k , m
                  q(jp+i) = q(jp+i) - alf*q(jp+k)
                  jp = jp + n
               enddo
            enddo
         endif
         l = l + k
         kp = kp + n
      enddo
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
               enddo
               alf = alf/q(kp+k)
               jp = kp
               do j = k , m
                  q(jp+i) = q(jp+i) - alf*q(jp+k)
                  jp = jp + n
               enddo
            enddo
            jp = kp
            do j = k , m
               q(jp+k) = -q(jp+k)
               jp = jp + n
            enddo
         endif
         q(kp+k) = q(kp+k) + one
      enddo
      end subroutine mxdrqf

!***********************************************************************
!> date: 92/12/01
!
! update of a qr decomposition. this qr decomposition is updated
! by the rule q*r:=q*r+alf*x*trans(y).
!
! parameters :
!  ii  n  number of columns of the matrix q.
!  ii  m  number of rows of the matrix q.
!  ru  q(m*n)  rectangular matrix stored rowwise in the
!         one-dimensional array (part of the orthogonal matrix).
!  ru  r(n*(n+1)/2)  upper triangular matrix stored in a packed form.
!  ri  alf  scalar parameter.
!  ri  x(m)  input vector.
!  ri  y(n)  input vector.
!  ra  z(n)  auxiliary vector.
!  io  inf  information. if inf=0 then x lies in the column space of q.
!         if inf=1 then x does not lie in the column space of q.
!
!### Method
! j.w.daniel, w.b.gragg, l.kaufman, g.w.steward : reorthogonalization
! and stable algorithms for updating the gram-schmidt qr factorization.
! mathematics of computation 30 (1976) 772-795.

      subroutine mxdrqu(n,m,q,r,alf,x,y,z,inf)
      implicit none
      integer :: n , m , inf
      double precision :: q(m*n) , r(n*(n+1)/2) , alf , x(m) , y(n) , z(n)
      double precision :: ck , cl , zk , zl !, mxvnor
      integer :: j , k , l , kj , kk , ier
      double precision :: one , con
      parameter (one=1.0d0,con=1.0d-10)
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
               enddo
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
            enddo
            kj = k
            do j = 1 , m
               call mxvrot(q(kj),q(kj+1),ck,cl,ier)
               kj = kj + n
            enddo
         endif
      enddo
      z(1) = alf*z(1)
      kj = 1
      do j = 1 , n
         r(kj) = r(kj) + z(1)*y(j)
         kj = kj + j
      enddo
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
            enddo
            kj = k
            do j = 1 , m
               call mxvrot(q(kj),q(kj+1),ck,cl,ier)
               kj = kj + n
            enddo
         endif
         kk = kk + l
      enddo
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
            enddo
         endif
      endif
      end subroutine mxdrqu

!***********************************************************************
!> date: 91/12/01
!
! a dense symmetric matrix a is augmented by the scaled unit matrix
! such that a:=a+alf*i (i is the unit matrix of order n).
!
! parameters :
!  ii  n  order of the matrix a.
!  ru  a(n*(n+1)/2)  dense symmetric matrix stored in the packed form.
!  ri  alf  scaling factor.

      subroutine mxdsda(n,a,alf)
      implicit none
      integer :: n
      double precision :: a(*) , alf
      integer :: i , j
      j = 0
      do i = 1 , n
         j = j + i
         a(j) = a(j) + alf
      enddo
      end subroutine mxdsda

!***********************************************************************
!> date: 91/12/01
!
! determination of the minimum diagonal element of a dense symmetric
! matrix.
!
! parameters :
!  ii  n  order of the matrix a.
!  ri  a(n*(n+1)/2)  dense symmetric matrix stored in the packed form.
!  io  inf  index of the minimum diagonal element of the matrix a.
!  rr  mxdsdl  minimum diagonal element of the matrix a.

      function mxdsdl(n,a,inf)
      implicit none
      integer :: n , inf
      double precision :: a(n*(n+1)/2) , mxdsdl
      double precision :: temp
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
      enddo
      mxdsdl = temp
      end function mxdsdl

!***********************************************************************
!> date: 91/12/01
!
! dense symmetric matrix augmented by the scaled dense symmetric matrix.
!
! parameters :
!  ii  n  order of the matrices.
!  ri  alf  scaling factor.
!  ri  a(n*(n+1)/2)  input matrix.
!  ri  b(n*(n+1)/2)  input matrix.
!  ro  c(n*(n+1)/2)  output matrix where c:=b+alf*a.

      subroutine mxdsma(n,alf,a,b,c)
      implicit none
      integer :: n
      double precision :: alf , a(n*(n+1)/2) , b(n*(n+1)/2) , c(n*(n+1)/2)
      integer :: i
      do i = 1 , n*(n+1)/2
         c(i) = b(i) + alf*a(i)
      enddo
      end subroutine mxdsma

!***********************************************************************
!> date: 91/12/01
!
! copying of a dense symmetric matrix.
!
! parameters :
!  ii  n  order of the matrices a and b.
!  ri  a(n*(n+1)/2)  dense symmetric matrix stored in the packed form.
!  ro  b(n*(n+1)/2)  dense symmetric matrix stored in the packed form
!         where b:=a.

      subroutine mxdsmc(n,a,b)
      implicit none
      integer :: n
      double precision :: a(n*(n+1)/2) , b(n*(n+1)/2)
      integer :: m , i
      m = n*(n+1)/2
      do i = 1 , m
         b(i) = a(i)
      enddo
      end subroutine mxdsmc

!***********************************************************************
!> date: 91/12/01
!
!  gershgorin bounds for eigenvalues of a dense symmetric matrix.
!  amin .le. any eigenvalue of a .le. amax.
!
! parameters :
!  ii  n  dimension of the matrix a.
!  ri  a(n*(n+1)/2)  dense symmetric matrix stored in the packed form.
!  ro  amin  lower bound for eigenvalues of a.
!  ro  amax  upper bound for eigenvalues of a.

      subroutine mxdsmg(n,a,amin,amax)
      implicit none
      integer :: n
      double precision :: a(n*(n+1)/2) , amin , amax
      double precision :: temp
      integer :: i , j , k , l
      double precision :: zero
      parameter (zero=0.0d0)
      amax = a(1)
      amin = a(1)
      k = 0
      do i = 1 , n
         temp = zero
         l = k
         do j = 1 , i - 1
            l = l + 1
            temp = temp + abs(a(l))
         enddo
         l = l + 1
         do j = i + 1 , n
            l = l + j - 1
            temp = temp + abs(a(l))
         enddo
         k = k + i
         amax = max(amax,a(k)+temp)
         amin = min(amin,a(k)-temp)
      enddo
      end subroutine mxdsmg

!***********************************************************************
!> date: 88/12/01
!
! dense symmetric matrix a is set to the unit matrix with the same
! order.
!
! parameters :
!  ii  n  order of the matrix a.
!  ro  a(n*(n+1)/2)  dense symmetric matrix stored in the packed form
!         which is set to the unit matrix (i.e. a:=i).

      subroutine mxdsmi(n,a)
      implicit none
      integer :: n
      double precision :: a(*)
      integer :: i , m
      m = n*(n+1)/2
      do i = 1 , m
         a(i) = 0.0d0
      enddo
      m = 0
      do i = 1 , n
         m = m + i
         a(m) = 1.0d0
      enddo
      end subroutine mxdsmi

!***********************************************************************
!> date: 89/12/01
!
! multiplication of a dense symmetric matrix a by a vector x.
!
! parameters :
!  ii  n  order of the matrix a.
!  ri  a(n*(n+1)/2)  dense symmetric matrix stored in the packed form.
!  ri  x(n)  input vector.
!  ro  y(n)  output vector equal to  a*x.

      subroutine mxdsmm(n,a,x,y)
      implicit none
      integer :: n
      double precision :: a(*) , x(*) , y(*)
      double precision :: temp
      integer :: i , j , k , l
      k = 0
      do i = 1 , n
         temp = 0.0d0
         l = k
         do j = 1 , i
            l = l + 1
            temp = temp + a(l)*x(j)
         enddo
         do j = i + 1 , n
            l = l + j - 1
            temp = temp + a(l)*x(j)
         enddo
         y(i) = temp
         k = k + i
      enddo
      end subroutine mxdsmm

!***********************************************************************
!> date: 91/12/01
!
! value of a quadratic form with a dense symmetric matrix a.
!
! parameters :
!  ii  n  order of the matrix a.
!  ri  a(n*(n+1)/2)  dense symmetric matrix stored in the packed form.
!  ri  x(n)  given vector.
!  ri  y(n)  given vector.
!  rr  mxdsmq  value of the quadratic form mxdsmq=trans(x)*a*y.

      double precision function mxdsmq(n,a,x,y)
      implicit none
      double precision :: zero
      parameter (zero=0.0d0)
      integer :: n
      double precision :: a(*) , x(*) , y(*)
      double precision :: temp , temp1 , temp2
      integer :: i , j , k
      temp = zero
      k = 0
      do i = 1 , n
         temp1 = zero
         temp2 = zero
         do j = 1 , i - 1
            k = k + 1
            temp1 = temp1 + a(k)*x(j)
            temp2 = temp2 + a(k)*y(j)
         enddo
         k = k + 1
         temp = temp + x(i)*(temp2+a(k)*y(i)) + y(i)*temp1
      enddo
      mxdsmq = temp
      end function mxdsmq

!***********************************************************************
!> date: 92/12/01
!
! plane rotation is applied to a dense symmetric matrix a. the case
! k=l+1 is required.
!
! parameters :
!  ii  n  order of the matrix a.
!  ru  a(n*(n+1)/2) dense symmetric matrix stored in the packed form.
!  ii  k  first index of plane rotation.
!  ii  l  second index of plane rotation.
!  ro  ck  diagonal element of the elementary orthogonal matrix.
!  ro  cl  off-diagonal element of the elementary orthogonal matrix.
!  io  ier  information on the transformation. ier<0-k or l out of
!         range. ier=0-plane rotation. ier=1-permutation.
!         ier=2-transformation suppressed.

      subroutine mxdsmr(n,a,k,l,ck,cl,ier)
      implicit none
      integer :: n , k , l , ier
      double precision :: a(*) , ck , cl
      double precision :: akk , akl , all , ckk , ckl , cll
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
      enddo
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
!
! parameters :
!  ii  n  order of the matrix a.
!  ru  a(n*(n+1)/2)  dense symmetric matrix stored in the packed form
!         which is scaled by the value alf (i.e. a:=alf*a).
!  ri  alf  scaling factor.

      subroutine mxdsms(n,a,alf)
      implicit none
      integer :: n
      double precision :: a(n*(n+1)/2) , alf
      integer :: i , m
      m = n*(n+1)/2
      do i = 1 , m
         a(i) = a(i)*alf
      enddo
      end subroutine mxdsms

!***********************************************************************
!> date: 89/12/01
!
! update of a dense symmetric matrix a. this update is defined as
! a:=a+alf*x*trans(x) where alf is a given scaling factor and x is
! a given vector.
!
! parameters :
!  ii  n  order of the matrix a.
!  ru  a(n*(n+1)/2)  dense symmetric matrix stored in the packed form.
!  ri  alf  scaling factor in the correction term.
!  ri  x(n)  vector in the correction term.

      subroutine mxdsmu(n,a,alf,x)
      implicit none
      integer :: n
      double precision :: a(n*(n+1)/2) , x(n) , alf
      double precision :: temp
      integer :: i , j , k
      k = 0
      do i = 1 , n
         temp = alf*x(i)
         do j = 1 , i
            k = k + 1
            a(k) = a(k) + temp*x(j)
         enddo
      enddo
      end subroutine mxdsmu

!***********************************************************************
!> date: 91/12/01
!
! k-th row of a dense symmetric matrix a is copied to the vector x.
!
! parameters :
!  ii  n  order of the matrix a.
!  ri  a(n*(n+1)/2)  dense symmetric matrix stored in the packed form.
!  ro  x(n)  output vector.
!  ii  k  index of copied row.

      subroutine mxdsmv(n,a,x,k)
      implicit none
      integer :: k , n
      double precision :: a(*) , x(*)
      integer :: i , l
      l = k*(k-1)/2
      do i = 1 , n
         if ( i<=k ) then
            l = l + 1
         else
            l = l + i - 1
         endif
         x(i) = a(l)
      enddo
      end subroutine mxdsmv

!***********************************************************************
!> date: 88/12/01
!
! copying of a vector.
!
! parameters :
!  ii  n  vector dimension.
!  ri  x(n)  input vector.
!  ro  y(n)  output vector where y:= x.

      subroutine mxvcop(n,x,y)
      implicit none
      integer :: n
      double precision :: x(*) , y(*)
      integer :: i
      do i = 1 , n
         y(i) = x(i)
      enddo
      end subroutine mxvcop

!***********************************************************************
!> date: 88/12/01
!
! vector difference.
!
! parameters :
!  ri  x(n)  input vector.
!  ri  y(n)  input vector.
!  ro  z(n)  output vector where z:= x - y.

      subroutine mxvdif(n,x,y,z)
      implicit none
      integer :: n
      double precision :: x(*) , y(*) , z(*)
      integer :: i
      do i = 1 , n
         z(i) = x(i) - y(i)
      enddo
      end subroutine mxvdif

!***********************************************************************
!> date: 91/12/01
!
! vector augmented by the scaled vector.
!
! parameters :
!  ii  n  vector dimension.
!  ri  a  scaling factor.
!  ri  x(n)  input vector.
!  ri  y(n)  input vector.
!  ro  z(n)  output vector where z:= y + a*x.

      subroutine mxvdir(n,a,x,y,z)
      implicit none
      double precision :: a
      integer :: n
      double precision :: x(*) , y(*) , z(*)
      integer :: i
      do i = 1 , n
         z(i) = y(i) + a*x(i)
      enddo
      end subroutine mxvdir

!***********************************************************************
!> date: 91/12/01
!
! dot product of two vectors.
!
! parameters :
!  ii  n  vector dimension.
!  ri  x(n)  input vector.
!  ri  y(n)  input vector.
!  rr  mxvdot  value of dot product mxvdot=trans(x)*y.
!
      function mxvdot(n,x,y)
      implicit none
      integer :: n
      double precision :: x(*) , y(*) , mxvdot
      double precision :: temp
      integer :: i
      temp = 0.0d0
      do i = 1 , n
         temp = temp + x(i)*y(i)
      enddo
      mxvdot = temp
      end function mxvdot

!***********************************************************************
!> date: 90/12/01
!
! elements of the integer :: vector are replaced by their absolute values.
!
! parameters :
!  ii  n dimension of the integer :: vector.
!  iu  ix(n)  integer :: vector which is updated so that ix(i):=abs(ix(i))
!         for all i.

      subroutine mxvina(n,ix)
      implicit none
      integer :: n
      integer :: ix(*)
      integer :: i
      do i = 1 , n
         ix(i) = abs(ix(i))
         if ( ix(i)>10 ) ix(i) = ix(i) - 10
      enddo
      end subroutine mxvina

!***********************************************************************
!> date: 91/12/01
!
! change of the integer :: vector element for the constraint addition.
!
! parameters :
!  iu  ix(n)  integer :: vector.
!  ii  i  index of the changed element.
!  ii job  change specification. is job.eq.0 then ix(i)=10-ix(i).

      subroutine mxvind(ix,i,job)
      implicit none
      integer :: ix(*) , i , job
      if ( job==0 ) ix(i) = 10 - ix(i)
      end subroutine mxvind

!***********************************************************************
!> date: 90/12/01
!
! initiation of the integer :: vector.
!
! parameters :
!  ii  n dimension of the integer :: vector.
!  ii  ip  integer :: parameter.
!  io  ix(n)  integer :: vector such that ix(i)=ip for all i.

      subroutine mxvins(n,ip,ix)
      implicit none
      integer :: n , ip , ix(n)
      integer :: i
      do i = 1 , n
         ix(i) = ip
      enddo
      end subroutine mxvins

!***********************************************************************
!> date: 91/12/01
!
! change of the integer :: vector element for the constraint addition.
!
! parameters :
!  ii  n  vector dimension.
!  iu  ix(n)  integer :: vector.
!  ii  i  index of the changed element.
!  ii  job  change specification

      subroutine mxvinv(ix,i,job)
      implicit none
      integer :: i , job
      integer :: ix(*)
      if ( (ix(i)==3 .or. ix(i)==5) .and. job<0 ) ix(i) = ix(i) + 1
      if ( (ix(i)==4 .or. ix(i)==6) .and. job>0 ) ix(i) = ix(i) - 1
      ix(i) = -ix(i)
      end subroutine mxvinv

!***********************************************************************
!> date: 91/12/01
!
! l-infinity norm of a vector.
!
! parameters :
!  ii  n  vector dimension.
!  ri  x(n)  input vector.
!  rr  mxvmax  l-infinity norm of the vector x.
!
      function mxvmax(n,x)
      implicit none
      integer :: n
      double precision :: x(*) , mxvmax
      integer :: i
      mxvmax = 0.0d0
      do i = 1 , n
         mxvmax = max(mxvmax,abs(x(i)))
      enddo
      end function mxvmax

!***********************************************************************
!> date: 88/12/01
!
! change the signs of vector elements.
!
! parameters :
!  ii  n  vector dimension.
!  ri  x(n)  input vector.
!  ro  y(n)  output vector where y:= - x.

      subroutine mxvneg(n,x,y)
      implicit none
      integer :: n
      double precision :: x(*) , y(*)
      integer :: i
      do i = 1 , n
         y(i) = -x(i)
      enddo
      end subroutine mxvneg

!***********************************************************************
!> date: 91/12/01
!
! euclidean norm of a vector.
!
! parameters :
!  ii  n  vector dimension.
!  ri  x(n)  input vector.
!  rr  mxvnor  euclidean norm of x.
!
      function mxvnor(n,x)
      implicit none
      integer :: n
      double precision :: x(*) , mxvnor
      double precision :: pom , den
      integer :: i
      double precision :: zero
      parameter (zero=0.0d0)
      den = zero
      do i = 1 , n
         den = max(den,abs(x(i)))
      enddo
      pom = zero
      if ( den>zero ) then
         do i = 1 , n
            pom = pom + (x(i)/den)**2
         enddo
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
      double precision :: ck , cl , xk , xl
      integer :: ier
      double precision :: den , pom
      if ( xl==0.0d0 ) then
         ier = 2
      elseif ( xk==0.0d0 ) then
         xk = xl
         xl = 0.0d0
         ier = 1
      else
         if ( abs(xk)>=abs(xl) ) then
            pom = xl/xk
            den = sqrt(1.0d0+pom*pom)
            ck = 1.0d0/den
            cl = pom/den
            xk = xk*den
         else
            pom = xk/xl
            den = sqrt(1.0d0+pom*pom)
            cl = 1.0d0/den
            ck = pom/den
            xk = xl*den
         endif
         xl = 0.0d0
         ier = 0
      endif
      end subroutine mxvort

!***********************************************************************
!> date: 91/12/01
!
! plane rotation is applied to two values.
!
! parameters :
!  ru  xk  first value for plane rotation.
!  ru  xl  second value for plane rotation.
!  ri  ck  diagonal element of the elementary orthogonal matrix.
!  ri  cl  off-diagonal element of the elementary orthogonal matrix.
!  ii  ier  information on the transformation. ier=0-general plane
!         rotation. ier=1-permutation. ier=2-transformation suppressed.

      subroutine mxvrot(xk,xl,ck,cl,ier)
      implicit none
      double precision :: ck , cl , xk , xl
      integer :: ier
      double precision :: yk , yl
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
! difference of two vectors returned in the substracted one.
!
! parameters :
!  ii  n  vector dimension.
!  ri  x(n)  input vector.
!  ru  y(n)  update vector where y:= x - y.

      subroutine mxvsav(n,x,y)
      implicit none
      integer :: n
      double precision :: x(*) , y(*)
      double precision :: temp
      integer :: i
      do i = 1 , n
         temp = y(i)
         y(i) = x(i) - y(i)
         x(i) = temp
      enddo
      end subroutine mxvsav

!***********************************************************************
!> date: 88/12/01
!
! scaling of a vector.
!
! parameters :
!  ii  n  vector dimension.
!  ri  x(n)  input vector.
!  ri  a  scaling factor.
!  ro  y(n)  output vector where y:= a*x.

      subroutine mxvscl(n,a,x,y)
      implicit none
      double precision :: a
      integer :: n
      double precision :: x(*) , y(*)
      integer :: i
      do i = 1 , n
         y(i) = a*x(i)
      enddo
      end subroutine mxvscl

!***********************************************************************
!> date: 88/12/01
!
! a scalar is set to all the elements of a vector.
!
! parameters :
!  ii  n  vector dimension.
!  ri  a  initial value.
!  ro  x(n)  output vector such that x(i)=a for all i.

      subroutine mxvset(n,a,x)
      implicit none
      double precision :: a
      integer :: n
      double precision :: x(*)
      integer :: i
      do i = 1 , n
         x(i) = a
      enddo
      end subroutine mxvset

!***********************************************************************
!> date: 88/12/01
!
! sum of two vectors.
!
! parameters :
!  ii  n  vector dimension.
!  ri  x(n)  input vector.
!  ri  y(n)  input vector.
!  ro  z(n)  output vector where z:= x + y.

      subroutine mxvsum(n,x,y,z)
      implicit none
      integer :: n
      double precision :: x(*) , y(*) , z(*)
      integer :: i
      do i = 1 , n
         z(i) = x(i) + y(i)
      enddo
      end subroutine mxvsum

    end module matrix_routines
