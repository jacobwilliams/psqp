
   module matrix_routines

      implicit NONE

      public

   contains

! SUBROUTINE MXDCMM               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! MULTIPLICATION OF A COLUMNWISE STORED DENSE RECTANGULAR MATRIX A
! BY A VECTOR X.
!
! PARAMETERS :
!  II  N  NUMBER OF ROWS OF THE MATRIX A.
!  II  M  NUMBER OF COLUMNS OF THE MATRIX A.
!  RI  A(N*M)  RECTANGULAR MATRIX STORED COLUMNWISE IN THE
!         ONE-DIMENSIONAL ARRAY.
!  RI  X(M)  INPUT VECTOR.
!  RO  Y(N)  OUTPUT VECTOR EQUAL TO A*X.
!
! SUBPROGRAMS USED :
!  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
!  S   MXVSET  INITIATION OF A VECTOR.
!
      SUBROUTINE MXDCMM(N,M,A,X,Y)
      IMPLICIT NONE
      INTEGER N , M
      DOUBLE PRECISION A(*) , X(*) , Y(*)
      INTEGER j , k
      CALL MXVSET(N,0.0D0,Y)
      k = 0
      DO j = 1 , M
         CALL MXVDIR(N,X(j),A(k+1),Y,Y)
         k = k + N
      ENDDO
      END SUBROUTINE MXDCMM

! SUBROUTINE MXDPGB                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! SOLUTION OF A SYSTEM OF LINEAR EQUATIONS WITH A DENSE SYMMETRIC
! POSITIVE DEFINITE MATRIX A+E USING THE FACTORIZATION A+E=L*D*TRANS(L)
! OBTAINED BY THE SUBROUTINE MXDPGF.
!
! PARAMETERS :
!  II  N ORDER OF THE MATRIX A.
!  RI  A(N*(N+1)/2) FACTORIZATION A+E=L*D*TRANS(L) OBTAINED BY THE
!         SUBROUTINE MXDPGF.
!  RU  X(N)  ON INPUT THE RIGHT HAND SIDE OF A SYSTEM OF LINEAR
!         EQUATIONS. ON OUTPUT THE SOLUTION OF A SYSTEM OF LINEAR
!         EQUATIONS.
!  II  JOB  OPTION. IF JOB=0 THEN X:=(A+E)**(-1)*X. IF JOB>0 THEN
!         X:=L**(-1)*X. IF JOB<0 THEN X:=TRANS(L)**(-1)*X.
!
! METHOD :
! BACK SUBSTITUTION
!
      SUBROUTINE MXDPGB(N,A,X,Job)
      IMPLICIT NONE
      INTEGER Job , N
      DOUBLE PRECISION A(*) , X(*)
      INTEGER i , ii , ij , j
      IF ( Job>=0 ) THEN
!
!     PHASE 1 : X:=L**(-1)*X
!
         ij = 0
         DO i = 1 , N
            DO j = 1 , i - 1
               ij = ij + 1
               X(i) = X(i) - A(ij)*X(j)
            ENDDO
            ij = ij + 1
         ENDDO
      ENDIF
      IF ( Job==0 ) THEN
!
!     PHASE 2 : X:=D**(-1)*X
!
         ii = 0
         DO i = 1 , N
            ii = ii + i
            X(i) = X(i)/A(ii)
         ENDDO
      ENDIF
      IF ( Job<=0 ) THEN
!
!     PHASE 3 : X:=TRANS(L)**(-1)*X
!
         ii = N*(N-1)/2
         DO i = N - 1 , 1 , -1
            ij = ii
            DO j = i + 1 , N
               ij = ij + j - 1
               X(i) = X(i) - A(ij)*X(j)
            ENDDO
            ii = ii - i
         ENDDO
      ENDIF
      END SUBROUTINE MXDPGB

! SUBROUTINE MXDPGD                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! COMPUTATION OF A DIRECTION OF NEGATIVE CURVATURE WITH RESPECT TO A
! DENSE SYMMETRIC MATRIX A USING THE FACTORIZATION A+E=L*D*TRANS(L)
!         OBTAINED BY THE SUBROUTINE MXDPGF.
!
! PARAMETERS :
!  II  N ORDER OF THE MATRIX A.
!  RI  A(N*(N+1)/2) FACTORIZATION A+E=L*D*TRANS(L) OBTAINED BY THE
!         SUBROUTINE MXDPGF.
!  RO  X(N)  COMPUTED DIRECTION OF NEGATIVE CURVATURE (I.E.
!         TRANS(X)*A*X<0) IF IT EXISTS.
!  II  INF  INFORMATION OBTAINED IN THE FACTORIZATION PROCESS. THE
!         DIRECTION OF NEGATIVE CURVATURE EXISTS ONLY IF INF>0.
!
! METHOD :
! P.E.GILL, W.MURRAY : NEWTON TYPE METHODS FOR UNCONSTRAINED AND
! LINEARLY CONSTRAINED OPTIMIZATION, MATH. PROGRAMMING 28 (1974)
! PP. 311-350.
!
      SUBROUTINE MXDPGD(N,A,X,Inf)
      IMPLICIT NONE
      INTEGER N , Inf
      DOUBLE PRECISION A(N*(N+1)/2) , X(N)
      INTEGER i , j , ii , ij
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!
!     RIGHT HAND SIDE FORMATION
!
      DO i = 1 , N
         X(i) = ZERO
      ENDDO
      IF ( Inf<=0 ) RETURN
      X(Inf) = ONE
!
!     BACK SUBSTITUTION
!
      ii = Inf*(Inf-1)/2
      DO i = Inf - 1 , 1 , -1
         ij = ii
         DO j = i + 1 , Inf
            ij = ij + j - 1
            X(i) = X(i) - A(ij)*X(j)
         ENDDO
         ii = ii - i
      ENDDO
      END SUBROUTINE MXDPGD

! SUBROUTINE MXDPGF                ALL SYSTEMS                89/12/01
! PORTABILITY : ALL SYSTEMS
! 89/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! FACTORIZATION A+E=L*D*TRANS(L) OF A DENSE SYMMETRIC POSITIVE DEFINITE
! MATRIX A+E WHERE D AND E ARE DIAGONAL POSITIVE DEFINITE MATRICES AND
! L IS A LOWER TRIANGULAR MATRIX. IF A IS SUFFICIENTLY POSITIVE
! DEFINITE THEN E=0.
!
! PARAMETERS :
!  II  N ORDER OF THE MATRIX A.
!  RU  A(N*(N+1)/2)  ON INPUT A GIVEN DENSE SYMMETRIC (USUALLY POSITIVE
!         DEFINITE) MATRIX A STORED IN THE PACKED FORM. ON OUTPUT THE
!         COMPUTED FACTORIZATION A+E=L*D*TRANS(L).
!  IO  INF  AN INFORMATION OBTAINED IN THE FACTORIZATION PROCESS. IF
!         INF=0 THEN A IS SUFFICIENTLY POSITIVE DEFINITE AND E=0. IF
!         INF<0 THEN A IS NOT SUFFICIENTLY POSITIVE DEFINITE AND E>0. IF
!         INF>0 THEN A IS INDEFINITE AND INF IS AN INDEX OF THE
!         MOST NEGATIVE DIAGONAL ELEMENT USED IN THE FACTORIZATION
!         PROCESS.
!  RU  ALF  ON INPUT A DESIRED TOLERANCE FOR POSITIVE DEFINITENESS. ON
!         OUTPUT THE MOST NEGATIVE DIAGONAL ELEMENT USED IN THE
!         FACTORIZATION PROCESS (IF INF>0).
!  RO  TAU  MAXIMUM DIAGONAL ELEMENT OF THE MATRIX E.
!
! METHOD :
! P.E.GILL, W.MURRAY : NEWTON TYPE METHODS FOR UNCONSTRAINED AND
! LINEARLY CONSTRAINED OPTIMIZATION, MATH. PROGRAMMING 28 (1974)
! PP. 311-350.
!
      SUBROUTINE MXDPGF(N,A,Inf,Alf,Tau)
      IMPLICIT NONE
      DOUBLE PRECISION Alf , Tau
      INTEGER Inf , N
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION bet , del , gam , rho , sig , tol
      INTEGER i , ij , ik , j , k , kj , kk , l
      l = 0
      Inf = 0
      tol = Alf
!
!     ESTIMATION OF THE MATRIX NORM
!
      Alf = 0.0D0
      bet = 0.0D0
      gam = 0.0D0
      Tau = 0.0D0
      kk = 0
      DO k = 1 , N
         kk = kk + k
         bet = MAX(bet,ABS(A(kk)))
         kj = kk
         DO j = k + 1 , N
            kj = kj + j - 1
            gam = MAX(gam,ABS(A(kj)))
         ENDDO
      ENDDO
      bet = MAX(tol,bet,gam/N)
!      DEL = TOL*BET
      del = tol*MAX(bet,1.0D0)
      kk = 0
      DO k = 1 , N
         kk = kk + k
!
!     DETERMINATION OF A DIAGONAL CORRECTION
!
         sig = A(kk)
         IF ( Alf>sig ) THEN
            Alf = sig
            l = k
         ENDIF
         gam = 0.0D0
         kj = kk
         DO j = k + 1 , N
            kj = kj + j - 1
            gam = MAX(gam,ABS(A(kj)))
         ENDDO
         gam = gam*gam
         rho = MAX(ABS(sig),gam/bet,del)
         IF ( Tau<rho-sig ) THEN
            Tau = rho - sig
            Inf = -1
         ENDIF
!
!     GAUSSIAN ELIMINATION
!
         A(kk) = rho
         kj = kk
         DO j = k + 1 , N
            kj = kj + j - 1
            gam = A(kj)
            A(kj) = gam/rho
            ik = kk
            ij = kj
            DO i = k + 1 , j
               ik = ik + i - 1
               ij = ij + 1
               A(ij) = A(ij) - A(ik)*gam
            ENDDO
         ENDDO
      ENDDO
      IF ( l>0 .AND. ABS(Alf)>del ) Inf = l
      END SUBROUTINE MXDPGF

! SUBROUTINE MXDPGN                ALL SYSTEMS               91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! ESTIMATION OF THE MINIMUM EIGENVALUE AND THE CORRESPONDING EIGENVECTOR
! OF A DENSE SYMMETRIC POSITIVE DEFINITE MATRIX A+E USING THE
! FACTORIZATION A+E=L*D*TRANS(L) OBTAINED BY THE SUBROUTINE MXDPGF.
!
! PARAMETERS :
!  II  N ORDER OF THE MATRIX A.
!  RI  A(N*(N+1)/2) FACTORIZATION A+E=L*D*TRANS(L) OBTAINED BY THE
!         SUBROUTINE MXDPGF.
!  RO  X(N)  ESTIMATED EIGENVECTOR.
!  RO  ALF  ESTIMATED EIGENVALUE.
!  II  JOB  OPTION. IF JOB=0 THEN ONLY ESTIMATED EIGENVALUE IS
!         COMPUTED. IS JOB>0 THEN BOTH ESTIMATED EIGENVALUE AND
!         ESTIMATED EIGENVECTOR ARE COMPUTED BY JOB ITERATIONS.
!
! SUBPROGRAMS USED :
!  S   MXDPGB  BACK SUBSTITUTION.
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!
! METHOD :
! A.K.CLINE, C.B.MOLER, G.W.STEWART, J.H.WILKINSON : AN ESTIMATE
! FOR THE CONDITION NUMBER OF A MATRIX. SIAM J. NUMER. ANAL. 16
! (1979) 368-373.
!
      SUBROUTINE MXDPGN(N,A,X,Alf,Job)
      IMPLICIT NONE
      INTEGER N , Job
      DOUBLE PRECISION A(N*(N+1)/2) , X(N) , Alf
      DOUBLE PRECISION xp , xm , fp , fm !, MXVDOT
      INTEGER i , k , ik , kk
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!
!     COMPUTATION OF THE VECTOR V WITH POSSIBLE MAXIMUM NORM SUCH
!     THAT  L*D**(1/2)*V=U  WHERE U HAS ELEMENTS +1 OR -1
!
      DO i = 1 , N
         X(i) = ZERO
      ENDDO
      kk = 0
      DO k = 1 , N
         kk = kk + k
         xp = -X(k) + ONE
         xm = -X(k) - ONE
         fp = ABS(xp)
         fm = ABS(xm)
         ik = kk
         DO i = k + 1 , N
            ik = ik + i - 1
            fp = fp + ABS(X(i)+A(ik)*xp)
            fm = fm + ABS(X(i)+A(ik)*xm)
         ENDDO
         IF ( fp>=fm ) THEN
            X(k) = xp
            ik = kk
            DO i = k + 1 , N
               ik = ik + i - 1
               X(i) = X(i) + A(ik)*xp
            ENDDO
         ELSE
            X(k) = xm
            ik = kk
            DO i = k + 1 , N
               ik = ik + i - 1
               X(i) = X(i) + A(ik)*xm
            ENDDO
         ENDIF
      ENDDO
!
!     COMPUTATION OF THE VECTOR X SUCH THAT
!     D**(1/2)*TRANS(L)*X=V
!
      fm = ZERO
      kk = 0
      DO k = 1 , N
         kk = kk + k
         IF ( Job<=0 ) THEN
            fp = SQRT(A(kk))
            X(k) = X(k)/fp
            fm = fm + X(k)*X(k)
         ELSE
            X(k) = X(k)/A(kk)
         ENDIF
      ENDDO
      fp = DBLE(N)
      IF ( Job<=0 ) THEN
!
!     FIRST ESTIMATION OF THE MINIMUM EIGENVALUE BY THE FORMULA
!     ALF=(TRANS(U)*U)/(TRANS(V)*V)
!
         Alf = fp/fm
         RETURN
      ENDIF
      fm = ZERO
      kk = N*(N+1)/2
      DO k = N , 1 , -1
         ik = kk
         DO i = k + 1 , N
            ik = ik + i - 1
            X(k) = X(k) - A(ik)*X(i)
         ENDDO
         fm = fm + X(k)*X(k)
         kk = kk - k
      ENDDO
      fm = SQRT(fm)
      IF ( Job<=1 ) THEN
!
!     SECOND ESTIMATION OF THE MINIMUM EIGENVALUE BY THE FORMULA
!     ALF=SQRT(TRANS(U)*U)/SQRT(TRANS(X)*X)
!
         Alf = SQRT(fp)/fm
      ELSE
!
!     INVERSE ITERATIONS
!
         DO k = 2 , Job
!
!     SCALING THE VECTOR X BY ITS NORM
!
            DO i = 1 , N
               X(i) = X(i)/fm
            ENDDO
            CALL MXDPGB(N,A,X,0)
            fm = SQRT(MXVDOT(N,X,X))
         ENDDO
         Alf = ONE/fm
      ENDIF
!
!     SCALING THE VECTOR X BY ITS NORM
!
      DO i = 1 , N
         X(i) = X(i)/fm
      ENDDO
      END SUBROUTINE MXDPGN

! FUNCTION MXDPGP                  ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! COMPUTATION OF THE NUMBER MXDPGP=TRANS(X)*D**(-1)*Y WHERE D IS A
! DIAGONAL MATRIX IN THE FACTORIZATION A+E=L*D*TRANS(L) OBTAINED BY THE
! SUBROUTINE MXDPGF.
!
! PARAMETERS :
!  II  N ORDER OF THE MATRIX A.
!  RI  A(N*(N+1)/2) FACTORIZATION A+E=L*D*TRANS(L) OBTAINED BY THE
!         SUBROUTINE MXDPGF.
!  RI  X(N)  INPUT VECTOR.
!  RI  Y(N)  INPUT VECTOR.
!  RR  MXDPGP  COMPUTED NUMBER MXDPGP=TRANS(X)*D**(-1)*Y.
!
      FUNCTION MXDPGP(N,A,X,Y)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION A(*) , X(*) , Y(*) , MXDPGP
      DOUBLE PRECISION temp
      INTEGER i , j
      j = 0
      temp = 0.0D0
      DO i = 1 , N
         j = j + i
         temp = temp + X(i)*Y(i)/A(j)
      ENDDO
      MXDPGP = temp
      END FUNCTION MXDPGP

! SUBROUTINE MXDPGS                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! SCALING OF A DENSE SYMMETRIC POSITIVE DEFINITE MATRIX A+E USING THE
! FACTORIZATION A+E=L*D*TRANS(L) OBTAINED BY THE SUBROUTINE MXDPGF.
!
! PARAMETERS :
!  II  N ORDER OF THE MATRIX A.
!  RU  A(N*(N+1)/2) FACTORIZATION A+E=L*D*TRANS(L) OBTAINED BY THE
!         SUBROUTINE MXDPGF.
!  RI  ALF  SCALING FACTOR.
!
      SUBROUTINE MXDPGS(N,A,Alf)
      IMPLICIT NONE
      DOUBLE PRECISION Alf
      INTEGER N
      DOUBLE PRECISION A(*)
      INTEGER i , j
      j = 0
      DO i = 1 , N
         j = j + i
         A(j) = A(j)*Alf
      ENDDO
      END SUBROUTINE MXDPGS

! SUBROUTINE MXDPGU                ALL SYSTEMS                89/12/01
! PORTABILITY : ALL SYSTEMS
! 89/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! CORRECTION OF A DENSE SYMMETRIC POSITIVE DEFINITE MATRIX A+E IN THE
! FACTORED FORM A+E=L*D*TRANS(L) OBTAINED BY THE SUBROUTINE MXDPGF.
! THE CORRECTION IS DEFINED AS A+E:=A+E+ALF*X*TRANS(X) WHERE ALF IS A
! GIVEN SCALING FACTOR AND X IS A GIVEN VECTOR.
!
! PARAMETERS :
!  II  N ORDER OF THE MATRIX A.
!  RU  A(N*(N+1)/2) FACTORIZATION A+E=L*D*TRANS(L) OBTAINED BY THE
!         SUBROUTINE MXDPGF.
!  RI  ALF  SCALING FACTOR IN THE CORRECTION TERM.
!  RI  X(N)  VECTOR IN THE CORRECTION TERM.
!  RA  Y(N) AUXILIARY VECTOR.
!
! METHOD :
! P.E.GILL, W.MURRAY, M.SAUNDERS: METHODS FOR COMPUTING AND MODIFYING
! THE LDV FACTORS OF A MATRIX, MATH. OF COMP. 29 (1974) PP. 1051-1077.
!
      SUBROUTINE MXDPGU(N,A,Alf,X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION ZERO , ONE , FOUR , CON
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,FOUR=4.0D0,CON=1.0D-8)
      DOUBLE PRECISION Alf , alfr
      INTEGER N
      DOUBLE PRECISION A(*) , X(*) , Y(*)
      DOUBLE PRECISION b , d , p , r , t , to
      INTEGER i , ii , ij , j
      IF ( Alf>=ZERO ) THEN
!
!     FORWARD CORRECTION IN CASE WHEN THE SCALING FACTOR IS NONNEGATIVE
!
         alfr = SQRT(Alf)
         CALL MXVSCL(N,alfr,X,Y)
         to = ONE
         ii = 0
         DO i = 1 , N
            ii = ii + i
            d = A(ii)
            p = Y(i)
            t = to + p*p/d
            r = to/t
            A(ii) = d/r
            b = p/(d*t)
            IF ( A(ii)<=FOUR*d ) THEN
!
!     AN EASY FORMULA FOR LIMITED DIAGONAL ELEMENT
!
               ij = ii
               DO j = i + 1 , N
                  ij = ij + j - 1
                  d = A(ij)
                  Y(j) = Y(j) - p*d
                  A(ij) = d + b*Y(j)
               ENDDO
            ELSE
!
!     A MORE COMPLICATE BUT NUMERICALLY STABLE FORMULA FOR UNLIMITED
!     DIAGONAL ELEMENT
!
               ij = ii
               DO j = i + 1 , N
                  ij = ij + j - 1
                  d = A(ij)
                  A(ij) = r*d + b*Y(j)
                  Y(j) = Y(j) - p*d
               ENDDO
            ENDIF
            to = t
         ENDDO
      ELSE
!
!     BACKWARD CORRECTION IN CASE WHEN THE SCALING FACTOR IS NEGATIVE
!
         alfr = SQRT(-Alf)
         CALL MXVSCL(N,alfr,X,Y)
         to = ONE
         ij = 0
         DO i = 1 , N
            d = Y(i)
            DO j = 1 , i - 1
               ij = ij + 1
               d = d - A(ij)*Y(j)
            ENDDO
            Y(i) = d
            ij = ij + 1
            to = to - d*d/A(ij)
         ENDDO
         IF ( to<=ZERO ) to = CON
         ii = N*(N+1)/2
         DO i = N , 1 , -1
            d = A(ii)
            p = Y(i)
            t = to + p*p/d
            A(ii) = d*to/t
            b = -p/(d*to)
            to = t
            ij = ii
            DO j = i + 1 , N
               ij = ij + j - 1
               d = A(ij)
               A(ij) = d + b*Y(j)
               Y(j) = Y(j) + p*d
            ENDDO
            ii = ii - i
         ENDDO
      ENDIF
      END SUBROUTINE MXDPGU

! SUBROUTINE MXDPRB                ALL SYSTEMS                89/12/01
! PORTABILITY : ALL SYSTEMS
! 89/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! SOLUTION OF A SYSTEM OF LINEAR EQUATIONS WITH A DENSE SYMMETRIC
! POSITIVE DEFINITE MATRIX A USING THE FACTORIZATION A=TRANS(R)*R.
!
! PARAMETERS :
!  II  N ORDER OF THE MATRIX A.
!  RI  A(N*(N+1)/2) FACTORIZATION A=TRANS(R)*R.
!  RU  X(N)  ON INPUT THE RIGHT HAND SIDE OF A SYSTEM OF LINEAR
!         EQUATIONS. ON OUTPUT THE SOLUTION OF A SYSTEM OF LINEAR
!         EQUATIONS.
!  II  JOB  OPTION. IF JOB=0 THEN X:=A**(-1)*X. IF JOB>0 THEN
!         X:=TRANS(R)**(-1)*X. IF JOB<0 THEN X:=R**(-1)*X.
!
! METHOD :
! BACK SUBSTITUTION
!
      SUBROUTINE MXDPRB(N,A,X,Job)
      IMPLICIT NONE
      INTEGER Job , N
      DOUBLE PRECISION A(*) , X(*)
      INTEGER i , ii , ij , j
      IF ( Job>=0 ) THEN
!
!     PHASE 1 : X:=TRANS(R)**(-1)*X
!
         ij = 0
         DO i = 1 , N
            DO j = 1 , i - 1
               ij = ij + 1
               X(i) = X(i) - A(ij)*X(j)
            ENDDO
            ij = ij + 1
            X(i) = X(i)/A(ij)
         ENDDO
      ENDIF
      IF ( Job<=0 ) THEN
!
!     PHASE 2 : X:=R**(-1)*X
!
         ii = N*(N+1)/2
         DO i = N , 1 , -1
            ij = ii
            DO j = i + 1 , N
               ij = ij + j - 1
               X(i) = X(i) - A(ij)*X(j)
            ENDDO
            X(i) = X(i)/A(ii)
            ii = ii - i
         ENDDO
      ENDIF
      END SUBROUTINE MXDPRB

! SUBROUTINE MXDPRC                ALL SYSTEMS                92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! CORRECTION OF A SINGULAR DENSE SYMMETRIC POSITIVE SEMIDEFINITE MATRIX
! A DECOMPOSED AS A=TRANS(R)*R.
!
! PARAMETERS :
!  II  N  ORDER OF THE MATRIX A.
!  RU  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
!  IO  INF  AN INFORMATION OBTAINED IN THE CORRECTION PROCESS. IF
!         INF=0 THEN A IS SUFFICIENTLY POSITIVE DEFINITE. IF
!         INF<0 THEN A IS NOT SUFFICIENTLY POSITIVE DEFINITE.
!         PROCESS.
!  RI  TOL  DESIRED TOLERANCE FOR POSITIVE DEFINITENESS.
!
      SUBROUTINE MXDPRC(N,A,Inf,Tol)
      IMPLICIT NONE
      INTEGER N , Inf
      DOUBLE PRECISION A(N*(N+1)/2) , Tol
      DOUBLE PRECISION tol1 , temp
      INTEGER l , i
      Inf = 0
      tol1 = SQRT(Tol)
      temp = tol1
      DO i = 1 , N*(N+1)/2
         temp = MAX(temp,ABS(A(i)))
      ENDDO
      temp = temp*tol1
      l = 0
      DO i = 1 , N
         l = l + i
         IF ( ABS(A(l))<=temp ) THEN
            A(l) = SIGN(temp,A(l))
            Inf = -1
         ENDIF
      ENDDO
      END SUBROUTINE MXDPRC

! SUBROUTINE MXDPRM                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! MULTIPLICATION OF A GIVEN VECTOR X BY A DENSE SYMMETRIC POSITIVE
! DEFINITE MATRIX A USING THE FACTORIZATION A=TRANS(R)*R.
!
! PARAMETERS :
!  II  N ORDER OF THE MATRIX A.
!  RI  A(N*(N+1)/2) FACTORIZATION A=TRANS(R)*R.
!  RU  X(N)  ON INPUT THE GIVEN VECTOR. ON OUTPUT THE RESULT OF
!         MULTIPLICATION.
!  II  JOB  OPTION. IF JOB=0 THEN X:=A*X. IF JOB>0 THEN X:=R*X.
!         IF JOB<0 THEN X:=TRANS(R)*X.
!
      SUBROUTINE MXDPRM(N,A,X,Job)
      IMPLICIT NONE
      INTEGER N , Job
      DOUBLE PRECISION A(N*(N+1)/2) , X(N)
      INTEGER i , j , ii , ij
      IF ( Job>=0 ) THEN
!
!     PHASE 1 : X:=R*X
!
         ii = 0
         DO i = 1 , N
            ii = ii + i
            X(i) = A(ii)*X(i)
            ij = ii
            DO j = i + 1 , N
               ij = ij + j - 1
               X(i) = X(i) + A(ij)*X(j)
            ENDDO
         ENDDO
      ENDIF
      IF ( Job<=0 ) THEN
!
!     PHASE 2 : X:=TRANS(R)*X
!
         ij = N*(N+1)/2
         DO i = N , 1 , -1
            X(i) = A(ij)*X(i)
            DO j = i - 1 , 1 , -1
               ij = ij - 1
               X(i) = X(i) + A(ij)*X(j)
            ENDDO
            ij = ij - 1
         ENDDO
      ENDIF
      END SUBROUTINE MXDPRM

! SUBROUTINE MXDRGR               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! PLANE ROTATION IS APPLIED TO A ROWWISE STORED DENSE RECTANGULAR
! MATRIX A.
!
! PARAMETERS :
!  II  N  NUMBER OF COLUMNS OF THE MATRIX A.
!  RU  A(M*N)  RECTANGULAR MATRIX STORED ROWWISE IN THE
!         ONE-DIMENSIONAL ARRAY.
!  II  K  FIRST INDEX OF THE PLANE ROTATION.
!  II  L  SECOND INDEX OF THE PLANE ROTATION.
!  RI  CK  DIAGONAL ELEMENT OF THE ELEMENTARY ORTHOGONAL MATRIX.
!  RI  CL  OFF-DIAGONAL ELEMENT OF THE ELEMENTARY ORTHOGONAL MATRIX.
!  II  IER  TYPE OF THE PLANE ROTATION. IER=0-GENERAL PLANE ROTATION.
!         IER=1-PERMUTATION. IER=2-TRANSFORMATION SUPPRESSED.
!
! SUBPROGRAMS USED :
!  S   MXVROT  PLANE ROTATION APPLIED TO TWO ELEMENTS.
!
      SUBROUTINE MXDRGR(N,A,K,L,Ck,Cl,Ier)
      IMPLICIT NONE
      INTEGER N , K , L , Ier
      DOUBLE PRECISION A(*) , Ck , Cl
      INTEGER i , ik , il
      IF ( Ier/=0 .AND. Ier/=1 ) RETURN
      ik = (K-1)*N
      il = (L-1)*N
      DO i = 1 , N
         ik = ik + 1
         il = il + 1
         CALL MXVROT(A(ik),A(il),Ck,Cl,Ier)
      ENDDO
      END SUBROUTINE MXDRGR

! SUBROUTINE MXDRMD               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! MULTIPLICATION OF A ROWWISE STORED DENSE RECTANGULAR MATRIX A BY
! A VECTOR X AND ADDITION OF A SCALED VECTOR ALF*Y.
!
! PARAMETERS :
!  II  N  NUMBER OF COLUMNS OF THE MATRIX A.
!  II  M  NUMBER OF ROWS OF THE MATRIX A.
!  RI  A(M*N)  RECTANGULAR MATRIX STORED ROWWISE IN THE
!         ONE-DIMENSIONAL ARRAY.
!  RI  X(N)  INPUT VECTOR.
!  RI  ALF  SCALING FACTOR.
!  RI  Y(M)  INPUT VECTOR.
!  RO  Z(M)  OUTPUT VECTOR EQUAL TO A*X+ALF*Y.
!
      SUBROUTINE MXDRMD(N,M,A,X,Alf,Y,Z)
      IMPLICIT NONE
      INTEGER N , M
      DOUBLE PRECISION A(M*N) , X(N) , Alf , Y(M) , Z(M)
      DOUBLE PRECISION temp
      INTEGER i , j , k
      k = 0
      DO j = 1 , M
         temp = Alf*Y(j)
         DO i = 1 , N
            temp = temp + A(k+i)*X(i)
         ENDDO
         Z(j) = temp
         k = k + N
      ENDDO
      END SUBROUTINE MXDRMD

! SUBROUTINE MXDRMI               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! ROWWISE STORED DENSE RECTANGULAR MATRIX A IS SET TO BE A PART OF THE
! UNIT MATRIX.
!
! PARAMETERS :
!  II  N  NUMBER OF COLUMNS OF THE MATRIX A.
!  II  M  NUMBER OF ROWS OF THE MATRIX A.
!  RO  A(M*N)  RECTANGULAR MATRIX STORED ROWWISE IN THE ONE-DIMENSIONAL
!          ARRAY. THIS MATRIX IS SET TO TRANS([I,0]).
!
      SUBROUTINE MXDRMI(N,M,A)
      IMPLICIT NONE
      INTEGER N , M
      DOUBLE PRECISION A(M*N)
      INTEGER i , j , k
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      k = 0
      DO j = 1 , M
         DO i = 1 , N
            A(i+k) = ZERO
            IF ( i==j ) A(i+k) = ONE
         ENDDO
         k = k + N
      ENDDO
      END SUBROUTINE MXDRMI

! SUBROUTINE MXDRMM               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! MULTIPLICATION OF A ROWWISE STORED DENSE RECTANGULAR MATRIX A BY
! A VECTOR X.
!
! PARAMETERS :
!  II  N  NUMBER OF COLUMNS OF THE MATRIX A.
!  II  M  NUMBER OF ROWS OF THE MATRIX A.
!  RI  A(M*N)  RECTANGULAR MATRIX STORED ROWWISE IN THE
!         ONE-DIMENSIONAL ARRAY.
!  RI  X(N)  INPUT VECTOR.
!  RO  Y(M)  OUTPUT VECTOR EQUAL TO A*X.
!
      SUBROUTINE MXDRMM(N,M,A,X,Y)
      IMPLICIT NONE
      INTEGER N , M
      DOUBLE PRECISION A(*) , X(*) , Y(*)
      DOUBLE PRECISION temp
      INTEGER i , j , k
      k = 0
      DO j = 1 , M
         temp = 0.0D0
         DO i = 1 , N
            temp = temp + A(k+i)*X(i)
         ENDDO
         Y(j) = temp
         k = k + N
      ENDDO
      END SUBROUTINE MXDRMM

! FUNCTION  MXDRMN               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! EUCLIDEAN NORM OF A PART OF THE I-TH COLUMN OF A ROWWISE STORED DENSE
! RECTANGULAR MATRIX A IS COMPUTED.
!
! PARAMETERS :
!  II  N  NUMBER OF COLUMNS OF THE MATRIX A.
!  II  M  NUMBER OF ROWS OF THE MATRIX A.
!  RI  A(M*N)  RECTANGULAR MATRIX STORED ROWWISE IN THE
!         ONE-DIMENSIONAL ARRAY.
!  II  I  INDEX OF THE COLUMN WHOSE NORM IS COMPUTED.
!  II  J  INDEX OF THE FIRST ELEMENT FROM WHICH THE NORM IS COMPUTED.
!
      FUNCTION MXDRMN(N,M,A,I,J)
      IMPLICIT NONE
      INTEGER N , M , I , J
      DOUBLE PRECISION A(M*N) , MXDRMN
      DOUBLE PRECISION pom , den
      INTEGER k , l
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      den = ZERO
      l = (J-1)*N
      DO k = J , M
         den = MAX(den,ABS(A(l+I)))
         l = l + N
      ENDDO
      pom = ZERO
      IF ( den>ZERO ) THEN
         l = (J-1)*N
         DO k = J , M
            pom = pom + (A(l+I)/den)**2
            l = l + N
         ENDDO
      ENDIF
      MXDRMN = den*SQRT(pom)
      END FUNCTION MXDRMN

! SUBROUTINE MXDRMV               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! K-TH COLUMN OF A ROWWISE STORED DENSE RECTANGULAR MATRIX A IS COPIED
! TO THE VECTOR X.
!
! PARAMETERS :
!  II  N  NUMBER OF COLUMNS OF THE MATRIX A.
!  II  M  NUMBER OF ROWS OF THE MATRIX A.
!  RI  A(M*N)  RECTANGULAR MATRIX STORED ROWWISE IN THE
!         ONE-DIMENSIONAL ARRAY.
!  RO  X(M)  OUTPUT VECTOR SUCH THAT X(J)=A(J,K) FOR ALL J.
!  II  K  INDEX OF THE ROW BEING COPIED TO THE OUTPUT VECTOR.
!
      SUBROUTINE MXDRMV(N,M,A,X,K)
      IMPLICIT NONE
      INTEGER N , M , K
      DOUBLE PRECISION A(*) , X(*)
      INTEGER i , j
      IF ( K<1 .OR. K>N ) RETURN
      i = K
      DO j = 1 , M
         X(j) = A(i)
         i = i + N
      ENDDO
      END SUBROUTINE MXDRMV

! SUBROUTINE MXDRQF               ALL SYSTEMS                92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! QR DECOMPOSITION OF ROWWISE STORED DENSE RECTANGULAR MATRIX Q USING
! HOUSEHOLDER TRANSFORMATIONS WITHOUT PIVOTING.
!
! PARAMETERS :
!  II  N  NUMBER OF COLUMNS OF THE MATRIX Q.
!  II  M  NUMBER OF ROWS OF THE MATRIX Q.
!  RU  Q(M*N)  RECTANGULAR MATRIX STORED ROWWISE IN THE
!         ONE-DIMENSIONAL ARRAY.
!  RO  R(N*(N+1)/2)  UPPER TRIANGULAR MATRIX STORED IN THE PACKED FORM.
!
! SUBPROGRAMS USED :
!  S   MXDRMN  EUCLIDEAN NORM OF A PART OF THE ROWWISE STORED
!         RECTANGULAR MATRIX COLUMN.
!
! METHOD :
! P.A.BUSSINGER, G.H.GOLUB : LINEAR LEAST SQUARES SOLUTION BY
! HOUSEHOLDER TRANSFORMATION. NUMER. MATH. 7 (1965) 269-276.
!
      SUBROUTINE MXDRQF(N,M,Q,R)
      IMPLICIT NONE
      INTEGER N , M
      DOUBLE PRECISION Q(M*N) , R(N*(N+1)/2)
      DOUBLE PRECISION alf , pom
      INTEGER i , j , k , l , jp , kp , nm
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      nm = MIN(N,M)
!
!     QR DECOMPOSITION
!
      l = 0
      kp = 0
      DO k = 1 , nm
         pom = MXDRMN(N,M,Q,k,k)
         IF ( Q(kp+k)/=ZERO ) pom = SIGN(pom,Q(kp+k))
         R(l+k) = -pom
         jp = 0
         DO j = 1 , k - 1
            R(l+j) = Q(jp+k)
            Q(jp+k) = ZERO
            jp = jp + N
         ENDDO
         IF ( pom/=ZERO ) THEN
!
!     HOUSEHOLDER TRANSFORMATION
!
            DO j = k , M
               Q(jp+k) = Q(jp+k)/pom
               jp = jp + N
            ENDDO
            Q(kp+k) = Q(kp+k) + ONE
            DO i = k + 1 , N
               alf = ZERO
               jp = kp
               DO j = k , M
                  alf = alf + Q(jp+k)*Q(jp+i)
                  jp = jp + N
               ENDDO
               alf = alf/Q(kp+k)
               jp = kp
               DO j = k , M
                  Q(jp+i) = Q(jp+i) - alf*Q(jp+k)
                  jp = jp + N
               ENDDO
            ENDDO
         ENDIF
         l = l + k
         kp = kp + N
      ENDDO
!
!     EXPLICIT FORMULATION OF THE ORTHOGONAL MATRIX
!
      kp = N*N
      DO k = N , 1 , -1
         kp = kp - N
         IF ( Q(kp+k)/=ZERO ) THEN
            DO i = k + 1 , N
               alf = ZERO
               jp = kp
               DO j = k , M
                  alf = alf + Q(jp+k)*Q(jp+i)
                  jp = jp + N
               ENDDO
               alf = alf/Q(kp+k)
               jp = kp
               DO j = k , M
                  Q(jp+i) = Q(jp+i) - alf*Q(jp+k)
                  jp = jp + N
               ENDDO
            ENDDO
            jp = kp
            DO j = k , M
               Q(jp+k) = -Q(jp+k)
               jp = jp + N
            ENDDO
         ENDIF
         Q(kp+k) = Q(kp+k) + ONE
      ENDDO
      END SUBROUTINE MXDRQF

! SUBROUTINE MXDRQU               ALL SYSTEMS                92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! UPDATE OF A QR DECOMPOSITION. THIS QR DECOMPOSITION IS UPDATED
! BY THE RULE Q*R:=Q*R+ALF*X*TRANS(Y).
!
! PARAMETERS :
!  II  N  NUMBER OF COLUMNS OF THE MATRIX Q.
!  II  M  NUMBER OF ROWS OF THE MATRIX Q.
!  RU  Q(M*N)  RECTANGULAR MATRIX STORED ROWWISE IN THE
!         ONE-DIMENSIONAL ARRAY (PART OF THE ORTHOGONAL MATRIX).
!  RU  R(N*(N+1)/2)  UPPER TRIANGULAR MATRIX STORED IN A PACKED FORM.
!  RI  ALF  SCALAR PARAMETER.
!  RI  X(M)  INPUT VECTOR.
!  RI  Y(N)  INPUT VECTOR.
!  RA  Z(N)  AUXILIARY VECTOR.
!  IO  INF  INFORMATION. IF INF=0 THEN X LIES IN THE COLUMN SPACE OF Q.
!         IF INF=1 THEN X DOES NOT LIE IN THE COLUMN SPACE OF Q.
!
! SUBPROGRAMS USED :
!  RF  MXVNOR  EUCLIDEAN NORM OF A VECTOR.
!  S   MXVORT  COMPUTATION OF THE PLANE ROTATION MATRIX.
!  S   MXVROT  PLANE ROTATION IS APPLIED TO TWO NUMBERS.
!
! METHOD :
! J.W.DANIEL, W.B.GRAGG, L.KAUFMAN, G.W.STEWARD : REORTHOGONALIZATION
! AND STABLE ALGORITHMS FOR UPDATING THE GRAM-SCHMIDT QR FACTORIZATION.
! MATHEMATICS OF COMPUTATION 30 (1976) 772-795.
      SUBROUTINE MXDRQU(N,M,Q,R,Alf,X,Y,Z,Inf)
      IMPLICIT NONE
      INTEGER N , M , Inf
      DOUBLE PRECISION Q(M*N) , R(N*(N+1)/2) , Alf , X(M) , Y(N) , Z(N)
      DOUBLE PRECISION ck , cl , zk , zl !, MXVNOR
      INTEGER j , k , l , kj , kk , ier
      DOUBLE PRECISION ONE , CON
      PARAMETER (ONE=1.0D0,CON=1.0D-10)
      Inf = 0
      kk = N*(N+1)/2
!
!     COMPUTATION OF THE VECTOR TRANS(Q)*X
!
      CALL MXDCMM(N,M,Q,X,Z)
      IF ( M>N ) THEN
!
!     IF X DOES NOT LIE IN THE COLUMN SPACE OF Q WE HAVE TO USE
!     A SUBPROBLEM WHOSE DIMENSION IS BY ONE GREATER (INF=1).
!
         zk = MXVNOR(M,X)
         CALL MXDRMD(N,M,Q,Z,-ONE,X,X)
         zl = MXVNOR(M,X)
         IF ( zl>CON*zk ) THEN
            Inf = 1
            CALL MXVSCL(M,-ONE/zl,X,X)
            CALL MXVORT(Z(N),zl,ck,cl,ier)
            IF ( ier==0 .OR. ier==1 ) THEN
               CALL MXVROT(R(kk),zl,ck,cl,ier)
               kj = N
               DO j = 1 , M
                  CALL MXVROT(Q(kj),X(j),ck,cl,ier)
                  kj = kj + N
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
!     APPLICATION OF PLANE ROTATIONS TO THE VECTOR Z SO THAT
!     TRANS(Q1)*Z=E1 WHERE Q1 IS AN ORTHOGONAL MATRIX (ACCUMULATION OF
!     THE PLANE ROTATIONS) AND E1 IS THE FIRST COLUMN OF THE UNIT
!     MATRIX. AT THE SAME TIME BOTH THE UPPER HESSENBERG MATRIX
!     TRANS(Q1)*R AND THE ORTHOGONAL MATRIX Q*Q1 ARE CONSTRUCTED SO THAT
!     Q*Q1*R1=Q*Q1*(TRANS(Q1)*R+ALF*E1*TRANS(Y)) WHERE R1 IS AN UPPER
!     HESSENBERG MATRIX.
!
      DO l = N , 2 , -1
         k = l - 1
         kk = kk - l
         CALL MXVORT(Z(k),Z(l),ck,cl,ier)
         IF ( ier==0 .OR. ier==1 ) THEN
            CALL MXVROT(R(kk),Z(l),ck,cl,ier)
            kj = kk
            DO j = l , N
               kj = kj + j - 1
               CALL MXVROT(R(kj),R(kj+1),ck,cl,ier)
            ENDDO
            kj = k
            DO j = 1 , M
               CALL MXVROT(Q(kj),Q(kj+1),ck,cl,ier)
               kj = kj + N
            ENDDO
         ENDIF
      ENDDO
      Z(1) = Alf*Z(1)
      kj = 1
      DO j = 1 , N
         R(kj) = R(kj) + Z(1)*Y(j)
         kj = kj + j
      ENDDO
!
!     APPLICATION OF PLANE ROTATIONS TO THE UPPER HESSENBERG MATRIX R1
!     GIVEN ABOVE SO THAT R2=TRANS(Q2)*R1 WHERE Q2 IS AN ORTHOGONAL
!     MATRIX (ACCUMULATION OF THE PLANE ROTATIONS) AND R2 IS AN UPPER
!     TRIANGULAR MATRIX. WE OBTAIN THE NEW QR DECOMPOSITION Q*Q1*Q2*R2.
!
      kk = 1
      DO l = 2 , N
         k = l - 1
         CALL MXVORT(R(kk),Z(l),ck,cl,ier)
         IF ( ier==0 .OR. ier==1 ) THEN
            kj = kk
            DO j = l , N
               kj = kj + j - 1
               CALL MXVROT(R(kj),R(kj+1),ck,cl,ier)
            ENDDO
            kj = k
            DO j = 1 , M
               CALL MXVROT(Q(kj),Q(kj+1),ck,cl,ier)
               kj = kj + N
            ENDDO
         ENDIF
         kk = kk + l
      ENDDO
!
!     BACK TRANSFORMATION OF THE GREATER SUBPROBLEM IF INF=1.
!
      IF ( Inf==1 ) THEN
         CALL MXVORT(R(kk),zl,ck,cl,ier)
         IF ( ier==0 .OR. ier==1 ) THEN
            kj = N
            DO j = 1 , M
               CALL MXVROT(Q(kj),X(j),ck,cl,ier)
               kj = kj + N
            ENDDO
         ENDIF
      ENDIF
      END SUBROUTINE MXDRQU

! SUBROUTINE MXDSDA                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! A DENSE SYMMETRIC MATRIX A IS AUGMENTED BY THE SCALED UNIT MATRIX
! SUCH THAT A:=A+ALF*I (I IS THE UNIT MATRIX OF ORDER N).
!
! PARAMETERS :
!  II  N  ORDER OF THE MATRIX A.
!  RU  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
!  RI  ALF  SCALING FACTOR.
!
      SUBROUTINE MXDSDA(N,A,Alf)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION A(*) , Alf
      INTEGER i , j
      j = 0
      DO i = 1 , N
         j = j + i
         A(j) = A(j) + Alf
      ENDDO
      END SUBROUTINE MXDSDA

! FUNCTION MXDSDL                  ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DETERMINATION OF THE MINIMUM DIAGONAL ELEMENT OF A DENSE SYMMETRIC
! MATRIX.
!
! PARAMETERS :
!  II  N  ORDER OF THE MATRIX A.
!  RI  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
!  IO  INF  INDEX OF THE MINIMUM DIAGONAL ELEMENT OF THE MATRIX A.
!  RR  MXDSDL  MINIMUM DIAGONAL ELEMENT OF THE MATRIX A.
!
      FUNCTION MXDSDL(N,A,Inf)
      IMPLICIT NONE
      INTEGER N , Inf
      DOUBLE PRECISION A(N*(N+1)/2) , MXDSDL
      DOUBLE PRECISION temp
      INTEGER i , j
      j = 1
      Inf = 1
      temp = A(1)
      DO i = 2 , N
         j = j + i
         IF ( temp>A(j) ) THEN
            Inf = j
            temp = A(j)
         ENDIF
      ENDDO
      MXDSDL = temp
      END FUNCTION MXDSDL

! SUBROUTINE MXDSMA             ALL SYSTEMS                 91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DENSE SYMMETRIC MATRIX AUGMENTED BY THE SCALED DENSE SYMMETRIC MATRIX.
!
! PARAMETERS :
!  II  N  ORDER OF THE MATRICES.
!  RI  ALF  SCALING FACTOR.
!  RI  A(N*(N+1)/2)  INPUT MATRIX.
!  RI  B(N*(N+1)/2)  INPUT MATRIX.
!  RO  C(N*(N+1)/2)  OUTPUT MATRIX WHERE C:=B+ALF*A.
!
      SUBROUTINE MXDSMA(N,Alf,A,B,C)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION Alf , A(N*(N+1)/2) , B(N*(N+1)/2) , C(N*(N+1)/2)
      INTEGER i
      DO i = 1 , N*(N+1)/2
         C(i) = B(i) + Alf*A(i)
      ENDDO
      END SUBROUTINE MXDSMA

! SUBROUTINE MXDSMC                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! COPYING OF A DENSE SYMMETRIC MATRIX.
!
! PARAMETERS :
!  II  N  ORDER OF THE MATRICES A AND B.
!  RI  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
!  RO  B(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM
!         WHERE B:=A.
!
      SUBROUTINE MXDSMC(N,A,B)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION A(N*(N+1)/2) , B(N*(N+1)/2)
      INTEGER m , i
      m = N*(N+1)/2
      DO i = 1 , m
         B(i) = A(i)
      ENDDO
      END SUBROUTINE MXDSMC

! SUBROUTINE MXDSMG                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
!  GERSHGORIN BOUNDS FOR EIGENVALUES OF A DENSE SYMMETRIC MATRIX.
!  AMIN .LE. ANY EIGENVALUE OF A .LE. AMAX.
!
! PARAMETERS :
!  II  N  DIMENSION OF THE MATRIX A.
!  RI  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
!  RO  AMIN  LOWER BOUND FOR EIGENVALUES OF A.
!  RO  AMAX  UPPER BOUND FOR EIGENVALUES OF A.
!
      SUBROUTINE MXDSMG(N,A,Amin,Amax)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION A(N*(N+1)/2) , Amin , Amax
      DOUBLE PRECISION temp
      INTEGER i , j , k , l
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      Amax = A(1)
      Amin = A(1)
      k = 0
      DO i = 1 , N
         temp = ZERO
         l = k
         DO j = 1 , i - 1
            l = l + 1
            temp = temp + ABS(A(l))
         ENDDO
         l = l + 1
         DO j = i + 1 , N
            l = l + j - 1
            temp = temp + ABS(A(l))
         ENDDO
         k = k + i
         Amax = MAX(Amax,A(k)+temp)
         Amin = MIN(Amin,A(k)-temp)
      ENDDO
      END SUBROUTINE MXDSMG

! SUBROUTINE MXDSMI                ALL SYSTEMS                88/12/01
! PORTABILITY : ALL SYSTEMS
! 88/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DENSE SYMMETRIC MATRIX A IS SET TO THE UNIT MATRIX WITH THE SAME
! ORDER.
!
! PARAMETERS :
!  II  N  ORDER OF THE MATRIX A.
!  RO  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM
!         WHICH IS SET TO THE UNIT MATRIX (I.E. A:=I).
!
      SUBROUTINE MXDSMI(N,A)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION A(*)
      INTEGER i , m
      m = N*(N+1)/2
      DO i = 1 , m
         A(i) = 0.0D0
      ENDDO
      m = 0
      DO i = 1 , N
         m = m + i
         A(m) = 1.0D0
      ENDDO
      END SUBROUTINE MXDSMI

! SUBROUTINE MXDSMM                ALL SYSTEMS                89/12/01
! PORTABILITY : ALL SYSTEMS
! 89/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! MULTIPLICATION OF A DENSE SYMMETRIC MATRIX A BY A VECTOR X.
!
! PARAMETERS :
!  II  N  ORDER OF THE MATRIX A.
!  RI  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
!  RI  X(N)  INPUT VECTOR.
!  RO  Y(N)  OUTPUT VECTOR EQUAL TO  A*X.
!
      SUBROUTINE MXDSMM(N,A,X,Y)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION A(*) , X(*) , Y(*)
      DOUBLE PRECISION temp
      INTEGER i , j , k , l
      k = 0
      DO i = 1 , N
         temp = 0.0D0
         l = k
         DO j = 1 , i
            l = l + 1
            temp = temp + A(l)*X(j)
         ENDDO
         DO j = i + 1 , N
            l = l + j - 1
            temp = temp + A(l)*X(j)
         ENDDO
         Y(i) = temp
         k = k + i
      ENDDO
      END SUBROUTINE MXDSMM

! FUNCTION MXDSMQ                  ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! VALUE OF A QUADRATIC FORM WITH A DENSE SYMMETRIC MATRIX A.
!
! PARAMETERS :
!  II  N  ORDER OF THE MATRIX A.
!  RI  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
!  RI  X(N)  GIVEN VECTOR.
!  RI  Y(N)  GIVEN VECTOR.
!  RR  MXDSMQ  VALUE OF THE QUADRATIC FORM MXDSMQ=TRANS(X)*A*Y.
!
      DOUBLE PRECISION FUNCTION MXDSMQ(N,A,X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER N
      DOUBLE PRECISION A(*) , X(*) , Y(*)
      DOUBLE PRECISION temp , temp1 , temp2
      INTEGER i , j , k
      temp = ZERO
      k = 0
      DO i = 1 , N
         temp1 = ZERO
         temp2 = ZERO
         DO j = 1 , i - 1
            k = k + 1
            temp1 = temp1 + A(k)*X(j)
            temp2 = temp2 + A(k)*Y(j)
         ENDDO
         k = k + 1
         temp = temp + X(i)*(temp2+A(k)*Y(i)) + Y(i)*temp1
      ENDDO
      MXDSMQ = temp
      END FUNCTION MXDSMQ

! SUBROUTINE MXDSMR               ALL SYSTEMS                92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! PLANE ROTATION IS APPLIED TO A DENSE SYMMETRIC MATRIX A. THE CASE
! K=L+1 IS REQUIRED.
!
! PARAMETERS :
!  II  N  ORDER OF THE MATRIX A.
!  RU  A(N*(N+1)/2) DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
!  II  K  FIRST INDEX OF PLANE ROTATION.
!  II  L  SECOND INDEX OF PLANE ROTATION.
!  RO  CK  DIAGONAL ELEMENT OF THE ELEMENTARY ORTHOGONAL MATRIX.
!  RO  CL  OFF-DIAGONAL ELEMENT OF THE ELEMENTARY ORTHOGONAL MATRIX.
!  IO  IER  INFORMATION ON THE TRANSFORMATION. IER<0-K OR L OUT OF
!         RANGE. IER=0-PLANE ROTATION. IER=1-PERMUTATION.
!         IER=2-TRANSFORMATION SUPPRESSED.
!
! SUBPROGRAMS USED :
!  S   MXVROT  PLANE ROTATION IS APPLIED TO TWO NUMBERS.
!
      SUBROUTINE MXDSMR(N,A,K,L,Ck,Cl,Ier)
      IMPLICIT NONE
      INTEGER N , K , L , Ier
      DOUBLE PRECISION A(*) , Ck , Cl
      DOUBLE PRECISION akk , akl , all , ckk , ckl , cll
      INTEGER j , kj , lj , kk , kl , ll
      IF ( Ier/=0 .AND. Ier/=1 ) RETURN
      IF ( K/=L+1 ) THEN
         Ier = -1
         RETURN
      ENDIF
      lj = L*(L-1)/2
      DO j = 1 , N
         IF ( j<=L ) THEN
            lj = lj + 1
            kj = lj + L
         ELSE
            lj = lj + j - 1
            kj = lj + 1
         ENDIF
         IF ( j/=K .AND. j/=L ) CALL MXVROT(A(kj),A(lj),Ck,Cl,Ier)
      ENDDO
      IF ( Ier==0 ) THEN
         ckk = Ck**2
         ckl = Ck*Cl
         cll = Cl**2
         ll = L*(L+1)/2
         kl = ll + L
         kk = ll + K
         akl = (ckl+ckl)*A(kl)
         akk = ckk*A(kk) + cll*A(ll) + akl
         all = cll*A(kk) + ckk*A(ll) - akl
         A(kl) = (cll-ckk)*A(kl) + ckl*(A(kk)-A(ll))
         A(kk) = akk
         A(ll) = all
      ELSE
         ll = L*(L+1)/2
         kk = ll + K
         akk = A(kk)
         A(kk) = A(ll)
         A(ll) = akk
      ENDIF
      END SUBROUTINE MXDSMR

! SUBROUTINE MXDSMS                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! SCALING OF A DENSE SYMMETRIC MATRIX.
!
! PARAMETERS :
!  II  N  ORDER OF THE MATRIX A.
!  RU  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM
!         WHICH IS SCALED BY THE VALUE ALF (I.E. A:=ALF*A).
!  RI  ALF  SCALING FACTOR.
!
      SUBROUTINE MXDSMS(N,A,Alf)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION A(N*(N+1)/2) , Alf
      INTEGER i , m
      m = N*(N+1)/2
      DO i = 1 , m
         A(i) = A(i)*Alf
      ENDDO
      END SUBROUTINE MXDSMS

! SUBROUTINE MXDSMU                ALL SYSTEMS                89/12/01
! PORTABILITY : ALL SYSTEMS
! 89/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! UPDATE OF A DENSE SYMMETRIC MATRIX A. THIS UPDATE IS DEFINED AS
! A:=A+ALF*X*TRANS(X) WHERE ALF IS A GIVEN SCALING FACTOR AND X IS
! A GIVEN VECTOR.
!
! PARAMETERS :
!  II  N  ORDER OF THE MATRIX A.
!  RU  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
!  RI  ALF  SCALING FACTOR IN THE CORRECTION TERM.
!  RI  X(N)  VECTOR IN THE CORRECTION TERM.
!
      SUBROUTINE MXDSMU(N,A,Alf,X)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION A(N*(N+1)/2) , X(N) , Alf
      DOUBLE PRECISION temp
      INTEGER i , j , k
      k = 0
      DO i = 1 , N
         temp = Alf*X(i)
         DO j = 1 , i
            k = k + 1
            A(k) = A(k) + temp*X(j)
         ENDDO
      ENDDO
      END SUBROUTINE MXDSMU

! SUBROUTINE MXDSMV                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! K-TH ROW OF A DENSE SYMMETRIC MATRIX A IS COPIED TO THE VECTOR X.
!
! PARAMETERS :
!  II  N  ORDER OF THE MATRIX A.
!  RI  A(N*(N+1)/2)  DENSE SYMMETRIC MATRIX STORED IN THE PACKED FORM.
!  RO  X(N)  OUTPUT VECTOR.
!  II  K  INDEX OF COPIED ROW.
!
      SUBROUTINE MXDSMV(N,A,X,K)
      IMPLICIT NONE
      INTEGER K , N
      DOUBLE PRECISION A(*) , X(*)
      INTEGER i , l
      l = K*(K-1)/2
      DO i = 1 , N
         IF ( i<=K ) THEN
            l = l + 1
         ELSE
            l = l + i - 1
         ENDIF
         X(i) = A(l)
      ENDDO
      END SUBROUTINE MXDSMV

! SUBROUTINE MXVCOP                ALL SYSTEMS                88/12/01
! PORTABILITY : ALL SYSTEMS
! 88/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! COPYING OF A VECTOR.
!
! PARAMETERS :
!  II  N  VECTOR DIMENSION.
!  RI  X(N)  INPUT VECTOR.
!  RO  Y(N)  OUTPUT VECTOR WHERE Y:= X.
!
      SUBROUTINE MXVCOP(N,X,Y)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(*) , Y(*)
      INTEGER i
      DO i = 1 , N
         Y(i) = X(i)
      ENDDO
      END SUBROUTINE MXVCOP

! SUBROUTINE MXVDIF                ALL SYSTEMS                88/12/01
! PORTABILITY : ALL SYSTEMS
! 88/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! VECTOR DIFFERENCE.
!
! PARAMETERS :
!  RI  X(N)  INPUT VECTOR.
!  RI  Y(N)  INPUT VECTOR.
!  RO  Z(N)  OUTPUT VECTOR WHERE Z:= X - Y.
!
      SUBROUTINE MXVDIF(N,X,Y,Z)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(*) , Y(*) , Z(*)
      INTEGER i
      DO i = 1 , N
         Z(i) = X(i) - Y(i)
      ENDDO
      END SUBROUTINE MXVDIF

! SUBROUTINE MXVDIR                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! VECTOR AUGMENTED BY THE SCALED VECTOR.
!
! PARAMETERS :
!  II  N  VECTOR DIMENSION.
!  RI  A  SCALING FACTOR.
!  RI  X(N)  INPUT VECTOR.
!  RI  Y(N)  INPUT VECTOR.
!  RO  Z(N)  OUTPUT VECTOR WHERE Z:= Y + A*X.
!
      SUBROUTINE MXVDIR(N,A,X,Y,Z)
      IMPLICIT NONE
      DOUBLE PRECISION A
      INTEGER N
      DOUBLE PRECISION X(*) , Y(*) , Z(*)
      INTEGER i
      DO i = 1 , N
         Z(i) = Y(i) + A*X(i)
      ENDDO
      END SUBROUTINE MXVDIR

! FUNCTION MXVDOT                  ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DOT PRODUCT OF TWO VECTORS.
!
! PARAMETERS :
!  II  N  VECTOR DIMENSION.
!  RI  X(N)  INPUT VECTOR.
!  RI  Y(N)  INPUT VECTOR.
!  RR  MXVDOT  VALUE OF DOT PRODUCT MXVDOT=TRANS(X)*Y.
!
      FUNCTION MXVDOT(N,X,Y)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(*) , Y(*) , MXVDOT
      DOUBLE PRECISION temp
      INTEGER i
      temp = 0.0D0
      DO i = 1 , N
         temp = temp + X(i)*Y(i)
      ENDDO
      MXVDOT = temp
      END FUNCTION MXVDOT

! SUBROUTINE MXVINA             ALL SYSTEMS                   90/12/01
! PORTABILITY : ALL SYSTEMS
! 90/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! ELEMENTS OF THE INTEGER VECTOR ARE REPLACED BY THEIR ABSOLUTE VALUES.
!
! PARAMETERS :
!  II  N DIMENSION OF THE INTEGER VECTOR.
!  IU  IX(N)  INTEGER VECTOR WHICH IS UPDATED SO THAT IX(I):=ABS(IX(I))
!         FOR ALL I.
!
      SUBROUTINE MXVINA(N,Ix)
      IMPLICIT NONE
      INTEGER N
      INTEGER Ix(*)
      INTEGER i
      DO i = 1 , N
         Ix(i) = ABS(Ix(i))
         IF ( Ix(i)>10 ) Ix(i) = Ix(i) - 10
      ENDDO
      END SUBROUTINE MXVINA

! SUBROUTINE MXVIND               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! CHANGE OF THE INTEGER VECTOR ELEMENT FOR THE CONSTRAINT ADDITION.
!
! PARAMETERS :
!  IU  IX(N)  INTEGER VECTOR.
!  II  I  INDEX OF THE CHANGED ELEMENT.
!  II JOB  CHANGE SPECIFICATION. IS JOB.EQ.0 THEN IX(I)=10-IX(I).
!
      SUBROUTINE MXVIND(Ix,I,Job)
      IMPLICIT NONE
      INTEGER Ix(*) , I , Job
      IF ( Job==0 ) Ix(I) = 10 - Ix(I)
      END SUBROUTINE MXVIND

! SUBROUTINE MXVINS             ALL SYSTEMS                   90/12/01
! PORTABILITY : ALL SYSTEMS
! 90/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! INITIATION OF THE INTEGER VECTOR.
!
! PARAMETERS :
!  II  N DIMENSION OF THE INTEGER VECTOR.
!  II  IP  INTEGER PARAMETER.
!  IO  IX(N)  INTEGER VECTOR SUCH THAT IX(I)=IP FOR ALL I.
!
      SUBROUTINE MXVINS(N,Ip,Ix)
      IMPLICIT NONE
      INTEGER N , Ip , Ix(N)
      INTEGER i
      DO i = 1 , N
         Ix(i) = Ip
      ENDDO
      END SUBROUTINE MXVINS

! SUBROUTINE MXVINV               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! CHANGE OF THE INTEGER VECTOR ELEMENT FOR THE CONSTRAINT ADDITION.
!
! PARAMETERS :
!  II  N  VECTOR DIMENSION.
!  IU  IX(N)  INTEGER VECTOR.
!  II  I  INDEX OF THE CHANGED ELEMENT.
!  II  JOB  CHANGE SPECIFICATION
!
      SUBROUTINE MXVINV(Ix,I,Job)
      IMPLICIT NONE
      INTEGER I , Job
      INTEGER Ix(*)
      IF ( (Ix(I)==3 .OR. Ix(I)==5) .AND. Job<0 ) Ix(I) = Ix(I) + 1
      IF ( (Ix(I)==4 .OR. Ix(I)==6) .AND. Job>0 ) Ix(I) = Ix(I) - 1
      Ix(I) = -Ix(I)
      END SUBROUTINE MXVINV

! FUNCTION MXVMAX               ALL SYSTEMS                   91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! L-INFINITY NORM OF A VECTOR.
!
! PARAMETERS :
!  II  N  VECTOR DIMENSION.
!  RI  X(N)  INPUT VECTOR.
!  RR  MXVMAX  L-INFINITY NORM OF THE VECTOR X.
!
      FUNCTION MXVMAX(N,X)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(*) , MXVMAX
      INTEGER i
      MXVMAX = 0.0D0
      DO i = 1 , N
         MXVMAX = MAX(MXVMAX,ABS(X(i)))
      ENDDO
      END FUNCTION MXVMAX

! SUBROUTINE MXVNEG                ALL SYSTEMS                88/12/01
! PORTABILITY : ALL SYSTEMS
! 88/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! CHANGE THE SIGNS OF VECTOR ELEMENTS.
!
! PARAMETERS :
!  II  N  VECTOR DIMENSION.
!  RI  X(N)  INPUT VECTOR.
!  RO  Y(N)  OUTPUT VECTOR WHERE Y:= - X.
!
      SUBROUTINE MXVNEG(N,X,Y)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(*) , Y(*)
      INTEGER i
      DO i = 1 , N
         Y(i) = -X(i)
      ENDDO
      END SUBROUTINE MXVNEG

! FUNCTION  MXVNOR               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! EUCLIDEAN NORM OF A VECTOR.
!
! PARAMETERS :
!  II  N  VECTOR DIMENSION.
!  RI  X(N)  INPUT VECTOR.
!  RR  MXVNOR  EUCLIDEAN NORM OF X.
!
      FUNCTION MXVNOR(N,X)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(*) , MXVNOR
      DOUBLE PRECISION pom , den
      INTEGER i
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      den = ZERO
      DO i = 1 , N
         den = MAX(den,ABS(X(i)))
      ENDDO
      pom = ZERO
      IF ( den>ZERO ) THEN
         DO i = 1 , N
            pom = pom + (X(i)/den)**2
         ENDDO
      ENDIF
      MXVNOR = den*SQRT(pom)
      END FUNCTION MXVNOR

! SUBROUTINE MXVORT               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DETERMINATION OF AN ELEMENTARY ORTHOGONAL MATRIX FOR PLANE ROTATION.
!
! PARAMETERS :
!  RU  XK  FIRST VALUE FOR PLANE ROTATION (XK IS TRANSFORMED TO
!         SQRT(XK**2+XL**2))
!  RU  XL  SECOND VALUE FOR PLANE ROTATION (XL IS TRANSFORMED TO
!         ZERO)
!  RO  CK  DIAGONAL ELEMENT OF THE ELEMENTARY ORTHOGONAL MATRIX.
!  RO  CL  OFF-DIAGONAL ELEMENT OF THE ELEMENTARY ORTHOGONAL MATRIX.
!  IO  IER  INFORMATION ON THE TRANSFORMATION. IER=0-GENERAL PLANE
!         ROTATION. IER=1-PERMUTATION. IER=2-TRANSFORMATION SUPPRESSED.
!
      SUBROUTINE MXVORT(Xk,Xl,Ck,Cl,Ier)
      IMPLICIT NONE
      DOUBLE PRECISION Ck , Cl , Xk , Xl
      INTEGER Ier
      DOUBLE PRECISION den , pom
      IF ( Xl==0.0D0 ) THEN
         Ier = 2
      ELSEIF ( Xk==0.0D0 ) THEN
         Xk = Xl
         Xl = 0.0D0
         Ier = 1
      ELSE
         IF ( ABS(Xk)>=ABS(Xl) ) THEN
            pom = Xl/Xk
            den = SQRT(1.0D0+pom*pom)
            Ck = 1.0D0/den
            Cl = pom/den
            Xk = Xk*den
         ELSE
            pom = Xk/Xl
            den = SQRT(1.0D0+pom*pom)
            Cl = 1.0D0/den
            Ck = pom/den
            Xk = Xl*den
         ENDIF
         Xl = 0.0D0
         Ier = 0
      ENDIF
      END SUBROUTINE MXVORT

! SUBROUTINE MXVROT               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! PLANE ROTATION IS APPLIED TO TWO VALUES.
!
! PARAMETERS :
!  RU  XK  FIRST VALUE FOR PLANE ROTATION.
!  RU  XL  SECOND VALUE FOR PLANE ROTATION.
!  RI  CK  DIAGONAL ELEMENT OF THE ELEMENTARY ORTHOGONAL MATRIX.
!  RI  CL  OFF-DIAGONAL ELEMENT OF THE ELEMENTARY ORTHOGONAL MATRIX.
!  II  IER  INFORMATION ON THE TRANSFORMATION. IER=0-GENERAL PLANE
!         ROTATION. IER=1-PERMUTATION. IER=2-TRANSFORMATION SUPPRESSED.
!
      SUBROUTINE MXVROT(Xk,Xl,Ck,Cl,Ier)
      IMPLICIT NONE
      DOUBLE PRECISION Ck , Cl , Xk , Xl
      INTEGER Ier
      DOUBLE PRECISION yk , yl
      IF ( Ier==0 ) THEN
         yk = Xk
         yl = Xl
         Xk = Ck*yk + Cl*yl
         Xl = Cl*yk - Ck*yl
      ELSEIF ( Ier==1 ) THEN
         yk = Xk
         Xk = Xl
         Xl = yk
      ENDIF
      END SUBROUTINE MXVROT

! SUBROUTINE MXVSAV                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DIFFERENCE OF TWO VECTORS RETURNED IN THE SUBSTRACTED ONE.
!
! PARAMETERS :
!  II  N  VECTOR DIMENSION.
!  RI  X(N)  INPUT VECTOR.
!  RU  Y(N)  UPDATE VECTOR WHERE Y:= X - Y.
!
      SUBROUTINE MXVSAV(N,X,Y)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(*) , Y(*)
      DOUBLE PRECISION temp
      INTEGER i
      DO i = 1 , N
         temp = Y(i)
         Y(i) = X(i) - Y(i)
         X(i) = temp
      ENDDO
      END SUBROUTINE MXVSAV

! SUBROUTINE MXVSCL                ALL SYSTEMS                88/12/01
! PORTABILITY : ALL SYSTEMS
! 88/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! SCALING OF A VECTOR.
!
! PARAMETERS :
!  II  N  VECTOR DIMENSION.
!  RI  X(N)  INPUT VECTOR.
!  RI  A  SCALING FACTOR.
!  RO  Y(N)  OUTPUT VECTOR WHERE Y:= A*X.
!
      SUBROUTINE MXVSCL(N,A,X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION A
      INTEGER N
      DOUBLE PRECISION X(*) , Y(*)
      INTEGER i
      DO i = 1 , N
         Y(i) = A*X(i)
      ENDDO
      END SUBROUTINE MXVSCL

! SUBROUTINE MXVSET                ALL SYSTEMS                88/12/01
! PORTABILITY : ALL SYSTEMS
! 88/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! A SCALAR IS SET TO ALL THE ELEMENTS OF A VECTOR.
!
! PARAMETERS :
!  II  N  VECTOR DIMENSION.
!  RI  A  INITIAL VALUE.
!  RO  X(N)  OUTPUT VECTOR SUCH THAT X(I)=A FOR ALL I.
!
      SUBROUTINE MXVSET(N,A,X)
      IMPLICIT NONE
      DOUBLE PRECISION A
      INTEGER N
      DOUBLE PRECISION X(*)
      INTEGER i
      DO i = 1 , N
         X(i) = A
      ENDDO
      END SUBROUTINE MXVSET

! SUBROUTINE MXVSUM                ALL SYSTEMS                88/12/01
! PORTABILITY : ALL SYSTEMS
! 88/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! SUM OF TWO VECTORS.
!
! PARAMETERS :
!  II  N  VECTOR DIMENSION.
!  RI  X(N)  INPUT VECTOR.
!  RI  Y(N)  INPUT VECTOR.
!  RO  Z(N)  OUTPUT VECTOR WHERE Z:= X + Y.
!
      SUBROUTINE MXVSUM(N,X,Y,Z)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(*) , Y(*) , Z(*)
      INTEGER i
      DO i = 1 , N
         Z(i) = X(i) + Y(i)
      ENDDO
      END SUBROUTINE MXVSUM

    end module matrix_routines
