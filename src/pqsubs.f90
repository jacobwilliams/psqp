
! SUBROUTINE PA0GS1             ALL SYSTEMS                 97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE:
! NUMERICAL COMPUTATION OF THE GRADIENT OF THE MODEL FUNCTION.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  II  KA  INDEF OF THE APPROXIMATED FUNCTION.
!  RI  X(N)  VECTOR OF VARIABLES.
!  RO  GA(N)  GRADIENT OF THE APPROXIMATED FUNCTION.
!  RI  FA  VALUE OF THE APPROXIMATED FUNCTION.
!  RI  ETA1  PRECISION OF THE COMPUTED VALUES.
!  IU  NAV  NUMBER OF APPROXIMATED FUNCTION EVALUATIONS.
!
      SUBROUTINE PA0GS1(N,Ka,X,Ga,Fa,Eta1,Nav)
      IMPLICIT NONE
      INTEGER N , Ka , Nav
      DOUBLE PRECISION X(*) , Ga(*) , Fa , Eta1
      DOUBLE PRECISION xstep , xtemp , ftemp , eta
      INTEGER ivar
      eta = SQRT(Eta1)
      ftemp = Fa
      DO ivar = 1 , N
!
!     STEP SELECTION
!
         xstep = 1.0D0
         xstep = eta*MAX(ABS(X(ivar)),xstep)*SIGN(1.0D0,X(ivar))
         xtemp = X(ivar)
         X(ivar) = X(ivar) + xstep
         xstep = X(ivar) - xtemp
         Nav = Nav + 1
         CALL FUN(N,Ka,X,Fa)
!
!     NUMERICAL DIFFERENTIATION
!
         Ga(ivar) = (Fa-ftemp)/xstep
         X(ivar) = xtemp
      ENDDO
      Fa = ftemp
      END

! SUBROUTINE PA1SQ1             ALL SYSTEMS                 97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! COMPUTATION OF THE VALUE AND THE GRADIENT OF THE OBJECTIVE FUNCTION
! WHICH IS DEFINED AS A SUM OF SQUARES.
!
! PARAMETERS:
!  II  N  NUMBER OF VARIABLES.
!  RI  X(N)  VECTOR OF VARIABLES.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!  RO  AF(N)  VALUES OF THE APPROXIMATED FUNCTIONS.
!  RI  GA(NF)  GRADIENT OF THE APPROXIMATED FUNCTION.
!  RI  AG(N*N)  RECTANGULAR MATRIX WHICH IS USED FOR THE DIRECTION
!         VECTOR DETERMINATION.
!  RO  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RI  ETA1  PRECISION OF THE COMPUTES FUNCTION VALUES.
!  II  KDA  DEGREE OF COMPUTED DERIVATIVES.
!  II  KD  DEGREE OF REQUIRED DERVATIVES.
!  IO  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!  IU  NFV  NUMBER OF OBJECTIVE FUNCTION VALUES COMPUTED.
!  IU  NFG  NUMBER OF OBJECTIVE FUNCTION GRADIENTS COMPUTED.
!
! SUBPROGRAMS USED :
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
!  S   MXVSET  INITIATION OF A VECTOR.
!
      SUBROUTINE PA1SQ1(N,X,F,Af,Ga,Ag,G,Eta1,Kda,Kd,Ld,Nfv,Nfg)
      IMPLICIT NONE
      INTEGER N , Kda , Kd , Ld , Nfv , Nfg
      DOUBLE PRECISION X(*) , F , Af(*) , Ga(*) , Ag(*) , G(*) , Eta1
      INTEGER ka , nav
      DOUBLE PRECISION fa
      IF ( Kd<=Ld ) RETURN
      IF ( Kd>=0 .AND. Ld<0 ) THEN
         F = 0.0D0
         Nfv = Nfv + 1
      ENDIF
      IF ( Kd>=1 .AND. Ld<1 ) THEN
         CALL MXVSET(N,0.0D0,G)
         IF ( Kda>0 ) Nfg = Nfg + 1
      ENDIF
      nav = 0
      DO ka = 1 , N
         IF ( Kd>=0 ) THEN
            IF ( Ld>=0 ) THEN
               fa = Af(ka)
               GOTO 20
            ELSE
               CALL FUN(N,ka,X,fa)
               Af(ka) = fa
            ENDIF
            IF ( Ld<0 ) F = F + fa*fa
 20         IF ( Kd>=1 ) THEN
               IF ( Kda>0 ) THEN
                  CALL DFUN(N,ka,X,Ga)
               ELSE
                  CALL PA0GS1(N,ka,X,Ga,fa,Eta1,nav)
               ENDIF
               CALL MXVDIR(N,fa,Ga,G,G)
               CALL MXVCOP(N,Ga,Ag((ka-1)*N+1))
            ENDIF
         ENDIF
      ENDDO
      Nfv = Nfv + nav/N
      IF ( Kd>=0 .AND. Ld<0 ) F = 0.5D0*F
      Ld = Kd
      END

! SUBROUTINE PA2SQ1             ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
!  COMPUTATION OF THE VALUE AND THE GRADIENT AND THE HESSIAN MATRIX
!  OF THE OBJECTIVE FUNCTION WHICH IS DEFINED AS A SUM OF SQUARES.
!
! PARAMETERS:
!  II  NF  NUMBER OF VARIABLES.
!  II  NA  NUMBER OF APPROXIMATED FUNCTIONS.
!  RI  GA(NF)  GRADIENT OF THE APPROXIMATED FUNCTION.
!  RO  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RI  HA(NF*(NF+1)/2)  HESSIAN MATRIX OF THE APPROXIMATED FUNCTION.
!  RO  H(NF*(NF+1)/2)  HESSIAN MATRIX OF THE OBJECTIVE FUNCTION.
!  RI  FA  VALUE OF THE SELECTED FUNCTION.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!
! COMMON DATA :
!  II  KDA  DEGREE OF ANALYTICALLY COMPUTED DERIVATIVES.
!  II  KD  DEGREE OF REQUIRED DERVATIVES.
!  IO  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!  IU  NAV  NUMBER OF APPROXIMATED FUNCTION VALUES COMPUTED.
!  IU  NAG  NUMBER OF APPROXIMATED FUNCTION GRADIENTS COMPUTED.
!  IU  NAH  NUMBER OF APPROXIMATED FUNCTION HESSIAN MATRICES COMPUTED.
!  IU  NFV  NUMBER OF OBJECTIVE FUNCTION VALUES COMPUTED.
!  IU  NFG  NUMBER OF OBJECTIVE FUNCTION GRADIENTS COMPUTED.
!  IU  NFH  NUMBER OF OBJECTIVE FUNCTION HESSIAN MATRICES COMPUTED.
!  IU  IDECF  DECOMPOSITION INDICATOR. IDECF=0-NO DECOMPOSITION.
!
! STATUS VARIABLES :
!  NS,ISP,TSS
!
! SUBPROGRAMS USED :
!  S  UYPRO1  PROLOGUE.
!  S  UYEPI1  EPILOGUE.
!  S  UYSET0  STATUS DEFINITION.
!  S  MXVSET  INITIATION OF A VECTOR.
!  S  MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
!  S  MXDSMO  INITIATION OF A DENSE SYMMETRIC MATRIX.
!  S  MXDSMA  DENSE SYMMETRIC MATRIX AUGMENTED BY THE SCALED DENSE
!         SYMMETRIC MATRIX.
!  S  MXDSMU  CORRECTION OF A DENSE SYMMETRIC MATRIX.
!
      SUBROUTINE PA2SQ1(Nf,Na,X,F,Af,Ga,G,H,Eta1,Kda,Kd,Ld,Nfv,Nfg)
      IMPLICIT NONE
      INTEGER Nf , Na , Kda , Kd , Ld , Nfv , Nfg
      DOUBLE PRECISION X(*) , F , Af(*) , Ga(*) , G(*) , H(*) , Eta1
      INTEGER ka , nav
      DOUBLE PRECISION fa
      IF ( Kd<=Ld ) RETURN
      IF ( Kd>=0 .AND. Ld<0 ) THEN
         F = 0.0D0
         Nfv = Nfv + 1
      ENDIF
      IF ( Kd>=1 .AND. Ld<1 ) THEN
         CALL MXVSET(Nf,0.0D0,G)
         IF ( Kda>0 ) Nfg = Nfg + 1
      ENDIF
      IF ( Kd>=2 .AND. Ld<2 ) CALL MXVSET(Nf*(Nf+1)/2,0.0D0,H)
      nav = 0
      DO ka = 1 , Na
         IF ( Kd>=0 ) THEN
            IF ( Ld>=0 ) THEN
               fa = Af(ka)
               GOTO 20
            ELSE
               CALL FUN(Nf,ka,X,fa)
               Af(ka) = fa
            ENDIF
            IF ( Ld<0 ) F = F + fa*fa
 20         IF ( Kd>=1 ) THEN
               IF ( Kda>0 ) THEN
                  CALL DFUN(Nf,ka,X,Ga)
               ELSE
                  CALL PA0GS1(Nf,ka,X,Ga,fa,Eta1,nav)
               ENDIF
               CALL MXVDIR(Nf,fa,Ga,G,G)
               IF ( Kd>=2 ) CALL MXDSMU(Nf,H,1.0D0,Ga)
            ENDIF
         ENDIF
      ENDDO
      Nfv = Nfv + nav/Na
      IF ( Kd>=0 .AND. Ld<0 ) F = 0.5D0*F
      Ld = Kd
      END

! SUBROUTINE PC1F01             ALL SYSTEMS                 97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! COMPUTATION OF THE VALUE AND THE GRADIENT OF THE CONSTRAINT FUNCTION.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  II  NC  NUMBER OF CONSTRAINTS.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  RI  FC  VALUE OF THE SELECTED CONSTRAINT FUNCTION.
!  RI  CF(NC)  VECTOR CONTAINING VALUES OF CONSTRAINT FUNCTIONS.
!  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  RI  GC(NF)  GRADIENT OF THE SELECTED CONSTRAINT FUNCTION.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE GRADIENTS OF CONSTRAINT
!         FUNCTIONS.
!  RO  CMAX  MAXIMUM CONSTRAINT VIOLATION.
!  II  KD  DEGREE OF REQUIRED DERVATIVES.
!  II  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!
! SUBPROGRAMS USED :
!  S   MXVCOP  COPYING OF A VECTOR.
!
      SUBROUTINE PC1F01(Nf,Nc,X,Fc,Cf,Cl,Cu,Ic,Gc,Cg,Cmax,Kd,Ld)
      IMPLICIT NONE
      DOUBLE PRECISION Fc , Cmax
      INTEGER Kd , Ld , Nc , Nf
      DOUBLE PRECISION Cf(*) , Cg(*) , Cl(*) , Cu(*) , Gc(*) , X(*)
      INTEGER Ic(*)
      DOUBLE PRECISION pom , temp
      INTEGER kc
      IF ( Kd<=Ld ) RETURN
      IF ( Ld<0 ) Cmax = 0.0D0
      DO kc = 1 , Nc
         IF ( Kd>=0 ) THEN
            IF ( Ld>=0 ) THEN
               Fc = Cf(kc)
               GOTO 20
            ELSE
               CALL CON(Nf,kc,X,Fc)
               Cf(kc) = Fc
            ENDIF
            IF ( Ic(kc)>0 ) THEN
               pom = 0.0D0
               temp = Cf(kc)
               IF ( Ic(kc)==1 .OR. Ic(kc)>=3 )                          &
                    pom = MIN(pom,temp-Cl(kc))
               IF ( Ic(kc)==2 .OR. Ic(kc)>=3 )                          &
                    pom = MIN(pom,Cu(kc)-temp)
               IF ( pom<0.0D0 ) Cmax = MAX(Cmax,-pom)
            ENDIF
 20         IF ( Kd>=1 ) THEN
               IF ( Ld>=1 ) THEN
                  CALL MXVCOP(Nf,Cg((kc-1)*Nf+1),Gc)
               ELSE
                  CALL DCON(Nf,kc,X,Gc)
                  CALL MXVCOP(Nf,Gc,Cg((kc-1)*Nf+1))
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      Ld = Kd
      END

! SUBROUTINE PF1F01                ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! COMPUTATION OF THE VALUE AND THE GRADIENT OF THE OBJECTIVE FUNCTION.
!
! PARAMETERS:
!  II  NF  NUMBER OF VARIABLES.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  RI  GF(NF)  GRADIENT OF THE MODEL FUNCTION.
!  RO  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RI  FF  VALUE OF THE MODEL FUNCTION.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!  II  KD  DEGREE OF REQUIRED DERIVATIVES.
!  II  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!  II  IEXT  TYPE OF EXTREMUM. IEXT=0-MINIMUM. IEXT=1-MAXIMUM.
!
! SUBPROGRAMS USED :
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVNEG  COPYING OF A VECTOR WITH CHANGE OF THE SIGN.
!
      SUBROUTINE PF1F01(Nf,X,Gf,G,Ff,F,Kd,Ld,Iext)
      IMPLICIT NONE
      DOUBLE PRECISION F , Ff
      INTEGER Iext , Kd , Ld , Nf
      DOUBLE PRECISION Gf(*) , G(*) , X(*)
      INTEGER NADd , NDEc , NFG , NFH , NFV , NIT , NREm , NREs
      COMMON /STAT  / NREs , NDEc , NREm , NADd , NIT , NFV , NFG , NFH
      IF ( Kd<=Ld ) RETURN
      IF ( Ld<0 ) THEN
         NFV = NFV + 1
         CALL OBJ(Nf,X,Ff)
         IF ( Iext<=0 ) THEN
            F = Ff
         ELSE
            F = -Ff
         ENDIF
      ENDIF
      IF ( Kd>=1 ) THEN
         IF ( Ld<1 ) THEN
            NFG = NFG + 1
            CALL DOBJ(Nf,X,Gf)
            IF ( Iext>0 ) CALL MXVNEG(Nf,Gf,G)
         ENDIF
      ENDIF
      Ld = Kd
      END

! SUBROUTINE PLADB0               ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! NEW LINEAR CONSTRAINT OR A NEW SIMPLE BOUND IS ADDED TO THE
! ACTIVE SET.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  IU  N  ACTUAL NUMBER OF VARIABLES.
!  IU  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RU  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RU  CZ(NF*NF)  MATRIX WHOSE COLUMNS ARE BASIC VECTORS FROM THE
!         CURRENT REDUCED SUBSPACE.
!  RA  S(NF)  AUXILIARY VECTOR.
!  RI  EPS7  TOLERANCE FOR LINEAR INDEPENDENCE OF CONSTRAINTS.
!  RO  GMAX  MAXIMUM ABSOLUTE VALUE OF A PARTIAL DERIVATIVE.
!  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF A NEGATIVE LAGRANGE MULTIPLIER.
!  II  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
!  IU  NADD  NUMBER OF CONSTRAINT ADDITIONS.
!  IO  IER  ERROR INDICATOR.
!
! SUBPROGRAMS USED :
!  S   PLADR0  CORRECTION OF KERNEL OF THE ORTHOGONAL PROJECTION
!         AFTER CONSTRAINT ADDITION.
!  S   MXDRMM  PREMULTIPLICATION OF A VECTOR BY A ROWWISE STORED DENSE
!         RECTANGULAR MATRIX.
!  S   MXDRMV  COPY OF THE SELECTED COLUMN OF A ROWWISE STORED DENSE
!         RECTANGULAR MATRIX.
!  S   MXDRGR  PLANE ROTATION OF A TRANSPOSED DENSE RECTANGULAR MATRIX.
!  S   MXVORT  DETERMINATION OF AN ELEMENTARY ORTHOGONAL MATRIX FOR
!         PLANE ROTATION.
!
      SUBROUTINE PLADB0(Nf,N,Ica,Cg,Cr,Cz,S,Eps7,Gmax,Umax,Inew,Nadd,   &
                        Ier)
      IMPLICIT NONE
      INTEGER Nf , N , Ica(*) , Inew , Nadd , Ier
      DOUBLE PRECISION Cg(*) , Cr(*) , Cz(*) , S(*) , Eps7 , Gmax , Umax
      DOUBLE PRECISION ck , cl
      INTEGER k , l , n1
      CALL PLADR0(Nf,N,Ica,Cg,Cr,S,Eps7,Gmax,Umax,Inew,Nadd,Ier)
      IF ( Ier/=0 ) RETURN
      IF ( N>0 ) THEN
         n1 = N + 1
         IF ( Inew>0 ) THEN
            CALL MXDRMM(Nf,n1,Cz,Cg((Inew-1)*Nf+1),S)
         ELSE
            CALL MXDRMV(Nf,n1,Cz,S,-Inew)
         ENDIF
         DO l = 1 , N
            k = l + 1
            CALL MXVORT(S(k),S(l),ck,cl,Ier)
            CALL MXDRGR(Nf,Cz,k,l,ck,cl,Ier)
            IF ( Ier<0 ) RETURN
         ENDDO
      ENDIF
      Ier = 0
      END

! SUBROUTINE PLADB4               ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! NEW LINEAR CONSTRAINT OR A NEW SIMPLE BOUND IS ADDED TO THE ACTIVE
! SET. TRANSFORMED HESSIAN MATRIX APPROXIMATION OR ITS INVERSION
! IS UPDATED.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  IU  N  ACTUAL NUMBER OF VARIABLES.
!  IU  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RU  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RU  CZ(NF*NF)  MATRIX WHOSE COLUMNS ARE BASIC VECTORS FROM THE
!         CURRENT REDUCED SUBSPACE.
!  RU  H(NF*(NF+1)/2)  TRANSFORMED HESSIAN MATRIX APPROXIMATION OR
!         ITS INVERSION.
!  RA  S(NF)  AUXILIARY VECTOR.
!  RI  EPS7  TOLERANCE FOR LINEAR INDEPENDENCE OF CONSTRAINTS.
!  RO  GMAX  MAXIMUM ABSOLUTE VALUE OF A PARTIAL DERIVATIVE.
!  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF A NEGATIVE LAGRANGE MULTIPLIER.
!  IU  IDECF  DECOMPOSITION INDICATOR. IDECF=0-NO DECOMPOSITION.
!         IDECF=1-GILL-MURRAY DECOMPOSITION. IDECF=9-INVERSION.
!         IDECF=10-DIAGONAL MATRIX.
!  II  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
!  IU  NADD  NUMBER OF CONSTRAINT ADDITIONS.
!  IO  IER  ERROR INDICATOR.
!
! SUBPROGRAMS USED :
!  S   PLADR0  CORRECTION OF KERNEL OF THE ORTHOGONAL PROJECTION
!         AFTER CONSTRAINT ADDITION.
!  S   MXDRMM  PREMULTIPLICATION OF A VECTOR BY A ROWWISE STORED DENSE
!         RECTANGULAR MATRIX.
!  S   MXDRMV  COPY OF THE SELECTED COLUMN OF A ROWWISE STORED DENSE
!         RECTANGULAR MATRIX.
!  S   MXDRGR  PLANE ROTATION OF A TRANSPOSED DENSE RECTANGULAR MATRIX.
!         RECTANGULAR MATRIX.
!  S   MXDSMR  PLANE ROTATION OF A DENSE SYMMETRIC MATRIX.
!  S   MXVORT  DETERMINATION OF AN ELEMENTARY ORTHOGONAL MATRIX FOR
!         PLANE ROTATION.
!
      SUBROUTINE PLADB4(Nf,N,Ica,Cg,Cr,Cz,H,S,Eps7,Gmax,Umax,Idecf,Inew,&
                        Nadd,Ier)
      IMPLICIT NONE
      INTEGER Nf , N , Ica(*) , Idecf , Inew , Nadd , Ier
      DOUBLE PRECISION Cg(*) , Cr(*) , Cz(*) , H(*) , S(*) , Eps7 ,     &
                       Gmax , Umax
      DOUBLE PRECISION ck , cl
      INTEGER i , j , k , l , n1
      IF ( Idecf/=0 .AND. Idecf/=9 ) THEN
         Ier = -2
         RETURN
      ENDIF
      CALL PLADR0(Nf,N,Ica,Cg,Cr,S,Eps7,Gmax,Umax,Inew,Nadd,Ier)
      IF ( Ier/=0 ) RETURN
      IF ( N>0 ) THEN
         n1 = N + 1
         IF ( Inew>0 ) THEN
            CALL MXDRMM(Nf,n1,Cz,Cg((Inew-1)*Nf+1),S)
         ELSE
            CALL MXDRMV(Nf,n1,Cz,S,-Inew)
         ENDIF
         DO l = 1 , N
            k = l + 1
            CALL MXVORT(S(k),S(l),ck,cl,Ier)
            CALL MXDRGR(Nf,Cz,k,l,ck,cl,Ier)
            CALL MXDSMR(n1,H,k,l,ck,cl,Ier)
            IF ( Ier<0 ) RETURN
         ENDDO
         IF ( Idecf==9 ) THEN
            l = N*(N+1)/2
            IF ( H(l+n1)/=0.0D0 ) THEN
               cl = 1.0D0/H(l+n1)
               k = 0
               DO i = 1 , N
                  ck = cl*H(l+i)
                  DO j = 1 , i
                     k = k + 1
                     H(k) = H(k) - ck*H(l+j)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      Ier = 0
      END

! SUBROUTINE PLADR0               ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! TRIANGULAR DECOMPOSITION OF KERNEL OF THE ORTHOGONAL PROJECTION
! IS UPDATED AFTER CONSTRAINT ADDITION.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  IU  N  ACTUAL NUMBER OF VARIABLES.
!  IU  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RU  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RA  S(NF)  AUXILIARY VECTOR.
!  RI  EPS7  TOLERANCE FOR LINEAR INDEPENDENCE OF CONSTRAINTS.
!  RO  GMAX  MAXIMUM ABSOLUTE VALUE OF A PARTIAL DERIVATIVE.
!  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF A NEGATIVE LAGRANGE MULTIPLIER.
!  II  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
!  IU  NADD  NUMBER OF CONSTRAINT ADDITIONS.
!  IO  IER  ERROR INDICATOR.
!
! SUBPROGRAMS USED :
!  S   MXSPRB  SPARSE BACK SUBSTITUTION.
!  S   MXVCOP  COPYING OF A VECTOR.
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!
      SUBROUTINE PLADR0(Nf,N,Ica,Cg,Cr,S,Eps7,Gmax,Umax,Inew,Nadd,Ier)
      IMPLICIT NONE
      INTEGER Nf , N , Ica(*) , Inew , Nadd , Ier
      DOUBLE PRECISION Cg(*) , Cr(*) , S(*) , Eps7 , Gmax , Umax
      DOUBLE PRECISION MXVDOT
      INTEGER nca , ncr , i , j , k , l
      Ier = 0
      IF ( N<=0 ) Ier = 2
      IF ( Inew==0 ) Ier = 3
      IF ( Ier/=0 ) RETURN
      nca = Nf - N
      ncr = nca*(nca+1)/2
      IF ( Inew>0 ) THEN
         CALL MXVCOP(Nf,Cg((Inew-1)*Nf+1),S)
         Gmax = MXVDOT(Nf,Cg((Inew-1)*Nf+1),S)
         DO j = 1 , nca
            l = Ica(j)
            IF ( l>0 ) THEN
               Cr(ncr+j) = MXVDOT(Nf,Cg((l-1)*Nf+1),S)
            ELSE
               i = -l
               Cr(ncr+j) = S(i)
            ENDIF
         ENDDO
      ELSE
         k = -Inew
         Gmax = 1.0D0
         DO j = 1 , nca
            l = Ica(j)
            IF ( l>0 ) THEN
               Cr(ncr+j) = Cg((l-1)*Nf+k)*Gmax
            ELSE
               Cr(ncr+j) = 0.0D0
            ENDIF
         ENDDO
      ENDIF
      IF ( nca==0 ) THEN
         Umax = Gmax
      ELSE
         CALL MXDPRB(nca,Cr,Cr(ncr+1),1)
         Umax = Gmax - MXVDOT(nca,Cr(ncr+1),Cr(ncr+1))
      ENDIF
      IF ( Umax<=Eps7*Gmax ) THEN
         Ier = 1
         RETURN
      ELSE
         N = N - 1
         nca = nca + 1
         ncr = ncr + nca
         Ica(nca) = Inew
         Cr(ncr) = SQRT(Umax)
         Nadd = Nadd + 1
      ENDIF
      END

! SUBROUTINE PLLPB1             ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DETERMINATION OF THE INITIAL FEASIBLE POINT AND THE LINEAR PROGRAMMING
! SUBROUTINE.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  II  NC  NUMBER OF LINEAR CONSTRAINTS.
!  RU  X(NF)  VECTOR OF VARIABLES.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  RO  XO(NF)  SAVED VECTOR OF VARIABLES.
!  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
!  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
!  RU  CF(NF)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
!  RA  CFD(NF)  VECTOR CONTAINING INCREMENTS OF THE CONSTRAINT
!         FUNCTIONS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  IO  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RO  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RO  CZ(NF*NF)  MATRIX WHOSE COLUMNS ARE BASIC VECTORS FROM THE
!         CURRENT REDUCED SUBSPACE.
!  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RO  GO(NF)  SAVED GRADIENT OF THE OBJECTIVE FUNCTION.
!  RA  S(NF)  DIRECTION VECTOR.
!  II  MFP  TYPE OF FEASIBLE POINT. MFP=1-ARBITRARY FEASIBLE POINT.
!         MFP=2-OPTIMUM FEASIBLE POINT. MFP=3-REPEATED SOLUTION.
!  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
!         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
!  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
!         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
!         SIDED LINEAR CONSTRAINTS.
!  RI  ETA9  MAXIMUM FOR REAL NUMBERS.
!  RI  EPS7  TOLERANCE FOR LINEAR INDEPENDENCE OF CONSTRAINTS.
!  RI  EPS9  TOLERANCE FOR ACTIVITY OF CONSTRAINTS.
!  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF A NEGATIVE LAGRANGE MULTIPLIER.
!  RO  GMAX  MAXIMUM ABSOLUTE VALUE OF A PARTIAL DERIVATIVE.
!  IO  N  DIMENSION OF THE MANIFOLD DEFINED BY ACTIVE CONSTRAINTS.
!  IO  ITERL  TYPE OF FEASIBLE POINT. ITERL=1-ARBITRARY FEASIBLE POINT.
!         ITERL=2-OPTIMUM FEASIBLE POINT. ITERL=-1 FEASIBLE POINT DOES
!         NOT EXISTS. ITERL=-2 OPTIMUM FEASIBLE POINT DOES NOT EXISTS.
!
! SUBPROGRAMS USED :
!  S   PLINIT  DETERMINATION OF INITIAL POINT SATISFYING SIMPLE BOUNDS.
!  S   PLMAXL  MAXIMUM STEPSIZE USING LINEAR CONSTRAINTS.
!  S   PLMAXS  MAXIMUM STEPSIZE USING SIMPLE BOUNDS.
!  S   PLMAXT  MAXIMUM STEPSIZE USING TRUST REGION BOUNDS.
!  S   PLNEWL  IDENTIFICATION OF ACTIVE LINEAR CONSTRAINTS.
!  S   PLNEWS  IDENTIFICATION OF ACTIVE SIMPLE BOUNDS.
!  S   PLNEWT  IDENTIFICATION OF ACTIVE TRUST REGION BOUNDS.
!  S   PLDIRL  NEW VALUES OF CONSTRAINT FUNCTIONS.
!  S   PLDIRS  NEW VALUES OF VARIABLES.
!  S   PLSETC  INITIAL VALUES OF CONSTRAINT FUNCTIONS.
!  S   PLSETG  DETERMINATION OF THE FIRST PHASE GRADIENT VECTOR.
!  S   PLTRBG  DETERMINATION OF LAGRANGE MULTIPLIERS AND COMPUTATION
!  S   PLADB0  CONSTRAINT ADDITION.
!  S   PLRMB0  CONSTRAINT DELETION.
!  S   MXDCMM  PREMULTIPLICATION OF A VECTOR BY A DENSE RECTANGULAR
!         MATRIX STORED BY COLUMNS.
!  S   MXDRMM  PREMULTIPLICATION OF A VECTOR BY A DENSE RECTANGULAR
!         MATRIX STORED BY ROWS.
!  S   MXDSMI  DETERMINATION OF THE INITIAL UNIT DENSE SYMMETRIC
!         MATRIX.
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVINA  ABSOLUTE VALUES OF ELEMENTS OF AN INTEGER VECTOR.
!  S   MXVINC  UPDATE OF AN INTEGER VECTOR.
!  S   MXVIND  CHANGE OF THE INTEGER VECTOR FOR CONSTRAINT ADDITION.
!  S   MXVINT  CHANGE OF THE INTEGER VECTOR FOR TRUST REGION BOUND
!         ADDITION.
!  S   MXVMUL  DIAGONAL PREMULTIPLICATION OF A VECTOR.
!  S   MXVNEG  COPYING OF A VECTOR WITH CHANGE OF THE SIGN.
!  S   MXVSET  INITIATION OF A VECTOR.
!
      SUBROUTINE PLLPB1(Nf,Nc,X,Ix,Xo,Xl,Xu,Cf,Cfd,Ic,Ica,Cl,Cu,Cg,Cr,  &
                        Cz,G,Go,S,Mfp,Kbf,Kbc,Eta9,Eps7,Eps9,Umax,Gmax, &
                        N,Iterl)
      IMPLICIT NONE
      INTEGER Nf , Nc , Ix(*) , Ic(*) , Ica(*) , Mfp , Kbf , Kbc , N ,  &
              Iterl
      DOUBLE PRECISION X(*) , Xo(*) , Xl(*) , Xu(*) , Cf(*) , Cfd(*) ,  &
                       Cl(*) , Cu(*) , Cg(*) , Cr(*) , Cz(*) , G(*) ,   &
                       Go(*) , S(*) , Eta9 , Eps7 , Eps9 , Umax , Gmax
      DOUBLE PRECISION pom , con , dmax
      INTEGER nca , ncr , ncz , ipom , i , k , iold , inew , ier ,      &
              krem , kc , nred
      INTEGER NREs , NDEc , NREm , NADd , NIT , NFV , NFG , NFH
      COMMON /STAT  / NREs , NDEc , NREm , NADd , NIT , NFV , NFG , NFH
      con = Eta9
!
!     INITIATION
!
      CALL MXVCOP(Nf,X,Xo)
      ipom = 0
      nred = 0
      krem = 0
      Iterl = 1
      dmax = 0.0D0
      IF ( Mfp==3 ) GOTO 200
      IF ( Kbf>0 ) CALL MXVINA(Nf,Ix)
!
!     SHIFT OF VARIABLES FOR SATISFYING SIMPLE BOUNDS
!
      CALL PLINIT(Nf,X,Ix,Xl,Xu,Eps9,Kbf,inew,Iterl)
      IF ( Iterl<0 ) GOTO 99999
      N = 0
      nca = 0
      ncz = 0
      DO i = 1 , Nf
         IF ( Kbf>0 .AND. Ix(i)<0 ) THEN
            nca = nca + 1
            Ica(nca) = -i
         ELSE
            N = N + 1
            CALL MXVSET(Nf,0.0D0,Cz(ncz+1))
            Cz(ncz+i) = 1.0D0
            ncz = ncz + Nf
         ENDIF
      ENDDO
      CALL MXDSMI(nca,Cr)
      IF ( Nc>0 ) THEN
         CALL MXDRMM(Nf,Nc,Cg,X,Cf)
!
!     ADDITION OF ACTIVE CONSTRAINTS AND INITIAL CHECK OF FEASIBILITY
!
         CALL MXVINA(Nc,Ic)
!      IF (NF.GT.N) CALL PLSETC(NF,NC,X,XO,CF,IC,CG,S)
         DO kc = 1 , Nc
            IF ( Ic(kc)/=0 ) THEN
               inew = 0
               CALL PLNEWL(kc,Cf,Ic,Cl,Cu,Eps9,inew)
               CALL PLADB0(Nf,N,Ica,Cg,Cr,Cz,S,Eps7,Gmax,Umax,inew,NADd,&
                           ier)
               CALL MXVIND(Ic,kc,ier)
               IF ( Ic(kc)<-10 ) ipom = 1
            ENDIF
         ENDDO
      ENDIF
 100  IF ( ipom==1 ) THEN
!
!     CHECK OF FEASIBILITY AND UPDATE OF THE FIRST PHASE OBJECTIVE
!     FUNCTION
!
         CALL PLSETG(Nf,Nc,Ic,Cg,Go,inew)
         IF ( inew==0 ) ipom = 0
      ENDIF
      IF ( ipom==0 .AND. Iterl==0 ) THEN
!
!     FEASIBILITY ACHIEVED
!
         Iterl = 1
         CALL MXVCOP(Nf,G,Go)
         IF ( Mfp==1 ) GOTO 99999
      ENDIF
!
!     LAGRANGE MULTIPLIERS AND REDUCED GRADIENT DETERMINATION
!
 200  CALL PLTRBG(Nf,N,Nc,Ix,Ic,Ica,Cg,Cr,Cz,Go,S,Eps7,Gmax,Umax,iold)
      inew = 0
      IF ( Gmax/=0.0D0 ) THEN
!
!     DIRECTION DETERMINATION
!
         nca = Nf - N
         ncr = nca*(nca+1)/2
         CALL MXDCMM(Nf,N,Cz,S,Cr(ncr+1))
         CALL MXVNEG(Nf,Cr(ncr+1),S)
!
!     STEPSIZE SELECTION
!
         pom = con
         CALL PLMAXL(Nf,Nc,Cf,Cfd,Ic,Cl,Cu,Cg,S,pom,Kbc,krem,inew)
         CALL PLMAXS(Nf,X,Ix,Xl,Xu,S,pom,Kbf,krem,inew)
         IF ( inew/=0 ) THEN
!
!     STEP REALIZATION
!
            CALL PLDIRS(Nf,X,Ix,S,pom,Kbf)
            CALL PLDIRL(Nc,Cf,Cfd,Ic,pom,Kbc)
!
!     CONSTRAINT ADDITION
!
            IF ( inew>0 ) THEN
               kc = inew
               inew = 0
               CALL PLNEWL(kc,Cf,Ic,Cl,Cu,Eps9,inew)
               CALL PLADB0(Nf,N,Ica,Cg,Cr,Cz,S,Eps7,Gmax,Umax,inew,NADd,&
                           ier)
               CALL MXVIND(Ic,kc,ier)
            ELSEIF ( inew+Nf>=0 ) THEN
               i = -inew
               inew = 0
               CALL PLNEWS(X,Ix,Xl,Xu,Eps9,i,inew)
               CALL PLADB0(Nf,N,Ica,Cg,Cr,Cz,S,Eps7,Gmax,Umax,inew,NADd,&
                           ier)
               CALL MXVIND(Ix,i,ier)
            ENDIF
            dmax = pom
            IF ( dmax>0.0D0 ) nred = nred + 1
            GOTO 100
         ELSEIF ( ipom==0 ) THEN
!
!     BOUNDED SOLUTION DOES NOT EXIST
!
            Iterl = -2
         ELSE
!
!     FEASIBLE SOLUTION DOES NOT EXIST
!
            Iterl = -3
         ENDIF
!
!     OPTIMUM ON A LINEAR MANIFOLD OBTAINED
!
      ELSEIF ( iold/=0 ) THEN
!
!     CONSTRAINT DELETION
!
         CALL PLRMB0(Nf,N,Ica,Cg,Cr,Cz,Go,S,iold,krem,NREm,ier)
         kc = Ica(Nf-N+1)
         IF ( kc>0 ) THEN
            Ic(kc) = -Ic(kc)
         ELSE
            k = -kc
            Ix(k) = -Ix(k)
         ENDIF
         dmax = 0.0D0
         GOTO 200
      ELSEIF ( ipom==0 ) THEN
!
!     OPTIMAL SOLUTION ACHIEVED
!
         Iterl = 2
      ELSE
         ipom = 0
         DO kc = 1 , Nc
            IF ( Ic(kc)<-10 ) THEN
               inew = 0
               CALL PLNEWL(kc,Cf,Ic,Cl,Cu,Eps9,inew)
               IF ( Ic(kc)<-10 ) ipom = 1
            ENDIF
         ENDDO
         IF ( ipom==0 ) THEN
!
!     OPTIMAL SOLUTION ACHIEVED
!
            CALL MXVCOP(Nf,Go,G)
            Iterl = 2
         ELSE
!
!     FEASIBLE SOLUTION DOES NOT EXIST
!
            CALL MXVCOP(Nf,Go,G)
            Iterl = -1
         ENDIF
      ENDIF
99999 END

! SUBROUTINE PLRMB0               ALL SYSTEMS                92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! OLD LINEAR CONSTRAINT OR AN OLD SIMPLE BOUND IS REMOVED FROM THE
! ACTIVE SET.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  II  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RU  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RU  CZ(NF*NF)  MATRIX WHOSE COLUMNS ARE BASIC VECTORS FROM THE
!         CURRENT REDUCED SUBSPACE.
!  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RU  GN(NF)  TRANSFORMED GRADIENT OF THE OBJECTIVE FUNCTION.
!  II  IOLD  INDEX OF THE OLD ACTIVE CONSTRAINT.
!  IO  KREM  AUXILIARY VARIABLE.
!  IU  NREM NUMBER OF CONSTRAINT DELETION.
!  IO  IER  ERROR INDICATOR.
!
! SUBPROGRAMS USED :
!  S   PLRMR0  CORRECTION OF KERNEL OF THE ORTHOGONAL PROJECTION
!         AFTER CONSTRAINT DELETION.
!  S   MXDPRB  BACK SUBSTITUTION.
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!  S   MXVMUL  DIAGONAL PREMULTIPLICATION OF A VECTOR.
!  S   MXVSET  INITIATION OF A VECTOR.
!
      SUBROUTINE PLRMB0(Nf,N,Ica,Cg,Cr,Cz,G,Gn,Iold,Krem,Nrem,Ier)
      IMPLICIT NONE
      INTEGER Nf , N , Ica(*) , Iold , Krem , Nrem , Ier
      DOUBLE PRECISION Cg(*) , Cr(*) , Cz(*) , G(*) , Gn(*)
      DOUBLE PRECISION MXVDOT
      INTEGER nca , ncr , ncz , i , j , kc
      Ier = 0
      IF ( N==Nf ) Ier = 2
      IF ( Iold==0 ) Ier = 3
      IF ( Ier/=0 ) RETURN
      nca = Nf - N
      ncr = nca*(nca-1)/2
      ncz = N*Nf
      CALL PLRMR0(Nf,Ica,Cr,Cz(ncz+1),N,Iold,Krem,Ier)
      CALL MXVSET(nca,0.0D0,Cz(ncz+1))
      Cz(ncz+nca) = 1.0D0
      CALL MXDPRB(nca,Cr,Cz(ncz+1),-1)
      CALL MXVCOP(nca,Cz(ncz+1),Cr(ncr+1))
      N = N + 1
      CALL MXVSET(Nf,0.0D0,Cz(ncz+1))
      DO j = 1 , nca
         kc = Ica(j)
         IF ( kc>0 ) THEN
            CALL MXVDIR(Nf,Cr(ncr+j),Cg((kc-1)*Nf+1),Cz(ncz+1),Cz(ncz+1)&
                        )
         ELSE
            i = -kc
            Cz(ncz+i) = Cz(ncz+i) + Cr(ncr+j)
         ENDIF
      ENDDO
      Gn(N) = MXVDOT(Nf,Cz(ncz+1),G)
      Nrem = Nrem + 1
      Ier = 0
      END

! SUBROUTINE PLQDB1             ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DUAL RANGE SPACE QUADRATIC PROGRAMMING METHOD.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  IO  N  DIMENSION OF THE MANIFOLD DEFINED BY ACTIVE CONSTRAINTS.
!  II  NC  NUMBER OF LINEAR CONSTRAINTS.
!  RU  X(NF)   VECTOR OF VARIABLES.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  II  IXA(NF)  VECTOR CONTAINING INFORMATION ON TRUST REGION ACTIVITY.
!  RI  XN(NF)  VECTOR OF SCALING FACTORS.
!  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
!  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
!  RI  CF(NF)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
!  RO  CFD(NC)  VECTOR CONTAINING INCREMENTS OF THE CONSTRAINT
!            FUNCTIONS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  II  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RI  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RO  CZ(NF)  VECTOR OF LAGRANGE MULTIPLIERS.
!  RO  G(NF)  GRADIENT OF THE LAGRANGIAN FUNCTION.
!  RO  GO(NF)  SAVED GRADIENT OF THE OBJECTIVE FUNCTION.
!  RU  H(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OR INVERSION OF THE
!         HESSIAN MATRIX APPROXIMATION.
!  RI  S(NF)  DIRECTION VECTOR.
!  RI  ETA2  TOLERANCE FOR POSITIVE DEFINITENESS OF THE HESSIAN MATRIX.
!  RI  ETA9  MAXIMUM FOR REAL NUMBERS.
!  RI  EPS7  TOLERANCE FOR LINEAR INDEPENDENCE OF CONSTRAINTS.
!  RI  EPS9  TOLERANCE FOR ACTIVITY OF CONSTRAINTS.
!  RU  XDEL  TRUST REGION BOUND.
!  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF A NEGATIVE LAGRANGE MULTIPLIER.
!  RO  GMAX  MAXIMUM ABSOLUTE VALUE OF A PARTIAL DERIVATIVE.
!  II  MFP  TYPE OF FEASIBLE POINT. MFP=1-ARBITRARY FEASIBLE POINT.
!         MFP=2-OPTIMUM FEASIBLE POINT. MFP=3-REPEATED SOLUTION.
!
! COMMON DATA :
!  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
!         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
!  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
!         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
!         SIDED LINEAR CONSTRAINTS.
!  II  NORMF  SCALING SPECIFICATION. NORMF=0-NO SCALING PERFORMED.
!         NORMF=1-SCALING FACTORS ARE DETERMINED AUTOMATICALLY.
!         NORMF=2-SCALING FACTORS ARE SUPPLIED BY USER.
!  IU  IDECF  DECOMPOSITION INDICATOR. IDECF=0-NO DECOMPOSITION.
!         IDECF=1-GILL-MURRAY DECOMPOSITION. IDECF=9-INVERSION.
!         IDECF=10-DIAGONAL MATRIX.
!  IU  NDECF  NUMBER OF DECOMPOSITIONS.
!  IO  ITERQ  TYPE OF FEASIBLE POINT. ITERQ=1-ARBITRARY FEASIBLE POINT.
!         ITERQ=2-OPTIMUM FEASIBLE POINT. ITERQ=-1 FEASIBLE POINT DOES
!         NOT EXISTS. ITERQ=-2 OPTIMUM FEASIBLE POINT DOES NOT EXISTS.
!
! SUBPROGRAMS USED :
!  S   PLMINS  DETERMINATION OF THE NEW ACTIVE SIMPLE BOUND.
!  S   PLMINL  DETERMINATION OF THE NEW ACTIVE LINEAR CONSTRAINT.
!  S   PLMINT  DETERMINATION OF THE NEW ACTIVE TRUST REGION BOUND.
!  S   PLADR1  ADDITION OF A NEW ACTIVE CONSTRAINT.
!  S   PLRMR0  CONSTRAIN DELETION.
!  S   PLSOB1  TRANSFORMATION OF THE LOCAL SOLUTION TO THE SOLUTION
!         OF THE ORIGINAL QP PROBLEM.
!  S   MXDPGF  GILL-MURRAY DECOMPOSITION OF A DENSE SYMMETRIC MATRIX.
!  S   MXDPGB  BACK SUBSTITUTION AFTER GILL-MURRAY DECOMPOSITION.
!  S   MXDPRB  BACK SUBSTITUTION.
!  S   MXDSMM  MATRIX VECTOR PRODUCT.
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
!  S   MXVINA  ABSOLUTE VALUES OF ELEMENTS OF AN INTEGER VECTOR.
!  S   MXVINV  CHANGE OF AN INTEGER VECTOR AFTER CONSTRAINT ADDITION.
!  S   MXVNEG  COPYING OF A VECTOR WITH CHANGE OF THE SIGN.
!
      SUBROUTINE PLQDB1(Nf,Nc,X,Ix,Xl,Xu,Cf,Cfd,Ic,Ica,Cl,Cu,Cg,Cr,Cz,G,&
                        Go,H,S,Mfp,Kbf,Kbc,Idecf,Eta2,Eta9,Eps7,Eps9,   &
                        Umax,Gmax,N,Iterq)
      IMPLICIT NONE
      INTEGER Nf , Nc , Ix(*) , Ic(*) , Ica(*) , Mfp , Kbf , Kbc ,      &
              Idecf , N , Iterq
      DOUBLE PRECISION X(*) , Xl(*) , Xu(*) , Cf(*) , Cfd(*) , Cl(*) ,  &
                       Cu(*) , Cg(*) , Cr(*) , Cz(*) , G(*) , Go(*) ,   &
                       H(*) , S(*) , Eta2 , Eta9 , Eps7 , Eps9 , Umax , &
                       Gmax
      DOUBLE PRECISION con , temp , step , step1 , step2 , dmax , par , &
                       snorm
      INTEGER nca , ncr , i , j , k , iold , jold , inew , jnew , knew ,&
              inf , ier , krem , kc , nred
      INTEGER NREs , NDEc , NREm , NADd , NIT , NFV , NFG , NFH
      COMMON /STAT  / NREs , NDEc , NREm , NADd , NIT , NFV , NFG , NFH
      con = Eta9
      IF ( Idecf<0 ) Idecf = 1
      IF ( Idecf==0 ) THEN
!
!     GILL-MURRAY DECOMPOSITION
!
         temp = Eta2
         CALL MXDPGF(Nf,H,inf,temp,step)
         NDEc = NDEc + 1
         Idecf = 1
      ENDIF
      IF ( Idecf>=2 .AND. Idecf<=8 ) THEN
         Iterq = -10
         RETURN
      ENDIF
!
!     INITIATION
!
      nred = 0
      jold = 0
      jnew = 0
      Iterq = 0
      dmax = 0.0D0
      IF ( Mfp/=3 ) THEN
         N = Nf
         nca = 0
         ncr = 0
         IF ( Kbf>0 ) CALL MXVINA(Nf,Ix)
         IF ( Kbc>0 ) CALL MXVINA(Nc,Ic)
      ENDIF
!
!     DIRECTION DETERMINATION
!
 100  CALL MXVNEG(Nf,Go,S)
      DO j = 1 , nca
         kc = Ica(j)
         IF ( kc>0 ) THEN
            CALL MXVDIR(Nf,Cz(j),Cg((kc-1)*Nf+1),S,S)
         ELSE
            k = -kc
            S(k) = S(k) + Cz(j)
         ENDIF
      ENDDO
      CALL MXVCOP(Nf,S,G)
      IF ( Idecf==1 ) THEN
         CALL MXDPGB(Nf,H,S,0)
      ELSE
         CALL MXDSMM(Nf,H,G,S)
      ENDIF
      IF ( Iterq/=3 ) THEN
!
!     CHECK OF FEASIBILITY
!
         inew = 0
         par = 0.0D0
         CALL PLMINN(Nf,Nc,Cf,Cfd,Ic,Cl,Cu,Cg,S,Eps9,par,Kbc,inew,knew)
         CALL PLMINS(Nf,Ix,X,Xl,Xu,S,Kbf,inew,knew,Eps9,par)
         IF ( inew==0 ) THEN
!
!     SOLUTION ACHIEVED
!
            CALL MXVNEG(Nf,G,G)
            Iterq = 2
            GOTO 99999
         ELSE
            snorm = 0.0D0
         ENDIF
 150     ier = 0
!
!     STEPSIZE DETERMINATION
!
         CALL PLADR1(Nf,N,Ica,Cg,Cr,H,S,G,Eps7,Gmax,Umax,Idecf,inew,    &
                     NADd,ier,1)
         CALL MXDPRB(nca,Cr,G,-1)
         IF ( knew<0 ) CALL MXVNEG(nca,G,G)
!
!     PRIMAL STEPSIZE
!
         IF ( ier/=0 ) THEN
            step1 = con
         ELSE
            step1 = -par/Umax
         ENDIF
!
!     DUAL STEPSIZE
!
         iold = 0
         step2 = con
         DO j = 1 , nca
            kc = Ica(j)
            IF ( kc>=0 ) THEN
               k = Ic(kc)
            ELSE
               i = -kc
               k = Ix(i)
            ENDIF
            IF ( k<=-5 ) THEN
            ELSEIF ( (k==-1 .OR. k==-3.) .AND. G(j)<=0.0D0 ) THEN
            ELSEIF ( .NOT.((k==-2 .OR. k==-4.) .AND. G(j)>=0.0D0) ) THEN
               temp = Cz(j)/G(j)
               IF ( step2>temp ) THEN
                  iold = j
                  step2 = temp
               ENDIF
            ENDIF
         ENDDO
!
!     FINAL STEPSIZE
!
         step = MIN(step1,step2)
         IF ( step>=con ) THEN
!
!     FEASIBLE SOLUTION DOES NOT EXIST
!
            Iterq = -1
            GOTO 99999
         ENDIF
!
!     NEW LAGRANGE MULTIPLIERS
!
         dmax = step
         CALL MXVDIR(nca,-step,G,Cz,Cz)
         snorm = snorm + SIGN(1,knew)*step
         par = par - (step/step1)*par
         IF ( step==step1 ) THEN
            IF ( N<=0 ) THEN
!
!     IMPOSSIBLE SITUATION
!
               Iterq = -5
               GOTO 99999
            ENDIF
!
!     CONSTRAINT ADDITION
!
            IF ( ier==0 ) THEN
               N = N - 1
               nca = nca + 1
               ncr = ncr + nca
               Cz(nca) = snorm
            ENDIF
            IF ( inew>0 ) THEN
               kc = inew
               CALL MXVINV(Ic,kc,knew)
            ELSEIF ( ABS(knew)==1 ) THEN
               i = -inew
               CALL MXVINV(Ix,i,knew)
            ELSE
               i = -inew
               IF ( knew>0 ) Ix(i) = -3
               IF ( knew<0 ) Ix(i) = -4
            ENDIF
            nred = nred + 1
            NADd = NADd + 1
            jnew = inew
            jold = 0
            GOTO 100
         ELSE
!
!     CONSTRAINT DELETION
!
            DO j = iold , nca - 1
               Cz(j) = Cz(j+1)
            ENDDO
            CALL PLRMF0(Nf,Nc,Ix,Ic,Ica,Cr,Ic,G,N,iold,krem,ier)
            ncr = ncr - nca
            nca = nca - 1
            jold = iold
            jnew = 0
            IF ( Kbc>0 ) CALL MXVINA(Nc,Ic)
            IF ( Kbf>0 ) CALL MXVINA(Nf,Ix)
            DO j = 1 , nca
               kc = Ica(j)
               IF ( kc>0 ) THEN
                  Ic(kc) = -Ic(kc)
               ELSE
                  kc = -kc
                  Ix(kc) = -Ix(kc)
               ENDIF
            ENDDO
            GOTO 150
         ENDIF
      ENDIF
99999 END

! SUBROUTINE PLADR1               ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! TRIANGULAR DECOMPOSITION OF KERNEL OF THE GENERAL PROJECTION
! IS UPDATED AFTER CONSTRAINT ADDITION.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  IU  N  ACTUAL NUMBER OF VARIABLES.
!  IU  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RI  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RU  H(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OR INVERSION OF THE
!         HESSIAN MATRIX APPROXIMATION.
!  RA  S(NF)  AUXILIARY VECTOR.
!  RO  G(NF)  VECTOR USED IN THE DUAL RANGE SPACE QUADRATIC PROGRAMMING
!         METHOD.
!  RI  EPS7  TOLERANCE FOR LINEAR INDEPENDENCE OF CONSTRAINTS.
!  RO  GMAX  MAXIMUM ABSOLUTE VALUE OF A PARTIAL DERIVATIVE.
!  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF A NEGATIVE LAGRANGE MULTIPLIER.
!  RO  E  AUXILIARY VARIABLE.
!  RI  T  AUXILIARY VARIABLE.
!  IU  IDECF  DECOMPOSITION INDICATOR. IDECF=0-NO DECOMPOSITION.
!         IDECF=1-GILL-MURRAY DECOMPOSITION. IDECF=9-INVERSION.
!         IDECF=10-DIAGONAL MATRIX.
!  II  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
!  IU  NADD  NUMBER OF CONSTRAINT ADDITIONS.
!  IO  IER  ERROR INDICATOR.
!  II  JOB  SPECIFICATION OF COMPUTATION. OUTPUT VECTOR G IS NOT OR IS
!         COMPUTED IN CASE WHEN N.LE.0 IF JOB=0 OR JOB=1 RESPECTIVELY.
!
! SUBPROGRAMS USED :
!  S   MXDPGB  BACK SUBSTITUTION.
!  S   MXDPRB  BACK SUBSTITUTION.
!  S   MXDSMM  MATRIX-VECTOR PRODUCT.
!  S   MXDSMV  COPYING OF A ROW OF DENSE SYMMETRIC MATRIX.
!  S   MXVCOP  COPYING OF A VECTOR.
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!
      SUBROUTINE PLADR1(Nf,N,Ica,Cg,Cr,H,S,G,Eps7,Gmax,Umax,Idecf,Inew, &
                        Nadd,Ier,Job)
      IMPLICIT NONE
      INTEGER Nf , N , Ica(*) , Idecf , Inew , Nadd , Ier , Job
      DOUBLE PRECISION Cg(*) , Cr(*) , H(*) , S(*) , G(*) , Eps7 ,      &
                       Gmax , Umax
      DOUBLE PRECISION MXVDOT
      INTEGER nca , ncr , jcg , j , k , l
      Ier = 0
      IF ( Job==0 .AND. N<=0 ) Ier = 2
      IF ( Inew==0 ) Ier = 3
      IF ( Idecf/=1 .AND. Idecf/=9 ) Ier = -2
      IF ( Ier/=0 ) RETURN
      nca = Nf - N
      ncr = nca*(nca+1)/2
      IF ( Inew>0 ) THEN
         jcg = (Inew-1)*Nf + 1
         IF ( Idecf==1 ) THEN
            CALL MXVCOP(Nf,Cg(jcg),S)
            CALL MXDPGB(Nf,H,S,0)
         ELSE
            CALL MXDSMM(Nf,H,Cg(jcg),S)
         ENDIF
         Gmax = MXVDOT(Nf,Cg(jcg),S)
      ELSE
         k = -Inew
         IF ( Idecf==1 ) THEN
            CALL MXVSET(Nf,0.0D0,S)
            S(k) = 1.0D0
            CALL MXDPGB(Nf,H,S,0)
         ELSE
            CALL MXDSMV(Nf,H,S,k)
         ENDIF
         Gmax = S(k)
      ENDIF
      DO j = 1 , nca
         l = Ica(j)
         IF ( l>0 ) THEN
            G(j) = MXVDOT(Nf,Cg((l-1)*Nf+1),S)
         ELSE
            l = -l
            G(j) = S(l)
         ENDIF
      ENDDO
      IF ( N==0 ) THEN
         CALL MXDPRB(nca,Cr,G,1)
         Umax = 0.0D0
         Ier = 2
         RETURN
      ELSEIF ( nca==0 ) THEN
         Umax = Gmax
      ELSE
         CALL MXDPRB(nca,Cr,G,1)
         Umax = Gmax - MXVDOT(nca,G,G)
         CALL MXVCOP(nca,G,Cr(ncr+1))
      ENDIF
      IF ( Umax<=Eps7*Gmax ) THEN
         Ier = 1
         RETURN
      ELSE
         nca = nca + 1
         ncr = ncr + nca
         Ica(nca) = Inew
         Cr(ncr) = SQRT(Umax)
         IF ( Job==0 ) THEN
            N = N - 1
            Nadd = Nadd + 1
         ENDIF
      ENDIF
      END

! SUBROUTINE PLDIRL               ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DETERMINATION OF THE NEW VALUES OF THE CONSTRAINT FUNCTIONS.
!
! PARAMETERS :
!  II  NC  NUMBER OF CONSTRAINTS.
!  RU  CF(NF)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
!  RI  CFD(NF)  VECTOR CONTAINING INCREMENTS OF THE CONSTRAINT
!         FUNCTIONS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  RI  STEP  CURRENT STEPSIZE.
!  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
!         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
!         SIDED LINEAR CONSTRAINTS.
!
      SUBROUTINE PLDIRL(Nc,Cf,Cfd,Ic,Step,Kbc)
      IMPLICIT NONE
      INTEGER Nc , Ic(*) , Kbc
      DOUBLE PRECISION Cf(*) , Cfd(*) , Step
      INTEGER kc
      IF ( Kbc>0 ) THEN
         DO kc = 1 , Nc
            IF ( Ic(kc)>=0 .AND. Ic(kc)<=10 ) THEN
               Cf(kc) = Cf(kc) + Step*Cfd(kc)
            ELSEIF ( Ic(kc)<-10 ) THEN
               Cf(kc) = Cf(kc) + Step*Cfd(kc)
            ENDIF
         ENDDO
      ENDIF
      END

! SUBROUTINE PLDIRS               ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DETERMINATION OF THE NEW VECTOR OF VARIABLES.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  RU  X(NF)  VECTOR OF VARIABLES.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  RI  S(NF)  DIRECTION VECTOR.
!  RI  STEP  CURRENT STEPSIZE.
!  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
!         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
!
      SUBROUTINE PLDIRS(Nf,X,Ix,S,Step,Kbf)
      IMPLICIT NONE
      INTEGER Nf , Ix(*) , Kbf
      DOUBLE PRECISION X(*) , S(*) , Step
      INTEGER i
      DO i = 1 , Nf
         IF ( Kbf<=0 ) THEN
            X(i) = X(i) + Step*S(i)
         ELSEIF ( Ix(i)>=0 .AND. Ix(i)<=10 ) THEN
            X(i) = X(i) + Step*S(i)
         ELSEIF ( Ix(i)<-10 ) THEN
            X(i) = X(i) + Step*S(i)
         ENDIF
      ENDDO
      END

! SUBROUTINE PLINIT             ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DETERMINATION OF THE INITIAL POINT WHICH SATISFIES SIMPLE BOUNDS.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  RU  X(NF)  VECTOR OF VARIABLES.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
!  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
!  RI  EPS9  TOLERANCE FOR ACTIVE CONSTRAINTS.
!  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
!  IO  IND  INDICATOR. IF IND.NE.0 THEN TRUST REGION BOUNDS CANNOT
!         BE SATISFIED.
!
! SUBPROGRAMS USED :
!  S   PLNEWS  TEST ON ACTIVITY OF A GIVEN SIMPLE BOUND.
!
      SUBROUTINE PLINIT(Nf,X,Ix,Xl,Xu,Eps9,Kbf,Inew,Ind)
      IMPLICIT NONE
      INTEGER Nf , Ix(*) , Kbf , Inew , Ind
      DOUBLE PRECISION X(*) , Xl(*) , Xu(*) , Eps9
      INTEGER i
      Ind = 0
      IF ( Kbf>0 ) THEN
         DO i = 1 , Nf
            CALL PLNEWS(X,Ix,Xl,Xu,Eps9,i,Inew)
            IF ( Ix(i)<5 ) THEN
            ELSEIF ( Ix(i)==5 ) THEN
               Ix(i) = -5
            ELSEIF ( Ix(i)==11 .OR. Ix(i)==13 ) THEN
               X(i) = Xl(i)
               Ix(i) = 10 - Ix(i)
            ELSEIF ( Ix(i)==12 .OR. Ix(i)==14 ) THEN
               X(i) = Xu(i)
               Ix(i) = 10 - Ix(i)
            ENDIF
         ENDDO
      ENDIF
      END

! SUBROUTINE PLMAXL               ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DETERMINATION OF THE MAXIMUM STEPSIZE USING LINEAR CONSTRAINTS.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  II  NC  NUMBER OF CURRENT LINEAR CONSTRAINTS.
!  RI  CF(NF)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCYIONS.
!  RO  CFD(NF)  VECTOR CONTAINING INCREMENTS OF THE CONSTRAINT
!         FUNCTIONS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RI  S(NF)  DIRECTION VECTOR.
!  RO  STEP  MAXIMUM STEPSIZE.
!  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
!         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
!         SIDED LINEAR CONSTRAINTS.
!  II  KREM  INDICATION OF LINEARLY DEPENDENT GRADIENTS.
!  IO  INEW  INDEX OF THE NEW ACTIVE FUNCTION.
!
! SUBPROGRAMS USED :
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!
      SUBROUTINE PLMAXL(Nf,Nc,Cf,Cfd,Ic,Cl,Cu,Cg,S,Step,Kbc,Krem,Inew)
      IMPLICIT NONE
      INTEGER Nf , Nc , Ic(*) , Kbc , Krem , Inew
      DOUBLE PRECISION Cf(*) , Cfd(*) , Cl(*) , Cu(*) , Cg(*) , S(*) ,  &
                       Step
      DOUBLE PRECISION temp , MXVDOT
      INTEGER jcg , kc
      IF ( Kbc>0 ) THEN
         jcg = 1
         DO kc = 1 , Nc
            IF ( Krem>0 .AND. Ic(kc)>10 ) Ic(kc) = Ic(kc) - 10
            IF ( Ic(kc)>0 .AND. Ic(kc)<=10 ) THEN
               temp = MXVDOT(Nf,Cg(jcg),S)
               Cfd(kc) = temp
               IF ( temp<0.0D0 ) THEN
                  IF ( Ic(kc)==1 .OR. Ic(kc)>=3 ) THEN
                     temp = (Cl(kc)-Cf(kc))/temp
                     IF ( temp<=Step ) THEN
                        Inew = kc
                        Step = temp
                     ENDIF
                  ENDIF
               ELSEIF ( temp>0.0D0 ) THEN
                  IF ( Ic(kc)==2 .OR. Ic(kc)>=3 ) THEN
                     temp = (Cu(kc)-Cf(kc))/temp
                     IF ( temp<=Step ) THEN
                        Inew = kc
                        Step = temp
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF ( Ic(kc)<-10 ) THEN
               temp = MXVDOT(Nf,Cg(jcg),S)
               Cfd(kc) = temp
               IF ( temp>0.0D0 ) THEN
                  IF ( Ic(kc)==-11 .OR. Ic(kc)==-13 .OR. Ic(kc)==-15 )  &
                       THEN
                     temp = (Cl(kc)-Cf(kc))/temp
                     IF ( temp<=Step ) THEN
                        Inew = kc
                        Step = temp
                     ENDIF
                  ENDIF
               ELSEIF ( temp<0.0D0 ) THEN
                  IF ( Ic(kc)==-12 .OR. Ic(kc)==-14 .OR. Ic(kc)==-16 )  &
                       THEN
                     temp = (Cu(kc)-Cf(kc))/temp
                     IF ( temp<=Step ) THEN
                        Inew = kc
                        Step = temp
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            jcg = jcg + Nf
         ENDDO
      ENDIF
      END

! SUBROUTINE PLMAXS               ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DETERMINATION OF THE MAXIMUM STEPSIZE USING THE SIMPLE BOUNDS
! FOR VARIABLES.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
!  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
!  RI  S(NF)  DIRECTION VECTOR.
!  RO  STEP  MAXIMUM STEPSIZE.
!  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
!         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
!  IO  KREM  INDICATION OF LINEARLY DEPENDENT GRADIENTS.
!  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
!
      SUBROUTINE PLMAXS(Nf,X,Ix,Xl,Xu,S,Step,Kbf,Krem,Inew)
      IMPLICIT NONE
      INTEGER Nf , Ix(*) , Kbf , Krem , Inew
      DOUBLE PRECISION X(*) , Xl(*) , Xu(*) , S(*) , Step
      DOUBLE PRECISION temp
      INTEGER i
      IF ( Kbf>0 ) THEN
         DO i = 1 , Nf
            IF ( Krem>0 .AND. Ix(i)>10 ) Ix(i) = Ix(i) - 10
            IF ( Ix(i)>0 .AND. Ix(i)<=10 ) THEN
               IF ( S(i)<0.0D0 ) THEN
                  IF ( Ix(i)==1 .OR. Ix(i)>=3 ) THEN
                     temp = (Xl(i)-X(i))/S(i)
                     IF ( temp<=Step ) THEN
                        Inew = -i
                        Step = temp
                     ENDIF
                  ENDIF
               ELSEIF ( S(i)>0.0D0 ) THEN
                  IF ( Ix(i)==2 .OR. Ix(i)>=3 ) THEN
                     temp = (Xu(i)-X(i))/S(i)
                     IF ( temp<=Step ) THEN
                        Inew = -i
                        Step = temp
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      Krem = 0
      END

! SUBROUTINE PLNEWL             ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! TEST ON ACTIVITY OF A GIVEN LINEAR CONSTRAINT.
!
! PARAMETERS :
!  II  KC  INDEX OF A GIVEN CONSTRAINT.
!  RI  CF(NC)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
!  IU  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  EPS9  TOLERANCE FOR ACTIVE CONSTRAINTS.
!  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
!
      SUBROUTINE PLNEWL(Kc,Cf,Ic,Cl,Cu,Eps9,Inew)
      IMPLICIT NONE
      INTEGER Kc , Ic(*) , Inew
      DOUBLE PRECISION Cf(*) , Cl(*) , Cu(*) , Eps9
      DOUBLE PRECISION temp
      IF ( Ic(Kc)<-10 ) Ic(Kc) = -Ic(Kc) - 10
      IF ( Ic(Kc)<=0 ) THEN
      ELSEIF ( Ic(Kc)==1 ) THEN
         temp = Eps9*MAX(ABS(Cl(Kc)),1.0D0)
         IF ( Cf(Kc)>Cl(Kc)+temp ) THEN
         ELSEIF ( Cf(Kc)>=Cl(Kc)-temp ) THEN
            Ic(Kc) = 11
            Inew = Kc
         ELSE
            Ic(Kc) = -11
         ENDIF
      ELSEIF ( Ic(Kc)==2 ) THEN
         temp = Eps9*MAX(ABS(Cu(Kc)),1.0D0)
         IF ( Cf(Kc)<Cu(Kc)-temp ) THEN
         ELSEIF ( Cf(Kc)<=Cu(Kc)+temp ) THEN
            Ic(Kc) = 12
            Inew = Kc
         ELSE
            Ic(Kc) = -12
         ENDIF
      ELSEIF ( Ic(Kc)==3 .OR. Ic(Kc)==4 ) THEN
         temp = Eps9*MAX(ABS(Cl(Kc)),1.0D0)
         IF ( Cf(Kc)>Cl(Kc)+temp ) THEN
            temp = Eps9*MAX(ABS(Cu(Kc)),1.0D0)
            IF ( Cf(Kc)<Cu(Kc)-temp ) THEN
            ELSEIF ( Cf(Kc)<=Cu(Kc)+temp ) THEN
               Ic(Kc) = 14
               Inew = Kc
            ELSE
               Ic(Kc) = -14
            ENDIF
         ELSEIF ( Cf(Kc)>=Cl(Kc)-temp ) THEN
            Ic(Kc) = 13
            Inew = Kc
         ELSE
            Ic(Kc) = -13
         ENDIF
      ELSEIF ( Ic(Kc)==5 .OR. Ic(Kc)==6 ) THEN
         temp = Eps9*MAX(ABS(Cl(Kc)),1.0D0)
         IF ( Cf(Kc)>Cl(Kc)+temp ) THEN
            temp = Eps9*MAX(ABS(Cu(Kc)),1.0D0)
            IF ( Cf(Kc)<Cu(Kc)-temp ) THEN
            ELSEIF ( Cf(Kc)<=Cu(Kc)+temp ) THEN
               Ic(Kc) = 16
               Inew = Kc
            ELSE
               Ic(Kc) = -16
            ENDIF
         ELSEIF ( Cf(Kc)>=Cl(Kc)-temp ) THEN
            Ic(Kc) = 15
            Inew = Kc
         ELSE
            Ic(Kc) = -15
         ENDIF
      ENDIF
      END

! SUBROUTINE PLMINN             ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DETERMINATION OF THE NEW ACTIVE LINEAR CONSTRAINT.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  II  NC  NUMBER OF CONSTRAINTS.
!  RI  CF(NC)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
!  RO  CFD(NC)  VECTOR CONTAINING INCREMENTS OF THE CONSTRAINT
!            FUNCTIONS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RI  S(NF)  DIRECTION VECTOR.
!  RI  EPS9  TOLERANCE FOR ACTIVE CONSTRAINTS.
!  RA  PAR  AUXILIARY VARIABLE.
!  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
!         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
!         SIDED LINEAR CONSTRAINTS.
!  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
!  IO  KNEW  SIGNUM OF THE NEW ACTIVE NORMAL.
!
! SUBPROGRAMS USED :
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!
      SUBROUTINE PLMINN(Nf,Nc,Cf,Cfd,Ic,Cl,Cu,Cg,S,Eps9,Par,Kbc,Inew,   &
                        Knew)
      IMPLICIT NONE
      INTEGER Nf , Nc , Ic(*) , Kbc , Inew , Knew
      DOUBLE PRECISION Cf(*) , Cfd(*) , Cl(*) , Cu(*) , Cg(*) , S(*) ,  &
                       Eps9 , Par
      DOUBLE PRECISION temp , pom , MXVDOT
      INTEGER jcg , kc
      IF ( Kbc>0 ) THEN
         jcg = 1
         DO kc = 1 , Nc
            IF ( Ic(kc)>0 ) THEN
               temp = MXVDOT(Nf,Cg(jcg),S)
               Cfd(kc) = temp
               temp = Cf(kc) + temp
               IF ( Ic(kc)==1 .OR. Ic(kc)>=3 ) THEN
                  pom = temp - Cl(kc)
                  IF ( pom<MIN(Par,-Eps9*MAX(ABS(Cl(kc)),1.0D0)) ) THEN
                     Inew = kc
                     Knew = 1
                     Par = pom
                  ENDIF
               ENDIF
               IF ( Ic(kc)==2 .OR. Ic(kc)>=3 ) THEN
                  pom = Cu(kc) - temp
                  IF ( pom<MIN(Par,-Eps9*MAX(ABS(Cu(kc)),1.0D0)) ) THEN
                     Inew = kc
                     Knew = -1
                     Par = pom
                  ENDIF
               ENDIF
            ENDIF
            jcg = jcg + Nf
         ENDDO
      ENDIF
      END

! SUBROUTINE PLMINS             ALL SYSTEMS                   91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DETERMINATION OF THE NEW ACTIVE SIMPLE BOUND.
!
! PARAMETERS :
!  II  NF DECLARED NUMBER OF VARIABLES.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  RI  XO(NF)  SAVED VECTOR OF VARIABLES.
!  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
!  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
!  RI  S(NF)  DIRECTION VECTOR.
!  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
!         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
!  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
!  IO  KNEW  SIGNUM OF THE NEW NORMAL.
!  RI  EPS9  TOLERANCE FOR ACTIVE CONSTRAINTS.
!  RA  PAR  AUXILIARY VARIABLE.
!
      SUBROUTINE PLMINS(Nf,Ix,Xo,Xl,Xu,S,Kbf,Inew,Knew,Eps9,Par)
      IMPLICIT NONE
      DOUBLE PRECISION Eps9 , Par
      INTEGER Inew , Kbf , Knew , Nf
      DOUBLE PRECISION S(*) , Xl(*) , Xo(*) , Xu(*)
      INTEGER Ix(*)
      DOUBLE PRECISION pom , temp
      INTEGER i
      IF ( Kbf>0 ) THEN
         DO i = 1 , Nf
            IF ( Ix(i)>0 ) THEN
               temp = 1.0D0
               IF ( Ix(i)==1 .OR. Ix(i)>=3 ) THEN
                  pom = Xo(i) + S(i)*temp - Xl(i)
                  IF ( pom<MIN(Par,-Eps9*MAX(ABS(Xl(i)),temp)) ) THEN
                     Inew = -i
                     Knew = 1
                     Par = pom
                  ENDIF
               ENDIF
               IF ( Ix(i)==2 .OR. Ix(i)>=3 ) THEN
                  pom = Xu(i) - S(i)*temp - Xo(i)
                  IF ( pom<MIN(Par,-Eps9*MAX(ABS(Xu(i)),temp)) ) THEN
                     Inew = -i
                     Knew = -1
                     Par = pom
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      END

! SUBROUTINE PLNEWS             ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! TEST ON ACTIVITY OF A GIVEN SIMPLE BOUND.
!
! PARAMETERS :
!  RI  X(NF)  VECTOR OF VARIABLES.
!  IU  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
!  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
!  RI  EPS9  TOLERANCE FOR ACTIVE CONSTRAINTS.
!  II  I  INDEX OF TESTED SIMPLE BOUND.
!  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
!
      SUBROUTINE PLNEWS(X,Ix,Xl,Xu,Eps9,I,Inew)
      IMPLICIT NONE
      INTEGER Ix(*) , I , Inew
      DOUBLE PRECISION X(*) , Xl(*) , Xu(*) , Eps9
      DOUBLE PRECISION temp
      temp = 1.0D0
      IF ( Ix(I)<=0 ) THEN
      ELSEIF ( Ix(I)==1 ) THEN
         IF ( X(I)<=Xl(I)+Eps9*MAX(ABS(Xl(I)),temp) ) THEN
            Ix(I) = 11
            Inew = -I
         ENDIF
      ELSEIF ( Ix(I)==2 ) THEN
         IF ( X(I)>=Xu(I)-Eps9*MAX(ABS(Xu(I)),temp) ) THEN
            Ix(I) = 12
            Inew = -I
         ENDIF
      ELSEIF ( Ix(I)==3 .OR. Ix(I)==4 ) THEN
         IF ( X(I)<=Xl(I)+Eps9*MAX(ABS(Xl(I)),temp) ) THEN
            Ix(I) = 13
            Inew = -I
         ENDIF
         IF ( X(I)>=Xu(I)-Eps9*MAX(ABS(Xu(I)),temp) ) THEN
            Ix(I) = 14
            Inew = -I
         ENDIF
      ENDIF
      END

! SUBROUTINE PLREDL               ALL SYSTEMS                   98/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! TRANSFORMATION OF THE INCOMPATIBLE QUADRATIC PROGRAMMING SUBPROBLEM.
!
! PARAMETERS :
!  II  NC  NUMBER OF CURRENT LINEAR CONSTRAINTS.
!  RI  CF(NF)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
!         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
!         SIDED LINEAR CONSTRAINTS.
!
      SUBROUTINE PLREDL(Nc,Cf,Ic,Cl,Cu,Kbc)
      IMPLICIT NONE
      INTEGER Nc , Ic(Nc) , Kbc , k
      DOUBLE PRECISION Cf(*) , Cl(*) , Cu(*)
      DOUBLE PRECISION temp
      INTEGER kc
      IF ( Kbc>0 ) THEN
         DO kc = 1 , Nc
            k = Ic(kc)
            IF ( ABS(k)==1 .OR. ABS(k)==3 .OR. ABS(k)==4 ) THEN
               temp = (Cf(kc)-Cl(kc))
               IF ( temp<0 ) Cf(kc) = Cl(kc) + 0.1D0*temp
            ENDIF
            IF ( ABS(k)==2 .OR. ABS(k)==3 .OR. ABS(k)==4 ) THEN
               temp = (Cf(kc)-Cu(kc))
               IF ( temp>0 ) Cf(kc) = Cu(kc) + 0.1D0*temp
            ENDIF
            IF ( ABS(k)==5 .OR. ABS(k)==6 ) THEN
               temp = (Cf(kc)-Cl(kc))
               Cf(kc) = Cl(kc) + 0.1D0*temp
            ENDIF
         ENDDO
      ENDIF
      END

! SUBROUTINE PLRMF0             ALL SYSTEMS                   91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! OPERATIONS AFTER CONSTRAINT DELETION.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  II  NC  NUMBER OF CONSTRAINTS.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  II  IA(NA)  VECTOR CONTAINING TYPES OF DEVIATIONS.
!  IU  IAA(NF+1)  VECTOR CONTAINING INDICES OF ACTIVE FUNCTIONS.
!  RU  AR((NF+1)*(NF+2)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  RA  S(NF+1)  AUXILIARY VECTOR.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  II  IOLD  INDEX OF THE OLD ACTIVE CONSTRAINT.
!  IO  KREM  AUXILIARY VARIABLE.
!  IO  IER  ERROR INDICATOR.
!
! SUBPROGRAMS USED :
!  S   PLRMR0  CORRECTION OF KERNEL OF THE ORTHOGONAL PROJECTION
!         AFTER CONSTRAINT DELETION.
!
      SUBROUTINE PLRMF0(Nf,Nc,Ix,Ia,Iaa,Ar,Ic,S,N,Iold,Krem,Ier)
      IMPLICIT NONE
      INTEGER Ier , Iold , Krem , N , Nc , Nf
      DOUBLE PRECISION Ar(*) , S(*)
      INTEGER Ia(*) , Iaa(*) , Ic(*) , Ix(*)
      INTEGER NADd , NDEc , NFG , NFH , NFV , NIT , NREm , NREs
      INTEGER l
      COMMON /STAT  / NREs , NDEc , NREm , NADd , NIT , NFV , NFG , NFH
      CALL PLRMR0(Nf,Iaa,Ar,S,N,Iold,Krem,Ier)
      N = N + 1
      NREm = NREm + 1
      l = Iaa(Nf-N+1)
      IF ( l>Nc ) THEN
         l = l - Nc
         Ia(l) = -Ia(l)
      ELSEIF ( l>0 ) THEN
         Ic(l) = -Ic(l)
      ELSE
         l = -l
         Ix(l) = -Ix(l)
      ENDIF
      END

! SUBROUTINE PLRMR0               ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! TRIANGULAR DECOMPOSITION OF KERNEL OF THE ORTHOGONAL PROJECTION IS
! UPDATED AFTER CONSTRAINT DELETION.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  IU  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RU  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RA  G(NF)  AUXILIARY VECTOR.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  II  IOLD  INDEX OF THE OLD ACTIVE CONSTRAINT.
!  IO  KREM  AUXILIARY VARIABLE.
!  IO  IER  ERROR INDICATOR.
!
! SUBPROGRAMS USED :
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVORT  DETERMINATION OF AN ELEMENTARY ORTHOGONAL MATRIX FOR
!         PLANE ROTATION.
!  S   MXVROT  PLANE ROTATION OF A VECTOR.
!  S   MXVSET  INITIATION OF A VECTOR.
!
      SUBROUTINE PLRMR0(Nf,Ica,Cr,G,N,Iold,Krem,Ier)
      IMPLICIT NONE
      INTEGER Ier , Iold , Krem , N , Nf
      DOUBLE PRECISION Cr(*) , G(*)
      INTEGER Ica(*)
      DOUBLE PRECISION ck , cl
      INTEGER i , j , k , kc , l , nca
      nca = Nf - N
      IF ( Iold<nca ) THEN
         k = Iold*(Iold-1)/2
         kc = Ica(Iold)
         CALL MXVCOP(Iold,Cr(k+1),G)
         CALL MXVSET(nca-Iold,0.0D0,G(Iold+1))
         k = k + Iold
         DO i = Iold + 1 , nca
            k = k + i
            CALL MXVORT(Cr(k-1),Cr(k),ck,cl,Ier)
            CALL MXVROT(G(i-1),G(i),ck,cl,Ier)
            l = k
            DO j = i , nca - 1
               l = l + j
               CALL MXVROT(Cr(l-1),Cr(l),ck,cl,Ier)
            ENDDO
         ENDDO
         k = Iold*(Iold-1)/2
         DO i = Iold , nca - 1
            l = k + i
            Ica(i) = Ica(i+1)
            CALL MXVCOP(i,Cr(l+1),Cr(k+1))
            k = l
         ENDDO
         Ica(nca) = kc
         CALL MXVCOP(nca,G,Cr(k+1))
      ENDIF
      Krem = 1
      END

! SUBROUTINE PLSETC             ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DETERMINATION OF INITIAL VALUES OF THE CONSTRAINT FUNCTIONS.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  II  NC  NUMBER OF CURRENT LINEAR CONSTRAINTS.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  RI  XO(NF)  SAVED VECTOR OF VARIABLES.
!  RU  CF(NF)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCYIONS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  RI  CG(NF*MCL)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RA  S(NF)  AUXILIARY VECTOR.
!
! SUBPROGRAMS USED :
!  S   MXVDIF  DIFFERENCE OF TWO VECTORS.
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!  S   MXVMUL  DIAGONAL PREMULTIPLICATION OF A VECTOR.
!
      SUBROUTINE PLSETC(Nf,Nc,X,Xo,Cf,Ic,Cg,S)
      IMPLICIT NONE
      INTEGER Nf , Nc , Ic(*)
      DOUBLE PRECISION X(*) , Xo(*) , Cf(*) , Cg(*) , S(*)
      DOUBLE PRECISION MXVDOT
      INTEGER jcg , kc
      CALL MXVDIF(Nf,X,Xo,S)
      jcg = 0
      DO kc = 1 , Nc
         IF ( Ic(kc)/=0 ) Cf(kc) = Cf(kc) + MXVDOT(Nf,S,Cg(jcg+1))
         jcg = jcg + Nf
      ENDDO
      END

! SUBROUTINE PLSETG             ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! GRADIENT DETERMINATION IN THE FIRST PHASE OF LP SUBROUTINE.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  II  NC  NUMBER OF CONSTRAINTS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RO  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  IO  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
!
! SUBPROGRAMS USED :
!  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
!  S   MXVSET  INITIATION OF A VECTOR.
!
      SUBROUTINE PLSETG(Nf,Nc,Ic,Cg,G,Inew)
      IMPLICIT NONE
      INTEGER Nf , Nc , Ic(*) , Inew
      DOUBLE PRECISION Cg(*) , G(*)
      INTEGER kc
      CALL MXVSET(Nf,0.0D0,G)
      Inew = 0
      DO kc = 1 , Nc
         IF ( Ic(kc)>=-10 ) THEN
         ELSEIF ( Ic(kc)==-11 .OR. Ic(kc)==-13 .OR. Ic(kc)==-15 ) THEN
            CALL MXVDIR(Nf,-1.0D0,Cg((kc-1)*Nf+1),G,G)
            Inew = 1
         ELSEIF ( Ic(kc)==-12 .OR. Ic(kc)==-14 .OR. Ic(kc)==-16 ) THEN
            CALL MXVDIR(Nf,1.0D0,Cg((kc-1)*Nf+1),G,G)
            Inew = 1
         ENDIF
      ENDDO
      END

! SUBROUTINE PLTLAG               ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! MAXIMUM ABSOLUTE VALUE OF THE NEGATIVE LAGRANGE MULTIPLIER IS
! COMPUTED.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  II  NC  NUMBER OF LINEARIZED CONSTRAINTS.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  II  IA(NA)  VECTOR CONTAINING TYPES OF DEVIATIONS.
!  II  IAA(NF+1)  VECTOR CONTAINING INDICES OF ACTIVE FUNCTIONS.
!  RI  AZ(NF+1)  VECTOR OF LAGRANGE MULTIPLIERS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  RI  EPS7  TOLERANCE FOR LINEAR AND QUADRATIC PROGRAMMING.
!  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF THE NEGATIVE LAGRANGE MULTIPLIER.
!  IO  IOLD  INDEX OF THE REMOVED CONSTRAINT.
!
      SUBROUTINE PLTLAG(Nf,N,Nc,Ix,Ia,Iaa,Az,Ic,Eps7,Umax,Iold)
      IMPLICIT NONE
      INTEGER Nf , N , Nc , Ix(*) , Ia(*) , Iaa(*) , Ic(*) , Iold
      DOUBLE PRECISION Az(*) , Eps7 , Umax
      DOUBLE PRECISION temp
      INTEGER naa , j , k , l
      Iold = 0
      Umax = 0.0D0
      naa = Nf - N
      DO j = 1 , naa
         temp = Az(j)
         l = Iaa(j)
         IF ( l>Nc ) THEN
            l = l - Nc
            k = Ia(l)
         ELSEIF ( l>0 ) THEN
            k = Ic(l)
         ELSE
            l = -l
            k = Ix(l)
         ENDIF
         IF ( k<=-5 ) THEN
         ELSEIF ( (k==-1 .OR. k==-3) .AND. Umax+temp>=0.0D0 ) THEN
         ELSEIF ( .NOT.((k==-2 .OR. k==-4) .AND. Umax-temp>=0.0D0) )    &
                  THEN
            Iold = j
            Umax = ABS(temp)
         ENDIF
      ENDDO
      IF ( Umax<=Eps7 ) Iold = 0
      END

! SUBROUTINE PLTRBG               ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! GRADIENT OF THE OBJECTIVE FUNCTION IS SCALED AND REDUCED. LAGRANGE
! MULTIPLIERS ARE DETERMINED. TEST VALUES GMAX AND UMAX ARE COMPUTED.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  II  NC  NUMBER OF CURRENT LINEAR CONSTRAINTS.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  II  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RI  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RU  CZ(NF*NF)  MATRIX WHOSE COLUMNS ARE BASIC VECTORS FROM THE
!         CURRENT REDUCED SUBSPACE. VECTOR CZ(1,NF) CONTAINS LAGRANGE
!         MULTIPLIERS BEING DETERMINED.
!  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RO  GN(NF)  TRANSFORMED GRADIENT OF THE OBJECTIVE FUNCTION IF IT IS
!         NONZERO.
!  RI  EPS7  TOLERANCE FOR LINEAR AND QUADRATIC PROGRAMMING.
!  RO  GMAX  NORM OF THE TRANSFORMED GRADIENT.
!  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF THE NEGATIVE LAGRANGE MULTIPLIER.
!  IO  IOLD  INDEX OF THE REMOVED CONSTRAINT.
!
! SUBPROGRAMS USED :
!  S   PLVLAG  GRADIENT IS PREMULTIPLIED BY THE MATRIX WHOSE COLUMNS
!         ARE NORMALS OF THE ACTIVE CONSTRAINTS.
!  S   PLTLAG  COMPUTATION OF THE MAXIMUM ABSOLUTE VALUE OF THE NEGATIVE
!         LAGRANGE MULTIPLIER.
!  S   MXDRMM  PREMULTIPLICATION OF A VECTOR BY A ROWWISE STORED DENSE
!         RECTANGULAR MATRIX.
!  S   MXDPRB  BACK SUBSTITUTION AFTER A CHOLESKI DECOMPOSITION.
!  RF  MXVMAX  L-INFINITY NORM OF A VECTOR.
!  S   MXVSET  INITIATION OF A VECTOR.
!
      SUBROUTINE PLTRBG(Nf,N,Nc,Ix,Ic,Ica,Cg,Cr,Cz,G,Gn,Eps7,Gmax,Umax, &
                        Iold)
      IMPLICIT NONE
      INTEGER Nf , N , Nc , Ix(*) , Ic(*) , Ica(*) , Iold
      DOUBLE PRECISION Cg(*) , Cr(*) , Cz(*) , G(*) , Gn(*) , Eps7 ,    &
                       Gmax , Umax
      DOUBLE PRECISION MXVMAX
      INTEGER nca , ncz
      Gmax = 0.0D0
      IF ( N>0 ) THEN
         CALL MXDRMM(Nf,N,Cz,G,Gn)
         Gmax = MXVMAX(N,Gn)
      ENDIF
      IF ( Nf>N .AND. Gmax<=Eps7 ) THEN
         nca = Nf - N
         ncz = N*Nf
         CALL PLVLAG(Nf,N,Nc,Ica,Cg,Cg,G,Cz(ncz+1))
         CALL MXDPRB(nca,Cr,Cz(ncz+1),0)
         CALL PLTLAG(Nf,N,Nc,Ix,Ic,Ica,Cz(ncz+1),Ic,Eps7,Umax,Iold)
         IF ( Umax<=Eps7 ) Iold = 0
         CALL MXVSET(N,0.0D0,Gn)
         Gmax = 0.0D0
      ELSE
         Iold = 0
         Umax = 0.0D0
      ENDIF
      END

! SUBROUTINE PLVLAG               ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! GRADIENT OF THE OBJECTIVE FUNCTION IS PREMULTIPLIED BY TRANSPOSE
! OF THE MATRIX WHOSE COLUMNS ARE NORMALS OF CURRENT ACTIVE CONSTRAINTS
! AND GRADIENTS OF CURRENT ACTIVE FUNCTIONS.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  II  NC  NUMBER OF LINEARIZED CONSTRAINTS.
!  II  IAA(NF+1)  VECTOR CONTAINING INDICES OF ACTIVE FUNCTIONS.
!  RI  AG(NF*NA)  VECTOR CONTAINING SCALING PARAMETERS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RO  GN(NF+1)  OUTPUT VECTOR.
!
! SUBPROGRAMS USED :
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!
      SUBROUTINE PLVLAG(Nf,N,Nc,Iaa,Ag,Cg,G,Gn)
      IMPLICIT NONE
      INTEGER Nf , N , Nc , Iaa(*)
      DOUBLE PRECISION Ag(*) , Cg(*) , G(*) , Gn(*)
      DOUBLE PRECISION MXVDOT
      INTEGER naa , j , l
      naa = Nf - N
      DO j = 1 , naa
         l = Iaa(j)
         IF ( l>Nc ) THEN
            l = l - Nc
            Gn(j) = MXVDOT(Nf,Ag((l-1)*Nf+1),G)
         ELSEIF ( l>0 ) THEN
            Gn(j) = MXVDOT(Nf,Cg((l-1)*Nf+1),G)
         ELSE
            l = -l
            Gn(j) = G(l)
         ENDIF
      ENDDO
      END

! SUBROUTINE PNINT1                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! EXTRAPOLATION OR INTERPOLATION FOR LINE SEARCH WITH DIRECTIONAL
! DERIVATIVES.
!
! PARAMETERS :
!  RI  RL  LOWER VALUE OF THE STEPSIZE PARAMETER.
!  RI  RU  UPPER VALUE OF THE STEPSIZE PARAMETER.
!  RI  FL  VALUE OF THE OBJECTIVE FUNCTION FOR R=RL.
!  RI  FU  VALUE OF THE OBJECTIVE FUNCTION FOR R=RU.
!  RI  PL  DIRECTIONAL DERIVATIVE FOR R=RL.
!  RI  PU  DIRECTIONAL DERIVATIVE FOR R=RU.
!  RO  R  VALUE OF THE STEPSIZE PARAMETER OBTAINED.
!  II  MODE  MODE OF LINE SEARCH.
!  II  MTYP  METHOD SELECTION. MTYP=1-BISECTION. MTYP=2-QUADRATIC
!         INTERPOLATION (WITH ONE DIRECTIONAL DERIVATIVE).
!         MTYP=3-QUADRATIC INTERPOLATION (WITH TWO DIRECTIONAL
!         DERIVATIVES). MTYP=4-CUBIC INTERPOLATION. MTYP=5-CONIC
!         INTERPOLATION.
!  IO  MERR  ERROR INDICATOR. MERR=0 FOR NORMAL RETURN.
!
! METHOD :
! EXTRAPOLATION OR INTERPOLATION WITH STANDARD MODEL FUNCTIONS.
!
      SUBROUTINE PNINT1(Rl,Ru,Fl,Fu,Pl,Pu,R,Mode,Mtyp,Merr)
      IMPLICIT NONE
      DOUBLE PRECISION Rl , Ru , Fl , Fu , Pl , Pu , R
      INTEGER Mode , Mtyp , Merr , ntyp
      DOUBLE PRECISION a , b , c , d , dis , den
      DOUBLE PRECISION C1L , C1U , C2L , C2U , C3L
      PARAMETER (C1L=1.1D0,C1U=1.0D3,C2L=1.0D-2,C2U=0.9D0,C3L=0.1D0)
      Merr = 0
      IF ( Mode<=0 ) RETURN
      IF ( Pl>=0.0D0 ) THEN
         Merr = 2
         RETURN
      ELSEIF ( Ru<=Rl ) THEN
         Merr = 3
         RETURN
      ENDIF
      DO ntyp = Mtyp , 1 , -1
         IF ( ntyp==1 ) THEN
!
!     BISECTION
!
            IF ( Mode==1 ) THEN
               R = 4.0D0*Ru
               RETURN
            ELSE
               R = 0.5D0*(Rl+Ru)
               RETURN
            ENDIF
         ELSEIF ( ntyp==Mtyp ) THEN
            a = (Fu-Fl)/(Pl*(Ru-Rl))
            b = Pu/Pl
         ENDIF
         IF ( ntyp==2 ) THEN
!
!     QUADRATIC EXTRAPOLATION OR INTERPOLATION WITH ONE DIRECTIONAL
!     DERIVATIVE
!
            den = 2.0D0*(1.0D0-a)
         ELSEIF ( ntyp==3 ) THEN
!
!     QUADRATIC EXTRAPOLATION OR INTERPOLATION WITH TWO DIRECTIONAL
!     DERIVATIVES
!
            den = 1.0D0 - b
         ELSEIF ( ntyp==4 ) THEN
!
!     CUBIC EXTRAPOLATION OR INTERPOLATION
!
            c = b - 2.0D0*a + 1.0D0
            d = b - 3.0D0*a + 2.0D0
            dis = d*d - 3.0D0*c
            IF ( dis<0.0D0 ) GOTO 100
            den = d + SQRT(dis)
         ELSEIF ( ntyp==5 ) THEN
!
!     CONIC EXTRAPOLATION OR INTERPOLATION
!
            dis = a*a - b
            IF ( dis<0.0D0 ) GOTO 100
            den = a + SQRT(dis)
            IF ( den<=0.0D0 ) GOTO 100
            den = 1.0D0 - b*(1.0D0/den)**3
         ENDIF
         IF ( Mode==1 .AND. den>0.0D0 .AND. den<1.0D0 ) THEN
!
!     EXTRAPOLATION ACCEPTED
!
            R = Rl + (Ru-Rl)/den
            R = MAX(R,C1L*Ru)
            R = MIN(R,C1U*Ru)
            RETURN
         ELSEIF ( Mode==2 .AND. den>1.0D0 ) THEN
!
!     INTERPOLATION ACCEPTED
!
            R = Rl + (Ru-Rl)/den
            IF ( Rl==0.0D0 ) THEN
               R = MAX(R,Rl+C2L*(Ru-Rl))
            ELSE
               R = MAX(R,Rl+C3L*(Ru-Rl))
            ENDIF
            R = MIN(R,Rl+C2U*(Ru-Rl))
            RETURN
         ENDIF
 100  ENDDO
      END

! SUBROUTINE PNINT3                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! EXTRAPOLATION OR INTERPOLATION FOR LINE SEARCH WITHOUT DIRECTIONAL
! DERIVATIVES.
!
! PARAMETERS :
!  RI  RO  INITIAL VALUE OF THE STEPSIZE PARAMETER.
!  RI  RL  LOWER VALUE OF THE STEPSIZE PARAMETER.
!  RI  RU  UPPER VALUE OF THE STEPSIZE PARAMETER.
!  RI  RI  INNER VALUE OF THE STEPSIZE PARAMETER.
!  RI  FO  VALUE OF THE OBJECTIVE FUNCTION FOR R=RO.
!  RI  FL  VALUE OF THE OBJECTIVE FUNCTION FOR R=RL.
!  RI  FU  VALUE OF THE OBJECTIVE FUNCTION FOR R=RU.
!  RI  FI  VALUE OF THE OBJECTIVE FUNCTION FOR R=RI.
!  RO  PO  INITIAL VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RO  R  VALUE OF THE STEPSIZE PARAMETER OBTAINED.
!  II  MODE  MODE OF LINE SEARCH.
!  II  MTYP  METHOD SELECTION. MTYP=1-BISECTION. MTYP=2-TWO POINT
!         QUADRATIC INTERPOLATION. MTYP=2-THREE POINT QUADRATIC
!         INTERPOLATION.
!  IO  MERR  ERROR INDICATOR. MERR=0 FOR NORMAL RETURN.
!
! METHOD :
! EXTRAPOLATION OR INTERPOLATION WITH STANDARD MODEL FUNCTIONS.
!
      SUBROUTINE PNINT3(Ro,Rl,Ru,Ri,Fo,Fl,Fu,Fi,Po,R,Mode,Mtyp,Merr)
      IMPLICIT NONE
      DOUBLE PRECISION ZERO , HALF , ONE , TWO , THREE , C1L , C1U ,    &
                       C2L , C2U , C3L
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0, &
                 C1L=1.1D0,C1U=1.0D3,C2L=1.0D-2,C2U=0.9D0,C3L=1.0D-1)
      DOUBLE PRECISION Fi , Fl , Fo , Fu , Po , R , Ri , Rl , Ro , Ru
      INTEGER Merr , Mode , Mtyp
      DOUBLE PRECISION ai , al , au , den , dis
      INTEGER ntyp
      LOGICAL l1 , l2
      Merr = 0
      IF ( Mode<=0 ) RETURN
      IF ( Po>=ZERO ) THEN
         Merr = 2
         RETURN

      ELSEIF ( Ru<=Rl ) THEN
         Merr = 3
         RETURN
      ENDIF
      l1 = Rl<=Ro
      l2 = Ri<=Rl
      DO ntyp = Mtyp , 1 , -1
         IF ( ntyp==1 ) THEN
!
!     BISECTION
!
            IF ( Mode==1 ) THEN
               R = TWO*Ru
               RETURN
            ELSEIF ( Ri-Rl<=Ru-Ri ) THEN
               R = HALF*(Ri+Ru)
               RETURN
            ELSE
               R = HALF*(Rl+Ri)
               RETURN
            ENDIF
         ELSEIF ( ntyp==Mtyp .AND. l1 ) THEN
            IF ( .NOT.l2 ) ai = (Fi-Fo)/(Ri*Po)
            au = (Fu-Fo)/(Ru*Po)
         ENDIF
         IF ( l1 .AND. (ntyp==2 .OR. l2) ) THEN
!
!     TWO POINT QUADRATIC EXTRAPOLATION OR INTERPOLATION
!
            IF ( au>=ONE ) GOTO 100
            R = HALF*Ru/(ONE-au)
         ELSEIF ( .NOT.l1 .OR. .NOT.l2 .AND. ntyp==3 ) THEN
!
!     THREE POINT QUADRATIC EXTRAPOLATION OR INTERPOLATION
!
            al = (Fi-Fl)/(Ri-Rl)
            au = (Fu-Fi)/(Ru-Ri)
            den = au - al
            IF ( den<=ZERO ) GOTO 100
            R = Ri - HALF*(au*(Ri-Rl)+al*(Ru-Ri))/den
         ELSEIF ( l1 .AND. .NOT.l2 .AND. ntyp==4 ) THEN
!
!     THREE POINT CUBIC EXTRAPOLATION OR INTERPOLATION
!
            dis = (ai-ONE)*(Ru/Ri)
            den = (au-ONE)*(Ri/Ru) - dis
            dis = au + ai - den - TWO*(ONE+dis)
            dis = den*den - THREE*dis
            IF ( dis<ZERO ) GOTO 100
            den = den + SQRT(dis)
            IF ( den==ZERO ) GOTO 100
            R = (Ru-Ri)/den
         ELSE
            GOTO 100
         ENDIF
         IF ( Mode==1 .AND. R>Ru ) THEN
!
!     EXTRAPOLATION ACCEPTED
!
            R = MAX(R,C1L*Ru)
            R = MIN(R,C1U*Ru)
            RETURN
         ELSEIF ( Mode==2 .AND. R>Rl .AND. R<Ru ) THEN
!
!     INTERPOLATION ACCEPTED
!
            IF ( Ri==ZERO .AND. ntyp/=4 ) THEN
               R = MAX(R,Rl+C2L*(Ru-Rl))
            ELSE
               R = MAX(R,Rl+C3L*(Ru-Rl))
            ENDIF
            R = MIN(R,Rl+C2U*(Ru-Rl))
            IF ( R/=Ri ) RETURN
         ENDIF
 100  ENDDO
      END

! SUBROUTINE PNSTEP                ALL SYSTEMS                89/12/01
! PORTABILITY : ALL SYSTEMS
! 89/01/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DETERMINATION OF A SCALING FACTOR FOR THE BOUNDARY STEP.
!
! PARAMETERS :
!  RI  DEL  MAXIMUM STEPSIZE.
!  RI  A  INPUT PARAMETER.
!  RI  B  INPUT PARAMETER.
!  RI  C  INPUT PARAMETER.
!  RO  ALF  SCALING FACTOR FOR THE BOUNDARY STEP SUCH THAT
!         A**2+2*B*ALF+C*ALF**2=DEL**2.
!
      SUBROUTINE PNSTEP(Del,A,B,C,Alf)
      IMPLICIT NONE
      DOUBLE PRECISION Del , A , B , C , Alf
      DOUBLE PRECISION den , dis
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      Alf = ZERO
      den = (Del+A)*(Del-A)
      IF ( den<=ZERO ) RETURN
      dis = B*B + C*den
      IF ( B>=ZERO ) THEN
         Alf = den/(SQRT(dis)+B)
      ELSE
         Alf = (SQRT(dis)-B)/C
      ENDIF
      END

! SUBROUTINE PP0AF8             ALL SYSTEMS                 97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! COMPUTATION OF VALUE OF THE AUGMENTED LAGRANGIAN FUNCTION.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  II  N  DIMENSION OF THE CONSTRAINT NULL SPACE.
!  II  NC  NUMBER OF CONSTRAINTS.
!  RI  CF(NC+1)  VECTOR CONTAINING VALUES OF THE CONSTRAINTS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  II  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CZ(NC)  VECTOR OF LAGRANGE MULTIPLIERS.
!  RI  RPF  PENALTY COEFFICIENT.
!  RO  FC  VALUE OF THE PENALTY TERM.
!  RO  F  VALUE OF THE PENALTY FUNCTION.
!
      SUBROUTINE PP0AF8(Nf,N,Nc,Cf,Ic,Ica,Cl,Cu,Cz,Rpf,Fc,F)
      IMPLICIT NONE
      INTEGER Nf , N , Nc , Ic(*) , Ica(*)
      DOUBLE PRECISION Cf(*) , Cl(*) , Cu(*) , Cz(*) , Rpf , Fc , F
      DOUBLE PRECISION pom , temp
      INTEGER j , kc
      Fc = 0.0D0
      DO kc = 1 , Nc
         IF ( Ic(kc)>0 ) THEN
            pom = 0.0D0
            temp = Cf(kc)
            IF ( Ic(kc)==1 .OR. Ic(kc)>=3 ) pom = MIN(pom,temp-Cl(kc))
            IF ( Ic(kc)==2 .OR. Ic(kc)>=3 ) pom = MIN(pom,Cu(kc)-temp)
            Fc = Fc + Rpf*ABS(pom)
         ENDIF
      ENDDO
      DO j = 1 , Nf - N
         kc = Ica(j)
         IF ( kc>0 ) THEN
            pom = 0.0D0
            temp = Cf(kc)
            IF ( Ic(kc)==1 .OR. Ic(kc)==3 .OR. Ic(kc)==5 )              &
                 pom = MIN(pom,temp-Cl(kc))
            IF ( Ic(kc)==2 .OR. Ic(kc)==4 .OR. Ic(kc)==6 )              &
                 pom = MAX(pom,temp-Cu(kc))
            Fc = Fc - Cz(j)*pom
         ENDIF
      ENDDO
      F = Cf(Nc+1) + Fc
      END

! SUBROUTINE PPSET2             ALL SYSTEMS                   97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! COMPUTATION OF THE NEW PENALTY PARAMETERS.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  II  NC  NUMBER OF CONSTRAINTS.
!  II  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CZ(NF)  VECTOR OF LAGRANGE MULTIPLIERS.
!  RI  CP(NC)  VECTOR CONTAINING PENALTY PARAMETERS.
!
      SUBROUTINE PPSET2(Nf,N,Nc,Ica,Cz,Cp)
      IMPLICIT NONE
      INTEGER Nf , N , Nc , Ica(*)
      DOUBLE PRECISION Cz(*) , Cp(*)
      DOUBLE PRECISION temp
      INTEGER j , l , kc
      DO kc = 1 , Nc
         Cp(kc) = 0.5D0*Cp(kc)
      ENDDO
      DO j = 1 , Nf - N
         l = Ica(j)
         IF ( l>0 ) THEN
            temp = ABS(Cz(j))
            Cp(l) = MAX(temp,Cp(l)+0.5D0*temp)
         ENDIF
      ENDDO
      END

! SUBROUTINE PS0G01                ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! SIMPLE SEARCH WITH TRUST REGION UPDATE.
!
! PARAMETERS :
!  RO  R  VALUE OF THE STEPSIZE PARAMETER.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  INITIAL VALUE OF THE OBJECTIVE FUNCTION.
!  RI  PO  INITIAL VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RI  PP  QUADRATIC PART OF THE PREDICTED FUNCTION VALUE.
!  RU  XDEL  TRUST REGION BOUND.
!  RO  XDELO  PREVIOUS TRUST REGION BOUND.
!  RI  XMAX MAXIMUM STEPSIZE.
!  RI  RMAX  MAXIMUM VALUE OF THE STEPSIZE PARAMETER.
!  RI  SNORM  EUCLIDEAN NORM OF THE DIRECTION VECTOR.
!  RI  BET1  LOWER BOUND FOR STEPSIZE REDUCTION.
!  RI  BET2  UPPER BOUND FOR STEPSIZE REDUCTION.
!  RI  GAM1  LOWER BOUND FOR STEPSIZE EXPANSION.
!  RI  GAM2  UPPER BOUND FOR STEPSIZE EXPANSION.
!  RI  EPS4  FIRST TOLERANCE FOR RATIO DF/DFPRED. STEP BOUND IS
!         DECREASED IF DF/DFPRED<EPS4.
!  RI  EPS5  SECOND TOLERANCE FOR RATIO DF/DFPRED. STEP BOUND IS
!         INCREASED IF IT IS ACTIVE AND DF/DFPRED>EPS5.
!  II  KD  DEGREE OF REQUIRED DERVATIVES.
!  IO  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!  IU  IDIR INDICATOR FOR DIRECTION DETERMINATION.
!         IDIR=0-BASIC DETERMINATION. IDIR=1-DETERMINATION
!         AFTER STEPSIZE REDUCTION. IDIR=2-DETERMINATION AFTER
!         STEPSIZE EXPANSION.
!  IO  ITERS  TERMINATION INDICATOR. ITERS=0-ZERO STEP. ITERS=1-STEP
!         BOUND WAS DECREASED. ITERS=2-STEP BOUND WAS UNCHANGED.
!         ITERS=3-STEP BOUND WAS INCREASED. ITERS=6-FIRST STEPSIZE.
!  II  ITERD TERMINATION INDICATOR. ITERD<0-BAD DECOMPOSITION.
!         ITERD=0-DESCENT DIRECTION. ITERD=1-NEWTON LIKE STEP.
!         ITERD=2-INEXACT NEWTON LIKE STEP. ITERD=3-BOUNDARY STEP.
!         ITERD=4-DIRECTION WITH THE NEGATIVE CURVATURE.
!         ITERD=5-MARQUARDT STEP.
!  IO  MAXST  MAXIMUM STEPSIZE INDICATOR. MAXST=0 OR MAXST=1 IF MAXIMUM
!         STEPSIZE WAS NOT OR WAS REACHED.
!  IO  NRED  ACTUAL NUMBER OF EXTRAPOLATIONS OR INTERPOLATIONS.
!  II  MRED  MAXIMUM NUMBER OF EXTRAPOLATIONS OR INTERPOLATIONS.
!  II  KTERS  TERMINATION SELECTION. KTERS=1-NORMAL TERMINATION.
!         KTERS=6-FIRST STEPSIZE.
!  II  MES1  SWITCH FOR EXTRAPOLATION. MES1=1-CONSTANT INCREASING OF
!         THE INTERVAL. MES1=2-EXTRAPOLATION SPECIFIED BY THE PARAMETER
!         MES. MES1=3 SUPPRESSED EXTRAPOLATION.
!  II  MES2  SWITCH FOR TERMINATION. MES2=1-NORMAL TERMINATION.
!         MES2=2-TERMINATION AFTER AT LEAST TWO STEPS (ASYMPTOTICALLY
!         PERFECT LINE SEARCH).
!  II  MES3  SAFEGUARD AGAINST ROUNDING ERRORS. MES3=0-SAFEGUARD
!         SUPPRESSED. MES3=1-FIRST LEVEL OF SAFEGUARD. MES3=2-SECOND
!         LEVEL OF SAFEGUARD.
!  IU  ISYS  CONTROL PARAMETER.
!
! COMMON DATA :
!
! METHOD :
! G.A.SCHULTZ, R.B.SCHNABEL, R.H.BYRD: A FAMILY OF TRUST-REGION-BASED
! ALGORITHMS FOR UNCONSTRAINED MINIMIZATION WITH STRONG GLOBAL
! CONVERGENCE PROPERTIES, SIAM J. NUMER.ANAL. 22 (1985) PP. 47-67.
!
      SUBROUTINE PS0G01(R,F,Fo,Po,Pp,Xdel,Xdelo,Xmax,Rmax,Snorm,Bet1,   &
                        Bet2,Gam1,Gam2,Eps4,Eps5,Kd,Ld,Idir,Iters,Iterd,&
                        Maxst,Nred,Mred,Kters,Mes1,Mes2,Mes3,Isys)
      IMPLICIT NONE
      INTEGER Kd , Ld , Idir , Iters , Iterd , Maxst , Nred , Mred ,    &
              Kters , Mes1 , Mes2 , Mes3 , Isys
      DOUBLE PRECISION R , F , Fo , Po , Pp , Xdel , Xdelo , Xmax ,     &
                       Rmax , Snorm , Bet1 , Bet2 , Gam1 , Gam2 , Eps4 ,&
                       Eps5
      DOUBLE PRECISION df , dfpred
      INTEGER nred1 , nred2
      SAVE nred1 , nred2
      IF ( Isys==1 ) THEN
         IF ( Kters<0 .OR. Kters>5 ) THEN
            Iters = 6
         ELSE
            df = Fo - F
            dfpred = -R*(Po+R*Pp)
            IF ( df<Eps4*dfpred ) THEN
!
!     STEP IS TOO LARGE, IT HAS TO BE REDUCED
!
               IF ( Mes1==1 ) THEN
                  Xdel = Bet2*Snorm
               ELSEIF ( Mes1==2 ) THEN
                  Xdel = Bet2*MIN(0.5D0*Xdel,Snorm)
               ELSE
                  Xdel = 0.5D0*Po*Snorm/(Po+df)
                  Xdel = MAX(Xdel,Bet1*Snorm)
                  Xdel = MIN(Xdel,Bet2*Snorm)
               ENDIF
               Iters = 1
               IF ( Mes3<=1 ) THEN
                  nred2 = nred2 + 1
               ELSE
                  IF ( Iterd>2 ) nred2 = nred2 + 1
               ENDIF
            ELSEIF ( df<=Eps5*dfpred ) THEN
!
!     STEP IS SUITABLE
!
               Iters = 2
            ELSE
!
!     STEP IS TOO SMALL, IT HAS TO BE ENLARGED
!
               IF ( Mes2==2 ) THEN
                  Xdel = MAX(Xdel,Gam1*Snorm)
               ELSEIF ( Iterd>2 ) THEN
                  Xdel = Gam1*Xdel
               ENDIF
               Iters = 3
            ENDIF
            Xdel = MIN(Xdel,Xmax,Gam2*Snorm)
            IF ( Fo<=F ) THEN
               IF ( nred1>=Mred ) THEN
                  Iters = -1
               ELSE
                  Idir = 1
                  Iters = 0
                  nred1 = nred1 + 1
               ENDIF
            ENDIF
         ENDIF
         Maxst = 0
         IF ( Xdel>=Xmax ) Maxst = 1
         IF ( Mes3==0 ) THEN
            Nred = nred1
         ELSE
            Nred = nred2
         ENDIF
         Isys = 0
         GOTO 99999
      ENDIF
!      GO TO (1,2) ISYS+1
      IF ( Idir==0 ) THEN
         nred1 = 0
         nred2 = 0
      ENDIF
      Idir = 0
      Xdelo = Xdel
!
!     COMPUTATION OF THE NEW FUNCTION VALUE
!
      R = MIN(1.0D0,Rmax)
      Kd = 0
      Ld = -1
      Isys = 1
      RETURN
99999 END

! SUBROUTINE PS0L02                ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
!  EXTENDED LINE SEARCH WITHOUT DIRECTIONAL DERIVATIVES.
!
! PARAMETERS :
!  RO  R  VALUE OF THE STEPSIZE PARAMETER.
!  RO  RO  INITIAL VALUE OF THE STEPSIZE PARAMETER.
!  RO  RP  PREVIOUS VALUE OF THE STEPSIZE PARAMETER.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  INITIAL VALUE OF THE OBJECTIVE FUNCTION.
!  RO  FP  PREVIOUS VALUE OF THE OBJECTIVE FUNCTION.
!  RI  PO  INITIAL VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RO  PP  PREVIOUS VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RI  FMIN  LOWER BOUND FOR VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FMAX  UPPER BOUND FOR VALUE OF THE OBJECTIVE FUNCTION.
!  RI  RMIN  MINIMUM VALUE OF THE STEPSIZE PARAMETER
!  RI  RMAX  MAXIMUM VALUE OF THE STEPSIZE PARAMETER
!  RI  TOLS  TERMINATION TOLERANCE FOR LINE SEARCH (IN TEST ON THE
!         CHANGE OF THE FUNCTION VALUE).
!  II  KD  DEGREE OF REQUIRED DERVATIVES.
!  IO  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!  II  NIT  ACTUAL NUMBER OF ITERATIONS.
!  II  KIT  NUMBER OF THE ITERATION AFTER LAST RESTART.
!  IO  NRED  ACTUAL NUMBER OF EXTRAPOLATIONS OR INTERPOLATIONS.
!  II  MRED  MAXIMUM NUMBER OF EXTRAPOLATIONS OR INTERPOLATIONS.
!  IO  MAXST  MAXIMUM STEPSIZE INDICATOR. MAXST=0 OR MAXST=1 IF MAXIMUM
!         STEPSIZE WAS NOT OR WAS REACHED.
!  II  IEST  LOWER BOUND SPECIFICATION. IEST=0 OR IEST=1 IF LOWER BOUND
!         IS NOT OR IS GIVEN.
!  II  INITS  CHOICE OF THE INITIAL STEPSIZE. INITS=0-INITIAL STEPSIZE
!         IS SPECIFIED IN THE CALLING PROGRAM. INITS=1-UNIT INITIAL
!         STEPSIZE. INITS=2-COMBINED UNIT AND QUADRATICALLY ESTIMATED
!         INITIAL STEPSIZE. INITS=3-QUADRATICALLY ESTIMATED INITIAL
!         STEPSIZE.
!  IO  ITERS  TERMINATION INDICATOR. ITERS=0-ZERO STEP. ITERS=1-PERFECT
!         LINE SEARCH. ITERS=2 GOLDSTEIN STEPSIZE. ITERS=3-CURRY
!         STEPSIZE. ITERS=4-EXTENDED CURRY STEPSIZE.
!         ITERS=5-ARMIJO STEPSIZE. ITERS=6-FIRST STEPSIZE.
!         ITERS=7-MAXIMUM STEPSIZE. ITERS=8-UNBOUNDED FUNCTION.
!         ITERS=-1-MRED REACHED. ITERS=-2-POSITIVE DIRECTIONAL
!         DERIVATIVE. ITERS=-3-ERROR IN INTERPOLATION.
!  II  KTERS  TERMINATION SELECTION. KTERS=1-PERFECT LINE SEARCH.
!         KTERS=2-GOLDSTEIN STEPSIZE. KTERS=3-CURRY STEPSIZE.
!         KTERS=4-EXTENDED CURRY STEPSIZE. KTERS=5-ARMIJO STEPSIZE.
!         KTERS=6-FIRST STEPSIZE.
!  II  MES  METHOD SELECTION. MES=1-BISECTION. MES=2-QUADRATIC
!         INTERPOLATION (WITH ONE DIRECTIONAL DERIVATIVE).
!         MES=3-QUADRATIC INTERPOLATION (WITH TWO DIRECTIONAL
!         DERIVATIVES). MES=4-CUBIC INTERPOLATION. MES=5-CONIC
!         INTERPOLATION.
!  IU  ISYS  CONTROL PARAMETER.
!
! SUBPROGRAM USED :
!  S   PNINT3  EXTRAPOLATION OR INTERPOLATION WITHOUT DIRECTIONAL
!         DERIVATIVES.
!
! METHOD :
! SAFEGUARDED EXTRAPOLATION AND INTERPOLATION WITH EXTENDED TERMINATION
! CRITERIA.
!
      SUBROUTINE PS0L02(R,Ro,Rp,F,Fo,Fp,Po,Pp,Fmin,Fmax,Rmin,Rmax,Tols, &
                        Kd,Ld,Nit,Kit,Nred,Mred,Maxst,Iest,Inits,Iters, &
                        Kters,Mes,Isys)
      IMPLICIT NONE
      INTEGER Kd , Ld , Nit , Kit , Nred , Mred , Maxst , Iest , Inits ,&
              Iters , Kters , Mes , Isys
      DOUBLE PRECISION R , Ro , Rp , F , Fo , Fp , Po , Pp , Fmin ,     &
                       Fmax , Rmin , Rmax , Tols
      DOUBLE PRECISION rl , fl , ru , fu , ri , fi , rtemp , TOL
      INTEGER mtyp , merr , mode , init1 , mes1 , mes2
      LOGICAL l1 , l2 , l3 , l4 , l6 , l7
      PARAMETER (TOL=1.0D-4)
      SAVE mtyp , mode , mes1 , mes2
      SAVE rl , fl , ru , fu , ri , fi
      IF ( Isys/=1 ) THEN
!      GO TO (1,3) ISYS+1
         mes1 = 2
         mes2 = 2
         Iters = 0
         IF ( Po>=0.0D0 ) THEN
            R = 0.0D0
            Iters = -2
            Isys = 0
            GOTO 99999
         ENDIF
         IF ( Rmax<=0.0D0 ) THEN
            Iters = 0
            Isys = 0
            GOTO 99999
         ENDIF
!
!     INITIAL STEPSIZE SELECTION
!
         IF ( Inits>0 ) THEN
            rtemp = Fmin - F
         ELSEIF ( Iest==0 ) THEN
            rtemp = F - Fp
         ELSE
            rtemp = MAX(F-Fp,1.0D1*(Fmin-F))
         ENDIF
         init1 = ABS(Inits)
         Rp = 0.0D0
         Fp = Fo
         Pp = Po
         IF ( init1==0 ) THEN
         ELSEIF ( init1==1 .OR. Inits>=1 .AND. Iest==0 ) THEN
            R = 1.0D0
         ELSEIF ( init1==2 ) THEN
            R = MIN(1.0D0,4.0D0*rtemp/Po)
         ELSEIF ( init1==3 ) THEN
            R = MIN(1.0D0,2.0D0*rtemp/Po)
         ELSEIF ( init1==4 ) THEN
            R = 2.0D0*rtemp/Po
         ENDIF
         rtemp = R
         R = MAX(R,Rmin)
         R = MIN(R,Rmax)
         mode = 0
         rl = 0.0D0
         fl = Fo
         ru = 0.0D0
         fu = Fo
         ri = 0.0D0
         fi = Fo
      ELSEIF ( Iters/=0 ) THEN
         Isys = 0
         GOTO 99999
      ELSE
         IF ( F<=Fmin ) THEN
            Iters = 7
            Isys = 0
            GOTO 99999
         ELSE
            l1 = R<=Rmin .AND. Nit/=Kit
            l2 = R>=Rmax
            l3 = F - Fo<=Tols*R*Po .OR. F - Fmin<=(Fo-Fmin)/1.0D1
            l4 = F - Fo>=(1.0D0-Tols)*R*Po .OR. mes2==2 .AND. mode==2
            l6 = ru - rl<=TOL*ru .AND. mode==2
            l7 = mes2<=2 .OR. mode/=0
            Maxst = 0
            IF ( l2 ) Maxst = 1
         ENDIF
!
!     TEST ON TERMINATION
!
         IF ( l1 .AND. .NOT.l3 ) THEN
            Iters = 0
            Isys = 0
            GOTO 99999
         ELSEIF ( l2 .AND. .NOT.F>=fu ) THEN
            Iters = 7
            Isys = 0
            GOTO 99999
         ELSEIF ( l6 ) THEN
            Iters = 1
            Isys = 0
            GOTO 99999
         ELSEIF ( l3 .AND. l7 .AND. Kters==5 ) THEN
            Iters = 5
            Isys = 0
            GOTO 99999
         ELSEIF ( l3 .AND. l4 .AND. l7 .AND.                            &
                  (Kters==2 .OR. Kters==3 .OR. Kters==4) ) THEN
            Iters = 2
            Isys = 0
            GOTO 99999
         ELSEIF ( Kters<0 .OR. Kters==6 .AND. l7 ) THEN
            Iters = 6
            Isys = 0
            GOTO 99999
         ELSEIF ( ABS(Nred)>=Mred ) THEN
            Iters = -1
            Isys = 0
            GOTO 99999
         ELSE
            Rp = R
            Fp = F
            mode = MAX(mode,1)
            mtyp = ABS(Mes)
            IF ( F>=Fmax ) mtyp = 1
         ENDIF
         IF ( mode==1 ) THEN
!
!     INTERVAL CHANGE AFTER EXTRAPOLATION
!
            rl = ri
            fl = fi
            ri = ru
            fi = fu
            ru = R
            fu = F
            IF ( F>=fi ) THEN
               Nred = 0
               mode = 2
            ELSEIF ( mes1==1 ) THEN
               mtyp = 1
            ENDIF
!
!     INTERVAL CHANGE AFTER INTERPOLATION
!
         ELSEIF ( R<=ri ) THEN
            IF ( F<=fi ) THEN
               ru = ri
               fu = fi
               ri = R
               fi = F
            ELSE
               rl = R
               fl = F
            ENDIF
         ELSEIF ( F<=fi ) THEN
            rl = ri
            fl = fi
            ri = R
            fi = F
         ELSE
            ru = R
            fu = F
         ENDIF
      ENDIF
!
!     NEW STEPSIZE SELECTION (EXTRAPOLATION OR INTERPOLATION)
!
      CALL PNINT3(Ro,rl,ru,ri,Fo,fl,fu,fi,Po,R,mode,mtyp,merr)
      IF ( merr>0 ) THEN
         Iters = -merr
         Isys = 0
         GOTO 99999
      ELSEIF ( mode==1 ) THEN
         Nred = Nred - 1
         R = MIN(R,Rmax)
      ELSEIF ( mode==2 ) THEN
         Nred = Nred + 1
      ENDIF
!
!     COMPUTATION OF THE NEW FUNCTION VALUE
!
      Kd = 0
      Ld = -1
      Isys = 1
      RETURN
99999 END

! SUBROUTINE PS1L01                ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
!  STANDARD LINE SEARCH WITH DIRECTIONAL DERIVATIVES.
!
! PARAMETERS :
!  RO  R  VALUE OF THE STEPSIZE PARAMETER.
!  RO  RP  PREVIOUS VALUE OF THE STEPSIZE PARAMETER.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  INITIAL VALUE OF THE OBJECTIVE FUNCTION.
!  RO  FP  PREVIOUS VALUE OF THE OBJECTIVE FUNCTION.
!  RO  P  VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RI  PO  INITIAL VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RO  PP  PREVIOUS VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RI  FMIN  LOWER BOUND FOR VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FMAX  UPPER BOUND FOR VALUE OF THE OBJECTIVE FUNCTION.
!  RI  RMIN  MINIMUM VALUE OF THE STEPSIZE PARAMETER
!  RI  RMAX  MAXIMUM VALUE OF THE STEPSIZE PARAMETER
!  RI  TOLS  TERMINATION TOLERANCE FOR LINE SEARCH (IN TEST ON THE
!         CHANGE OF THE FUNCTION VALUE).
!  RI  TOLP  TERMINATION TOLERANCE FOR LINE SEARCH (IN TEST ON THE
!         CHANGE OF THE DIRECTIONAL DERIVATIVE).
!  RO  PAR1  PARAMETER FOR CONTROLLED SCALING OF VARIABLE METRIC
!         UPDATES.
!  RO  PAR2  PARAMETER FOR CONTROLLED SCALING OF VARIABLE METRIC
!         UPDATES.
!  II  KD  DEGREE OF REQUIRED DERVATIVES.
!  IO  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!  II  NIT  ACTUAL NUMBER OF ITERATIONS.
!  II  KIT  NUMBER OF THE ITERATION AFTER LAST RESTART.
!  IO  NRED  ACTUAL NUMBER OF EXTRAPOLATIONS OR INTERPOLATIONS.
!  II  MRED  MAXIMUM NUMBER OF EXTRAPOLATIONS OR INTERPOLATIONS.
!  IO  MAXST  MAXIMUM STEPSIZE INDICATOR. MAXST=0 OR MAXST=1 IF MAXIMUM
!         STEPSIZE WAS NOT OR WAS REACHED.
!  II  IEST  LOWER BOUND SPECIFICATION. IEST=0 OR IEST=1 IF LOWER BOUND
!         IS NOT OR IS GIVEN.
!  II  INITS  CHOICE OF THE INITIAL STEPSIZE. INITS=0-INITIAL STEPSIZE
!         IS SPECIFIED IN THE CALLING PROGRAM. INITS=1-UNIT INITIAL
!         STEPSIZE. INITS=2-COMBINED UNIT AND QUADRATICALLY ESTIMATED
!         INITIAL STEPSIZE. INITS=3-QUADRATICALLY ESTIMATED INITIAL
!         STEPSIZE.
!  IO  ITERS  TERMINATION INDICATOR. ITERS=0-ZERO STEP. ITERS=1-PERFECT
!         LINE SEARCH. ITERS=2 GOLDSTEIN STEPSIZE. ITERS=3-CURRY
!         STEPSIZE. ITERS=4-EXTENDED CURRY STEPSIZE.
!         ITERS=5-ARMIJO STEPSIZE. ITERS=6-FIRST STEPSIZE.
!         ITERS=7-MAXIMUM STEPSIZE. ITERS=8-UNBOUNDED FUNCTION.
!         ITERS=-1-MRED REACHED. ITERS=-2-POSITIVE DIRECTIONAL
!         DERIVATIVE. ITERS=-3-ERROR IN INTERPOLATION.
!  II  KTERS  TERMINATION SELECTION. KTERS=1-PERFECT LINE SEARCH.
!         KTERS=2-GOLDSTEIN STEPSIZE. KTERS=3-CURRY STEPSIZE.
!         KTERS=4-EXTENDED CURRY STEPSIZE. KTERS=5-ARMIJO STEPSIZE.
!         KTERS=6-FIRST STEPSIZE.
!  II  MES  METHOD SELECTION. MES=1-BISECTION. MES=2-QUADRATIC
!         INTERPOLATION (WITH ONE DIRECTIONAL DERIVATIVE).
!         MES=3-QUADRATIC INTERPOLATION (WITH TWO DIRECTIONAL
!         DERIVATIVES). MES=4-CUBIC INTERPOLATION. MES=5-CONIC
!         INTERPOLATION.
!  IU  ISYS  CONTROL PARAMETER.
!
! SUBPROGRAM USED :
!  S   PNINT1  EXTRAPOLATION OR INTERPOLATION WITH DIRECTIONAL
!         DERIVATIVES.
!
! METHOD :
! SAFEGUARDED EXTRAPOLATION AND INTERPOLATION WITH STANDARD TERMINATION
! CRITERIA.
!
      SUBROUTINE PS1L01(R,Rp,F,Fo,Fp,P,Po,Pp,Fmin,Fmax,Rmin,Rmax,Tols,  &
                        Tolp,Par1,Par2,Kd,Ld,Nit,Kit,Nred,Mred,Maxst,   &
                        Iest,Inits,Iters,Kters,Mes,Isys)
      IMPLICIT NONE
      INTEGER Kd , Ld , Nit , Kit , Nred , Mred , Maxst , Iest , Inits ,&
              Iters , Kters , Mes , Isys
      DOUBLE PRECISION R , Rp , F , Fo , Fp , P , Po , Pp , Fmin ,      &
                       Fmax , Rmin , Rmax , Tols , Tolp , Par1 , Par2
      DOUBLE PRECISION rl , fl , pl , ru , fu , pu , rtemp
      INTEGER mtyp , merr , mode , init1 , mes1 , mes2 , mes3
      LOGICAL l1 , l2 , l3 , l5 , l7 , m1 , m2 , m3
      DOUBLE PRECISION CON , CON1
      PARAMETER (CON=1.0D-2,CON1=1.0D-13)
      SAVE mtyp , mode , mes1 , mes2 , mes3
      SAVE rl , fl , pl , ru , fu , pu
      IF ( Isys==1 ) THEN
         IF ( mode==0 ) THEN
            Par1 = P/Po
            Par2 = F - Fo
         ENDIF
         IF ( Iters/=0 ) THEN
            Isys = 0
            GOTO 99999
         ELSE
            IF ( F<=Fmin ) THEN
               Iters = 7
               Isys = 0
               GOTO 99999
            ELSE
               l1 = R<=Rmin .AND. Nit/=Kit
               l2 = R>=Rmax
               l3 = F - Fo<=Tols*R*Po
               l5 = P>=Tolp*Po .OR. mes2==2 .AND. mode==2
               l7 = mes2<=2 .OR. mode/=0
               m1 = .FALSE.
               m2 = .FALSE.
               m3 = l3
               IF ( mes3>=1 ) THEN
                  m1 = ABS(P)<=CON*ABS(Po) .AND. Fo - F>=(CON1/CON)     &
                       *ABS(Fo)
                  l3 = l3 .OR. m1
               ENDIF
               IF ( mes3>=2 ) THEN
                  m2 = ABS(P)<=0.5D0*ABS(Po) .AND. ABS(Fo-F)            &
                       <=2.0D0*CON1*ABS(Fo)
                  l3 = l3 .OR. m2
               ENDIF
               Maxst = 0
               IF ( l2 ) Maxst = 1
            ENDIF
!
!     TEST ON TERMINATION
!
            IF ( l1 .AND. .NOT.l3 ) THEN
               Iters = 0
               Isys = 0
               GOTO 99999
            ELSEIF ( l2 .AND. l3 .AND. .NOT.l5 ) THEN
               Iters = 7
               Isys = 0
               GOTO 99999
            ELSEIF ( m3 .AND. mes1==3 ) THEN
               Iters = 5
               Isys = 0
               GOTO 99999
            ELSEIF ( l3 .AND. l5 .AND. l7 ) THEN
               Iters = 4
               Isys = 0
               GOTO 99999
            ELSEIF ( Kters<0 .OR. Kters==6 .AND. l7 ) THEN
               Iters = 6
               Isys = 0
               GOTO 99999
            ELSEIF ( ABS(Nred)>=Mred ) THEN
               Iters = -1
               Isys = 0
               GOTO 99999
            ELSE
               Rp = R
               Fp = F
               Pp = P
               mode = MAX(mode,1)
               mtyp = ABS(Mes)
               IF ( F>=Fmax ) mtyp = 1
            ENDIF
            IF ( mode==1 ) THEN
!
!     INTERVAL CHANGE AFTER EXTRAPOLATION
!
               rl = ru
               fl = fu
               pl = pu
               ru = R
               fu = F
               pu = P
               IF ( .NOT.l3 ) THEN
                  Nred = 0
                  mode = 2
               ELSEIF ( mes1==1 ) THEN
                  mtyp = 1
               ENDIF
!
!     INTERVAL CHANGE AFTER INTERPOLATION
!
            ELSEIF ( .NOT.l3 ) THEN
               ru = R
               fu = F
               pu = P
            ELSE
               rl = R
               fl = F
               pl = P
            ENDIF
         ENDIF
      ELSE
!      GO TO (1,3) ISYS+1
         mes1 = 2
         mes2 = 2
         mes3 = 2
         Iters = 0
         IF ( Po>=0.0D0 ) THEN
            R = 0.0D0
            Iters = -2
            Isys = 0
            GOTO 99999
         ENDIF
         IF ( Rmax<=0.0D0 ) THEN
            Iters = 0
            Isys = 0
            GOTO 99999
         ENDIF
!
!     INITIAL STEPSIZE SELECTION
!
         IF ( Inits>0 ) THEN
            rtemp = Fmin - F
         ELSEIF ( Iest==0 ) THEN
            rtemp = F - Fp
         ELSE
            rtemp = MAX(F-Fp,1.0D1*(Fmin-F))
         ENDIF
         init1 = ABS(Inits)
         Rp = 0.0D0
         Fp = Fo
         Pp = Po
         IF ( init1==0 ) THEN
         ELSEIF ( init1==1 .OR. Inits>=1 .AND. Iest==0 ) THEN
            R = 1.0D0
         ELSEIF ( init1==2 ) THEN
            R = MIN(1.0D0,4.0D0*rtemp/Po)
         ELSEIF ( init1==3 ) THEN
            R = MIN(1.0D0,2.0D0*rtemp/Po)
         ELSEIF ( init1==4 ) THEN
            R = 2.0D0*rtemp/Po
         ENDIF
         R = MAX(R,Rmin)
         R = MIN(R,Rmax)
         mode = 0
         ru = 0.0D0
         fu = Fo
         pu = Po
      ENDIF
!
!     NEW STEPSIZE SELECTION (EXTRAPOLATION OR INTERPOLATION)
!
      CALL PNINT1(rl,ru,fl,fu,pl,pu,R,mode,mtyp,merr)
      IF ( merr>0 ) THEN
         Iters = -merr
         Isys = 0
         GOTO 99999
      ELSEIF ( mode==1 ) THEN
         Nred = Nred - 1
         R = MIN(R,Rmax)
      ELSEIF ( mode==2 ) THEN
         Nred = Nred + 1
      ENDIF
!
!     COMPUTATION OF THE NEW FUNCTION VALUE AND THE NEW DIRECTIONAL
!     DERIVATIVE
!
      Kd = 1
      Ld = -1
      Isys = 1
      RETURN
99999 END

! SUBROUTINE PUDBG1                ALL SYSTEMS                92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! VARIABLE METRIC UPDATE OF A DENSE SYMMETRIC POSITIVE DEFINITE MATRIX
! USING THE FACTORIZATION B=L*D*TRANS(L).
!
! PARAMETERS :
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  RU  H(M)  FACTORIZATION B=L*D*TRANS(L) OF A POSITIVE
!         DEFINITE APPROXIMATION OF THE HESSIAN MATRIX.
!  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RA  S(NF)  AUXILIARY VECTOR.
!  RU  XO(NF)  VECTORS OF VARIABLES DIFFERENCE.
!  RI  GO(NF)  GRADIENTS DIFFERENCE.
!  RI  R  VALUE OF THE STEPSIZE PARAMETER.
!  RI  PO  OLD VALUE OF THE DIRECTIONAL DERIVATIVE.
!  II  NIT  ACTUAL NUMBER OF ITERATIONS.
!  II  KIT  NUMBER OF THE ITERATION AFTER LAST RESTART.
!  IO  ITERH  TERMINATION INDICATOR. ITERH<0-BAD DECOMPOSITION.
!         ITERH=0-SUCCESSFUL UPDATE. ITERH>0-NONPOSITIVE PARAMETERS.
!  II  MET1  SELECTION OF SELF SCALING. MET1=1-SELF SCALING SUPPRESSED.
!         MET1=2 INITIAL SELF SCALING.
!  II  MEC  CORRECTION IF THE NEGATIVE CURVATURE OCCURS.
!         MEC=1-CORRECTION SUPPRESSED. MEC=2-POWELL'S CORRECTION.
!
! SUBPROGRAMS USED :
!  S   MXDPGU  CORRECTION OF A DENSE SYMMETRIC POSITIVE DEFINITE
!         MATRIX IN THE FACTORED FORM B=L*D*TRANS(L).
!  S   MXDPGS  SCALING OF A DENSE SYMMETRIC POSITIVE DEFINITE MATRIX
!         IN THE FACTORED FORM B=L*D*TRANS(L).
!  S   MXVDIF  DIFFERENCE OF TWO VECTORS.
!  RF  MXVDOT  DOT PRODUCT OF VECTORS.
!  S   MXVSCL  SCALING OF A VECTOR.
!
! METHOD :
! BFGS VARIABLE METRIC METHOD.
!
      SUBROUTINE PUDBG1(N,H,G,S,Xo,Go,R,Po,Nit,Kit,Iterh,Met,Met1,Mec)
      IMPLICIT NONE
      DOUBLE PRECISION Po , R
      INTEGER Iterh , Kit , Met , Met1 , Mec , N , Nit
      DOUBLE PRECISION G(*) , Go(*) , H(*) , S(*) , Xo(*)
      DOUBLE PRECISION a , b , c , gam , par , den , dis
      LOGICAL l1 , l3
      DOUBLE PRECISION MXVDOT , MXDPGP
      l1 = Met1>=3 .OR. Met1==2 .AND. Nit==Kit
      l3 = .NOT.l1
!
!     DETERMINATION OF THE PARAMETERS B, C
!
      b = MXVDOT(N,Xo,Go)
      a = 0.0D0
      IF ( l1 ) THEN
         CALL MXVCOP(N,Go,S)
         CALL MXDPGB(N,H,S,1)
         a = MXDPGP(N,H,S,S)
         IF ( a<=0.0D0 ) THEN
            Iterh = 1
            RETURN
         ENDIF
      ENDIF
      CALL MXVDIF(N,Go,G,S)
      CALL MXVSCL(N,R,S,S)
      c = -R*Po
      IF ( c<=0.0D0 ) THEN
         Iterh = 3
         RETURN
      ENDIF
      IF ( Mec>1 ) THEN
         IF ( b<=1.0D-4*c ) THEN
!
!     POWELL'S CORRECTION
!
            dis = (1.0D0-0.1D0)*c/(c-b)
            CALL MXVDIF(N,Go,S,Go)
            CALL MXVDIR(N,dis,Go,S,Go)
            b = c + dis*(b-c)
            IF ( l1 ) a = c + 2.0D0*(1.0D0-dis)*(b-c) + dis*dis*(a-c)
         ENDIF
      ELSEIF ( b<=1.0D-4*c ) THEN
         Iterh = 2
         RETURN
      ENDIF
      IF ( l1 ) THEN
!
!     DETERMINATION OF THE PARAMETER GAM (SELF SCALING)
!
         IF ( Met==1 ) THEN
            par = c/b
         ELSEIF ( a<=0.0D0 ) THEN
            par = c/b
         ELSE
            par = SQRT(c/a)
         ENDIF
         gam = par
         IF ( Met1>1 ) THEN
            IF ( Nit/=Kit ) l3 = gam<0.5D0 .OR. gam>4.0D0
         ENDIF
      ENDIF
      IF ( l3 ) THEN
         gam = 1.0D0
         par = gam
      ENDIF
      IF ( Met==1 ) THEN
!
!     BFGS UPDATE
!
         CALL MXDPGU(N,H,par/b,Go,Xo)
         CALL MXDPGU(N,H,-1.0D0/c,S,Xo)
      ELSE
!
!     HOSHINO UPDATE
!
         den = par*b + c
         dis = 0.5D0*b
         CALL MXVDIR(N,par,Go,S,S)
         CALL MXDPGU(N,H,par/dis,Go,Xo)
         CALL MXDPGU(N,H,-1.0D0/den,S,Xo)
      ENDIF
      Iterh = 0
      IF ( gam==1.0D0 ) RETURN
      CALL MXDPGS(N,H,1.0D0/gam)
      END

! SUBROUTINE PUDBI1                ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! VARIABLE METRIC UPDATES.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  RU  H(N*(N+1)/2)  UPDATED APPROXIMATION OF THE INVERSE HESSIAN
!         MATRIX.
!  RA  S(N)  AUXILIARY VECTOR.
!  RI  XO(N)  VECTOR OF VARIABLES DIFFERENCE.
!  RI  GO(N)  GRADIENT DIFFERENCE.
!  RO  R  VALUE OF THE STEPSIZE PARAMETER.
!  RI  PO  INITIAL VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RO  PAR1  PARAMETER FOR CONTROL SCALING.
!  RO  PAR2  PARAMETER FOR CONTROL SCALING.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  INITIAL VALUE OF THE OBJECTIVE FUNCTION.
!  RI  P  CURRENT VALUE OF THE DIRECTIONAL DERIVATIVE.
!  II  NIT  NUMBER OF ITERATIONS.
!  II  KIT  INDEX OF THE ITERATION WITH THE LAST RESTART.
!  II  MET  VARIABLE METRIC UPDATE.
!  II  MET1  SCALING STRATEGY.
!  II  MET2  CORRECTION RULE.
!  IU  IDECF  DECOMPOSITION INDICATOR.
!  II  ITERD  TERMINATION INDICATOR. ITERD<0-BAD DECOMPOSITION.
!         ITERD=0-DESCENT DIRECTION. ITERD=1-NEWTON LIKE STEP.
!         ITERD=2-INEXACT NEWTON LIKE STEP. ITERD=3-BOUNDARY STEP.
!         ITERD=4-DIRECTION WITH THE NEGATIVE CURVATURE.
!         ITERD=5-MARQUARDT STEP.
!  IO  ITERH  UPDATE INDICATOR. ITERH=0-SUCCESSFUL UPDATE.
!         ITERH>0-UNSUCCESSFUL UPDATE.
!
! METHOD :
! VARIOUS VARIABLE METRIC UPDATES INCLUDING BFGS UPDATE.
!
      SUBROUTINE PUDBI1(N,H,S,Xo,Go,R,Po,Par1,Par2,F,Fo,P,Nit,Kit,Met,  &
                        Met1,Met2,Idecf,Iterd,Iterh)
      IMPLICIT NONE
      INTEGER N , Nit , Kit , Met , Met1 , Met2 , Idecf , Iterd , Iterh
      DOUBLE PRECISION H(*) , S(*) , Xo(*) , Go(*) , R , Po
      DOUBLE PRECISION Par1 , Par2
      DOUBLE PRECISION F , Fo , P
      DOUBLE PRECISION aa , cc
      DOUBLE PRECISION MXVDOT
      DOUBLE PRECISION dis , pom , pom3 , pom4 , a , b , c , gam , rho ,&
                       par
      DOUBLE PRECISION den
      LOGICAL l1 , l2 , l3
      IF ( Met>0 ) THEN
         IF ( Idecf/=9 ) THEN
            Iterh = -1
            GOTO 100
         ENDIF
         l1 = ABS(Met1)>=3 .OR. ABS(Met1)==2 .AND. Nit==Kit
         l3 = .NOT.l1
!
!     DETERMINATION OF THE PARAMETERS A, B, C
!
         b = MXVDOT(N,Xo,Go)
         IF ( b<=0.0D0 ) THEN
            Iterh = 2
            GOTO 100
         ENDIF
         CALL MXDSMM(N,H,Go,S)
         a = MXVDOT(N,Go,S)
         IF ( a<=0.0D0 ) THEN
            Iterh = 1
            GOTO 100
         ENDIF
         IF ( .NOT.(Met/=1 .OR. l1) ) THEN
            c = 0.0D0
         ELSEIF ( Iterd/=1 ) THEN
            c = 0.0D0
         ELSE
            c = -R*Po
            IF ( c<=0.0D0 ) THEN
               Iterh = 3
               GOTO 100
            ENDIF
         ENDIF
!
!     DETERMINATION OF THE PARAMETER RHO (NONQUADRATIC PROPERTIES)
!
         IF ( Met2==1 ) THEN
            rho = 1.0D0
         ELSEIF ( Fo-F+P==0 ) THEN
            rho = 1.0D0
         ELSE
            rho = 0.5D0*b/(Fo-F+P)
         ENDIF
         IF ( rho<=1.0D-2 ) rho = 1.0D0
         IF ( rho>=1.0D2 ) rho = 1.0D0
         aa = a/b
         cc = c/b
         IF ( l1 ) THEN
!
!     DETERMINATION OF THE PARAMETER GAM (SELF SCALING)
!
            par = a/b
            pom3 = 0.7D0
            pom4 = 6.0D0
            gam = rho/par
            IF ( Nit/=Kit ) THEN
               IF ( Met1==3 ) THEN
                  l2 = Par2<=0.0D0
                  l3 = l2 .AND. ABS(Par1)<=0.2D0
                  l3 = l3 .OR. (.NOT.l2 .AND. gam>1.0D0)
                  l3 = l3 .OR. (l2 .AND. Par1<0.0D0 .AND. gam>1.0D0)
                  l3 = l3 .OR. (l2 .AND. Par1>0.0D0 .AND. gam<1.0D0)
                  l3 = l3 .OR. gam<pom3
                  l3 = l3 .OR. gam>pom4
               ELSEIF ( Met1==4 ) THEN
                  l3 = gam<pom3 .OR. gam>pom4
               ENDIF
            ENDIF
         ENDIF
         IF ( l3 ) THEN
            gam = 1.0D0
            par = rho/gam
         ENDIF
         IF ( Met/=1 ) THEN
!
!     NEW UPDATE
!
            pom = 1.0D0/(aa*cc)
            IF ( pom<1.0D0 ) THEN
               pom = MAX(1.0D-15,(SQRT(c/a)-pom)/(1.0D0-pom))
!
!     GENERAL UPDATE
!
               den = par + pom*aa
               dis = pom/den
               CALL MXDSMU(N,H,(par*dis-1.0D0)/a,S)
               CALL MXVDIR(N,-dis,S,Xo,S)
               CALL MXDSMU(N,H,den/b,S)
               GOTO 50
            ENDIF
         ENDIF
!
!     BFGS UPDATE
!
         pom = 1.0D0
         dis = par + aa
         CALL MXVDIR(N,-dis,Xo,S,Xo)
         dis = 1.0D0/(b*dis)
         CALL MXDSMU(N,H,dis,Xo)
         CALL MXDSMU(N,H,-dis,S)
!
!     SCALING
!
 50      IF ( gam/=1.0D0 ) CALL MXDSMS(N,H,gam)
      ENDIF
 100  Iterh = 0
      END

! SUBROUTINE PUDBM2                ALL SYSTEMS                92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! VARIABLE METRIC UPDATE OF A DENSE SYMMETRIC POSITIVE DEFINITE MATRIX.
!
! PARAMETERS :
!  RU  H(M)  POSITIVE DEFINITE APPROXIMATION OF THE HESSIAN
!         MATRIX.
!  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RA  S(NF)  AUXILIARY VECTOR.
!  RU  XO(NF)  VECTORS OF VARIABLES DIFFERENCE.
!  RI  GO(NF)  GRADIENTS DIFFERENCE.
!
! COMMON DATA :
!  II  NF DECLARED NUMBER OF VARIABLES.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  II  M  NUMBER OF NONZERO ELEMENTS OF THE MATRIX.
!  II  MET  METHOD SELECTION. MET=1-BFGS UPDATE. MET=2-DFP UPDATE.
!         MET=3-HOSHINO UPDATE.
!  II  MET1  SELECTION OF SELF SCALING.  MET1=1-SELF SCALING SUPPRESSED.
!         MET1=2-INITIAL SELF SCALING.
!  II  MET2  SELECTION OF THE LINE SEARCH MODEL. MET2=1-QUADRATIC MODEL.
!         MET2=2 USE OF TAYLOR EXPANSION.
!  II  MET3  METHOD CORRECTION. MET3=1-NO CORRECTION.
!         MET3=2-POWELL'S CORRECTION.
!  II  ITERD  TERMINATION INDICATOR. ITERD<0-BAD DECOMPOSITION.
!         ITERD=0-DESCENT DIRECTION. ITERD=1-NEWTON LIKE STEP.
!         ITERD=2-INEXACT NEWTON LIKE STEP. ITERD=3-BOUNDARY STEP.
!         ITERD=4-DIRECTION WITH THE NEGATIVE CURVATURE.
!         ITERD=5-MARQUARDT STEP.
!  IO  ITERH  TERMINATION INDICATOR. ITERH<0-BAD DECOMPOSITION.
!         ITERH=0-SUCCESSFUL UPDATE. ITERH>0-NONPOSITIVE PARAMETERS.
!  II  IDECF  DECOMPOSITION INDICATOR. IDECF=0-NO DECOMPOSITION.
!         IDECF=1-GILL-MURRAY DECOMPOSITION. IDECF=2-BUNCH-PARLETT
!         DECOMPOSITION. IDECF=3-INVERSION.
!  II  ITRAN  TRANSFORMATION INDICATOR. ITRAN=0 OR ITRAN=1 IF
!         TRANSFORMATION IS NOT OR IS USED.
!  II  NIT  ACTUAL NUMBER OF ITERATIONS.
!  II  KIT  NUMBER OF THE ITERATION AFTER LAST RESTART.
!  RI  R  VALUE OF THE STEPSIZE PARAMETER.
!  RI  F  NEW VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  OLD VALUE OF THE OBJECTIVE FUNCTION.
!  RI  P  NEW VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RI  PO  OLD VALUE OF THE DIRECTIONAL DERIVATIVE.
!  TO  TUXX  TEXT INFORMATION ON THE CORRECTION USED.
!
! SUBPROGRAMS USED :
!  S   MXDSMM  MATRIX-VECTOR PRODUCT.
!  S   MXDSMU  CORRECTION OF A DENSE SYMMETRIC MATRIX.
!  S   MXDSMS  SCALING OF A DENSE SYMMETRIC MATRIX.
!  S   MXVDIF  DIFFERENCE OF TWO VECTORS.
!  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
!  RF  MXVDOT  DOT PRODUCT OF VECTORS.
!  S   MXVNEG  COPYING OF A VECTOR WITH THE CHANGE OF THE SIGN.
!  S   MXVSCL  SCALING OF A VECTOR.
!  S   UOU1D1  PRINT OF ENTRY TO VARIABLE METRIC UPDATE.
!  S   UOU1D2  PRINT OF EXIT FROM VARIABLE METRIC UPDATE.
!
! METHOD :
! BASIC VARIABLE METRIC METHODS.
!
      SUBROUTINE PUDBM2(Nf,N,H,Hh,S,Xo,Go,So,Fo,Par,Met1,Met3,Idecf,    &
                        Iterh)
      IMPLICIT NONE
      INTEGER Nf , N , Met1 , Met3 , Idecf , Iterh
      DOUBLE PRECISION H(Nf*(Nf+1)/2) , Hh(Nf*(Nf+1)/2) , S(Nf) , Xo(Nf)&
                       , Go(Nf) , So(Nf) , Fo , Par
      DOUBLE PRECISION den , a , b , c , gam , pom , MXVDOT
      LOGICAL l1
      DOUBLE PRECISION CON
      PARAMETER (CON=1.0D-8)
      IF ( Idecf/=0 ) THEN
         Iterh = -1
         RETURN
      ENDIF
      l1 = Met1>=2
!
!     DETERMINATION OF THE PARAMETERS B, C
!
      CALL MXDSMM(N,H,Xo,S)
      CALL MXVDIF(N,Go,S,So)
      IF ( Met3==2 ) CALL MXVSCL(N,1.0D0/SQRT(Fo),So,So)
      b = MXVDOT(N,Xo,So)
      IF ( b<=0.0D0 ) l1 = .FALSE.
      a = 0.0D0
      CALL MXDSMM(N,Hh,Xo,S)
      c = MXVDOT(N,Xo,S)
      IF ( c<=0.0D0 ) l1 = .FALSE.
      IF ( l1 ) THEN
!
!     DETERMINATION OF THE PARAMETER GAM (SELF SCALING)
!
         gam = c/b
      ELSE
         gam = 1.0D0
      ENDIF
      Par = gam
!
!     RANK ONE UPDATE
!
      den = Par*b - c
      IF ( ABS(den)<=CON*MAX(CON,ABS(Par*b),ABS(c)) ) THEN
         IF ( b>0.0D0 .AND. c>0.0D0 ) THEN
!
!     BFGS UPDATE
!
            pom = 0.0D0
            CALL MXDSMU(N,Hh,Par/b,So)
            IF ( c>0.0D0 ) CALL MXDSMU(N,Hh,-1.0D0/c,S)
            GOTO 100
         ELSE
            Iterh = 4
            RETURN
         ENDIF
      ENDIF
      pom = Par*b/den
      CALL MXVDIR(N,-Par,So,S,S)
      CALL MXDSMU(N,Hh,1.0D0/den,S)
 100  Iterh = 0
      IF ( gam/=1.0D0 ) CALL MXDSMS(N,Hh,1.0D0/gam)
      END

! SUBROUTINE PUDBQ1                ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! BROYDEN GOOD UPDATE OF A RECTANGULAR MATRIX AFTER THE QR
! DECOMPOSITION.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  II  NA  NUMBER OF EQUATIONS.
!  RU  H(N*(N+1)/2)  UPDATED UPPER TRIANGULAR MATRIX.
!  RI  ETA2  PARAMETER WHICH CONTROLS A NONSINGULARITY
!  RU  AG(N*NA)  UPDATED RECTANGULAR MATRIX.
!  RA  S(N)  AUXILIARY VECTOR.
!  RI  XO(N)  VECTOR OF VARIABLES DIFFERENCE.
!  RI  AFO(NA)  RIGHT HAND SIDES DIFFERENCE.
!  II  MET  VARIABLE METRIC UPDATE.
!  IO  ITERH  UPDATE INDICATOR. ITERH=0-SUCCESSFUL UPDATE.
!         ITERH>0-UNSUCCESSFUL UPDATE.
!  IU  IDECA  DECOMPOSITION INDICATOR.
!  II  NDECA  NUMBER OF DECOMPOSITIONS.
!
! METHOD :
! VARIOUS VARIABLE METRIC UPDATES INCLUDING BFGS UPDATE.
!
      SUBROUTINE PUDBQ1(N,Na,H,Eta2,Ag,S,Xo,Afo,Met,Iterh,Ideca,Ndeca)
      IMPLICIT NONE
      INTEGER N , Na , Met , inf , Iterh , Ideca , Ndeca
      DOUBLE PRECISION H(*) , Eta2 , Ag(*) , S(*) , Xo(*) , Afo(*)
      DOUBLE PRECISION den , MXVDOT
      IF ( Met<=0 ) RETURN
      IF ( Ideca==0 ) THEN
!
!     QR DECOMPOSITION
!
         den = Eta2
         CALL MXDRQF(N,Na,Ag,H)
         CALL MXDPRC(N,H,inf,den)
         Ndeca = Ndeca + 1
         Ideca = 1
      ELSEIF ( Ideca/=1 ) THEN
         Iterh = -1
         GOTO 99999
      ENDIF
!
!     THE GOOD BROYDEN UPDATE
!
      den = MXVDOT(N,Xo,Xo)
      IF ( den<=0.0D0 ) THEN
         Iterh = 2
         GOTO 99999
      ENDIF
      CALL MXVCOP(N,Xo,S)
      CALL MXVNEG(N,Xo,Xo)
      CALL MXDPRM(N,H,Xo,1)
      CALL MXDRMD(N,Na,Ag,Xo,1.0D0,Afo,Afo)
      CALL MXDRQU(N,Na,Ag,H,1.0D0/den,Afo,S,Xo,inf)
      Iterh = 0
99999 END

! SUBROUTINE PUDFM1                ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! VARIABLE METRIC UPDATE OF A DENSE SYMMETRIC POSITIVE DEFINITE MATRIX.
!
! PARAMETERS :
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  RU  B(M)  POSITIVE DEFINITE APPROXIMATION OF THE HESSIAN MATRIX.
!  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RA  S(NF)  AUXILIARY VECTOR.
!  RU  XO(NF)  VECTORS OF VARIABLES DIFFERENCE.
!  RI  GO(NF)  GRADIENTS DIFFERENCE.
!  RI  F  CURRENT VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  PREVIOUS VALUE OF THE OBJECTIVE FUNCTION.
!  RI  ETA5  TOLERANCE FOR A HYBRID METHOD.
!  RI  ETA9  MAXIMUM FOR REAL NUMBERS.
!  II  IPOM1  METHOD INDICATOR.
!  IO  IPOM2  INDICATOR FOR SCALING.
!  II  MET  METHOD SELECTION. MET=0-NO UPDATE. MET=1-BFGS UPDATE.
!  II  MET1  SELECTION OF SELF SCALING.  MET1=1-SELF SCALING SUPPRESSED.
!         MET1=2 SELF SCALING IN THE FIRST ITERATION AFTER RESTART.
!         MET1=3-SELF SCALING IN EACH ITERATION.
!
! COMMON DATA :
!  IO  ITERH  TERMINATION INDICATOR. ITERH<0-BAD DECOMPOSITION.
!         ITERH=0-SUCCESSFUL UPDATE. ITERH>0-NONPOSITIVE PARAMETERS.
!  II  IDECF  DECOMPOSITION INDICATOR. IDECF=0-NO DECOMPOSITION.
!  TO  TUXX  TEXT INFORMATION ON THE CORRECTION USED.
!
! SUBPROGRAMS USED :
!  S   MXDSMM  MATRIX-VECTOR PRODUCT.
!  S   MXDSMU  CORRECTION OF A DENSE SYMMETRIC MATRIX.
!  S   MXDSMS  SCALING OF A DENSE SYMMETRIC MATRIX.
!  RF  MXVDOT  DOT PRODUCT OF VECTORS.
!  S   UOERR1  ERROR MESAGES.
!  S   UOU1D1  PRINT OF ENTRY TO VARIABLE METRIC UPDATE.
!  S   UOU1D2  PRINT OF EXIT FROM VARIABLE METRIC UPDATE.
!
! METHOD :
!  FLETCHER'S COMBINATION OF THE GAUSS-NEWTON AND THE BFGS METHODS.
!
      SUBROUTINE PUDFM1(N,B,S,Xo,Go,F,Fo,Eta5,Ipom1,Ipom2,Met1,Idecf,   &
                        Iterh)
      IMPLICIT NONE
      INTEGER N , Ipom1 , Ipom2 , Met1 , Idecf , Iterh
      DOUBLE PRECISION B(N*(N+1)/2) , S(N) , Xo(N) , Go(N) , F , Fo ,   &
                       Eta5
      DOUBLE PRECISION MXVDOT
      DOUBLE PRECISION ab , bb , cb , gam , par
      LOGICAL l1
      IF ( Idecf/=0 ) THEN
         Iterh = -1
         RETURN
      ENDIF
      par = 0.0D0
!
!     DETERMINATION OF THE PARAMETERS A,B,C
!
      bb = MXVDOT(N,Xo,Go)
      IF ( bb<=0.0D0 ) THEN
         Iterh = 2
         Ipom1 = 0
         RETURN
      ENDIF
      ab = 0.0D0
      CALL MXDSMM(N,B,Xo,S)
      cb = MXVDOT(N,Xo,S)
      IF ( cb<=0.0D0 ) THEN
         Iterh = 3
         RETURN
      ENDIF
      l1 = Met1==4 .OR. Met1==3 .AND. Ipom2>=1 .OR. Met1==2 .AND.       &
           Ipom2==1
      IF ( Fo-F>=Eta5*Fo ) THEN
         Ipom1 = 0
      ELSE
         Ipom1 = 1
      ENDIF
      IF ( l1 ) THEN
!
!     DETERMINATION OF THE PARAMETER GAM (SELF SCALING)
!
         gam = cb/bb
      ELSE
         gam = 1.0D0
      ENDIF
!
!     BFGS UPDATE
!
      CALL MXDSMU(N,B,gam/bb,Go)
      CALL MXDSMU(N,B,-1.0D0/cb,S)
      Iterh = 0
      Ipom2 = 0
      IF ( l1 ) CALL MXDSMS(N,B,1.0D0/gam)
      END

! SUBROUTINE PUDRV1                ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DRIVER FOR HYBRID QUASI-NEWTON UPDATES.
!
! PARAMETERS:
!  RI  R  VALUE OF THE STEPSIZE PARAMETER.
!  RI  FO  PREVIOUS VALUE OF THE OBJECTIVE FUNCTION.
!  RI  F  CURRENT VALUE OF THE OBJECTIVE FUNCTION.
!  RI  PO  PREVIOUS VALUE OF THE DIRECTIONAL DERIVATIVE.
!  II  IPOM1  UPDATE SELECTION.
!  II  IPOM2  METHOD SELECTION.
!  IO  NRED  ACTUAL NUMBER OF EXTRAPOLATIONS OR INTERPOLATIONS.
!  II  IREST  RESTART SPECIFICATION. IF IREST=0 DOES NOT HOLD THEN A
!         RESTART IS PERFORMED.
!
      SUBROUTINE PUDRV1(R,Fo,F,Po,Ipom1,Ipom2,Nred,Irest)
      IMPLICIT NONE
      INTEGER Ipom1 , Ipom2 , Nred , Irest
      DOUBLE PRECISION R , Fo , F , Po
      DOUBLE PRECISION pom
      DOUBLE PRECISION CON2
      PARAMETER (CON2=1.0D-2)
      pom = (Fo-F)/Fo
      SELECT CASE (Ipom2)
      CASE (2)
         Irest = 1
         IF ( pom>=CON2 ) THEN
            Ipom1 = 0
         ELSEIF ( F-Fo<=R*Po ) THEN
            Ipom1 = 0
         ELSE
            Ipom1 = 1
            Irest = 0
         ENDIF
      CASE (3)
         Irest = 1
         IF ( Nred<=0 ) THEN
            IF ( Ipom1/=1 ) THEN
               Ipom1 = 2
               Irest = 0
            ELSE
               Ipom1 = 0
            ENDIF
         ELSEIF ( pom>=CON2 ) THEN
            Ipom1 = 0
         ELSEIF ( Ipom1/=2 ) THEN
            Ipom1 = 1
            Irest = 0
         ELSE
            Ipom1 = 0
         ENDIF
      CASE (4)
         Irest = 1
         Ipom1 = 0
      CASE DEFAULT
         Irest = 1
         IF ( Nred<=0 ) THEN
            Ipom1 = 2
            Irest = 0
         ELSE
            Ipom1 = 0
         ENDIF
      END SELECT
      END

! SUBROUTINE PUDSD2                ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! INITIATION OF A DENSE SYMMETRIC POSITIVE DEFINITE MATRIX
!
! PARAMETERS :
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  RU  H(N*(N+1)/2)  POSITIVE DEFINITE APPROXIMATION OF THE HESSIAN
!         MATRIX
!  RU  B(N*(N+1)/2)  POSITIVE DEFINITE APPROXIMATION OF THE HESSIAN
!         MATRIX
!  RI  F  CURRENT VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  PREVIOUS VALUE OF THE OBJECTIVE FUNCTION.
!  RI  ETA5  TOLERANCE FOR A HYBRID METHOD.
!  II  MET3  TYPE OF STRUCTURED UPDATE. MET3=1-STANDARD STRUCTURED
!         UPDATE. MET3=2-TOTALLY STRUCTURED UPDATE.
!
! COMMON DATA :
!  RU  RAN  RANDOM NUMBER.
!  II  IDECF  DECOMPOSITION INDICATOR. IDECF=0-NO DECOMPOSITION.
!  II  IREST  RESTART SPECIFICATION. IF IREST=0 DOES NOT HOLD THEN A
!         RESTART IS PERFORMED.
!  II  IDIR INDICATOR OF DIRECTION DETERMINATION. IDIR=0-BASIC
!         DETERMINATION. IDIR=1-DETERMINATION AFTER STEPSIZE
!         REDUCTION. IDIR=2-DETERMINATION AFTER STEPSIZE EXPANSION.
!  IU  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!  TO  TUXX   TEXT INFORMATION ON THE RESTART USED.
!
! SUBPROGRAMS USED :
!  S   MXDSMI  GENERATION OF THE UNIT MATRIX.
!  S   MXDSDO  INITIATION OF A DIAGONAL MATRIX.
!  S   MXDSMA  DENSE SYMMETRIC MATRIX AUGMENTED
!              BY THE SCALED DENSE SYMMETRIC MATRIX.
!  S   MXDSMO  A SCALAR IS SET TO ALL ELEMENTS
!              OF A DENSE SYMMETRIC MATRIX.
!  S   MXDSMS  SCALING OF A DENSE SYMMETRIC MATRIX.
!  S   UYSET3  DEFINITION OF THE RESTART VARIABLES.
!
      SUBROUTINE PUDSD2(N,H,B,F,Fo,Eta5,Met3,Ld,Idir,Idecf,Irest,Ind)
      IMPLICIT NONE
      INTEGER N , Met3 , Ld , Idir , Idecf , Irest , Ind
      DOUBLE PRECISION H(N*(N+1)/2) , B(N*(N+1)/2) , F , Fo , Eta5
      INTEGER iudsd
      SAVE iudsd
      Ind = 0
      IF ( Irest<0 ) THEN
         CALL MXDSMI(N,B)
         IF ( F<1.0D0 ) CALL MXDSMS(N,B,SQRT(F))
         iudsd = 1
      ELSEIF ( Irest==0 ) THEN
         IF ( Idir<=0 ) THEN
            IF ( Fo-F<=Eta5*Fo ) THEN
               IF ( Met3==2 ) THEN
                  CALL MXDSMA(N,SQRT(F),B,H,H)
               ELSE
                  CALL MXDSMA(N,1.0D0,B,H,H)
               ENDIF
               Ld = MIN(Ld,1)
               iudsd = 0
            ENDIF
         ENDIF
      ELSEIF ( iudsd==0 ) THEN
         IF ( Met3==2 ) THEN
            CALL MXDSMA(N,-SQRT(F),B,H,H)
         ELSE
            CALL MXDSMA(N,-1.0D0,B,H,H)
         ENDIF
         CALL MXDSMI(N,B)
         IF ( F<1.0D0 ) CALL MXDSMS(N,B,SQRT(F))
         iudsd = 1
      ELSE
         CALL MXDSMI(N,H)
         Ld = MIN(Ld,1)
         iudsd = 1
         Ind = 1
      ENDIF
      Idecf = 0
      END

! SUBROUTINE PUDSD3                ALL SYSTEMS                97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! INITIATION OF A DENSE SYMMETRIC POSITIVE DEFINITE MATRIX
!
! PARAMETERS :
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  RU  H(N*(N+1)/2)  FACTORIZATION H=L*D*TRANS(L) OF A POSITIVE
!         SEMIDEFINITE HESSIAN MATRIX.
!  RU  B(N*(N+1)/2)  FACTORIZATION B=L*D*TRANS(L) OF A POSITIVE
!         DEFINITE APPROXIMATION OF THE HESSIAN MATRIX.
!  IU  IPOM1  METHOD INDICATOR.
!  IU  IPOM2  INDICATOR FOR SCALING.
!
! COMMON DATA :
!  RU  RAN  RANDOM NUMBER.
!  II  IDECF  DECOMPOSITION INDICATOR. IDECF=0-NO DECOMPOSITION.
!  II  IREST  RESTART SPECIFICATION. IF IREST=0 DOES NOT HOLD THEN A
!         RESTART IS PERFORMED.
!  II  ITERS  TERMINATION INDICATOR. ITERS=0-ZERO STEP.
!  IU  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!
! SUBPROGRAMS USED :
!  S   MXDSMI  GENERATION OF THE UNIT MATRIX.
!  S   MXDSDO  INITIATION OF A DIAGONAL MATRIX.
!  S   UUDSMC  COPYING OF A DENSE SYMMETRIC MATRIX.
!  RF  UNRAN1  RANDOM NUMBER GENERATOR.
!  S   UYSET3  DEFINITION OF THE RESTART VARIABLES.
!
      SUBROUTINE PUDSD3(N,H,B,Ipom1,Ipom2,Ld,Idecf,Iters,Irest,Ind)
      IMPLICIT NONE
      INTEGER N , Ipom1 , Ipom2 , Ld , Idecf , Iters , Irest , Ind
      DOUBLE PRECISION H(N*(N+1)/2) , B(N*(N+1)/2)
      INTEGER kdecf
      SAVE kdecf
      Ind = 0
      IF ( .NOT.(Irest==0 .OR. Irest>0 .AND. Ipom1==0) ) THEN
         CALL MXDSMI(N,B)
         Ipom2 = 1
         kdecf = -1
      ENDIF
      IF ( Ipom1==1 ) THEN
         IF ( Iters>0 .OR. Irest>0 ) THEN
            CALL MXDSMC(N,B,H)
            Ld = MIN(Ld,1)
            IF ( Ipom2==1 ) Idecf = kdecf
            IF ( Irest>0 ) Ind = 1
         ENDIF
      ELSEIF ( Irest>0 ) THEN
         Ipom1 = 1
         CALL MXDSMC(N,B,H)
         Ld = MIN(Ld,1)
         IF ( Ipom2==1 ) Idecf = kdecf
      ENDIF
      END

! SUBROUTINE PYADB4             ALL SYSTEMS                   98/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! NEW LINEAR CONSTRAINTS OR NEW SIMPLE BOUNDS ARE ADDED TO THE ACTIVE
! SET. GILL-MURRAY FACTORIZATION OF THE TRANSFORMED HESSIAN MATRIX
! APPROXIMATION IS UPDATED.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  IU  N  ACTUAL NUMBER OF VARIABLES.
!  II  NC  NUMBER OF LINEARIZED CONSTRAINTS.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  IU  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
!  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
!  RI  CF(NC)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
!  RI  CFD(NC) VECTOR CONTAINING INCREMENTS OF THE CONSTRAINT FUNCTIONS.
!  IU  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  IU  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RU  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RU  CZ(NF*NF)  MATRIX WHOSE COLUMNS ARE BASIC VECTORS FROM THE
!         CURRENT REDUCED SUBSPACE.
!  RU  H(NF*(NF+1)/2)  GILL-MURRAY FACTORIZATION OF THE TRANSFORMED
!         HESSIAN MATRIX APPROXIMATION.
!  RA  S(NF)  AUXILIARY VECTOR.
!  RI  R  VALUE OF THE STEPSIZE PARAMETER.
!  RI  EPS7  TOLERANCE FOR LINEAR INDEPENDENCE OF CONSTRAINTS.
!  RI  EPS9  TOLERANCE FOR ACTIVE CONSTRAINTS.
!  RO  GMAX  MAXIMUM ABSOLUTE VALUE OF A PARTIAL DERIVATIVE.
!  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF A NEGATIVE LAGRANGE MULTIPLIER.
!  II  KBF  TYPE OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS. KBF=1-ONE
!         SIDED SIMPLE BOUNDS. KBF=2-TWO SIDED SIMPLE BOUNDS.
!  II  KBC  TYPE OF CONSTRAINTS. KBC=0-NO CONSTRAINTS. KBC=1-CONSTRAINTS
!         WITH ONE SIDED BOUNDS. KBC=2-CONSTRAINTS WITH TWO SIDED
!         BOUNDS.
!  IU  INEW  INDEX OF THE NEW ACTIVE CONSTRAINT.
!  IO  IER  ERROR INDICATOR.
!  IO  ITERM  TERMINATION INDICATOR.
!
! COMMON DATA :
!  IU  NADD  NUMBER OF CONSTRAINT ADDITIONS.
!
! SUBPROGRAMS USED :
!  S   PLADB4  ADDITION OF A NEW ACTIVE CONSTRAINT.
!  S   PLNEWS  IDENTIFICATION OF ACTIVE UPPER BOUNDS.
!  S   PLNEWL  IDENTIFICATION OF ACTIVE LINEAR CONSTRAINRS.
!  S   PLDIRL  NEW VALUES OF CONSTRAINT FUNCTIONS.
!  S   MXVIND  CHANGE OF THE INTEGER VECTOR FOR CONSTRAINT ADDITION.
!
      SUBROUTINE PYADB4(Nf,N,Nc,X,Ix,Xl,Xu,Cf,Cfd,Ic,Ica,Cl,Cu,Cg,Cr,Cz,&
                        H,S,R,Eps7,Eps9,Gmax,Umax,Kbf,Kbc,Inew,Ier,     &
                        Iterm)
      IMPLICIT NONE
      INTEGER Nf , N , Nc , Ix(*) , Ic(*) , Ica(*) , Kbf , Kbc , Inew , &
              Ier , Iterm
      DOUBLE PRECISION X(*) , Xl(*) , Xu(*) , Cf(*) , Cfd(*) , Cl(*) ,  &
                       Cu(*) , Cg(*) , Cr(*) , Cz(*) , H(*) , S(*) , R ,&
                       Eps7 , Eps9 , Gmax , Umax
      INTEGER i , j , k , l , ij , ik , kc , kj , kk , ll
      DOUBLE PRECISION den , temp
      INTEGER NREs , NDEc , NREm , NADd , NIT , NFV , NFG , NFH
      COMMON /STAT  / NREs , NDEc , NREm , NADd , NIT , NFV , NFG , NFH
      IF ( Kbc>0 ) THEN
         IF ( R/=0.0D0 ) CALL PLDIRL(Nc,Cf,Cfd,Ic,R,Kbc)
         IF ( Inew/=0 ) THEN
            IF ( Kbf>0 ) THEN
               DO i = 1 , Nf
                  Inew = 0
                  CALL PLNEWS(X,Ix,Xl,Xu,Eps9,i,Inew)
                  CALL PLADB4(Nf,N,Ica,Cg,Cr,Cz,H,S,Eps7,Gmax,Umax,9,   &
                              Inew,NADd,Ier)
                  CALL MXVIND(Ix,i,Ier)
                  IF ( Ier<0 ) THEN
                     Iterm = -15
                     RETURN
                  ENDIF
               ENDDO
            ENDIF
            DO kc = 1 , Nc
               Inew = 0
               CALL PLNEWL(kc,Cf,Ic,Cl,Cu,Eps9,Inew)
               CALL PLADB4(Nf,N,Ica,Cg,Cr,Cz,H,S,Eps7,Gmax,Umax,9,Inew, &
                           NADd,Ier)
               CALL MXVIND(Ic,kc,Ier)
               IF ( Ier<0 ) THEN
                  Iterm = -15
                  RETURN
               ENDIF
            ENDDO
         ENDIF
      ELSEIF ( Kbf>0 ) THEN
         k = 0
         DO l = 1 , Nf
            IF ( Ix(l)>=0 ) k = k + 1
            Inew = 0
            CALL PLNEWS(X,Ix,Xl,Xu,Eps9,l,Inew)
            IF ( Inew/=0 ) THEN
               Ix(l) = 10 - Ix(l)
               kk = k*(k-1)/2
               den = H(kk+k)
               IF ( den/=0.0D0 ) THEN
                  ij = 0
                  kj = kk
                  DO j = 1 , N
                     IF ( j<=k ) THEN
                        kj = kj + 1
                     ELSE
                        kj = kj + j - 1
                     ENDIF
                     IF ( j/=k ) temp = H(kj)/den
                     ik = kk
                     DO i = 1 , j
                        IF ( i<=k ) THEN
                           ik = ik + 1
                        ELSE
                           ik = ik + i - 1
                        ENDIF
                        ij = ij + 1
                        IF ( i/=k .AND. j/=k ) H(ij) = H(ij)            &
                             + temp*H(ik)
                     ENDDO
                  ENDDO
               ENDIF
               ll = kk + k
               DO i = k + 1 , N
                  DO j = 1 , i
                     ll = ll + 1
                     IF ( j/=k ) THEN
                        kk = kk + 1
                        H(kk) = H(ll)
                     ENDIF
                  ENDDO
               ENDDO
               N = N - 1
            ENDIF
         ENDDO
      ENDIF
      END

! SUBROUTINE PYFUT1                ALL SYSTEMS                98/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! TERMINATION CRITERIA AND TEST ON RESTART.
!
! PARAMETERS :
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  RI  F  NEW VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  OLD VALUE OF THE OBJECTIVE FUNCTION.
!  RI  UMAX  MAXIMUN ABSOLUTE VALUE OF THE NEGATIVE LAGRANGE MULTIPLIER.
!  RO  GMAX  NORM OF THE TRANSFORMED GRADIENT.
!  RI  DMAX  MAXIMUM RELATIVE DIFFERENCE OF VARIABLES.
!  RI  TOLX  LOWER BOUND FOR STEPLENGTH.
!  RI  TOLF  LOWER BOUND FOR FUNCTION DECREASE.
!  RI  TOLB  LOWER BOUND FOR FUNCTION VALUE.
!  RI  TOLG  LOWER BOUND FOR GRADIENT.
!  II  KD  DEGREE OF REQUIRED DERIVATIVES.
!  IU  NIT  ACTUAL NUMBER OF ITERATIONS.
!  II  KIT  NUMBER OF THE ITERATION AFTER RESTART.
!  II  MIT  MAXIMUM NUMBER OF ITERATIONS.
!  IU  NFV  ACTUAL NUMBER OF COMPUTED FUNCTION VALUES.
!  II  MFV  MAXIMUM NUMBER OF COMPUTED FUNCTION VALUES.
!  IU  NFG  ACTUAL NUMBER OF COMPUTED GRADIENT VALUES.
!  II  MFG  MAXIMUM NUMBER OF COMPUTED GRADIENT VALUES.
!  IU  NTESX  ACTUAL NUMBER OF TESTS ON STEPLENGTH.
!  II  MTESX  MAXIMUM NUMBER OF TESTS ON STEPLENGTH.
!  IU  NTESF  ACTUAL NUMBER OF TESTS ON FUNCTION DECREASE.
!  II  MTESF  MAXIMUM NUMBER OF TESTS ON FUNCTION DECREASE.
!  II  ITES  SYSTEM VARIBLE WHICH SPECIFIES TERMINATION. IF ITES=0
!         THEN TERMINATION IS SUPPRESSED.
!  II  IRES1  RESTART SPECIFICATION. RESTART IS PERFORMED AFTER
!         IRES1*N+IRES2 ITERATIONS.
!  II  IRES2  RESTART SPECIFICATION. RESTART IS PERFORMED AFTER
!         IRES1*N+IRES2 ITERATIONS.
!  IU  IREST  RESTART INDICATOR. RESTART IS PERFORMED IF IREST>0.
!  II  ITERS  TERMINATION INDICATOR FOR STEPLENGTH DETERMINATION.
!         ITERS=0 FOR ZERO STEP.
!  IO  ITERM  TERMINATION INDICATOR. ITERM=1-TERMINATION AFTER MTESX
!         UNSUFFICIENT STEPLENGTHS. ITERM=2-TERMINATION AFTER MTESF
!         UNSUFFICIENT FUNCTION DECREASES. ITERM=3-TERMINATION ON LOWER
!         BOUND FOR FUNCTION VALUE. ITERM=4-TERMINATION ON LOWER BOUND
!         FOR GRADIENT. ITERM=11-TERMINATION AFTER MAXIMUM NUMBER OF
!         ITERATIONS. ITERM=12-TERMINATION AFTER MAXIMUM NUMBER OF
!         COMPUTED FUNCTION VALUES.
!
      SUBROUTINE PYFUT1(N,F,Fo,Umax,Gmax,Dmax,Tolx,Tolf,Tolb,Tolg,Kd,   &
                        Nit,Kit,Mit,Nfv,Mfv,Nfg,Mfg,Ntesx,Mtesx,Ntesf,  &
                        Mtesf,Ites,Ires1,Ires2,Irest,Iters,Iterm)
      IMPLICIT NONE
      INTEGER N , Kd , Nit , Kit , Mit , Nfv , Mfv , Nfg , Mfg , Ntesx ,&
              Mtesx , Ntesf , Mtesf , Ites , Ires1 , Ires2 , Irest ,    &
              Iters , Iterm
      DOUBLE PRECISION F , Fo , Umax , Gmax , Dmax , Tolx , Tolf ,      &
                       Tolg , Tolb
      DOUBLE PRECISION temp
      IF ( Iterm<0 ) RETURN
      IF ( Ites>0 ) THEN
         IF ( Iters/=0 ) THEN
            IF ( Nit<=0 ) Fo = F + MIN(SQRT(ABS(F)),ABS(F)/1.0D1)
            IF ( F<=Tolb ) THEN
               Iterm = 3
               RETURN
            ENDIF
            IF ( Kd>0 ) THEN
               IF ( Gmax<=Tolg .AND. Umax<=Tolg ) THEN
                  Iterm = 4
                  RETURN
               ENDIF
            ENDIF
            IF ( Nit<=0 ) THEN
               Ntesx = 0
               Ntesf = 0
            ENDIF
            IF ( Dmax<=Tolx ) THEN
               Iterm = 1
               Ntesx = Ntesx + 1
               IF ( Ntesx>=Mtesx ) RETURN
            ELSE
               Ntesx = 0
            ENDIF
            temp = ABS(Fo-F)/MAX(ABS(F),1.0D0)
            IF ( temp<=Tolf ) THEN
               Iterm = 2
               Ntesf = Ntesf + 1
               IF ( Ntesf>=Mtesf ) RETURN
            ELSE
               Ntesf = 0
            ENDIF
         ENDIF
         IF ( Nit>=Mit ) THEN
            Iterm = 11
            RETURN
         ENDIF
         IF ( Nfv>=Mfv ) THEN
            Iterm = 12
            RETURN
         ENDIF
         IF ( Nfg>=Mfg ) THEN
            Iterm = 13
            RETURN
         ENDIF
      ENDIF
      Iterm = 0
      IF ( N>0 .AND. Nit-Kit>=Ires1*N+Ires2 ) Irest = MAX(Irest,1)
      Nit = Nit + 1
      END

! SUBROUTINE PYRMB1               ALL SYSTEMS                98/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! OLD LINEAR CONSTRAINT OR AN OLD SIMPLE BOUND IS REMOVED FROM THE
! ACTIVE SET. TRANSFORMED GRADIENT OF THE OBJECTIVE FUNCTION AND
! TRANSFORMED HESSIAN MATRIX APPROXIMATION ARE UPDATED.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  IU  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  IU  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  IU  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RU  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RU  CZ(NF*NF)  MATRIX WHOSE COLUMNS ARE BASIC VECTORS FROM THE
!         CURRENT REDUCED SUBSPACE.
!  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RU  GN(NF)  TRANSFORMED GRADIENT OF THE OBJECTIVE FUNCTION.
!  RU  H(NF*(NF+1)/2)  TRANSFORMED HESSIAN MATRIX APPROXIMATION.
!  RI  EPS8  TOLERANCE FOR CONSTRAINT TO BE REMOVED.
!  RI  UMAX  MAXIMUN ABSOLUTE VALUE OF THE NEGATIVE LAGRANGE MULTIPLIER.
!  RI  GMAX  NORM OF THE TRANSFORMED GRADIENT.
!  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
!         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
!  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
!         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
!         SIDED LINEAR CONSTRAINTS.
!  II  IOLD  INDEX OF THE REMOVED CONSTRAINT.
!  IA  KOLD  AUXILIARY VARIABLE.
!  IA  KREM  AUXILIARY VARIABLE.
!  IO  IER  ERROR INDICATOR.
!  IO  ITERM  TERMINATION INDICATOR.
!
! COMMON DATA :
!  IU  NREM  NUMBER OF CONSTRAINT DELETIONS.
!
! SUBPROGRAMS USED :
!  S   PLRMB0  CONSTRAINT DELETION.
!  S   MXVSET  INITIATION OF A VECTOR.
!
      SUBROUTINE PYRMB1(Nf,N,Ix,Ic,Ica,Cg,Cr,Cz,G,Gn,H,Eps8,Umax,Gmax,  &
                        Kbf,Kbc,Iold,Kold,Krem,Ier,Iterm)
      IMPLICIT NONE
      INTEGER Nf , N , Ix(*) , Ic(*) , Ica(*) , Kbf , Kbc , Iold ,      &
              Kold , Krem , Ier , Iterm
      DOUBLE PRECISION Cg(*) , Cr(*) , Cz(*) , G(*) , Gn(*) , H(*) ,    &
                       Eps8 , Umax , Gmax
      INTEGER i , j , k , kc , l
      INTEGER NREs , NDEc , NREm , NADd , NIT , NFV , NFG , NFH
      COMMON /STAT  / NREs , NDEc , NREm , NADd , NIT , NFV , NFG , NFH
      IF ( Kbc>0 ) THEN
         IF ( Umax>Eps8*Gmax ) THEN
            CALL PLRMB0(Nf,N,Ica,Cg,Cr,Cz,G,Gn,Iold,Krem,NREm,Ier)
            IF ( Ier<0 ) THEN
               Iterm = -16
            ELSEIF ( Ier>0 ) THEN
               Iold = 0
            ELSE
               k = N*(N-1)/2
               CALL MXVSET(N,0.0D0,H(k+1))
               H(k+N) = 1.0D0
               kc = Ica(Nf-N+1)
               IF ( kc>0 ) THEN
                  Ic(kc) = -Ic(kc)
               ELSE
                  k = -kc
                  Ix(k) = -Ix(k)
               ENDIF
            ENDIF
         ELSE
            Iold = 0
         ENDIF
      ELSEIF ( Kbf>0 ) THEN
         IF ( Umax>Eps8*Gmax ) THEN
            Ix(Iold) = MIN(ABS(Ix(Iold)),3)
            DO i = N , Kold , -1
               Gn(i+1) = Gn(i)
            ENDDO
            Gn(Kold) = G(Iold)
            N = N + 1
            k = N*(N-1)/2
            l = k + N
            DO i = N , Kold , -1
               DO j = i , 1 , -1
                  IF ( i/=Kold .AND. j/=Kold ) THEN
                     H(l) = H(k)
                     k = k - 1
                     l = l - 1
                  ELSEIF ( i==Kold .AND. j==Kold ) THEN
                     H(l) = 1.0D0
                     l = l - 1
                  ELSE
                     H(l) = 0.0D0
                     l = l - 1
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            Iold = 0
            Kold = 0
         ENDIF
      ENDIF
      END

! SUBROUTINE PYTRBD             ALL SYSTEMS                   98/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! VECTORS OF VARIABLES DIFFERENCE AND GRADIENTS DIFFERENCE ARE COMPUTED
! AND TRANSFORMED. TEST VALUE DMAX IS DETERMINED.
!
! PARAMETERS :
!  II  NF DECLARED NUMBER OF VARIABLES.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  RU  XO(NF)  VECTORS OF VARIABLES DIFFERENCE.
!  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RU  GO(NF)  GRADIENTS DIFFERENCE.
!  RI  CZ(NF*NF)  MATRIX WHOSE COLUMNS ARE BASIC VECTORS FROM CURRENT
!         REDUCED SUBSPACE.
!  RU  SN(NF)  TRANSFORMED DIRECTION VECTOR.
!  RI  R  VALUE OF THE STEPSIZE PARAMETER.
!  RU  F  NEW VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  OLD VALUE OF THE OBJECTIVE FUNCTION.
!  RU  P  NEW VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RU  PO  OLD VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RO  DMAX  MAXIMUM RELATIVE DIFFERENCE OF VARIABLES.
!  II  ITERS  TERMINATION INDICATOR FOR STEPLENGTH DETERMINATION.
!         ITERS=0 FOR ZERO STEP.
!  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
!         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
!  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
!         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
!         SIDED LINEAR CONSTRAINTS.
!
! SUBPROGRAMS USED :
!  S   MXDRMM  PREMULTIPLICATION OF A VECTOR BY TRANSPOSE OF A DENSE
!         RECTANGULAR MATRIX.
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVDIF  DIFFERENCE OF TWO VECTORS.
!  S   MXVMUL  DIAGONAL PREMULTIPLICATION OF A VECTOR.
!  S   MXVSAV  DIFFERENCE OF TWO VECTORS WITH COPYING AND SAVING THE
!         SUBSTRACTED ONE.
!  S   MXVSCL  SCALING OF A VECTOR.
!
      SUBROUTINE PYTRBD(Nf,N,X,Ix,Xo,G,Go,Cz,Sn,R,F,Fo,P,Po,Dmax,Iters, &
                        Kbf,Kbc)
      IMPLICIT NONE
      INTEGER Nf , N , Ix(*) , Iters , Kbf , Kbc
      DOUBLE PRECISION X(*) , Xo(*) , G(*) , Go(*) , Cz(*) , Sn(*) , R ,&
                       F , Fo , P , Po , Dmax
      INTEGER i , k
      IF ( Iters>0 ) THEN
         CALL MXVDIF(Nf,X,Xo,Xo)
         CALL MXVDIF(Nf,G,Go,Go)
         Po = R*Po
         P = R*P
      ELSE
         F = Fo
         P = Po
         CALL MXVSAV(Nf,X,Xo)
         CALL MXVSAV(Nf,G,Go)
      ENDIF
      Dmax = 0.0D0
      IF ( Kbc>0 ) THEN
         DO i = 1 , Nf
            Dmax = MAX(Dmax,ABS(Xo(i))/MAX(ABS(X(i)),1.0D0))
         ENDDO
         IF ( N>0 ) THEN
            CALL MXVSCL(N,R,Sn,Xo)
            CALL MXVCOP(Nf,Go,Sn)
            CALL MXDRMM(Nf,N,Cz,Sn,Go)
         ENDIF
      ELSEIF ( Kbf>0 ) THEN
         k = 0
         DO i = 1 , Nf
            IF ( Ix(i)>=0 ) THEN
               k = k + 1
               Dmax = MAX(Dmax,ABS(Xo(i))/MAX(ABS(X(i)),1.0D0))
               Xo(k) = Xo(i)
               Go(k) = Go(i)
            ENDIF
         ENDDO
      ELSE
         DO i = 1 , Nf
            Dmax = MAX(Dmax,ABS(Xo(i))/MAX(ABS(X(i)),1.0D0))
         ENDDO
      ENDIF
      END

! SUBROUTINE PYTRBG               ALL SYSTEMS                98/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! GRADIENT OF THE OBJECTIVE FUNCTION IS SCALED AND REDUCED.
! TEST VALUES GMAX AND UMAX ARE COMPUTED.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  II  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RI  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RU  CZ(NF*NF)  MATRIX WHOSE COLUMNS ARE BASIC VECTORS FROM THE
!         CURRENT REDUCED SUBSPACE.
!  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RO  GN(NF)  TRANSFORMED GRADIENT OF THE OBJECTIVE FUNCTION.
!  RI  EPS7  TOLERANCE FOR LINEAR INDEPENDENCE OF CONSTRAINTS.
!  RO  UMAX  MAXIMUM ABSOLUTE VALUE OF THE NEGATIVE LAGRANGE MULTIPLIER.
!  RO  GMAX  NORM OF THE TRANSFORMED GRADIENT.
!  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
!         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
!  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
!         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
!         SIDED LINEAR CONSTRAINTS.
!  II  IOLD  INDEX OF THE REMOVED CONSTRAINT.
!  IA  KOLD  AUXILIARY VARIABLE.
!
! SUBPROGRAMS USED :
!  S   MXDRMM  PREMULTIPLICATION OF A VECTOR BY A ROWWISE STORED DENSE
!         RECTANGULAR MATRIX.
!  S   MXDPRB  BACK SUBSTITUTION.
!  S   MXVCOP  COPYING OF A VECTOR.
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!  RF  MXVMAX  L-INFINITY NORM OF A VECTOR.
!  S   MXVMUL  DIAGONAL PREMULTIPLICATION OF A VECTOR.
!
      SUBROUTINE PYTRBG(Nf,N,Ix,Ic,Ica,Cg,Cr,Cz,G,Gn,Umax,Gmax,Kbf,Kbc, &
                        Iold,Kold)
      IMPLICIT NONE
      INTEGER Nf , N , Ix(*) , Ic(*) , Ica(*) , Kbf , Kbc , Iold , Kold
      DOUBLE PRECISION Cg(*) , Cr(*) , Cz(*) , G(*) , Gn(*) , Umax ,    &
                       Gmax
      DOUBLE PRECISION temp , MXVMAX , MXVDOT
      INTEGER nca , ncz , i , j , k , kc
      Iold = 0
      Kold = 0
      Umax = 0.0D0
      Gmax = 0.0D0
      IF ( Kbc>0 ) THEN
         IF ( Nf>N ) THEN
            nca = Nf - N
            ncz = N*Nf
            CALL MXVCOP(Nf,G,Gn)
            DO j = 1 , nca
               k = Ica(j)
               IF ( k>0 ) THEN
                  Cz(ncz+j) = MXVDOT(Nf,Cg((k-1)*Nf+1),Gn)
               ELSE
                  i = -k
                  Cz(ncz+j) = Gn(i)
               ENDIF
            ENDDO
            CALL MXDPRB(nca,Cr,Cz(ncz+1),0)
            DO j = 1 , nca
               temp = Cz(ncz+j)
               kc = Ica(j)
               IF ( kc>0 ) THEN
                  k = Ic(kc)
               ELSE
                  i = -kc
                  k = Ix(i)
               ENDIF
               IF ( k<=-5 ) THEN
               ELSEIF ( (k==-1 .OR. k==-3) .AND. Umax+temp>=0.0D0 ) THEN
               ELSEIF ( .NOT.((k==-2 .OR. k==-4) .AND. Umax-temp>=0.0D0)&
                        ) THEN
                  Iold = j
                  Umax = ABS(temp)
               ENDIF
            ENDDO
         ENDIF
         IF ( N>0 ) THEN
            CALL MXDRMM(Nf,N,Cz,G,Gn)
            Gmax = MXVMAX(N,Gn)
         ENDIF
      ELSEIF ( Kbf>0 ) THEN
         j = 0
         Iold = 0
         Kold = 0
         DO i = 1 , Nf
            temp = G(i)
            k = Ix(i)
            IF ( k>=0 ) THEN
               j = j + 1
               Gn(j) = temp
               Gmax = MAX(Gmax,ABS(temp))
            ELSEIF ( k<=-5 ) THEN
            ELSEIF ( (k==-1 .OR. k==-3) .AND. Umax+temp>=0.0D0 ) THEN
            ELSEIF ( .NOT.((k==-2 .OR. k==-4) .AND. Umax-temp>=0.0D0) ) &
                     THEN
               Iold = i
               Kold = j + 1
               Umax = ABS(temp)
            ENDIF
         ENDDO
         N = j
      ELSE
         DO i = 1 , Nf
            temp = G(i)
            Gmax = MAX(Gmax,ABS(temp))
         ENDDO
         N = Nf
      ENDIF
      END

! SUBROUTINE PYTRBH               ALL SYSTEMS                98/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! HESSIAN MATRIX OF THE OBJECTIVE FUNCTION OR ITS APPROXIMATION IS
! SCALED AND REDUCED.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  RI  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RI  CZ(NF*NF)  MATRIX WHOSE COLUMNS ARE BASIC VECTORS FROM THE
!         CURRENT REDUCED SUBSPACE.
!  RI  H(NF*(NF+1)/2)  HESSIAN MATRIX OR ITS APPROXIMATION.
!  RA  S(NF)  AUXILIARY VECTOR.
!  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
!         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
!  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
!         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
!         SIDED LINEAR CONSTRAINTS.
!  II  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!  II  ITERS  TERMINATION INDICATOR FOR STEPLENGTH DETERMINATION.
!
! SUBPROGRAMS USED :
!  S   MXDSMM  MATRIX VECTOR PRODUCT.
!  S   MXVCOP  COPYING OF A VECTOR.
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!
      SUBROUTINE PYTRBH(Nf,N,Ix,Cr,Cz,H,S,Kbf,Kbc,Ld,Iters)
      IMPLICIT NONE
      INTEGER Nf , N , Ix(*) , Kbf , Kbc , Ld , Iters
      DOUBLE PRECISION Cr(*) , Cz(*) , H(*) , S(*)
      DOUBLE PRECISION MXVDOT
      INTEGER nca , ncr , icz , jcz , i , j , k , l
      IF ( Ld/=2 .OR. Iters==0 ) RETURN
      IF ( Kbc>0 ) THEN
         IF ( N<=0 ) RETURN
         nca = Nf - N
         ncr = nca*(nca+1)/2
         k = ncr
         jcz = 1
         DO j = 1 , N
            CALL MXDSMM(Nf,H,Cz(jcz),S)
            icz = 1
            DO i = 1 , j
               k = k + 1
               Cr(k) = MXVDOT(Nf,Cz(icz),S)
               icz = icz + Nf
            ENDDO
            jcz = jcz + Nf
         ENDDO
         CALL MXVCOP(N*(N+1)/2,Cr(ncr+1),H)
      ELSEIF ( Kbf>0 ) THEN
         k = 0
         l = 0
         DO i = 1 , Nf
            DO j = 1 , i
               k = k + 1
               IF ( Ix(i)>=0 .AND. Ix(j)>=0 ) THEN
                  l = l + 1
                  H(l) = H(k)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      END

! SUBROUTINE PYTRBS               ALL SYSTEMS                98/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! SCALED AND REDUCED DIRECTION VECTOR IS BACK TRANSFORMED.
! VECTORS X,G AND VALUES F,P ARE SAVED.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  IU  N  ACTUAL NUMBER OF VARIABLES.
!  II  NC  NUMBER OF LINEARIZED CONSTRAINTS.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
!  RO  XO(NF)  SAVED VECTOR OF VARIABLES.
!  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
!  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
!  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RO  GO(NF)  SAVED GRADIENT OF THE OBJECTIVE FUNCTION.
!  RI  CF(NF)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCYIONS.
!  RO  CFD(NF)  VECTOR CONTAINING INCREMENTS OF THE CONSTRAINT
!         FUNCTIONS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RI  CZ(NF*NF)  MATRIX WHOSE COLUMNS ARE BASIC VECTORS FROM THE
!         CURRENT REDUCED SUBSPACE.
!  RI  SN(NF)  TRANSFORMED DIRECTION VECTOR.
!  RO  S(NF)  DIRECTION VECTOR.
!  RO  RO  SAVED VALUE OF THE STEPSIZE PARAMETER.
!  RO  FP  PREVIOUS VALUE OF THE OBJECTIVE FUNCTION.
!  RU  FO  SAVED VALUE OF THE OBJECTIVE FUNCTION.
!  RI  F  VALUE OF THE OBJECTIVE FUNCTION.
!  RO  PO  SAVED VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RI  P  VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RU  RMAX  MAXIMUM VALUE OF THE STEPSIZE PARAMETER.
!  II  KBF  SPECIFICATION OF SIMPLE BOUNDS. KBF=0-NO SIMPLE BOUNDS.
!         KBF=1-ONE SIDED SIMPLE BOUNDS. KBF=2=TWO SIDED SIMPLE BOUNDS.
!  II  KBC  SPECIFICATION OF LINEAR CONSTRAINTS. KBC=0-NO LINEAR
!         CONSTRAINTS. KBC=1-ONE SIDED LINEAR CONSTRAINTS. KBC=2=TWO
!         SIDED LINEAR CONSTRAINTS.
!  IO  KREM  INDICATION OF LINEARLY DEPENDENT GRADIENTS.
!  IO  INEW  INDEX OF THE NEW ACTIVE FUNCTION.
!
! SUBPROGRAMS USED :
!  S   PLMAXS  DETERMINATION OF THE MAXIMUM STEPSIZE USING SIMPLE
!         BOUNDS.
!  S   PLMAXL  DETERMINATION OF THE MAXIMUM STEPSIZE USING LINEAR
!         CONSTRAINTS.
!  S   MXDCMM  MATRIX VECTOR PRODUCT.
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVSET  INITIATION OF A VECTOR.
!
      SUBROUTINE PYTRBS(Nf,N,Nc,X,Ix,Xo,Xl,Xu,G,Go,Cf,Cfd,Ic,Cl,Cu,Cg,  &
                        Cz,Sn,S,Ro,Fp,Fo,F,Po,P,Rmax,Kbf,Kbc,Krem,Inew)
      IMPLICIT NONE
      INTEGER Nf , N , Nc , Ix(*) , Ic(*) , Kbf , Kbc , Krem , Inew
      DOUBLE PRECISION X(*) , Xo(*) , Xl(*) , Xu(*) , G(*) , Go(*) ,    &
                       Cf(*) , Cfd(*) , Cl(*) , Cu(*) , Cg(*) , Cz(*) , &
                       Sn(*) , S(*) , Ro , Fp , Fo , F , Po , P , Rmax
      INTEGER i , k
      Fp = Fo
      Ro = 0.0D0
      Fo = F
      Po = P
      CALL MXVCOP(Nf,X,Xo)
      CALL MXVCOP(Nf,G,Go)
      IF ( Kbc>0 ) THEN
         IF ( N>0 ) THEN
            CALL MXDCMM(Nf,N,Cz,Sn,S)
            Inew = 0
            CALL PLMAXL(Nf,Nc,Cf,Cfd,Ic,Cl,Cu,Cg,S,Rmax,Kbc,Krem,Inew)
            CALL PLMAXS(Nf,X,Ix,Xl,Xu,S,Rmax,Kbf,Krem,Inew)
         ELSE
            CALL MXVSET(Nf,0.0D0,S)
         ENDIF
      ELSEIF ( Kbf>0 ) THEN
         k = N + 1
         DO i = Nf , 1 , -1
            IF ( Ix(i)<0 ) THEN
               S(i) = 0.0D0
            ELSE
               k = k - 1
               S(i) = Sn(k)
            ENDIF
         ENDDO
         Inew = 0
         CALL PLMAXS(Nf,X,Ix,Xl,Xu,S,Rmax,Kbf,Krem,Inew)
      ENDIF
      END

! SUBROUTINE PYTRFD             ALL SYSTEMS                   90/12/01
! PORTABILITY : ALL SYSTEMS
! 90/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! PREPARATION OF VARIABLE METRIC UPDATE.
!
! PARAMETERS :
!  II  NF  DECLARED NUMBER OF VARIABLES.
!  II  NC  NUMBER OF CONSTRAINTS.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  RU  XO(NF)  SAVED VECTOR OF VARIABLES.
!  II  IAA(NF+1)  VECTOR CONTAINING INDICES OF ACTIVE FUNCTIONS.
!  RI  AG(NF*NA)  MATRIX WHOSE COLUMNS ARE GRADIENTS OF THE LINEAR
!          APPROXIMATED FUNCTIONS.
!  RI  AZ(NF+1)  VECTOR OF LAGRANGE MULTIPLIERS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RI  G(NF)  GRADIENT OF THE LAGRANGIAN FUNCTION.
!  RU  GO(NF)  SAVED GRADIENT OF THE LAGRANGIAN FUNCTION.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  II  KD  DEGREE OF REQUIRED DERVATIVES.
!  IU  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!  RU  R  VALUE OF THE STEPSIZE PARAMETER.
!  RU  F  VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  SAVED VALUE OF THE OBJECTIVE FUNCTION.
!  RU  P  VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RU  PO  SAVED VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RO  DMAX  RELATIVE STEPSIZE.
!  IO  ITERS  TERMINATION INDICATOR. ITERS=0-ZERO STEP. ITERS=1-PERFECT
!         LINE SEARCH. ITERS=2 GOLDSTEIN STEPSIZE. ITERS=3-CURRY
!         STEPSIZE. ITERS=4-EXTENDED CURRY STEPSIZE.
!         ITERS=5-ARMIJO STEPSIZE. ITERS=6-FIRST STEPSIZE.
!         ITERS=7-MAXIMUM STEPSIZE. ITERS=8-UNBOUNDED FUNCTION.
!         ITERS=-1-MRED REACHED. ITERS=-2-POSITIVE DIRECTIONAL
!         DERIVATIVE. ITERS=-3-ERROR IN INTERPOLATION.
!
! SUBPROGRAMS USED :
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVDIF  DIFFERENCE OF TWO VECTORS.
!  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
!  S   MXVSET  INITIATION OF A VECTOR.
!  S   MXVSAV  DIFFERENCE OF TWO VECTORS WITH COPYING AND SAVING THE
!         SUBSTRACTED ONE.
!
      SUBROUTINE PYTRFD(Nf,Nc,X,Xo,Iaa,Ag,Az,Cg,G,Go,N,Kd,Ld,R,F,Fo,P,  &
                        Po,Dmax,Iters)
      IMPLICIT NONE
      DOUBLE PRECISION Dmax , F , Fo , P , Po , R
      INTEGER Iters , Kd , Ld , N , Nc , Nf
      DOUBLE PRECISION Ag(*) , Az(*) , Cg(*) , G(*) , Go(*) , X(*) ,    &
                       Xo(*)
      INTEGER Iaa(*)
      INTEGER i , j , l
      CALL MXVSET(Nf,0.0D0,G)
      DO j = 1 , Nf - N
         l = Iaa(j)
         IF ( l>Nc ) THEN
            l = l - Nc
            CALL MXVDIR(Nf,-Az(j),Ag((l-1)*Nf+1),G,G)
         ELSEIF ( l>0 ) THEN
            CALL MXVDIR(Nf,-Az(j),Cg((l-1)*Nf+1),G,G)
         ELSE
            l = -l
            G(l) = G(l) - Az(j)
         ENDIF
      ENDDO
      IF ( Iters>0 ) THEN
         CALL MXVDIF(Nf,X,Xo,Xo)
         CALL MXVDIF(Nf,G,Go,Go)
         Po = R*Po
         P = R*P
      ELSE
         R = 0.0D0
         F = Fo
         P = Po
         CALL MXVSAV(Nf,X,Xo)
         CALL MXVSAV(Nf,G,Go)
         Ld = Kd
      ENDIF
      Dmax = 0.0D0
      DO i = 1 , Nf
         Dmax = MAX(Dmax,ABS(Xo(i))/MAX(ABS(X(i)),1.0D0))
      ENDDO
      N = Nf
      END

! SUBROUTINE PYTRND             ALL SYSTEMS                   91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! DUAL RANGE SPACE QUADRATIC PROGRAMMING METHOD FOR MINIMAX
! APPROXIMATION.
!
! PARAMETERS :
!  II  NF DECLARED NUMBER OF VARIABLES.
!  II  N  ACTUAL NUMBER OF VARIABLES.
!  II  NC NUMBER OF CONSTRAINTS.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  RI  XN(NF)  VECTOR OF SCALING FACTORS.
!  RO  XO(NF)  SAVED VECTOR OF VARIABLES.
!  II  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RO  CZ(NF)  VECTOR OF LAGRANGE MULTIPLIERS.
!  RO  CZS(NF)  SAVED VECTOR OF LAGRANGE MULTIPLIERS.
!  RI  G(NF)  GRADIENT OF THE LAGRANGIAN FUNCTION.
!  RI  GO(NF)  SAVED GRADIENT OF THE LAGRANGIAN FUNCTION.
!  RO  R  VALUE OF THE STEPSIZE PARAMETER.
!  RO  F  NEW VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  OLD VALUE OF THE OBJECTIVE FUNCTION.
!  RO  P  NEW VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RI  PO  OLD VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RI  CMAX  VALUE OF THE CONSTRAINT VIOLATION.
!  RO  CMAXO  SAVED VALUE OF THE CONSTRAINT VIOLATION.
!  RO  DMAX  MAXIMUM RELATIVE DIFFERENCE OF VARIABLES.
!
! COMMON DATA :
!  II  NORMF  SCALING SPECIFICATION. NORMF=0-NO SCALING PERFORMED.
!         NORMF=1-SCALING FACTORS ARE DETERMINED AUTOMATICALLY.
!         NORMF=2-SCALING FACTORS ARE SUPPLIED BY USER.
!  II  ITERS  TERMINATION INDICATOR FOR STEPLENGTH DETERMINATION.
!         ITERS=0 FOR ZERO STEP.
!
! SUBPROGRAMS USED :
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
!  S   MXVSET  INITIATION OF A VECTOR.
!  S   MXVSAV  DIFFERENCE OF TWO VECTORS WITH COPYING AND SAVING THE
!         SUBSTRACTED ONE.
!
      SUBROUTINE PYTRND(Nf,N,X,Xo,Ica,Cg,Cz,G,Go,R,F,Fo,P,Po,Cmax,Cmaxo,&
                        Dmax,Kd,Ld,Iters)
      IMPLICIT NONE
      INTEGER Nf , N , Kd , Ld , Iters
      INTEGER Ica(*)
      DOUBLE PRECISION X(*) , Xo(*) , Cg(*) , Cz(*) , G(*) , Go(*) , R ,&
                       F , Fo , P , Po , Cmax , Cmaxo , Dmax
      INTEGER i , j , l
      DO j = 1 , Nf - N
         l = Ica(j)
         IF ( l>0 ) THEN
            CALL MXVDIR(Nf,-Cz(j),Cg((l-1)*Nf+1),G,G)
         ELSE
            l = -l
            G(l) = G(l) - Cz(j)
         ENDIF
      ENDDO
      IF ( Iters>0 ) THEN
         CALL MXVDIF(Nf,X,Xo,Xo)
         CALL MXVDIF(Nf,G,Go,Go)
         Po = R*Po
         P = R*P
      ELSE
         F = Fo
         P = Po
         Cmax = Cmaxo
         CALL MXVSAV(Nf,X,Xo)
         CALL MXVSAV(Nf,G,Go)
         Ld = Kd
      ENDIF
      Dmax = 0.0D0
      DO i = 1 , Nf
         Dmax = MAX(Dmax,ABS(Xo(i))/MAX(ABS(X(i)),1.0D0))
      ENDDO
      N = Nf
      END

! SUBROUTINE PYTRUD             ALL SYSTEMS                   98/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! VECTORS OF VARIABLES DIFFERENCE AND GRADIENTS DIFFERENCE ARE COMPUTED
! AND SCALED. TEST VALUE DMAX IS DETERMINED.
!
! PARAMETERS :
!  II  NF DECLARED NUMBER OF VARIABLES.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  RU  XO(NF)  VECTORS OF VARIABLES DIFFERENCE.
!  RI  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RU  GO(NF)  GRADIENTS DIFFERENCE.
!  RO  R  VALUE OF THE STEPSIZE PARAMETER.
!  RO  F  NEW VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  OLD VALUE OF THE OBJECTIVE FUNCTION.
!  RO  P  NEW VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RI  PO  OLD VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RO  DMAX  MAXIMUM RELATIVE DIFFERENCE OF VARIABLES.
!  II  KD  DEGREE OF REQUIRED DERVATIVES.
!  IO  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!  II  ITERS  TERMINATION INDICATOR FOR STEPLENGTH DETERMINATION.
!         ITERS=0 FOR ZERO STEP.
!
! SUBPROGRAMS USED :
!  S   PYSET1  DEGREE DEFINITION OF THE COMPUTED DERIVATIVES.
!  S   MXVDIF  DIFFERENCE OF TWO VECTORS.
!  S   MXVSAV  DIFFERENCE OF TWO VECTORS WITH COPYING AND SAVING THE
!         SUBSTRACTED ONE.
!
      SUBROUTINE PYTRUD(Nf,X,Xo,G,Go,R,F,Fo,P,Po,Dmax,Kd,Ld,Iters)
      IMPLICIT NONE
      INTEGER Nf , Kd , Ld , Iters
      DOUBLE PRECISION X(*) , Xo(*) , G(*) , Go(*) , R , F , Fo , P ,   &
                       Po , Dmax
      INTEGER i
      IF ( Iters>0 ) THEN
         CALL MXVDIF(Nf,X,Xo,Xo)
         CALL MXVDIF(Nf,G,Go,Go)
         Po = R*Po
         P = R*P
      ELSE
         F = Fo
         P = Po
         CALL MXVSAV(Nf,X,Xo)
         CALL MXVSAV(Nf,G,Go)
         Ld = Kd
      ENDIF
      Dmax = 0.0D0
      DO i = 1 , Nf
         Dmax = MAX(Dmax,ABS(Xo(i))/MAX(ABS(X(i)),1.0D0))
      ENDDO
      END

! SUBROUTINE PYTRUF             ALL SYSTEMS                   98/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! VECTORS OF VARIABLES DIFFERENCE AND RIGHT HAND SIDES DIFFERENCE ARE
! COMPUTED AND SCALED. TEST VALUE DMAX IS DETERMINED.
!
! PARAMETERS :
!  II  NF DECLARED NUMBER OF VARIABLES.
!  II  NA NUMBER OF APPROXIMATED FUNCTIONS.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  RU  XO(NF)  VECTORS OF VARIABLES DIFFERENCE.
!  RI  AF(NA)  VECTOR OF RIGHT HAND SIDES.
!  RI  AFO(NA)  VECTOR OF RIGHT HAND SIDES DIFFERENCE.
!  RO  R  VALUE OF THE STEPSIZE PARAMETER.
!  RO  F  NEW VALUE OF THE OBJECTIVE FUNCTION.
!  RI  FO  OLD VALUE OF THE OBJECTIVE FUNCTION.
!  RO  P  NEW VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RI  PO  OLD VALUE OF THE DIRECTIONAL DERIVATIVE.
!  RO  DMAX  MAXIMUM RELATIVE DIFFERENCE OF VARIABLES.
!  II  KD  DEGREE OF REQUIRED DERVATIVES.
!  IO  LD  DEGREE OF PREVIOUSLY COMPUTED DERIVATIVES.
!  II  ITERS  TERMINATION INDICATOR FOR STEPLENGTH DETERMINATION.
!         ITERS=0 FOR ZERO STEP.
!
! SUBPROGRAMS USED :
!  S   PYSET1  DEGREE DEFINITION OF THE COMPUTED DERIVATIVES.
!  S   MXVDIF  DIFFERENCE OF TWO VECTORS.
!  S   MXVSAV  DIFFERENCE OF TWO VECTORS WITH COPYING AND SAVING THE
!         SUBSTRACTED ONE.
!
      SUBROUTINE PYTRUF(Nf,Na,X,Xo,Af,Afo,R,F,Fo,P,Po,Dmax,Kd,Ld,Iters)
      IMPLICIT NONE
      INTEGER Nf , Na , Kd , Ld , Iters
      DOUBLE PRECISION X(*) , Xo(*) , Af(*) , Afo(*) , R , F , Fo , P , &
                       Po , Dmax
      INTEGER i
      IF ( Iters>0 ) THEN
         CALL MXVDIF(Nf,X,Xo,Xo)
         CALL MXVDIF(Na,Af,Afo,Afo)
         Po = R*Po
         P = R*P
      ELSE
         R = 0.0D0
         F = Fo
         P = Po
         CALL MXVSAV(Nf,X,Xo)
         CALL MXVSAV(Na,Af,Afo)
         Ld = Kd
      ENDIF
      Dmax = 0.0D0
      DO i = 1 , Nf
         Dmax = MAX(Dmax,ABS(Xo(i))/MAX(ABS(X(i)),1.0D0))
      ENDDO
      END
