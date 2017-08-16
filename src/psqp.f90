
!***********************************************************************
! SUBROUTINE PSQPN              ALL SYSTEMS                   97/01/22
! PURPOSE :
! EASY TO USE SUBROUTINE FOR GENERAL NONLINEAR PROGRAMMING PROBLEMS.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  II  NB  CHOICE OF SIMPLE BOUNDS. NB=0-SIMPLE BOUNDS SUPPRESSED.
!         NB>0-SIMPLE BOUNDS ACCEPTED.
!  II  NC  NUMBER OF LINEAR CONSTRAINTS.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS. IX(I)=0-VARIABLE
!         X(I) IS UNBOUNDED. IX(I)=1-LOVER BOUND XL(I).LE.X(I).
!         IX(I)=2-UPPER BOUND X(I).LE.XU(I). IX(I)=3-TWO SIDE BOUND
!         XL(I).LE.X(I).LE.XU(I). IX(I)=5-VARIABLE X(I) IS FIXED.
!  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
!  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
!  RI  CF(NC+1)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!         IC(KC)=0-CONSTRAINT CF(KC) IS NOT USED. IC(KC)=1-LOVER
!         CONSTRAINT CL(KC).LE.CF(KC). IC(KC)=2-UPPER CONSTRAINT
!         CF(KC).LE.CU(KC). IC(KC)=3-TWO SIDE CONSTRAINT
!         CL(KC).LE.CF(KC).LE.CU(KC). IC(KC)=5-EQUALITY CONSTRAINT
!         CF(KC).EQ.CL(KC).
!  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  II  IPAR(6)  INTEGER PAREMETERS:
!      IPAR(1)  MAXIMUM NUMBER OF ITERATIONS.
!      IPAR(2)  MAXIMUM NUMBER OF FUNCTION EVALUATIONS.
!      IPAR(3)  THIS PARAMETER IS NOT USED IN THE SUBROUTINE PSQP.
!      IPAR(4)  THIS PARAMETER IS NOT USED IN THE SUBROUTINE PSQP.
!      IPAR(5)  VARIABLE METRIC UPDATE USED. IPAR(5)=1-THE BFGS UPDATE.
!         IPAR(5)-THE HOSHINO UPDATE.
!      IPAR(6)  CORRECTION OF THE VARIABLE METRIC UPDATE IF A NEGATIVE
!         CURVATURE OCCURS. IPAR(6)=1-NO CORRECTION. IPAR(6)=2-POWELL'S
!         CORRECTION.
!  RI  RPAR(5)  REAL PARAMETERS:
!      RPAR(1)  MAXIMUM STEPSIZE.
!      RPAR(2)  TOLERANCE FOR CHANGE OF VARIABLES.
!      RPAR(3)  TOLERANCE FOR CONSTRAINT VIOLATIONS.
!      RPAR(4)  TOLERANCE FOR THE GRADIENT OF THE LAGRANGIAN FUNCTION.
!      RPAR(5)  PENALTY COEFFICIENT.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!  RO  GMAX  MAXIMUM PARTIAL DERIVATIVE OF THE LAGRANGIAN FUNCTION.
!  RO  CMAX  MAXIMUM CONSTRAINT VIOLATION.
!  II  IPRNT  PRINT SPECIFICATION. IPRNT=0-NO PRINT.
!         ABS(IPRNT)=1-PRINT OF FINAL RESULTS.
!         ABS(IPRNT)=2-PRINT OF FINAL RESULTS AND ITERATIONS.
!         IPRNT>0-BASIC FINAL RESULTS. IPRNT<0-EXTENDED FINAL
!         RESULTS.
!  IO  ITERM  VARIABLE THAT INDICATES THE CAUSE OF TERMINATION.
!         ITERM=1-IF ABS(X-XO) WAS LESS THAN OR EQUAL TO TOLX IN
!                   MTESX (USUALLY TWO) SUBSEQUEBT ITERATIONS.
!         ITERM=2-IF ABS(F-FO) WAS LESS THAN OR EQUAL TO TOLF IN
!                   MTESF (USUALLY TWO) SUBSEQUEBT ITERATIONS.
!         ITERM=3-IF F IS LESS THAN OR EQUAL TO TOLB.
!         ITERM=4-IF GMAX IS LESS THAN OR EQUAL TO TOLG.
!         ITERM=11-IF NIT EXCEEDED MIT. ITERM=12-IF NFV EXCEEDED MFV.
!         ITERM=13-IF NFG EXCEEDED MFG. ITERM<0-IF THE METHOD FAILED.
!         IF ITERM=-6, THEN THE TERMINATION CRITERION HAS NOT BEEN
!         SATISFIED, BUT THE POINT OBTAINED IF USUALLY ACCEPTABLE.
!
! VARIABLES IN COMMON /STAT/ (STATISTICS) :
!  IO  NRES  NUMBER OF RESTARTS.
!  IO  NDEC  NUMBER OF MATRIX DECOMPOSITION.
!  IO  NREM  NUMBER OF CONSTRAINT DELETIONS.
!  IO  NADD  NUMBER OF CONSTRAINT ADDITIONS.
!  IO  NIT  NUMBER OF ITERATIONS.
!  IO  NFV  NUMBER OF FUNCTION EVALUATIONS.
!  IO  NFG  NUMBER OF GRADIENT EVALUATIONS.
!  IO  NFH  NUMBER OF HESSIAN EVALUATIONS.
!
! SUBPROGRAMS USED :
!  S   PSQP  RECURSIVE QUADRATIC PROGRAMMING METHOD WITH THE BFGS
!         VARIABLE METRIC UPDATE.
!
! EXTERNAL SUBROUTINES :
!  SE  OBJ  COMPUTATION OF THE VALUE OF THE OBJECTIVE FUNCTION.
!         CALLING SEQUENCE: CALL OBJ(NF,X,FF) WHERE NF IS THE NUMBER
!         OF VARIABLES, X(NF) IS A VECTOR OF VARIABLES AND FF IS THE
!         VALUE OF THE OBJECTIVE FUNCTION.
!  SE  DOBJ  COMPUTATION OF THE GRADIENT OF THE OBJECTIVE FUNCTION.
!         CALLING SEQUENCE: CALL DOBJ(NF,X,GF) WHERE NF IS THE NUMBER
!         OF VARIABLES, X(NF) IS A VECTOR OF VARIABLES AND GC(NF) IS
!         THE GRADIENT OF THE OBJECTIVE FUNCTION.
!  SE  CON  COMPUTATION OF THE VALUE OF THE CONSTRAINT FUNCTION.
!         CALLING SEQUENCE: CALL CON(NF,KC,X,FC) WHERE NF IS THE
!         NUMBER OF VARIABLES, KC IS THE INDEX OF THE CONSTRAINT
!         FUNCTION, X(NF) IS A VECTOR OF VARIABLES AND FC IS THE
!         VALUE OF THE CONSTRAINT FUNCTION.
!  SE  DCON  COMPUTATION OF THE GRADIENT OF THE CONSTRAINT FUNCTION.
!         CALLING SEQUENCE: CALL DCON(NF,KC,X,GC) WHERE NF IS THE
!         NUMBER OF VARIABLES, KC IS THE INDEX OF THE CONSTRAINT
!         FUNCTION, X(NF) IS A VECTOR OF VARIABLES AND GC(NF) IS THE
!         GRADIENT OF THE CONSTRAINT FUNCTION.
!
      SUBROUTINE PSQPN(Nf,Nb,Nc,X,Ix,Xl,Xu,Cf,Ic,Cl,Cu,Ipar,Rpar,F,Gmax,&
                     & Cmax,Iprnt,Iterm)
      IMPLICIT NONE
!
!     POINTERS FOR AUXILIARY ARRAYS
!
      DOUBLE PRECISION F , Cmax , Gmax
      INTEGER Iprnt , Iterm , Nb , Nc , Nf
      DOUBLE PRECISION Cf(*) , Cl(*) , Cu(*) , Rpar(5) , X(*) , Xl(*) , &
                     & Xu(*)
      INTEGER Ic(*) , Ipar(6) , Ix(*)
      INTEGER NADd , NDEc , NFG , NFH , NFV , NIT , NREm , NREs
      INTEGER lcfd , lcfo , lcg , lcp , lcr , lcz , lg , lgc , lgf ,    &
            & lgo , lh , lia , ls , lxo
      COMMON /STAT  / NREs , NDEc , NREm , NADd , NIT , NFV , NFG , NFH
      INTEGER ia(:)
      DOUBLE PRECISION ra(:)
      ALLOCATABLE ia , ra
      ALLOCATE (ia(Nf),ra((Nf+Nc+8)*Nf+3*Nc+1))
      lcg = 1
      lcfo = lcg + Nf*Nc
      lcfd = lcfo + Nc + 1
      lgc = lcfd + Nc
      lcr = lgc + Nf
      lcz = lcr + Nf*(Nf+1)/2
      lcp = lcz + Nf
      lgf = lcp + Nc
      lg = lgf + Nf
      lh = lg + Nf
      ls = lh + Nf*(Nf+1)/2
      lxo = ls + Nf
      lgo = lxo + Nf
      lia = 1
      CALL PSQP(Nf,Nb,Nc,X,Ix,Xl,Xu,Cf,Ic,Cl,Cu,ra,ra(lcfo),ra(lcfd),   &
              & ra(lgc),ia,ra(lcr),ra(lcz),ra(lcp),ra(lgf),ra(lg),ra(lh)&
              & ,ra(ls),ra(lxo),ra(lgo),Rpar(1),Rpar(2),Rpar(3),Rpar(4),&
              & Rpar(5),Cmax,Gmax,F,Ipar(1),Ipar(2),Ipar(5),Ipar(6),    &
              & Iprnt,Iterm)
      DEALLOCATE (ia,ra)
      END

!***********************************************************************
! SUBROUTINE PSQP               ALL SYSTEMS                   97/01/22
! PURPOSE :
! RECURSIVE QUADRATIC PROGRAMMING METHOD WITH THE BFGS VARIABLE METRIC
! UPDATE FOR GENERAL NONLINEAR PROGRAMMING PROBLEMS.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  II  NB  CHOICE OF SIMPLE BOUNDS. NB=0-SIMPLE BOUNDS SUPPRESSED.
!         NB>0-SIMPLE BOUNDS ACCEPTED.
!  II  NC  NUMBER OF LINEAR CONSTRAINTS.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  II  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS. IX(I)=0-VARIABLE
!         X(I) IS UNBOUNDED. IX(I)=1-LOVER BOUND XL(I).LE.X(I).
!         IX(I)=2-UPPER BOUND X(I).LE.XU(I). IX(I)=3-TWO SIDE BOUND
!         XL(I).LE.X(I).LE.XU(I). IX(I)=5-VARIABLE X(I) IS FIXED.
!  RI  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
!  RI  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
!  RO  CF(NC+1)  VECTOR CONTAINING VALUES OF THE CONSTRAINT FUNCTIONS.
!  II  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
!         IC(KC)=0-CONSTRAINT CF(KC) IS NOT USED. IC(KC)=1-LOVER
!         CONSTRAINT CL(KC).LE.CF(KC). IC(KC)=2-UPPER CONSTRAINT
!         CF(KC).LE.CU(KC). IC(KC)=3-TWO SIDE CONSTRAINT
!         CL(KC).LE.CF(KC).LE.CU(KC). IC(KC)=5-EQUALITY CONSTRAINT
!         CF(KC).EQ.CL(KC).
!  RI  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
!  RI  CG(NF*NC)  MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
!         CONSTRAINTS.
!  RA  CFO(NC)  VECTOR CONTAINING SAVED VALUES OF THE CONSTRAINT
!         FUNCTIONS.
!  RA  CFD(NC)  VECTOR CONTAINING INCREMENTS OF THE CONSTRAINT
!         FUNCTIONS.
!  RA  GC(NF)  GRADIENT OF THE SELECTED CONSTRAINT FUNCTION.
!  IO  ICA(NF)  VECTOR CONTAINING INDICES OF ACTIVE CONSTRAINTS.
!  RO  CR(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OF KERNEL OF THE
!         ORTHOGONAL PROJECTION.
!  RO  CZ(NF)  VECTOR OF LAGRANGE MULTIPLIERS.
!  RO  GF(NF)  GRADIENT OF THE MODEL FUNCTION.
!  RO  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RU  H(NF*(NF+1)/2)  TRIANGULAR DECOMPOSITION OR INVERSION OF THE
!         HESSIAN MATRIX APPROXIMATION.
!  RO  S(NF)  DIRECTION VECTOR.
!  RU  XO(NF)  VECTORS OF VARIABLES DIFFERENCE.
!  RI  GO(NF)  GRADIENTS DIFFERENCE.
!  RI  XMAX  MAXIMUM STEPSIZE.
!  RI  TOLX  TOLERANCE FOR CHANGE OF VARIABLES.
!  RI  TOLC  TOLERANCE FOR CONSTRAINT VIOLATIONS.
!  RI  TOLG  TOLERANCE FOR THE GRADIENT OF THE LAGRANGIAN FUNCTION.
!  RI  RPF  PENALTY COEFFICIENT.
!  RO  CMAX  MAXIMUM CONSTRAINT VIOLATION.
!  RO  GMAX  MAXIMUM PARTIAL DERIVATIVE OF THE LAGRANGIAN FUNCTION.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!         FUNCTIONS.
!  II  MIT  MAXIMUN NUMBER OF ITERATIONS.
!  II  MFV  MAXIMUN NUMBER OF FUNCTION EVALUATIONS.
!  II  MET  VARIABLE METRIC UPDATE USED. MET=1-THE BFGS UPDATE.
!         MET=2-THE HOSHINO UPDATE.
!  II  MEC  CORRECTION IF THE NEGATIVE CURVATURE OCCURS.
!         MEC=1-CORRECTION SUPPRESSED. MEC=2-POWELL'S CORRECTION.
!  II  IPRNT  PRINT SPECIFICATION. IPRNT=0-NO PRINT.
!         ABS(IPRNT)=1-PRINT OF FINAL RESULTS.
!         ABS(IPRNT)=2-PRINT OF FINAL RESULTS AND ITERATIONS.
!         IPRNT>0-BASIC FINAL RESULTS. IPRNT<0-EXTENDED FINAL
!         RESULTS.
!  IO  ITERM  VARIABLE THAT INDICATES THE CAUSE OF TERMINATION.
!         ITERM=1-IF ABS(X-XO) WAS LESS THAN OR EQUAL TO TOLX IN
!                   MTESX (USUALLY TWO) SUBSEQUEBT ITERATIONS.
!         ITERM=2-IF ABS(F-FO) WAS LESS THAN OR EQUAL TO TOLF IN
!                   MTESF (USUALLY TWO) SUBSEQUEBT ITERATIONS.
!         ITERM=3-IF F IS LESS THAN OR EQUAL TO TOLB.
!         ITERM=4-IF GMAX IS LESS THAN OR EQUAL TO TOLG.
!         ITERM=11-IF NIT EXCEEDED MIT. ITERM=12-IF NFV EXCEEDED MFV.
!         ITERM=13-IF NFG EXCEEDED MFG. ITERM<0-IF THE METHOD FAILED.
!         IF ITERM=-6, THEN THE TERMINATION CRITERION HAS NOT BEEN
!         SATISFIED, BUT THE POINT OBTAINED IF USUALLY ACCEPTABLE.
!
! VARIABLES IN COMMON /STAT/ (STATISTICS) :
!  IO  NRES  NUMBER OF RESTARTS.
!  IO  NDEC  NUMBER OF MATRIX DECOMPOSITION.
!  IO  NREM  NUMBER OF CONSTRAINT DELETIONS.
!  IO  NADD  NUMBER OF CONSTRAINT ADDITIONS.
!  IO  NIT  NUMBER OF ITERATIONS.
!  IO  NFV  NUMBER OF FUNCTION EVALUATIONS.
!  IO  NFG  NUMBER OF GRADIENT EVALUATIONS.
!  IO  NFH  NUMBER OF HESSIAN EVALUATIONS.
!
! SUBPROGRAMS USED :
!  S   PC1F01  COMPUTATION OF THE VALUE AND THE GRADIENT OF THE
!         CONSTRAINT FUNCTION.
!  S   PF1F01  COMPUTATION OF THE VALUE AND THE GRADIENT OF THE
!         OBJECTIVE FUNCTION.
!  S   PLQDB1  GENERAL QUADRATIC PROGRAMMING SUBROUTINE BASED ON THE
!         GOLDFARB-IDNANI DUAL METHOD.
!  S   PLNEWS  IDENTIFICATION OF ACTIVE SIMPLE BOUNDS.
!  S   PLREDL  TRANSFORMATION OF THE INCOMPATIBLE QUADRATIC PROGRAMMING
!         SUBPROBLEM.
!  S   PP0AF8  COMPUTATION OF VALUE OF THE AUGMENTED LAGRANGIAN
!         FUNCTION.
!  S   PPSET2  COMPUTATION OF THE NEW PENALTY PARAMETERS.
!  S   PS0L02  LINE SEARCH USING ONLY FUNCTION VALUES.
!  S   PYTRND  DETERMINATION OF DIFFERENCES FOR VARIABLE METRIC
!         UPDATES.
!  S   PUDBG1  VARIABLE METRIC UPDATE AFTER GILL-MURRAY DECOMPOSITION.
!  S   MXDSMI  SYMMETRIC MATRIX IS REPLACED BY THE UNIT MATRIX.
!  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVINA  ABSOLUTE VALUES OF ELEMENTS OF AN INTEGER VECTOR.
!  RF  MXVMAX  L-INFINITY NORM OF A VECTOR.
!  S   MXVSET  INITIATION OF A VECTOR.
!
! EXTERNAL SUBROUTINES :
!  SE  OBJ  COMPUTATION OF THE VALUE OF THE OBJECTIVE FUNCTION.
!         CALLING SEQUENCE: CALL OBJ(NF,X,FF) WHERE NF IS THE NUMBER
!         OF VARIABLES, X(NF) IS A VECTOR OF VARIABLES AND FF IS THE
!         VALUE OF THE OBJECTIVE FUNCTION.
!  SE  DOBJ  COMPUTATION OF THE GRADIENT OF THE OBJECTIVE FUNCTION.
!         CALLING SEQUENCE: CALL DOBJ(NF,X,GF) WHERE NF IS THE NUMBER
!         OF VARIABLES, X(NF) IS A VECTOR OF VARIABLES AND GC(NF) IS
!         THE GRADIENT OF THE OBJECTIVE FUNCTION.
!  SE  CON  COMPUTATION OF THE VALUE OF THE CONSTRAINT FUNCTION.
!         CALLING SEQUENCE: CALL CON(NF,KC,X,FC) WHERE NF IS THE
!         NUMBER OF VARIABLES, KC IS THE INDEX OF THE CONSTRAINT
!         FUNCTION, X(NF) IS A VECTOR OF VARIABLES AND FC IS THE
!         VALUE OF THE CONSTRAINT FUNCTION.
!  SE  DCON  COMPUTATION OF THE GRADIENT OF THE CONSTRAINT FUNCTION.
!         CALLING SEQUENCE: CALL DCON(NF,KC,X,GC) WHERE NF IS THE
!         NUMBER OF VARIABLES, KC IS THE INDEX OF THE CONSTRAINT
!         FUNCTION, X(NF) IS A VECTOR OF VARIABLES AND GC(NF) IS THE
!         GRADIENT OF THE CONSTRAINT FUNCTION.
!
! METHOD :
! RECURSIVE QUADRATIC PROGRAMMING METHOD WITH THE BFGS VARIABLE METRIC
! UPDATE.
!
      SUBROUTINE PSQP(Nf,Nb,Nc,X,Ix,Xl,Xu,Cf,Ic,Cl,Cu,Cg,Cfo,Cfd,Gc,Ica,&
                    & Cr,Cz,Cp,Gf,G,H,S,Xo,Go,Xmax,Tolx,Tolc,Tolg,Rpf,  &
                    & Cmax,Gmax,F,Mit,Mfv,Met,Mec,Iprnt,Iterm)
      IMPLICIT NONE
      DOUBLE PRECISION F , Cmax , Gmax , Rpf , Tolc , told , Tolg ,     &
                     & tols , Tolx , Xmax
      INTEGER Iprnt , Iterm , Met , met1 , Mec , mes , Mfv , Mit , Nb , &
            & Nc , Nf
      DOUBLE PRECISION Cf(*) , Cfd(*) , Cfo(*) , Cg(*) , Cl(*) , Cp(*) ,&
                     & Cr(*) , Cz(*) , Cu(*) , G(*) , Gc(*) , Gf(*) ,   &
                     & Go(*) , H(*) , S(*) , X(*) , Xl(*) , Xo(*) ,     &
                     & Xu(*)
      INTEGER Ic(*) , Ica(*) , Ix(*)
      INTEGER NADd , NDEc , NFG , NFH , NFV , NIT , NREm , NREs
      DOUBLE PRECISION alf1 , alf2 , cmaxo , dmax , eps7 , eps9 , eta0 ,&
                     & eta2 , eta9 , fmax , fmin , fo , gnorm , p , po ,&
                     & r , rmax , rmin , ro , snorm , tolb , tolf ,     &
                     & umax , rp , fp , pp , ff , fc
      INTEGER i , idecf , iext , irest , iterd , iterl , iterh , iterq ,&
            & iters , kbc , kbf , kc , kd , kit , ld , mred , mtesf ,   &
            & mtesx , n , k , ntesx , iest , inits , kters , maxst ,    &
            & isys , mfp , nred , ipom , lds
      DOUBLE PRECISION MXVDOT , MXVMAX
      COMMON /STAT  / NREs , NDEc , NREm , NADd , NIT , NFV , NFG , NFH
      IF ( ABS(Iprnt)>1 ) WRITE (6,'(1X,''ENTRY TO PSQP :'')')
!
!     INITIATION
!
      kbf = 0
      kbc = 0
      IF ( Nb>0 ) kbf = 2
      IF ( Nc>0 ) kbc = 2
      NREs = 0
      NDEc = 0
      NREm = 0
      NADd = 0
      NIT = 0
      NFV = 0
      NFG = 0
      NFH = 0
      isys = 0
      iest = 0
      iext = 0
      mtesx = 2
      mtesf = 2
      inits = 1
      Iterm = 0
      iters = 0
      iterd = 0
      iterq = 0
      mred = 20
      irest = 1
      iters = 2
      kters = 5
      idecf = 1
      eta0 = 1.0D-15
      eta2 = 1.0D-15
      eta9 = 1.0D60
      eps7 = 1.0D-15
      eps9 = 1.0D-8
      alf1 = 1.0D-10
      alf2 = 1.0D10
      fmax = 1.0D60
      fmin = -fmax
      tolb = -fmax
      dmax = eta9
      tolf = 1.0D-16
      IF ( Xmax<=0.0D0 ) Xmax = 1.0D+16
      IF ( Tolx<=0.0D0 ) Tolx = 1.0D-16
      IF ( Tolg<=0.0D0 ) Tolg = 1.0D-6
      IF ( Tolc<=0.0D0 ) Tolc = 1.0D-6
      told = 1.0D-8
      tols = 1.0D-4
      IF ( Rpf<=0.0D0 ) Rpf = 1.0D-4
      IF ( Met<=0 ) Met = 1
      met1 = 2
      IF ( Mec<=0 ) Mec = 2
      mes = 1
      IF ( Mit<=0 ) Mit = 1000
      IF ( Mfv<=0 ) Mfv = 2000
      kd = 1
      ld = -1
      kit = 0
      CALL MXVSET(Nc,0.0D0,Cp)
!
!     INITIAL OPERATIONS WITH SIMPLE BOUNDS
!
      IF ( kbf>0 ) THEN
         DO i = 1 , Nf
            IF ( (Ix(i)==3 .OR. Ix(i)==4) .AND. Xu(i)<=Xl(i) ) THEN
               Xu(i) = Xl(i)
               Ix(i) = 5
            ELSEIF ( Ix(i)==5 .OR. Ix(i)==6 ) THEN
               Xl(i) = X(i)
               Xu(i) = X(i)
               Ix(i) = 5
            ENDIF
            IF ( Ix(i)==1 .OR. Ix(i)==3 ) X(i) = MAX(X(i),Xl(i))
            IF ( Ix(i)==2 .OR. Ix(i)==3 ) X(i) = MIN(X(i),Xu(i))
         ENDDO
      ENDIF
!     INITIAL OPERATIONS WITH GENERAL CONSTRAINTS
!
      IF ( kbc>0 ) THEN
         k = 0
         DO kc = 1 , Nc
            IF ( (Ic(kc)==3 .OR. Ic(kc)==4) .AND. Cu(kc)<=Cl(kc) ) THEN
               Cu(kc) = Cl(kc)
               Ic(kc) = 5
            ELSEIF ( Ic(kc)==5 .OR. Ic(kc)==6 ) THEN
               Cu(kc) = Cl(kc)
               Ic(kc) = 5
            ENDIF
            k = k + Nf
         ENDDO
      ENDIF
      IF ( kbf>0 ) THEN
         DO i = 1 , Nf
            IF ( Ix(i)>=5 ) Ix(i) = -Ix(i)
            IF ( Ix(i)<=0 ) THEN
            ELSEIF ( (Ix(i)==1 .OR. Ix(i)==3) .AND. X(i)<=Xl(i) ) THEN
               X(i) = Xl(i)
            ELSEIF ( (Ix(i)==2 .OR. Ix(i)==3) .AND. X(i)>=Xu(i) ) THEN
               X(i) = Xu(i)
            ENDIF
            CALL PLNEWS(X,Ix,Xl,Xu,eps9,i,iterl)
            IF ( Ix(i)>10 ) Ix(i) = 10 - Ix(i)
         ENDDO
      ENDIF
      fo = fmin
      Gmax = eta9
      dmax = eta9
 100  lds = ld
      CALL PF1F01(Nf,X,Gf,Gf,ff,F,kd,ld,iext)
      ld = lds
      CALL PC1F01(Nf,Nc,X,fc,Cf,Cl,Cu,Ic,Gc,Cg,Cmax,kd,ld)
      Cf(Nc+1) = F
      IF ( ABS(Iprnt)>1 ) WRITE (6,                                     &
     &'(1X,''NIT='',I4,2X,''NFV='',I4,2X,''NFG='',I4,2X,       ''F='',G1&
     &2.6,2X,''C='',E7.1,2X,''G='',E7.1)') NIT , NFV , NFG , F , Cmax , &
    & Gmax
!
!     START OF THE ITERATION WITH TESTS FOR TERMINATION.
!
      IF ( Iterm<0 ) GOTO 500
      IF ( iters/=0 ) THEN
         IF ( F<=tolb ) THEN
            Iterm = 3
            GOTO 500
         ENDIF
         IF ( dmax<=Tolx ) THEN
            Iterm = 1
            ntesx = ntesx + 1
            IF ( ntesx>=mtesx ) GOTO 500
         ELSE
            ntesx = 0
         ENDIF
      ENDIF
      IF ( NIT>=Mit ) THEN
         Iterm = 11
         GOTO 500
      ENDIF
      IF ( NFV>=Mfv ) THEN
         Iterm = 12
         GOTO 500
      ENDIF
      Iterm = 0
      NIT = NIT + 1
!
!     RESTART
!
 200  n = Nf
      IF ( irest>0 ) THEN
         CALL MXDSMI(n,H)
         ld = MIN(ld,1)
         idecf = 1
         IF ( kit<NIT ) THEN
            NREs = NREs + 1
            kit = NIT
         ELSE
            Iterm = -10
            IF ( iters<0 ) Iterm = iters - 5
            GOTO 500
         ENDIF
      ENDIF
!
!     DIRECTION DETERMINATION USING A QUADRATIC PROGRAMMING PROCEDURE
!
      CALL MXVCOP(Nc+1,Cf,Cfo)
      mfp = 2
      ipom = 0
 300  CALL PLQDB1(Nf,Nc,X,Ix,Xl,Xu,Cf,Cfd,Ic,Ica,Cl,Cu,Cg,Cr,Cz,G,Gf,H, &
                & S,mfp,kbf,kbc,idecf,eta2,eta9,eps7,eps9,umax,Gmax,n,  &
                & iterq)
      IF ( iterq<0 ) THEN
         IF ( ipom<10 ) THEN
            ipom = ipom + 1
            CALL PLREDL(Nc,Cf,Ic,Cl,Cu,kbc)
            GOTO 300
         ENDIF
         iterd = iterq - 10
         GOTO 400
      ENDIF
      ipom = 0
      iterd = 1
      Gmax = MXVMAX(Nf,G)
      gnorm = SQRT(MXVDOT(Nf,G,G))
      snorm = SQRT(MXVDOT(Nf,S,S))
 400  IF ( iterd<0 ) Iterm = iterd
      IF ( Iterm==0 ) THEN
         CALL MXVCOP(Nc+1,Cfo,Cf)
!
!     TEST FOR SUFFICIENT DESCENT
!
         p = MXVDOT(Nf,G,S)
         irest = 1
         IF ( snorm<=0.0D0 ) THEN
         ELSEIF ( p+told*gnorm*snorm<=0.0D0 ) THEN
            irest = 0
         ENDIF
         IF ( irest/=0 ) GOTO 200
         nred = 0
         rmin = alf1*gnorm/snorm
         rmax = MIN(alf2*gnorm/snorm,Xmax/snorm)
         IF ( Gmax<=Tolg .AND. Cmax<=Tolc ) THEN
            Iterm = 4
            GOTO 500
         ENDIF
         CALL PPSET2(Nf,n,Nc,Ica,Cz,Cp)
         CALL MXVINA(Nc,Ic)
         CALL PP0AF8(Nf,n,Nc,Cf,Ic,Ica,Cl,Cu,Cz,Rpf,fc,F)
!
!     PREPARATION OF LINE SEARCH
!
         ro = 0.0D0
         fo = F
         po = p
         cmaxo = Cmax
         CALL MXVCOP(Nf,X,Xo)
         CALL MXVCOP(Nf,G,Go)
         CALL MXVCOP(Nf,Gf,Cr)
         CALL MXVCOP(Nc+1,Cf,Cfo)
!
!     LINE SEARCH WITHOUT DIRECTIONAL DERIVATIVES
!
 450     CALL PS0L02(r,ro,rp,F,fo,fp,po,pp,fmin,fmax,rmin,rmax,tols,kd, &
                   & ld,NIT,kit,nred,mred,maxst,iest,inits,iters,kters, &
                   & mes,isys)
         IF ( isys==0 ) THEN
            kd = 1
!
!     DECISION AFTER UNSUCCESSFUL LINE SEARCH
!
            IF ( iters<=0 ) THEN
               r = 0.0D0
               F = fo
               p = po
               CALL MXVCOP(Nf,Xo,X)
               CALL MXVCOP(Nf,Cr,Gf)
               CALL MXVCOP(Nc+1,Cfo,Cf)
               irest = 1
               ld = kd
               GOTO 200
            ENDIF
!
!     COMPUTATION OF THE VALUE AND THE GRADIENT OF THE OBJECTIVE
!     FUNCTION TOGETHER WITH THE VALUES AND THE GRADIENTS OF THE
!     APPROXIMATED FUNCTIONS
!
            IF ( kd>ld ) THEN
               lds = ld
               CALL PF1F01(Nf,X,Gf,Gf,ff,F,kd,ld,iext)
               ld = lds
               CALL PC1F01(Nf,Nc,X,fc,Cf,Cl,Cu,Ic,Gc,Cg,Cmax,kd,ld)
            ENDIF
!
!     PREPARATION OF VARIABLE METRIC UPDATE
!
            CALL MXVCOP(Nf,Gf,G)
            CALL PYTRND(Nf,n,X,Xo,Ica,Cg,Cz,G,Go,r,F,fo,p,po,Cmax,cmaxo,&
                      & dmax,kd,ld,iters)
!
!     VARIABLE METRIC UPDATE
!
            CALL PUDBG1(n,H,G,S,Xo,Go,r,po,NIT,kit,iterh,Met,met1,Mec)
!      IF (MER.GT.0.AND.ITERH.GT.0) IREST=1
!
!     END OF THE ITERATION
!
            GOTO 100
         ELSE
!      GO TO (11174,11172) ISYS+1
            CALL MXVDIR(Nf,r,S,Xo,X)
            lds = ld
            CALL PF1F01(Nf,X,Gf,G,ff,F,kd,ld,iext)
            ld = lds
            CALL PC1F01(Nf,Nc,X,fc,Cf,Cl,Cu,Ic,Gc,Cg,Cmax,kd,ld)
            Cf(Nc+1) = F
            CALL PP0AF8(Nf,n,Nc,Cf,Ic,Ica,Cl,Cu,Cz,Rpf,fc,F)
            GOTO 450
         ENDIF
      ENDIF

 500  IF ( Iprnt>1 .OR. Iprnt<0 ) WRITE (6,'(1X,''EXIT FROM PSQP :'')')
      IF ( Iprnt/=0 ) WRITE (6,                                         &
     &'(1X,''NIT='',I4,2X,''NFV='',I4,2X,''NFG='',I4,2X,       ''F='',G1&
     &2.6,2X,''C='',E7.1,2X,''G='',E7.1,2X,''ITERM='',I3)') NIT , NFV , &
    & NFG , F , Cmax , Gmax , Iterm
      IF ( Iprnt<0 ) WRITE (6,                                          &
                           &'(1X,''X='',5(G14.7,1X):/(3X,5(G14.7,1X)))')&
                          & (X(i),i=1,Nf)
      END
