************************************************************************
*     PROGRAM TSQPN
*
*     TEST PROGRAM FOR THE SUBROUTINE PMINU
*
      INTEGER          NADD,NDEC,NEXT,NFG,NFH,NFV,NIT,NREM,NRES
      DOUBLE PRECISION F,FMIN,CMAX,GMAX
      INTEGER          I,IERR,IPRNT,ITERM,ITIME,NB,NC,NF
      DOUBLE PRECISION CF(200),CL(200),CU(200),RPAR(5),X(40),
     +                 XL(40),XU(40)
      INTEGER          IC(200),IX(40),IPAR(6)
      INTEGER          NITER,NFVAL,NSUCC
      COMMON           /PROB/NEXT
      COMMON           /STAT/NRES,NDEC,NREM,NADD,NIT,NFV,NFG,NFH
      NITER=0
      NFVAL=0
      NSUCC=0
      CALL TYTIM1(ITIME)
*
*     LOOP FOR 34 TEST PROBLEMS
*
      DO 30 NEXT = 1,34
*
*     CHOICE OF INTEGER AND REAL PARAMETERS
*
          DO 10 I = 1,6
              IPAR(I) = 0
   10     CONTINUE
          DO 20 I = 1,5
              RPAR(I) = 0.0D0
   20     CONTINUE
          IPAR(5)=2
          IPRNT = 1
*
*     PROBLEM DIMENSION
*
          NF = 20
          NB = 20
          NC = 30
*
*     INITIATION OF X AND CHOICE OF RPAR(1)
*
          CALL TIND07(NF,NC,X,IX,XL,XU,IC,CL,CU,FMIN,RPAR(1),NEXT,IERR)
          IF (IERR.NE. 0) GO TO 30
          RPAR(1)=0.0D 0
*
*     SOLUTION
*
          CALL PSQPN(NF,NB,NC,X,IX,XL,XU,CF,IC,CL,CU,IPAR,RPAR,
     +                 F,CMAX,GMAX,IPRNT,ITERM)
      NITER=NITER+NIT
      NFVAL=NFVAL+NFV
      IF (ITERM.GT.0.AND.ITERM.LT.9) NSUCC=NSUCC+1
   30 CONTINUE
      WRITE(6,40) NITER,NFVAL,NSUCC
   40 FORMAT(' NITER =',I5,3X,' NFVAL =',I5,3X,' NSUCC =',I5)
      CALL TYTIM2(ITIME)
      STOP
      END
*
*     USER SUPPLIED SUBROUTINE (CALCULATION OF FF)
*
      SUBROUTINE OBJ(NF,X,FF)
*
*     FUNCTION EVALUATION
*
      DOUBLE PRECISION FF
      INTEGER          NF
      DOUBLE PRECISION X(*)
      INTEGER          NEXT
      EXTERNAL         TFFU07
      COMMON           /PROB/NEXT
      CALL TFFU07(NF,X,FF,NEXT)
      RETURN
      END
*
*     USER SUPPLIED SUBROUTINE (CALCULATION OF GF)
*
      SUBROUTINE DOBJ(NF,X,GF)
*
*     GRADIENT EVALUATION
*
      INTEGER          NF
      DOUBLE PRECISION GF(*),X(*)
      INTEGER          NEXT
      EXTERNAL         TFGU07
      COMMON           /PROB/NEXT
      CALL TFGU07(NF,X,GF,NEXT)
      RETURN
      END
*
*     USER SUPPLIED SUBROUTINE (CALCULATION OF FC)
*
      SUBROUTINE CON(NF,KC,X,FC)
*
*     FUNCTION EVALUATION
*
      DOUBLE PRECISION FC
      INTEGER          KC,NF
      DOUBLE PRECISION X(*)
      INTEGER          NEXT
      EXTERNAL         TCFU07
      COMMON           /PROB/NEXT
      CALL TCFU07(NF,KC,X,FC,NEXT)
      RETURN
      END
*
*     USER SUPPLIED SUBROUTINE (CALCULATION OF GC)
*
      SUBROUTINE DCON(NF,KC,X,GC)
*
*     GRADIENT EVALUATION
*
      INTEGER          KC,NF
      DOUBLE PRECISION GC(*),X(*)
      INTEGER          NEXT
      EXTERNAL         TCGU07
      COMMON           /PROB/NEXT
      CALL TCGU07(NF,KC,X,GC,NEXT)
      RETURN
      END
*
*     EMPTY SUBROUTINES
*
      SUBROUTINE FUN(NF,KA,X,FA)
      INTEGER KA,NF
      DOUBLE PRECISION FA,X(*)
      KA=NF
      FA=X(1)
      RETURN
      END
      SUBROUTINE DFUN(NF,KA,X,GA)
      INTEGER KA,NF
      DOUBLE PRECISION GA(*),X(*)
      KA=NF
      GA(1)=X(1)
      RETURN
      END
