* SUBROUTINE TILD03                ALL SYSTEMS                92/11/01
C PORTABILITY : ALL SYSTEMS
C 92/11/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  INITIAL VALUES OF THE VARIABLES FOR CONSTRAINED MINIMIZATION.
*  LINEARLY CONSTRAINED AND DENSE VERSION.
*
* PARAMETERS :
*  IO  NF  NUMBER OF VARIABLES.
*  IO  NC  NUMBER OF CONSTRAINTS.
*  II  NCL  NUMBER OF LINEAR CONSTRAINTS.
*  RO  X(NF)  VECTOR OF VARIABLES.
*  IO  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  RO  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
*  RO  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
*  IO  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  RO  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RO  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RO  CG(NF*NC) MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*  IO  IERR  ERROR INDICATOR.
*
      SUBROUTINE TILD03(NF,NC,X,IX,XL,XU,IC,CL,CU,CG,FMIN,XMAX,NEXT,
     * IERR)
      INTEGER NF,NC,IX(NF),IC(NC),NEXT,IEXT,IERR
      DOUBLE PRECISION X(NF),XL(NF),XU(NF),CL(NC),CU(NC),CG(NF*NC),
     * FMIN,XMAX
      INTEGER I,J,NCL
      DOUBLE PRECISION Y(235)
      COMMON /EMPR03/ Y
      DOUBLE PRECISION ETA9
      PARAMETER (ETA9=1.0D 60)
      FMIN=0.0D 0
      XMAX=1.0D 3
      IEXT=0
      IERR=0
      GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150) NEXT
   10 NF=2
      NC=1
      NCL=1
      X(1)=0.0D 0
      X(2)=0.0D 0
      IX(1)=0
      IX(2)=0
      IC(1)=5
      CL(1)=0.0D 0
      CU(1)=0.0D 0
      CG(1)=4.0D 0
      CG(2)=-3.0D 0
      FMIN=-ETA9
      RETURN
   20 NF=2
      NC=2
      NCL=2
      X(1)=1.0D 0
      X(2)=0.5D 0
      IX(1)=1
      IX(2)=1
      XL(1)=0.0D 0
      XL(2)=0.0D 0
      IC(1)=1
      IC(2)=3
      CL(1)=0.0D 0
      CL(2)=0.0D 0
      CU(2)=6.0D 0
      CG(1)=1.0D 0/SQRT(3.0D 0)
      CG(2)=-1.0D 0
      CG(3)=1.0D 0
      CG(4)=SQRT(3.0D 0)
      FMIN=-ETA9
      RETURN
   30 NF=3
      NC=1
      NCL=1
      X(1)=0.7D 0
      X(2)=0.2D 0
      X(3)=0.1D 0
      DO 31 I=1,NF
      IX(I)=3
      XL(I)=0.0D 0
      XU(I)=1.0D 0
      CG(I)=1.0D 0
   31 CONTINUE
      IC(1)=5
      CL(1)=1.0D 0
      CU(1)=1.0D 0
      FMIN=-ETA9
      RETURN
   40 NF=3
      NC=1
      NCL=1
      DO 41 I=1,NF
      X(I)=1.0D 1
      IX(I)=3
      XL(I)=0.0D 0
      XU(I)=4.2D 1
   41 CONTINUE
      IC(1)=3
      CL(1)=0.0D 0
      CU(1)=7.2D 1
      CG(1)=1.0D 0
      CG(2)=2.0D 0
      CG(3)=2.0D 0
      FMIN=-ETA9
      RETURN
   50 NF=4
      NC=1
      NCL=1
      DO 51 I=1,NF
      X(I)=2.0D 0
      IX(I)=3
      XL(I)=0.0D 0
      XU(I)=1.0D 0
   51 CONTINUE
      XU(4)=2.0D 0
      IC(1)=5
      CL(1)=0.0D 0
      CU(1)=0.0D 0
      CG(1)=1.0D 0
      CG(2)=2.0D 0
      CG(3)=2.0D 0
      CG(4)=-1.0D 0
      RETURN
   60 NF=4
      NC=1
      NCL=1
      X(1)=3.0D 0
      X(2)=-1.0D 0
      X(3)=0.0D 0
      X(4)=1.0D 0
      DO 61 I=1,NF
      IX(I)=2
      XU(I)=2.0D 1
      CG(I)=1.0D 0
   61 CONTINUE
      IC(1)=1
      CL(1)=1.0D 0
      RETURN
   70 NF=5
      NC=2
      NCL=2
      X(1)=1.0D 1
      X(2)=7.0D 0
      X(3)=2.0D 0
      X(4)=-3.0D 0
      X(5)=0.8D 0
      DO 71 I=1,NF
      IX(I)=0
      CG(I)=1.0D 0
      CG(I+NF)=0.0D 0
   71 CONTINUE
      IC(1)=5
      IC(2)=5
      CL(1)=7.0D 0
      CU(1)=7.0D 0
      CL(2)=6.0D 0
      CU(2)=6.0D 0
      CG(4)=4.0D 0
      CG(5)=0.0D 0
      CG(8)=1.0D 0
      CG(10)=5.0D 0
      RETURN
   80 NF=5
      NC=3
      NCL=3
      X(1)=3.5D 1
      X(2)=-3.1D 1
      X(3)=1.1D 1
      X(4)=5.0D 0
      X(5)=-5.0D 0
      DO 82 I=1,NF
      IX(I)=0
      DO 81 J=1,NC
      CG((J-1)*NF+I)=0.0D 0
   81 CONTINUE
   82 CONTINUE
      DO 83 J=1,NC
      IC(J)=5
      CL(J)=6.0D 0
      CU(J)=6.0D 0
      CG((J-1)*NF+J)=1.0D 0
      CG((J-1)*NF+J+1)=2.0D 0
      CG((J-1)*NF+J+2)=3.0D 0
   83 CONTINUE
      RETURN
   90 NF=5
      NC=3
      NCL=3
      DO 91 J=1,5
      X(J)=2.0D 0
      IX(J)=0
   91 CONTINUE
      DO 92 J=1,3
      IC(J)=5
   92 CONTINUE
      DO 93 J=1,3
      CL(J)=0.0D 0
      CU(J)=0.0D 0
   93 CONTINUE
      DO 94 J=1,15
      CG(J)=0.0D 0
   94 CONTINUE
      CG(1)=1.0D 0
      CG(2)=3.0D 0
      CG(8)=1.0D 0
      CG(9)=1.0D 0
      CG(10)=-2.0D 0
      CG(12)=1.0D 0
      CG(15)=-1.0D 0
      RETURN
  100 NF=5
      NC=10
      NCL=10
      DO 101 J=1,5
      X(J)=0.0D 0
      IX(J)=1
      XL(J)=0.0D 0
  101 CONTINUE
      X(5)=1.0D 0
      DO 102 J=1,10
      IC(J)=1
  102 CONTINUE
      DO 103 J=1,50
      CG(J)=0.0D 0
  103 CONTINUE
      CG(1)=-16.0D 0
      CG(2)=2.0D 0
      CG(4)=1.0D 0
      CG(7)=-2.0D 0
      CG(9)=4.0D 0
      CG(10)=2.0D 0
      CG(11)=-3.5D 0
      CG(13)=2.0D 0
      CG(17)=-2.0D 0
      CG(19)=-4.0D 0
      CG(20)=-1.0D 0
      CG(22)=-9.0D 0
      CG(23)=-2.0D 0
      CG(24)=1.0D 0
      CG(25)=-2.8D 0
      CG(26)=2.0D 0
      CG(28)=-4.0D 0
      DO 104 J=1,5
      CG(30+J)=-1.0D 0
  104 CONTINUE
      CG(36)=-1.0D 0
      CG(37)=-2.0D 0
      CG(38)=-3.0D 0
      CG(39)=-2.0D 0
      CG(40)=-1.0D 0
      CG(41)=1.0D 0
      CG(42)=2.0D 0
      CG(43)=3.0D 0
      CG(44)=4.0D 0
      CG(45)=5.0D 0
      DO 105 J=1,5
      CG(45+J)=1.0D 0
  105 CONTINUE
      CL(1)=-40.0D 0
      CL(2)=-2.0D 0
      CL(3)=-0.25D 0
      CL(4)=-4.0D 0
      CL(5)=-4.0D 0
      CL(6)=-1.0D 0
      CL(7)=-40.0D 0
      CL(8)=-60.0D 0
      CL(9)=5.0D 0
      CL(10)=1.0D 0
      FMIN=-ETA9
      RETURN
  110 NF=6
      NC=1
      NCL=1
      X(1)=6.0D 3
      X(2)=1.5D 0
      X(3)=4.0D 6
      X(4)=2.0D 0
      X(5)=3.0D-3
      X(6)=5.0D 7
      DO 111 J=1,6
      IX(J)=3
  111 CONTINUE
      DO 112 J=1,6
      XL(J)=0.0D 0
  112 CONTINUE
      XL(2)=-10.0D 0
      XL(5)=-1.0D 0
      XU(1)=2.0D 4
      XU(2)=10.0D 0
      XU(3)=1.0D 7
      XU(4)=20.0D 0
      XU(5)=1.0D 0
      XU(6)=2.0D 8
      IC(1)=5
      CL(1)=1.76D 4
      CU(1)=CL(1)
      DO 113 J=1,6
      CG(J)=0.0D 0
  113 CONTINUE
      CG(1)=1.0D 0
      CG(2)=4.0D 3
      FMIN=-ETA9
      RETURN
  120 NF=6
      NC=6
      NCL=6
      DO 121 J=1,6
      X(J)=0.0D 0
      IX(J)=1
      XL(J)=0.0D 0
  121 CONTINUE
      IX(1)=3
      IX(4)=3
      XU(1)=1.0D 0
      XU(4)=1.0D 0
      X(1)=1.0D 0
      X(2)=2.0D 0
      X(6)=2.0D 0
      DO 122 J=1,6
      IC(J)=5
  122 CONTINUE
      CL(1)=6.0D 0
      CL(2)=3.0D 0
      CL(3)=2.0D 0
      CL(4)=1.0D 0
      CL(5)=2.0D 0
      CL(6)=2.0D 0
      DO 123 J=1,6
      CU(J)=CL(J)
  123 CONTINUE
      DO 124 J=1,36
      CG(J)=0.0D 0
  124 CONTINUE
      CG(1)=1.0D 0
      CG(2)=2.0D 0
      CG(5)=5.0D 0
      CG(7)=1.0D 0
      CG(8)=1.0D 0
      CG(9)=1.0D 0
      CG(16)=1.0D 0
      CG(17)=1.0D 0
      CG(18)=1.0D 0
      CG(19)=1.0D 0
      CG(22)=1.0D 0
      CG(26)=1.0D 0
      CG(29)=1.0D 0
      CG(33)=1.0D 0
      CG(36)=1.0D 0
      RETURN
  130 NF=8
      NC=1
      NCL=1
      X(1)=0.1D 0
      X(2)=0.2D 0
      X(3)=100.0D 0
      X(4)=125.0D 0
      X(5)=175.0D 0
      X(6)=11.2D 0
      X(7)=13.2D 0
      X(8)=15.8D 0
      DO 131 J=1,8
      IX(J)=3
  131 CONTINUE
      XL(1)=0.001D 0
      XU(1)=0.499D 0
      XL(2)=XL(1)
      XU(2)=XU(1)
      XL(3)=100.0D 0
      XU(3)=180.0D 0
      XL(4)=130.0D 0
      XU(4)=210.0D 0
      XL(5)=170.0D 0
      XU(5)=240.0D 0
      DO 132 J=6,8
      XL(J)=5.0D 0
      XU(J)=25.0D 0
  132 CONTINUE
      IC(1)=2
      CU(1)=1.0D 0
      DO 133 J=3,8
      CG(J)=0.0D 0
  133 CONTINUE
      CG(1)=1.0D 0
      CG(2)=1.0D 0
      Y(1)=95.0D 0
      Y(2)=105.0D 0
      DO 1301  J=3,6
      Y(J)=110.0D 0
 1301 CONTINUE
      DO 1302 J=7,10
      Y(J)=115.0D 0
 1302 CONTINUE
      DO 1303 J=11,25
      Y(J)=120.0D 0
 1303 CONTINUE
      DO 1304 J=26,40
      Y(J)=125.0D 0
 1304 CONTINUE
      DO 1305 J=41,55
      Y(J)=130.0D 0
 1305 CONTINUE
      DO 1306 J=56,68
 1306 Y(J)=135.0D 0
      CONTINUE
      DO 1307 J=69,89
      Y(J)=140.0D 0
 1307 CONTINUE
      DO 1308 J=90,101
      Y(J)=145.0D 0
 1308 CONTINUE
      DO 1309 J=102,118
      Y(J)=150.0D 0
 1309 CONTINUE
      DO 1310 J=119,122
      Y(J)=155.0D 0
 1310 CONTINUE
      DO 1311 J=123,142
      Y(J)=160.0D 0
 1311 CONTINUE
      DO 1312 J=143,150
      Y(J)=165.0D 0
 1312 CONTINUE
      DO 1313 J=151,167
      Y(J)=170.0D 0
 1313 CONTINUE
      DO 1314 J=168,175
      Y(J)=175.0D 0
 1314 CONTINUE
      DO 1315 J=176,181
      Y(J)=180.0D 0
 1315 CONTINUE
      DO 1316 J=182,187
      Y(J)=185.0D 0
 1316 CONTINUE
      DO 1317 J=188,194
      Y(J)=190.0D 0
 1317 CONTINUE
      DO 1318 J=195,198
      Y(J)=195.0D 0
 1318 CONTINUE
      DO 1319 J=199,201
      Y(J)=200.0D 0
 1319 CONTINUE
      DO 1320 J=202,204
      Y(J)=205.0D 0
 1320 CONTINUE
      DO 1321 J=205,212
      Y(J)=210.0D 0
 1321 CONTINUE
      Y(213)=215.0D 0
      DO 1322 J=214,219
      Y(J)=220.0D 0
 1322 CONTINUE
      DO 1323 J=220,224
      Y(J)=230.0D 0
 1323 CONTINUE
      Y(225)=235.0D 0
      DO 1324 J=226,232
      Y(J)=240.0D 0
 1324 CONTINUE
      Y(233)=245.0D 0
      DO 1325 J=234,235
      Y(J)=250.0D 0
 1325 CONTINUE
      RETURN
  140 NF=10
      NC=3
      NCL=3
      DO 141 J=1,10
      X(J)=0.1D 0
      IX(J)=1
      XL(J)=1.0D -6
  141 CONTINUE
      IC(1)=5
      IC(2)=5
      IC(3)=5
      CL(1)=2.0D 0
      CL(2)=1.0D 0
      CL(3)=1.0D 0
      DO 142 J=1,3
      CU(J)=CL(J)
  142 CONTINUE
      DO 143 J=1,30
      CG(J)=0.0D 0
  143 CONTINUE
      CG(1)=1.0D 0
      CG(2)=2.0D 0
      CG(3)=2.0D 0
      CG(6)=1.0D 0
      CG(10)=1.0D 0
      CG(14)=1.0D 0
      CG(15)=2.0D 0
      CG(16)=1.0D 0
      CG(17)=1.0D 0
      CG(23)=1.0D 0
      CG(27)=1.0D 0
      CG(28)=1.0D 0
      CG(29)=2.0D 0
      CG(30)=1.0D 0
      FMIN=-ETA9
      RETURN
  150 NF=16
      NC=8
      NCL=8
      DO 151 J=1,16
      X(J)=10.0D 0
      IX(J)=3
      XL(J)=0.0D 0
      XU(J)=5.0D 0
  151 CONTINUE
      DO 152 J=1,8
      IC(J)=5
  152 CONTINUE
      DO 153 J=1,128
      CG(J)=0.0D 0
  153 CONTINUE
      CG(1)=0.22D 0
      CG(2)=0.20D 0
      CG(3)=0.19D 0
      CG(4)=0.25D 0
      CG(5)=0.15D 0
      CG(6)=0.11D 0
      CG(7)=0.12D 0
      CG(8)=0.13D 0
      CG(9)=1.0D 0
      CG(17)=-1.46D 0
      CG(19)=-1.30D 0
      CG(20)=1.82D 0
      CG(21)=-1.15D 0
      CG(23)=0.80D 0
      CG(26)=1.0D 0
      CG(33)=1.29D 0
      CG(34)=-0.89D 0
      CG(37)=-1.16D 0
      CG(38)=-0.96D 0
      CG(40)=-0.49D 0
      CG(43)=1.0D 0
      CG(49)=-1.10D 0
      CG(50)=-1.06D 0
      CG(51)=0.95D 0
      CG(52)=-0.54D 0
      CG(54)=-1.78D 0
      CG(55)=-0.41D 0
      CG(60)=1.0D 0
      CG(68)=-1.43D 0
      CG(69)=1.51D 0
      CG(70)=0.59D 0
      CG(71)=-0.33D 0
      CG(72)=-0.43D 0
      CG(77)=1.0D 0
      CG(82)=-1.72D 0
      CG(83)=-0.33D 0
      CG(85)=1.62D 0
      CG(86)=1.24D 0
      CG(87)=0.21D 0
      CG(88)=-0.26D 0
      CG(94)=1.0D 0
      CG(97)=1.12D 0
      CG(100)=0.31D 0
      CG(103)=1.12D 0
      CG(105)=-0.36D 0
      CG(111)=1.0D 0
      CG(114)=0.45D 0
      CG(115)=0.26D 0
      CG(116)=-1.10D 0
      CG(117)=0.58D 0
      CG(119)=-1.03D 0
      CG(120)=0.10D 0
      CG(128)=1.0D 0
      CL(1)=2.5D 0
      CL(2)=1.1D 0
      CL(3)=-3.1D 0
      CL(4)=-3.5D 0
      CL(5)=1.3D 0
      CL(6)=2.1D 0
      CL(7)=2.3D 0
      CL(8)=-1.5D 0
      DO 154 J=1,8
      CU(J)=CL(J)
  154 CONTINUE
      RETURN
      END
* SUBROUTINE TFFU03                ALL SYSTEMS                92/11/01
C PORTABILITY : ALL SYSTEMS
C 92/11/01 LU : ORIGINAL V0ERSION
*
* PURPOSE :
*  VALUES OF MODEL FUNCTIONS FOR CONSTRAINED MINIMIZATION.
*  UNIVERSAL VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  F  VALUE OF THE MODEL FUNCTION.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TFFU03(N,X,F,NEXT)
      INTEGER N,NEXT
      DOUBLE PRECISION X(N),F
      DOUBLE PRECISION A,B,C,H,SOUC,SQ,HK,PI
      DATA PI /3.14159265358979323D 0/
      INTEGER I,J
      DOUBLE PRECISION CC(5,5),AA(16,16),EV(5),DV(5),Y(235),CV(10)
      COMMON /EMPR03/ Y
      GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150) NEXT
   10 A=PI
      F=SIN(A*X(1)/1.2D 1)*COS(A*X(2)/1.6D 1)
      RETURN
   20 F=((X(1)-3.0D 0)**2-9.0D 0)*X(2)**3/(2.7D 1*SQRT(3.0D 0))
      RETURN
   30 A=X(3)+3.0D-2
      B=A+X(2)
      C=B+X(1)
      F=-32.174D 0*(255.0D 0*LOG(C/(B+9.0D-2*X(1)))+
     10280.0D 0*LOG(B/(A+7.0D-2*X(2)))+
     10290.0D 0*LOG(A/(3.0D-2+1.3D-1*X(3))))
      RETURN
   40 F=-X(1)*X(2)*X(3)
      RETURN
   50 F=2.0D 0-X(1)*X(2)*X(3)
      RETURN
   60 F=(X(1)+1.0D 1*X(2))**2+5.0D 0*(X(3)-X(4))**2+
     &  (X(2)-2.0D 0*X(3))**4+1.0D 1*(X(1)-X(4))**4
      RETURN
   70 F=(X(1)-X(2))**2+(X(3)-1.0D 0)**2+(X(4)-1.0D 0)**4+
     &  (X(5)-1.0D 0)**6
      RETURN
   80 F=(X(1)-X(2))**2+(X(2)-X(3))**2+(X(3)-X(4))**4+(X(4)-X(5))**2
      RETURN
   90 F=(X(1)-X(2))**2+(X(2)+X(3)-2.0D 0)**2+
     &   (X(4)-1.0D 0)**2+(X(5)-1.0D 0)**2
      RETURN
  100 EV(1)=-15.0D 0
      EV(2)=-27.0D 0
      EV(3)=-36.0D 0
      EV(4)=-18.0D 0
      EV(5)=-12.0D 0
      CC(1,1)=30.0D 0
      CC(1,2)=-20.0D 0
      CC(1,3)=-10.0D 0
      CC(1,4)=32.0D 0
      CC(1,5)=-10.0D 0
      CC(2,1)=-20.0D 0
      CC(2,2)=39.0D 0
      CC(2,3)=-6.0D 0
      CC(2,4)=-31.0D 0
      CC(2,5)=32.0D 0
      CC(3,1)=-10.0D 0
      CC(3,2)=-6.0D 0
      CC(3,3)=10.0D 0
      CC(3,4)=-6.0D 0
      CC(3,5)=-10.0D 0
      CC(4,1)=32.0D 0
      CC(4,2)=-31.0D 0
      CC(4,3)=-6.0D 0
      CC(4,4)=39.0D 0
      CC(4,5)=-20.0D 0
      CC(5,1)=-10.0D 0
      CC(5,2)=32.0D 0
      CC(5,3)=-10.0D 0
      CC(5,4)=-20.0D 0
      CC(5,5)=30.0D 0
      DV(1)=4.0D 0
      DV(2)=8.0D 0
      DV(3)=10.0D 0
      DV(4)=6.0D 0
      DV(5)=2.0D 0
      F=0.0D 0
      DO 101 J=1,5
      DO 102 I=1,5
      F=F+CC(I,J)*X(I)*X(J)
  102 CONTINUE
      F=F+EV(J)*X(J)+DV(J)*X(J)**3
  101 CONTINUE
      RETURN
  110 HK=0.96D 0*4.9D 13
      H=((X(1)-1.0D 4)**2/6.4D 7+(X(1)-1.0D 4)*(X(2)-1.0D 0)/2.0D 4+
     &  (X(2)-1.0D 0)**2)*(X(3)-2.0D 6)**2/HK+
     &  (X(4)-10.0D 0)**2/2.5D 3+(X(5)-1.0D-3)**2/2.5D-3+
     &  (X(6)-1.0D 8)**2/2.5D 17
      F=-EXP(-H/2.0D 0)
      RETURN
  120 F=X(1)+2.0D 0*X(2)+4.0D 0*X(5)+EXP(X(1)*X(4))
      RETURN
  130 F=0.0D 0
      SQ=SQRT(2.0D 0*PI)
      DO 131 I=1,235
      A=X(1)/X(6)*EXP(-(Y(I)-X(3))**2/(2.0D 0*X(6)**2))
      B=X(2)/X(7)*EXP(-(Y(I)-X(4))**2/(2.0D 0*X(7)**2))
      C=(1.0D 0-X(2)-X(1))/X(8)*
     &     EXP(-(Y(I)-X(5))**2/(2.0D 0*X(8)**2))
      F=F-LOG((A+B+C)/SQ)
  131 CONTINUE
      RETURN
  140 F=0.0D 0
      CV(1)=-6.089D 0
      CV(2)=-17.164D 0
      CV(3)=-34.054D 0
      CV(4)=-5.914D 0
      CV(5)=-24.721D 0
      CV(6)=-14.986D 0
      CV(7)=-24.100D 0
      CV(8)=-10.708D 0
      CV(9)=-26.662D 0
      CV(10)=-22.179D 0
      SOUC=0.0D 0
      DO 141 J=1,10
      SOUC=SOUC+X(J)
  141 CONTINUE
      DO 142 J=1,10
      F=F+X(J)*(CV(J)+LOG(X(J)/SOUC))
  142 CONTINUE
      RETURN
  150 F=0.0D 0
      DO 151 I=1,16
      DO 152 J=1,16
      AA(I,J)=0.0D 0
  152 CONTINUE
  151 CONTINUE
      DO 153 I=1,16
      AA(I,I)=1.0D 0
  153 CONTINUE
      AA(1,4)=1.0D 0
      AA(1,7)=1.0D 0
      AA(1,8)=1.0D 0
      AA(1,16)=1.0D 0
      AA(2,3)=1.0D 0
      AA(2,7)=1.0D 0
      AA(2,10)=1.0D 0
      AA(3,7)=1.0D 0
      AA(3,9)=1.0D 0
      AA(3,10)=1.0D 0
      AA(3,14)=1.0D 0
      AA(4,7)=1.0D 0
      AA(4,11)=1.0D 0
      AA(4,15)=1.0D 0
      AA(5,6)=1.0D 0
      AA(5,10)=1.0D 0
      AA(5,12)=1.0D 0
      AA(5,16)=1.0D 0
      AA(6,8)=1.0D 0
      AA(6,15)=1.0D 0
      AA(7,11)=1.0D 0
      AA(7,13)=1.0D 0
      AA(8,10)=1.0D 0
      AA(8,15)=1.0D 0
      AA(9,12)=1.0D 0
      AA(9,16)=1.0D 0
      AA(10,14)=1.0D 0
      AA(11,13)=1.0D 0
      AA(12,14)=1.0D 0
      AA(13,14)=1.0D 0
      DO 154 I=1,16
      DO 155 J=1,16
      F=F+AA(I,J)*(X(I)**2+X(I)+1.0D 0)*(X(J)**2+X(J)+1.0D 0)
  155 CONTINUE
  154 CONTINUE
      RETURN
      END
* SUBROUTINE TFGU03                ALL SYSTEMS                92/11/01
C PORTABILITY : ALL SYSTEMS
C 92/11/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  GRADIENTS OF MODEL FUNCTIONS FOR CONSTRAINED MINIMIZATION.
*  UNIVERSAL VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  G(N)  GRADIENT OF THE MODEL FUNCTION.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TFGU03(N,X,G,NEXT)
      INTEGER N,NEXT
      DOUBLE PRECISION X(N),G(N)
      DOUBLE PRECISION A,B,C,D,E,F,SOU,AEXP,BEXP,CEXP,HI,SQ,H,HK,PI
      DATA PI /3.14159265358979323D 0/
      INTEGER I,J
      DOUBLE PRECISION CC(5,5),EV(5),DV(5),CV(10),Y(235)
      DOUBLE PRECISION AA(16,16)
      COMMON /EMPR03/ Y
      GO TO (10,20,30,40,40,60,70,80,90,100,110,120,130,140,150) NEXT
   10 A=PI
      B=A/1.2D 1
      C=A/1.6D 1
      G(1)=B*COS(B*X(1))*COS(C*X(2))
      G(2)=-C*SIN(B*X(1))*SIN(C*X(2))
      RETURN
   20 A=2.7D 1*SQRT(3.0D 0)
      G(1)=2.0D 0*(X(1)-3.0D 0)*X(2)**3/A
      G(2)=((X(1)-3.0D 0)**2-9.0D 0)*3.0D 0*X(2)**2/A
      RETURN
   30 A=X(3)+3.0D-2
      B=A+X(2)
      C=B+X(1)
      D=1.3D-1*X(3)+3.0D-2
      E=A+7.0D-2*X(2)
      F=B+9.0D-2*X(1)
      G(1)=-32.174D 0*255.0D 0*(F-9.0D-2*C)/(F*C)
      G(2)=-32.174D 0*(255.0D 0*(F-C)/(F*C)+
     10280.0D 0*(E-7.0D-2*B)/(E*B))
      G(3)=-32.174D 0*(255.0D 0*(F-C)/(F*C)+
     &  280.0D 0*(E-B)/(E*B)+290.0D 0*(D-1.3D-1*A)/(D*A))
      RETURN
   40 G(1)=-X(2)*X(3)
      G(2)=-X(1)*X(3)
      G(3)=-X(1)*X(2)
      RETURN
   60 A=X(1)+1.0D 1*X(2)
      B=X(3)-X(4)
      C=(X(2)-2.0D 0*X(3))**3
      D=(X(1)-X(4))**3
      G(1)=2.0D 0*A+4.0D 1*D
      G(2)=2.0D 1*A+4.0D 0*C
      G(3)=1.0D 1*B-8.0D 0*C
      G(4)=-1.0D 1*B-4.0D 1*D
      RETURN
   70 A=2.0D 0*(X(1)-X(2))
      B=2.0D 0*(X(3)-1.0D 0)
      C=4.0D 0*(X(4)-1.0D 0)**3
      D=6.0D 0*(X(5)-1.0D 0)**5
      G(1)=A
      G(2)=-A
      G(3)=B
      G(4)=C
      G(5)=D
      RETURN
   80 A=2.0D 0*(X(1)-X(2))
      B=2.0D 0*(X(2)-X(3))
      C=4.0D 0*(X(3)-X(4))**3
      D=2.0D 0*(X(4)-X(5))
      G(1)=A
      G(2)=B-A
      G(3)=C-B
      G(4)=D-C
      G(5)=-D
      RETURN
   90 A=2.0D 0*(X(1)-X(2))
      B=2.0D 0*(X(2)+X(3)-2.0D 0)
      G(1)=A
      G(2)=-A+B
      G(3)=B
      G(4)=2.0D 0*(X(4)-1.0D 0)
      G(5)=2.0D 0*(X(5)-1.0D 0)
      RETURN
  100 EV(1)=-15.0D 0
      EV(2)=-27.0D 0
      EV(3)=-36.0D 0
      EV(4)=-18.0D 0
      EV(5)=-12.0D 0
      CC(1,1)=30.0D 0
      CC(1,2)=-20.0D 0
      CC(1,3)=-10.0D 0
      CC(1,4)=32.0D 0
      CC(1,5)=-10.0D 0
      CC(2,1)=-20.0D 0
      CC(2,2)=39.0D 0
      CC(2,3)=-6.0D 0
      CC(2,4)=-31.0D 0
      CC(2,5)=32.0D 0
      CC(3,1)=-10.0D 0
      CC(3,2)=-6.0D 0
      CC(3,3)=10.0D 0
      CC(3,4)=-6.0D 0
      CC(3,5)=-10.0D 0
      CC(4,1)=32.0D 0
      CC(4,2)=-31.0D 0
      CC(4,3)=-6.0D 0
      CC(4,4)=39.0D 0
      CC(4,5)=-20.0D 0
      CC(5,1)=-10.0D 0
      CC(5,2)=32.0D 0
      CC(5,3)=-10.0D 0
      CC(5,4)=-20.0D 0
      CC(5,5)=30.0D 0
      DV(1)=4.0D 0
      DV(2)=8.0D 0
      DV(3)=10.0D 0
      DV(4)=6.0D 0
      DV(5)=2.0D 0
      DO 105 J=1,5
      G(J)=0.0D 0
  105 CONTINUE
      DO 106 J=1,5
      DO 107 I=1,5
      G(J)=G(J)+(CC(I,J)+CC(J,I))*X(I)
  107 CONTINUE
      G(J)=G(J)+EV(J)+3.0D 0*DV(J)*X(J)**2
  106 CONTINUE
      RETURN
  110 HK=0.96D 0*4.9D 13
      H=((X(1)-1.0D 4)**2/6.4D 7+(X(1)-1.0D 4)*(X(2)-1.0D 0)/2.0D 4+
     &  (X(2)-1.0D 0)**2)*(X(3)-2.0D 6)**2/HK+
     &  (X(4)-10.0D 0)**2/2.5D 3+(X(5)-1.0D-3)**2/2.5D-3+
     &  (X(6)-1.0D 8)**2/2.5D 17
      D=EXP(-H/2.0D 0)/2.0D 0
      G(1)=D*(2.0D 0*(X(1)-1.0D 4)/6.4D 7+(X(2)-1.0D 0)/2.0D 4)*
     &     (X(3)-2.0D 6)**2/HK
      G(2)=D*((X(1)-1.0D 4)/2.0D 4+
     &     2.0D 0*(X(2)-1.0D 0))*(X(3)-2.0D 6)**2/HK
      G(3)=D*((X(1)-1.0D 4)**2/6.4D 7+
     &     (X(1)-1.0D 4)*(X(2)-1.0D 0)/2.0D 4+
     &     (X(2)-1.0D 0)**2)*2.0D 0*(X(3)-2.0D 6)/HK
      G(4)=D*2.0D 0*(X(4)-10.0D 0)/2.5D 3
      G(5)=D*2.0D 0*(X(5)-1.0D-3)/2.5D-3
      G(6)=D*2.0D 0*(X(6)-1.0D 8)/2.5D 17
      RETURN
  120 G(1)=1.0D 0+X(4)*EXP(X(1)*X(4))
      G(2)=2.0D 0
      G(3)=0.0D 0
      G(4)=X(1)*EXP(X(1)*X(4))
      G(5)=4.0D 0
      G(6)=0.0D 0
      RETURN
  130 DO 131 J=1,8
      G(J)=0.0D 0
  131 CONTINUE
      SQ=SQRT(2.0D 0*PI)
      DO 132 I=1,235
      AEXP=EXP(-(Y(I)-X(3))**2/(2.0D 0*X(6)**2))
      A=X(1)/X(6)*AEXP
      BEXP=EXP(-(Y(I)-X(4))**2/(2.0D 0*X(7)**2))
      B=X(2)/X(7)*BEXP
      CEXP=EXP(-(Y(I)-X(5))**2/(2.0D 0*X(8)**2))
      C=(1.0D 0-X(2)-X(1))/X(8)*CEXP
      HI=-SQ/(A+B+C)
      G(1)=G(1)+HI*(AEXP/X(6)-CEXP/X(8))
      G(2)=G(2)+HI*(BEXP/X(7)-CEXP/X(8))
      G(3)=G(3)+HI*X(1)/X(6)*AEXP*(Y(I)-X(3))/(X(6)**2)
      G(4)=G(4)+HI*X(2)/X(7)*BEXP*(Y(I)-X(4))/(X(7)**2)
      G(5)=G(5)+HI*(1.0D 0-X(2)-X(1))/X(8)*CEXP*(Y(I)-X(5))/(X(8)**2)
      G(6)=G(6)+HI*((-X(1)/(X(6)**2))*AEXP+
     &   X(1)/X(6)*AEXP*(Y(I)-X(3))**2/(X(6)**3))
      G(7)=G(7)+HI*((-X(2)/(X(7)**2))*BEXP+
     &   X(2)/X(7)*BEXP*(Y(I)-X(4))**2/(X(7)**3))
      G(8)=G(8)+HI*((X(1)+X(2)-1.0D 0)/(X(8)**2)*CEXP+
     &   (1.0D 0-X(2)-X(1))/X(8)*CEXP*(Y(I)-X(5))**2/(X(8)**3))
  132 CONTINUE
      RETURN
  140 CV(1)=-6.089D 0
      CV(2)=-17.164D 0
      CV(3)=-34.054D 0
      CV(4)=-5.914D 0
      CV(5)=-24.721D 0
      CV(6)=-14.986D 0
      CV(7)=-24.100D 0
      CV(8)=-10.708D 0
      CV(9)=-26.662D 0
      CV(10)=-22.179D 0
      SOU=0.0D 0
      DO 142 J=1,10
      SOU=SOU+X(J)
  142 CONTINUE
      DO 143 J=1,10
      G(J)=CV(J)+LOG(X(J)/SOU)
  143 CONTINUE
      RETURN
  150 DO 151 I=1,16
      DO 152 J=1,16
      AA(I,J)=0.0D 0
  152 CONTINUE
  151 CONTINUE
      DO 153 I=1,16
      AA(I,I)=1.0D 0
  153 CONTINUE
      AA(1,4)=1.0D 0
      AA(1,7)=1.0D 0
      AA(1,8)=1.0D 0
      AA(1,16)=1.0D 0
      AA(2,3)=1.0D 0
      AA(2,7)=1.0D 0
      AA(2,10)=1.0D 0
      AA(3,7)=1.0D 0
      AA(3,9)=1.0D 0
      AA(3,10)=1.0D 0
      AA(3,14)=1.0D 0
      AA(4,7)=1.0D 0
      AA(4,11)=1.0D 0
      AA(4,15)=1.0D 0
      AA(5,6)=1.0D 0
      AA(5,10)=1.0D 0
      AA(5,12)=1.0D 0
      AA(5,16)=1.0D 0
      AA(6,8)=1.0D 0
      AA(6,15)=1.0D 0
      AA(7,11)=1.0D 0
      AA(7,13)=1.0D 0
      AA(8,10)=1.0D 0
      AA(8,15)=1.0D 0
      AA(9,12)=1.0D 0
      AA(9,16)=1.0D 0
      AA(10,14)=1.0D 0
      AA(11,13)=1.0D 0
      AA(12,14)=1.0D 0
      AA(13,14)=1.0D 0
      DO 154 J=1,16
      G(J)=0.0D 0
  154 CONTINUE
      DO 155 J=1,16
      DO 156 I=1,16
      G(J)=G(J)+
     &  (AA(I,J)+AA(J,I))*(X(I)**2+X(I)+1.0D 0)*(2.0D 0*X(J)+1.0D 0)
  156 CONTINUE
  155 CONTINUE
      RETURN
      END
*
*     TEST SUBROUTINES FOR CONSTRAINED OPTIMIZATION
*
* SUBROUTINE TIND07             ALL SYSTEMS                 90/12/01
C PORTABILITY : ALL SYSTEMS
C 90/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  INITIAL VALUES OF THE VARIABLES FOR NONLINEAR PROGRAMMING.
*  DENSE VERSION.
*
* PARAMETERS :
*  IO  N  NUMBER OF VARIABLES.
*  IO  NC  NUMBER OF CONSTRAINTS.
*  IO  NCL NUMBER OF LINEAR CONSTRAINTS.
*  RO  X(N)  VECTOR OF VARIABLES.
*  IO  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  RO  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
*  RO  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
*  IO  IC(NC)  TYPES OF CONSTRAINTS.
*  RO  CL(NC)  LOWER BOUNDS OF CONSTRAINTS.
*  RO  CU(NC)  UPPER BOUNDS OF CONSTRAINTS.
*  RO  FMIN  LOWER BOUND FOR VALUE OF THE OBJECTIVE FUNCTION.
*  RO  XMAX  MAXIMUM STEPSIZE.
*  IO  NEXT  NUMBER OF THE TEST PROBLEM.
*  IO  IERR  ERROR INDICATOR.
*
      SUBROUTINE TIND07(N,NC,X,IX,XL,XU,IC,CL,CU,FMIN,XMAX,NEXT,
     * IERR)
      INTEGER I,N,NC,NCL,IX(N),IC(NC),NEXT,IERR
      DOUBLE PRECISION X(N),XL(N),XU(N),CL(NC),CU(NC),FMIN,XMAX
      DOUBLE PRECISION Y(128)
      COMMON /EMPR07/ Y
      NCL=0
      FMIN=-1.0D 60
      XMAX=1.0D 3
      IERR=0
      DO 1 I=1,N
      IX(I)=0
      XL(I)=0.0D 0
      XU(I)=0.0D 0
    1 CONTINUE
      DO 2 I=1,NC
      IC(I)=1
      CL(I)=0.0D 0
      CU(I)=0.0D 0
    2 CONTINUE
      GO TO(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,
     * 180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,
     * 340) NEXT
   10 IF(N.GE.2.AND.NC.GE.1) THEN
      N=2
      NC=1
      X(1)=2.0D 0
      X(2)=0.0D 0
      IC(1)=2
      CU(1)=1.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
   20 IF(N.GE.2.AND.NC.GE.1) THEN
      N=2
      NC=1
      X(1)=-10.0D 0
      X(2)= 10.0D 0
      IC(1)=2
      CU(1)=1.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
   30 IF(N.GE.2.AND.NC.GE.1) THEN
      N=2
      NC=1
      X(1)=0.0D 0
      X(2)=0.0D 0
      IC(1)=2
      CU(1)=25.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
   40 IF(N.GE.2.AND.NC.GE.1) THEN
      N=2
      NC=1
      X(1)=-2.0D 0
      X(2)=-2.0D 0
      IX(1)=1
      IX(2)=1
      ELSE
      IERR=1
      ENDIF
      RETURN
   50 IF(N.GE.2.AND.NC.GE.1) THEN
      N=2
      NC=1
      X(1)=4.2D-1
      X(2)=5.0D 0
      IX(1)=1
      IX(2)=1
      XL(1)= 0.4D 0
      XL(2)=-4.0D 0
      CL(1)=0.09D 0
      Y(1)=8.0D 0
      Y(2)=8.0D 0
      Y(3)=10.0D 0
      Y(4)=10.0D 0
      Y(5)=10.0D 0
      Y(6)=10.0D 0
      Y(7)=12.0D 0
      Y(8)=12.0D 0
      Y(9)=12.0D 0
      Y(10)=12.0D 0
      Y(11)=14.0D 0
      Y(12)=14.0D 0
      Y(13)=14.0D 0
      Y(14)=16.0D 0
      Y(15)=16.0D 0
      Y(16)=16.0D 0
      Y(17)=18.0D 0
      Y(18)=18.0D 0
      Y(19)=20.0D 0
      Y(20)=20.0D 0
      Y(21)=20.0D 0
      Y(22)=22.0D 0
      Y(23)=22.0D 0
      Y(24)=22.0D 0
      Y(25)=24.0D 0
      Y(26)=24.0D 0
      Y(27)=24.0D 0
      Y(28)=26.0D 0
      Y(29)=26.0D 0
      Y(30)=26.0D 0
      Y(31)=28.0D 0
      Y(32)=28.0D 0
      Y(33)=30.0D 0
      Y(34)=30.0D 0
      Y(35)=30.0D 0
      Y(36)=32.0D 0
      Y(37)=32.0D 0
      Y(38)=34.0D 0
      Y(39)=36.0D 0
      Y(40)=36.0D 0
      Y(41)=38.0D 0
      Y(42)=38.0D 0
      Y(43)=40.0D 0
      Y(44)=42.0D 0
      Y(45)=0.49D 0
      Y(46)=0.49D 0
      Y(47)=0.48D 0
      Y(48)=0.47D 0
      Y(49)=0.48D 0
      Y(50)=0.47D 0
      Y(51)=0.46D 0
      Y(52)=0.46D 0
      Y(53)=0.45D 0
      Y(54)=0.43D 0
      Y(55)=0.45D 0
      Y(56)=0.43D 0
      Y(57)=0.43D 0
      Y(58)=0.44D 0
      Y(59)=0.43D 0
      Y(60)=0.43D 0
      Y(61)=0.46D 0
      Y(62)=0.45D 0
      Y(63)=0.42D 0
      Y(64)=0.42D 0
      Y(65)=0.43D 0
      Y(66)=0.41D 0
      Y(67)=0.41D 0
      Y(68)=0.40D 0
      Y(69)=0.42D 0
      Y(70)=0.40D 0
      Y(71)=0.40D 0
      Y(72)=0.41D 0
      Y(73)=0.40D 0
      Y(74)=0.41D 0
      Y(75)=0.41D 0
      Y(76)=0.40D 0
      Y(77)=0.40D 0
      Y(78)=0.40D 0
      Y(79)=0.38D 0
      Y(80)=0.41D 0
      Y(81)=0.40D 0
      Y(82)=0.40D 0
      Y(83)=0.41D 0
      Y(84)=0.38D 0
      Y(85)=0.40D 0
      Y(86)=0.40D 0
      Y(87)=0.39D 0
      Y(88)=0.39D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
   60 IF(N.GE.2.AND.NC.GE.2) THEN
      N=2
      NC=2
      X(1)=-2.0D 0
      X(2)= 1.0D 0
      IX(1)=3
      IX(2)=2
      XL(1)=-0.5D 0
      XU(1)= 0.5D 0
      XU(2)= 1.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
   70 IF(N.GE.2.AND.NC.GE.2) THEN
      N=2
      NC=2
      X(1)=20.1D 0
      X(2)=5.84D 0
      IX(1)=3
      IX(2)=3
      XL(1)=1.3D 1
      XU(1)=1.0D 2
      XL(2)=0.0D 0
      XU(2)=1.0D 2
      IC(2)=2
      CL(1)=1.0D 2
      CU(2)=82.81D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
   80 IF(N.GE.2.AND.NC.GE.3) THEN
      N=2
      NC=3
      X(1)=-2.0D 0
      X(2)= 1.0D 0
      IX(1)=3
      XL(1)=-0.5D 0
      XU(1)= 0.5D 0
      CL(3)= 1.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
   90 IF(N.GE.2.AND.NC.GE.3) THEN
      N=2
      NC=3
      X(1)=90.0D 0
      X(2)=10.0D 0
      IX(1)=3
      IX(2)=3
      XL(1)=0.0D 0
      XU(1)=75.0D 0
      XL(2)=0.0D 0
      XU(2)=65.0D 0
      CL(1)=700.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  100 IF(N.GE.3.AND.NC.GE.1) THEN
      N=3
      NC=1
      X(1)=1.0D 0
      X(2)=1.0D 0
      X(3)=1.0D 0
      IC(1)=2
      CU(1)=48.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  110 IF(N.GE.3.AND.NC.GE.1) THEN
      N=3
      NC=1
      X(1)=1.0D 0
      X(2)=1.0D 0
      X(3)=1.0D 0
      IX(1)=3
      IX(2)=3
      IX(3)=3
      XL(1)= 1.0D 0
      XU(1)= 1.0D 1
      XL(2)=-1.0D 1
      XU(2)= 1.0D 1
      XL(3)=-1.0D 1
      XU(3)= 1.0D 1
      CL(1)= 1.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  120 IF(N.GE.3.AND.NC.GE.1) THEN
      N=3
      NC=1
      X(1)=1.0D 0
      X(2)=1.0D 0
      X(3)=1.0D 0
      IX(1)=3
      IX(2)=3
      IX(3)=3
      XL(1)=-1.0D 1
      XU(1)= 1.0D 1
      XL(2)= 1.0D 0
      XU(2)= 1.0D 1
      XL(3)=-1.0D 1
      XU(3)= 1.0D 0
      CL(1)= 1.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  130 IF(N.GE.3.AND.NC.GE.1) THEN
      N=3
      NC=1
      X(1)=1.0D 0
      X(2)=1.0D 0
      X(3)=1.0D 0
      IX(1)=1
      IX(2)=1
      IX(3)=1
      XL(1)=1.0D-5
      XL(2)=1.0D-5
      XL(3)=1.0D-5
      IC(1)=2
      CU(1)=1.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  140 IF(N.GE.3.AND.NC.GE.2) THEN
      N=3
      NC=2
      X(1)=0.0D 0
      X(2)=1.05D 0
      X(3)=2.9D 0
      IX(1)=3
      IX(2)=3
      IX(3)=3
      XL(1)=0.0D 0
      XU(1)=1.0D 2
      XL(2)=0.0D 0
      XU(2)=1.0D 2
      XL(3)=0.0D 0
      XU(3)=1.0D 1
      IC(1)=2
      IC(2)=2
      ELSE
      IERR=1
      ENDIF
      RETURN
  150 IF(N.GE.3.AND.NC.GE.7) THEN
      N=3
      NC=7
      X(1)=1745.0D 0
      X(2)=12000.0D 0
      X(3)=110.0D 0
      IX(1)=3
      IX(2)=3
      IX(3)=3
      XL(1)=1.0D-5
      XU(1)=2.0D 3
      XL(2)=1.0D-5
      XU(2)=1.6D 4
      XL(3)=1.0D-5
      XU(3)=1.2D 2
      DO 151 I=1,7
      IC(I)=3
  151 CONTINUE
      CL(1)=0.0D 0
      CL(2)=0.0D 0
      CL(3)=8.5D 1
      CL(4)=9.0D 1
      CL(5)=3.0D 0
      CL(6)=1.0D-2
      CL(7)=1.45D 2
      CU(1)=5.0D 3
      CU(2)=2.0D 3
      CU(3)=9.3D 1
      CU(4)=9.5D 1
      CU(5)=1.2D 1
      CU(6)=4.0D 0
      CU(7)=1.62D 2
      ELSE
      IERR=1
      ENDIF
      RETURN
  160 IF(N.GE.4.AND.NC.GE.1) THEN
      N=4
      NC=1
      X(1)=2.0D 0
      X(2)=4.0D 0
      X(3)=4.0D-2
      X(4)=2.0D 0
      DO 161 I=1,4
      IX(I)=3
      XL(I)=1.0D-5
      XU(I)=1.0D 2
  161 CONTINUE
      XU(3)=1.0D 0
      DO 162 I=2,19
      Y(I)=DBLE(I-1)
      Y(I+38)=LOG(Y(I))
  162 CONTINUE
      Y(1)=0.1D 0
      Y(1+38)=LOG(Y(1))
      Y(20)=0.00189D 0
      Y(21)=0.1038D 0
      Y(22)=0.268D 0
      Y(23)=0.506D 0
      Y(24)=0.577D 0
      Y(25)=0.604D 0
      Y(26)=0.725D 0
      Y(27)=0.898D 0
      Y(28)=0.947D 0
      Y(29)=0.845D 0
      Y(30)=0.702D 0
      Y(31)=0.528D 0
      Y(32)=0.385D 0
      Y(33)=0.257D 0
      Y(34)=0.159D 0
      Y(35)=0.0869D 0
      Y(36)=0.0453D 0
      Y(37)=0.01509D 0
      Y(38)=0.00189D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  170 IF(N.GE.4.AND.NC.GE.2) THEN
      N=4
      NC=2
      X(1)=1.0D 0
      X(2)=1.0D 0
      X(3)=1.0D 0
      X(4)=1.0D 0
      DO 171 I=1,4
      IX(I)=3
      XL(I)=1.0D-3
      XU(I)=DBLE(5-I)*1.0D 5
  171 CONTINUE
      IC(1)=2
      IC(2)=2
      CU(1)=0.0401D 0
      CU(2)=0.010085D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  180 IF(N.GE.4.AND.NC.GE.3) THEN
      N=4
      NC=3
      DO 181 I=1,N
      X(I)=0.0D 0
  181 CONTINUE
      IC(1)=2
      IC(2)=2
      IC(3)=2
      CU(1)=8.0D 0
      CU(2)=1.0D 1
      CU(3)=5.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  190 IF(N.GE.5.AND.NC.GE.3) THEN
      N=5
      NC=3
      X(1)=1.0D 0
      X(2)=1.0D 0
      X(3)=1.0D 0
      X(4)=1.0D 0
      X(5)=1.0D 0
      IC(1)=2
      CU(1)= 2.0D 1
      CL(2)=-2.0D 0
      CL(3)= 5.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  200 IF(N.GE.5.AND.NC.GE.3) THEN
      N=5
      NC=3
      X(1)=78.0D 0
      X(2)=33.0D 0
      X(3)=27.0D 0
      X(4)=27.0D 0
      X(5)=27.0D 0
      DO 201 I=1,5
      IX(I)=3
      XL(I)=27.0D 0
      XU(I)=45.0D 0
  201 CONTINUE
      XL(1)=78.0D 0
      XL(2)=33.0D 0
      XU(1)=10.2D 1
      IC(1)=3
      IC(2)=3
      IC(3)=3
      CL(1)=0.0D 0
      CU(1)=9.2D 1
      CL(2)=9.0D 1
      CU(2)=1.1D 2
      CL(3)=2.0D 1
      CU(3)=2.5D 1
      ELSE
      IERR=1
      ENDIF
      RETURN
  210 IF(N.GE.5.AND.NC.GE.3) THEN
      N=5
      NC=3
      X(1)=2.52D 0
      X(2)=2.0D 0
      X(3)=37.5D 0
      X(4)=9.25D 0
      X(5)=6.8D 0
      DO 211 I=1,5
      IX(I)=3
  211 CONTINUE
      XL(1)=0.0D 0
      XU(1)=1.0D 3
      XL(2)=1.2D 0
      XU(2)=2.4D 0
      XL(3)=2.0D 1
      XU(3)=6.0D 1
      XL(4)=9.0D 0
      XU(4)=9.3D 0
      XL(5)=6.5D 0
      XU(5)=7.0D 0
      IC(1)=3
      IC(2)=3
      IC(3)=3
      CL(1)=0.0D 0
      CU(1)=294.0D 3
      CL(2)=0.0D 0
      CU(2)=294.0D 3
      CL(3)=0.0D 0
      CU(3)=277.2D 3
      ELSE
      IERR=1
      ENDIF
      RETURN
  220 IF(N.GE.5.AND.NC.GE.21) THEN
      N=5
      NC=21
      X(1)=900.0D 0
      X(2)=80.0D 0
      X(3)=115.0D 0
      X(4)=267.0D 0
      X(5)=27.0D 0
      DO 221 I=1,5
      IX(I)=3
  221 CONTINUE
      XL(1)=704.4148D 0
      XU(1)=906.3855D 0
      XL(2)=68.6D 0
      XU(2)=288.88D 0
      XL(3)=0.0D 0
      XU(3)=134.75D 0
      XL(4)=193.0D 0
      XU(4)=287.0966D 0
      XL(5)=25.0D 0
      XU(5)=84.1988D 0
      IC(3)=2
      CU(3)=21.0D 0
      CL(4)=110.6D 0
      DO 222 I=5,21
      IC(I)=3
  222 CONTINUE
      CL(5)= 213.1D 0
      CL(6)= 17.505D 0
      CL(7)= 11.275D 0
      CL(8)= 214.228D 0
      CL(9)= 7.458D 0
      CL(10)= 0.961D 0
      CL(11)= 1.612D 0
      CL(12)= 0.146D 0
      CL(13)= 107.99D 0
      CL(14)= 922.693D 0
      CL(15)= 926.832D 0
      CL(16)= 18.766D 0
      CL(17)= 1072.163D 0
      CL(18)= 8961.448D 0
      CL(19)= 0.063D 0
      CL(20)= 71084.33D 0
      CL(21)= 2.802713D 6
      CU(5)= 405.23D 0
      CU(6)= 1053.6667D 0
      CU(7)= 35.03D 0
      CU(8)= 665.585D 0
      CU(9)= 584.463D 0
      CU(10)= 265.916D 0
      CU(11)= 7.046D 0
      CU(12)= 0.222D 0
      CU(13)= 273.366D 0
      CU(14)= 1286.105D 0
      CU(15)= 1444.046D 0
      CU(16)= 537.141D 0
      CU(17)= 3247.039D 0
      CU(18)= 26844.086D 0
      CU(19)= 0.386D 0
      CU(20)= 1.4D 5
      CU(21)= 1.2146108D 7
      ELSE
      IERR=1
      ENDIF
      RETURN
  230 IF(N.GE.6.AND.NC.GE.2) THEN
      N=6
      NC=2
      X(1)=1.0D 0
      X(2)=1.0D 0
      X(3)=1.0D 0
      X(4)=3.0D 0
      X(5)=0.0D 0
      X(6)=0.5D 0
      IX(6)=3
      XL(6)=4.0D 0
      XU(6)=8.0D 0
      IC(1)=2
      IC(2)=2
      CU(1)=5.0D 0
      CU(2)=1.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  240 IF(N.GE.6.AND.NC.GE.2) THEN
      N=6
      NC=2
      X(1)=5.54D 0
      X(2)=4.4D 0
      X(3)=12.02D 0
      X(4)=11.82D 0
      X(5)=0.702D 0
      X(6)=0.852D 0
      DO 241 I=1,6
      IX(I)=1
  241 CONTINUE
      IC(2)=2
      CL(1)=2.07D 0
      CU(2)=1.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  250 IF(N.GE.6.AND.NC.GE.4) THEN
      N=6
      NC=4
      DO 251 I=1,6
      X(I)=0.0D 0
      IX(I)=3
  251 CONTINUE
      XU(1)=0.31D 0
      XU(2)=0.046D 0
      XU(3)=0.068D 0
      XU(4)=0.042D 0
      XU(5)=0.028D 0
      XU(6)=0.134D-1
      CL(1)=32.97D 0
      CL(2)=25.12D 0
      CL(3)=-124.08D 0
      CL(4)=-173.02D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  260 IF(N.GE.7.AND.NC.GE.4) THEN
      N=7
      NC=4
      X(1)=1.0D 0
      X(2)=2.0D 0
      X(3)=0.0D 0
      X(4)=4.0D 0
      X(5)=0.0D 0
      X(6)=1.0D 0
      X(7)=1.0D 0
      IC(1)=2
      IC(2)=2
      IC(3)=2
      IC(4)=2
      CU(1)=127.0D 0
      CU(2)=282.0D 0
      CU(3)=196.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  270 IF(N.GE.7.AND.NC.GE.5) THEN
      N=7
      NC=5
      DO 271 I=1,7
      X(I)=6.0D 0
      IX(I)=3
      XL(I)=1.0D-1
      XU(I)=1.0D 1
  271 CONTINUE
      XL(7)=1.0D-2
      DO 272 I=1,4
      IC(I)=2
      CU(I)=1.0D 0
  272 CONTINUE
      IC(5)=3
      CL(5)=1.0D 2
      CU(5)=3.0D 3
      ELSE
      IERR=1
      ENDIF
      RETURN
  280 IF(N.GE.8.AND.NC.GE.5) THEN
      N=8
      NC=5
      X(1)=6.0D 0
      X(2)=3.0D 0
      X(3)=0.4D 0
      X(4)=0.2D 0
      X(5)=6.0D 0
      X(6)=6.0D 0
      X(7)=1.0D 0
      X(8)=0.5D 0
      DO 281 I=1,8
      IX(I)=3
      XL(I)=1.0D-1
      XU(I)=1.0D 1
  281 CONTINUE
      DO 282 I=1,4
      IC(I)=2
      CU(I)=1.0D 0
  282 CONTINUE
      IC(5)=3
      CL(5)=1.0D 0
      CU(5)=4.2D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  290 IF(N.GE.8.AND.NC.GE.6) THEN
      N=8
      NC=6
      DO 291 I=1,8
      X(I)=5.0D 3
      IX(I)=3
      XL(I)=1.0D 1
      XU(I)=1.0D 3
  291 CONTINUE
      X(4)=200.0D 0
      X(5)=350.0D 0
      X(6)=150.0D 0
      X(7)=225.0D 0
      X(8)=425.0D 0
      XL(1)=1.0D 2
      XL(2)=1.0D 3
      XL(3)=1.0D 3
      XU(1)=1.0D 4
      XU(2)=1.0D 4
      XU(3)=1.0D 4
      DO 292 I=1,3
      IC(I)=2
      CU(I)=1.0D 0
  292 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  300 IF(N.GE.9.AND.NC.GE.13) THEN
      N=9
      NC=13
      DO 301 I=1,9
      X(I)=1.0D 0
      IX(I)=1
      XL(I)=0.0D 0
      IC(I)=2
      CU(I)=1.0D 0
  301 CONTINUE
      IX(9)=1
      XL(9)=0.0D 0
      IC(13)=2
      ELSE
      IERR=1
      ENDIF
      RETURN
  310 IF(N.GE.10.AND.NC.GE.8) THEN
      N=10
      NC=8
      X(1)=2.0D 0
      X(2)=3.0D 0
      X(3)=5.0D 0
      X(4)=5.0D 0
      X(5)=1.0D 0
      X(6)=2.0D 0
      X(7)=7.0D 0
      X(8)=3.0D 0
      X(9)=6.0D 0
      X(10)=1.0D 1
      DO 311 I=1,8
      IC(I)=2
  311 CONTINUE
      CU(1)=120.0D 0
      CU(2)=40.0D 0
      CU(3)=30.0D 0
      CU(4)=0.0D 0
      CU(5)=105.0D 0
      CU(6)=0.0D 0
      CU(7)=0.0D 0
      CU(8)=12.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  320 IF(N.GE.15.AND.NC.GE.5) THEN
      N=15
      NC=5
      DO 321 I=1,15
      X(I)=0.001D 0
      IX(I)=1
      XL(I)=0.0D 0
  321 CONTINUE
      X(7)=60.0D 0
      CL(1)=15.0D 0
      CL(2)=27.0D 0
      CL(3)=36.0D 0
      CL(4)=18.0D 0
      CL(5)=12.0D 0
      DO 322 I=1,90
      Y(I)=0.0D 0
  322 CONTINUE
      Y(1)=-16.0D 0
      Y(3)=-3.5D 0
      Y(6)= 2.0D 0
      Y(7)=-1.0D 0
      Y(8)=-1.0D 0
      Y(9)= 1.0D 0
      Y(10)= 1.0D 0
      Y(11)= 2.0D 0
      Y(12)=-2.0D 0
      Y(14)=-2.0D 0
      Y(15)=-9.0D 0
      Y(17)=-1.0D 0
      Y(18)=-2.0D 0
      Y(19)= 2.0D 0
      Y(20)= 1.0D 0
      Y(23)= 2.0D 0
      Y(25)=-2.0D 0
      Y(26)=-4.0D 0
      Y(27)=-1.0D 0
      Y(28)=-3.0D 0
      Y(29)= 3.0D 0
      Y(30)= 1.0D 0
      Y(31)= 1.0D 0
      Y(32)= 4.0D 0
      Y(34)=-4.0D 0
      Y(35)= 1.0D 0
      Y(37)=-1.0D 0
      Y(38)=-2.0D 0
      Y(39)= 4.0D 0
      Y(40)= 1.0D 0
      Y(42)= 2.0D 0
      Y(44)=-1.0D 0
      Y(45)=-2.8D 0
      Y(47)=-1.0D 0
      Y(48)=-1.0D 0
      Y(49)= 5.0D 0
      Y(50)= 1.0D 0
      Y(51)=-4.0D 1
      Y(52)=-2.0D 0
      Y(53)=-2.5D-1
      Y(54)=-4.0D 0
      Y(55)=-4.0D 0
      Y(56)=-1.0D 0
      Y(57)=-4.0D 1
      Y(58)=-6.0D 1
      Y(59)= 5.0D 0
      Y(60)= 1.0D 0
      Y(61)= 3.0D 1
      Y(62)=-2.0D 1
      Y(63)=-1.0D 1
      Y(64)= 3.2D 1
      Y(65)=-1.0D 1
      Y(66)=-2.0D 1
      Y(67)= 3.9D 1
      Y(68)=-6.0D 0
      Y(69)=-3.1D 1
      Y(70)= 3.2D 1
      Y(71)=-1.0D 1
      Y(72)=-6.0D 0
      Y(73)= 1.0D 1
      Y(74)=-6.0D 0
      Y(75)=-1.0D 1
      Y(76)= 3.2D 1
      Y(77)=-3.1D 1
      Y(78)=-6.0D 0
      Y(79)= 3.9D 1
      Y(80)=-2.0D 1
      Y(81)=-1.0D 1
      Y(82)= 3.2D 1
      Y(83)=-1.0D 1
      Y(84)=-2.0D 1
      Y(85)= 3.0D 1
      Y(86)= 4.0D 0
      Y(87)= 8.0D 0
      Y(88)= 1.0D 1
      Y(89)= 6.0D 0
      Y(90)= 2.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
  330 IF(N.GE.16.AND.NC.GE.19) THEN
      N=16
      NC=19
      X(1)=0.8D 0
      X(2)=0.83D 0
      X(3)=0.85D 0
      X(4)=0.87D 0
      X(5)=0.9D 0
      X(6)=0.1D 0
      X(7)=0.12D 0
      X(8)=0.19D 0
      X(9)=0.25D 0
      X(10)=0.29D 0
      X(11)=512.0D 0
      X(12)=13.1D 0
      X(13)=71.8D 0
      X(14)=640.0D 0
      X(15)=650.0D 0
      X(16)=5.7D 0
      DO 331 I=1,16
      IX(I)=3
  331 CONTINUE
      DO 332 I=1,10
      XL(I)=0.1D 0
      XU(I)=0.9D 0
  332 CONTINUE
      XL(5)=0.9D 0
      XU(5)=1.0D 0
      XL(6)=1.0D-4
      XU(6)=0.1D 0
      XL(11)=1.0D 0
      XU(11)=1.0D 3
      XL(12)=1.0D-6
      XU(12)=5.0D 2
      XL(13)=1.0D 0
      XU(13)=5.0D 2
      DO 333 I=14,15
      XL(I)=5.0D 2
      XU(I)=1.0D 3
  333 CONTINUE
      XL(16)=1.0D-6
      XU(16)=5.0D 2
      DO 334 I=1,19
      IC(I)=2
      CU(I)=1.0D 0
  334 CONTINUE
      DO 335 I=1,5
      Y(I)=1.262626D 0
      Y(I+5)=-1.231060D 0
      Y(I+25)=1.0D 0
      Y(I+34)=0.002D 0
      Y(I+43)=1.D 0
      Y(I+51)=1.D 0
      Y(I+56)=1.D 0
  335 CONTINUE
      DO 336 I=11,23,3
      Y(I)=0.03475D 0
      Y(I+1)=0.975D 0
      Y(I+2)=-0.00975D 0
  336 CONTINUE
      Y(28)=-Y(28)
      Y(30)=0.002D 0
      Y(31)=Y(30)
      Y(32)=-Y(30)
      Y(33)=-Y(30)
      Y(34)=1.0D 0
      Y(37)=1.0D 0
      Y(38)=-Y(38)
      Y(39)=-Y(39)
      Y(40)=1.0D 0
      Y(41)=Y(40)
      Y(42)=5.0D 2
      Y(43)=-Y(42)
      Y(44)=-Y(44)
      Y(47)=Y(42)
      Y(48)=-Y(48)
      Y(49)=Y(43)
      Y(50)=0.9D 0
      Y(51)=0.002D 0
      Y(52)=-Y(51)
      Y(53)=Y(51)
      Y(54)=Y(52)
      ELSE
      IERR=1
      ENDIF
      RETURN
  340 IF(N.GE.20.AND.NC.GE.17) THEN
      N=20
      NC=17
      X(1)=2.0D 0
      X(2)=3.0D 0
      X(3)=5.0D 0
      X(4)=5.0D 0
      X(5)=1.0D 0
      X(6)=2.0D 0
      X(7)=7.0D 0
      X(8)=3.0D 0
      X(9)=6.0D 0
      X(10)=1.0D 1
      X(11)=2.0D 0
      X(12)=2.0D 0
      X(13)=6.0D 0
      X(14)=1.5D 1
      X(15)=1.0D 0
      X(16)=2.0D 0
      X(17)=1.0D 0
      X(18)=2.0D 0
      X(19)=1.0D 0
      X(20)=3.0D 0
      DO 341 I=1,17
      IC(I)=2
  341 CONTINUE
      CU(1)=120.0D 0
      CU(2)=40.0D 0
      CU(3)=30.0D 0
      CU(4)=0.0D 0
      CU(5)=105.0D 0
      CU(6)=0.0D 0
      CU(7)=0.0D 0
      CU(8)=12.0D 0
      CU(9)=0.0D 0
      CU(10)=28.0D 0
      CU(11)=87.0D 0
      CU(12)=10.0D 0
      CU(13)=92.0D 0
      CU(14)=54.0D 0
      CU(15)=68.0D 0
      CU(16)=-19.0D 0
      CU(17)=0.0D 0
      ELSE
      IERR=1
      ENDIF
      RETURN
      END
* SUBROUTINE TCFU07             ALL SYSTEMS                 91/12/01
C PORTABILITY : ALL SYSTEMS
C 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  VALUES OF TEST FUNCTIONS FOR NONLINEAR APPROXIMATION.
*  UNIVERSAL VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  KC  INDEX OF THE CONSTRAINT.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  FC VALUE OF THE CONSTRAINT AT THE SELECTED POINT.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TCFU07(N,KC,X,FC,NEXT)
      INTEGER N,KC,NEXT
      DOUBLE PRECISION X(N),FC
      DOUBLE PRECISION X1,X2,X3,X4,X5,X6,X7,X8
      INTEGER I,J,K
      DOUBLE PRECISION Y(128)
      COMMON /EMPR07/ Y
      GO TO(10,20,30,40,50,60,70,80,90,100,10,91,130,140,150,160,170,
     * 180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,
     * 340),NEXT
   10 FC=X(1)**2+X(2)**2
      RETURN
   20 FC=3.0D 0*X(1)**2-2.0D 0*X(1)*X(2)+X(2)**2
      RETURN
   30 FC=4.0D 0*X(1)**2+X(2)**2
      RETURN
   40 FC=(1.0D 0-X(1))**3-X(2)
      RETURN
   50 FC=0.49D 0*X(2)-X(1)*X(2)
      RETURN
   60 GO TO (61,62),KC
   61 FC=X(2)**2-X(1)
      RETURN
   62 FC=X(1)**2-X(2)
      RETURN
   70 GO TO (71,72),KC
   71 FC=(X(1)-5.0D 0)**2+(X(2)-5.0D 0)**2
      RETURN
   72 FC=(X(1)-6.0D 0)**2+(X(2)-5.0D 0)**2
      RETURN
   80 GO TO (81,82,10),KC
   81 FC=X(1)+X(2)**2
      RETURN
   82 FC=X(1)**2+X(2)
      RETURN
   90 GO TO (91,92,93),KC
   91 FC=X(1)*X(2)
      RETURN
   92 FC=X(2)-X(1)**2/125.0D 0
      RETURN
   93 FC=(X(2)-5.0D 1)**2-5.0D 0*(X(1)-5.5D 1)
      RETURN
  100 FC=X(1)**2+2.0D 0*X(2)**2+4.0D 0*X(3)**2
      RETURN
  130 FC=4.0D 0/X(1)+32.0D 0/X(2)+120.0D 0/X(3)
      RETURN
  140 GO TO (141,142),KC
  141 FC=EXP(X(1))-X(2)
      RETURN
  142 FC=EXP(X(2))-X(3)
      RETURN
  150 X2=1.6D 0*X(1)
  151 X3=1.22D 0*X2-X(1)
      X6=(X(2)+X3)/X(1)
      X1=0.01D 0*X(1)*(112.0D 0+13.167D 0*X6-0.6667D 0*X6**2)
      IF (ABS(X1-X2).GT.0.001D 0) THEN
      X2=X1
      GO TO 151
      ENDIF
      X4=93.0D 0
  152 X5=86.35D 0+1.098D 0*X6-0.038D 0*X6**2+0.325D 0*(X4-89.0D 0)
      X8=3.0D 0*X5-133.0D 0
      X7=35.82D 0-0.222D 0*X8
      X1=98000.0D 0*X(3)/(X2*X7+1000.0D 0*X(3))
      IF (ABS(X1-X4).GT.0.001D 0) THEN
      X4=X1
      GO TO 152
      ENDIF
      GO TO (153,154,155,156,157,158,159),KC
  153 FC=X2
      RETURN
  154 FC=X3
      RETURN
  155 FC=X4
      RETURN
  156 FC=X5
      RETURN
  157 FC=X6
      RETURN
  158 FC=X7
      RETURN
  159 FC=X8
      RETURN
  160 FC=X(3)+(1.0D 0-X(3))*X(4)
      RETURN
  170 GO TO (171,172),KC
  171 FC=4.0D 0/X(1)+2.25D 0/X(2)+1.0D 0/X(3)+0.25D 0/X(4)
      RETURN
  172 FC=0.16D 0/X(1)+0.36D 0/X(2)+0.64D 0/X(3)+0.64D 0/X(4)
      RETURN
  180 X1=X(1)*X(1)
      X2=X(2)*X(2)
      X3=X(3)*X(3)
      X4=X(4)*X(4)
      X5=X(1)+X(1)
      X6=X(2)+X(2)
      X7=X(3)+X(3)
      X8=X(4)+X(4)
      GO TO (181,182,183),KC
  181 FC=X1+X2+X3+X4+X(1)-X(2)+X(3)-X(4)
      RETURN
  182 FC=X1+X2+X3+X4+X4-X(1)-X(4)
      RETURN
  183 FC=X1+X1+X2+X3+X5-X(2)-X(4)
      RETURN
  190 GO TO (191,192,193),KC
  191 FC=X(1)**2+X(2)**2+X(3)**2+X(4)**2+X(5)**2
      RETURN
  192 FC=X(1)**2*X(3)+X(4)*X(5)
      RETURN
  193 FC=X(2)**2*X(4)+10.0D 0*X(1)*X(5)
      RETURN
  200 GO TO (201,202,203),KC
  201 FC=85.334407D 0+0.0056868D 0*X(2)*X(5)+0.0006262D 0*X(1)*X(4)-
     * 0.0022053D 0*X(3)*X(5)
      RETURN
  202 FC=80.51249D 0+0.0071317D 0*X(2)*X(5)+0.0029955D 0*X(1)*X(2)+
     * 0.0021813D 0*X(3)**2
      RETURN
  203 FC=9.300961D 0+0.0047026D 0*X(3)*X(5)+0.0012547D 0*X(1)*X(3)+
     * 0.0019085D 0*X(3)*X(4)
      RETURN
  210 GO TO (211,212,213),KC
  211 FC=-145421.4020D 0*X(1)+2931.15060D 0*X(1)*X(2)-
     * 40.4279320D 0*X(1)*X(3)+5106.1920D 0*X(1)*X(4)+
     * 15711.36D 0*X(1)*X(5)
      RETURN
  212 FC=-155011.1084D 0*X(1)+4360.53352D 0*X(1)*X(2)+
     * 12.9492344D 0*X(1)*X(3)+10236.884D 0*X(1)*X(4)+
     * 13176.786D 0*X(1)*X(5)
      RETURN
  213 FC=-326669.5104D 0*X(1)+7390.68412D 0*X(1)*X(2)-
     * 27.8986976D 0*X(1)*X(3)+16643.076D 0*X(1)*X(4)+
     * 30988.146D 0*X(1)*X(5)
      RETURN
  220 CONTINUE
      X1=146.312D 3
      Y(1)=X(2)+X(3)+41.6D 0
      Y(18)=0.024D 0*X(4)-4.62D 0
      Y(2)=12.5D 0/Y(18)+12.0D 0
      Y(19)=0.0003535D 0*X(1)**2+0.5311D 0*X(1)+0.08705D 0*Y(2)*X(1)
      Y(20)=0.052D 0*X(1)+78.0D 0+0.002377D 0*Y(2)*X(1)
      Y(3)=Y(19)/Y(20)
      Y(4)=19.0D 0*Y(3)
      Y(21)=0.04782D 0*(X(1)-Y(3))+0.1956D 0*(X(1)-Y(3))**2/X(2)+
     * 0.6376D 0*Y(4)+1.594D 0*Y(3)
      Y(22)=100.0D 0*X(2)
      Y(23)=X(1)-Y(3)-Y(4)
      Y(24)=0.95D 0-Y(21)/Y(22)
      Y(5)=Y(23)*Y(24)
      Y(6)=X(1)-Y(5)-Y(4)-Y(3)
      Y(25)=(Y(5)+Y(4))*0.995D 0
      Y(7)=Y(25)/Y(1)
      Y(8)=Y(25)/3798.0D 0
      Y(26)=Y(7)-0.0663D 0*Y(7)/Y(8)-0.3153D 0
      Y(9)=96.82D 0/Y(26)+0.321D 0*Y(1)
      Y(10)=1.29D 0*Y(5)+1.258D 0*Y(4)+2.29D 0*Y(3)+1.71D 0*Y(6)
      Y(11)=1.71D 0*X(1)-0.452D 0*Y(4)+0.58D 0*Y(3)
      Y(27)=12.3D 0/752.3D 0
      Y(28)=1.75D 0*Y(2)*0.995D 0*X(1)
      Y(29)=0.995D 0*Y(10)+1998.0D 0
      Y(12)=Y(27)*X(1)+Y(28)/Y(29)
      Y(13)=Y(29)-1.75D 0*Y(2)
      Y(14)=3623.0D 0+64.4D 0*X(2)+58.4D 0*X(3)+X1/(Y(9)+
     * X(5))
      Y(30)=0.995D 0*Y(10)+60.8D 0*X(2)+48.0D 0*X(4)-0.1121D 0*Y(14)-
     * 5095.0D 0
      Y(15)=Y(13)/Y(30)
      Y(16)=148.0D 3-331.0D 3*Y(15)+40.0D 0*Y(13)-61.0D 0*Y(15)*Y(13)
      Y(31)=2324.0D 0*Y(10)-2.874D 7*Y(2)
      Y(17)=1.413D 7-1328.0D 0*Y(10)-531.0D 0*Y(11)+Y(31)/Y(29)
      Y(32)=Y(13)/Y(15)-Y(13)/0.52D 0
      Y(33)=1.104D 0-0.72D 0*Y(15)
      Y(34)=Y(9)+X(5)
      IF(KC.GT.4) GO TO 225
      GO TO (221,222,223,224),KC
  221 FC=Y(4)-0.28D 0*Y(5)/0.72D 0
      RETURN
  222 FC=1.5D 0*X(2)-X(3)
      RETURN
  223 FC=3496.0D 0*Y(2)/Y(29)
      RETURN
  224 FC=62212.0D 0/Y(34)-Y(1)
      RETURN
  225 FC=Y(KC-4)
      RETURN
  230 GO TO (231,232),KC
  231 FC=X(1)**2+X(2)**2+X(3)**2
      RETURN
  232 FC=X(5)**2+(X(4)-3.0D 0)**2
      RETURN
  240 GO TO (241,242),KC
  241 FC=0.001D 0*X(1)*X(2)*X(3)*X(4)*X(5)*X(6)
      RETURN
  242 FC=0.00062D 0*X(1)*X(4)*X(5)**2*(X(1)+X(2)+X(3))+
     * 0.00058D 0*X(2)*X(3)*X(6)**2*(X(1)+1.57D 0*X(2)+X(4))
      RETURN
  250 GO TO (251,252,253,254),KC
  251 FC=17.1D 0*X(1)+38.2D 0*X(2)+204.2D 0*X(3)+212.3D 0*X(4)+
     * 623.4D 0*X(5)+1495.5D 0*X(6)-169.0D 0*X(1)*X(3)-
     * 3580.0D 0*X(3)*X(5)-3810.0D 0*X(4)*X(5)-18500.0D 0*X(4)*X(6)-
     * 24300.0D 0*X(5)*X(6)
      RETURN
  252 FC=17.9D 0*X(1)+36.8D 0*X(2)+113.9D 0*X(3)+169.7D 0*X(4)+
     * 337.8D 0*X(5)+1385.2D 0*X(6)-139.0D 0*X(1)*X(3)-
     * 2450.0D 0*X(4)*X(5)-16600.0D 0*X(4)*X(6)-17200.0D 0*X(5)*X(6)
      RETURN
  253 FC=-273.0D 0*X(2)-70.0D 0*X(4)-819.0D 0*X(5)+
     * 26000.0D 0*X(4)*X(5)
      RETURN
  254 FC=159.9D 0*X(1)-311.0D 0*X(2)+587.0D 0*X(4)+391.0D 0*X(5)+
     * 2198.0D 0*X(6)-14000.0D 0*X(1)*X(6)
      RETURN
  260 GO TO (261,262,263,264),KC
  261 FC=2.0D 0*X(1)**2+3.0D 0*X(2)**4+X(3)+4.0D 0*X(4)**2+
     &5.0D 0*X(5)
      RETURN
  262 FC=7.0D 0*X(1)+3.0D 0*X(2)+1.0D 1*X(3)**2+X(4)-X(5)
      RETURN
  263 FC=2.3D 1*X(1)+X(2)**2+6.0D 0*X(6)**2-8.0D 0*X(7)
      RETURN
  264 FC=4.0D 0*X(1)**2+X(2)**2-3.0D 0*X(1)*X(2)+2.0D 0*X(3)**2+
     &5.0D 0*X(6)-1.1D 1*X(7)
      RETURN
  270 GO TO (271,272,273,274,275),KC
  271 FC=0.5D 0*X(1)**0.5D 0*X(7)/(X(3)*X(6)**2)+
     * 0.7D 0*X(1)**3*X(2)*X(6)*X(7)**0.5D 0/X(3)**2+
     * 0.2D 0*X(3)*X(6)**(2.0D 0/3.0D 0)*X(7)**0.25D 0/
     * (X(2)*X(4)**0.5D 0)
      RETURN
  272 FC=1.3D 0*X(2)*X(6)/(X(1)**0.5D 0*X(3)*X(5))+
     * 0.8D 0*X(3)*X(6)**2/(X(4)*X(5))+
     * 3.1D 0*X(2)**0.5D 0*X(6)**(1.0D 0/3.0D 0)/
     * (X(1)*X(4)**2*X(5))
      RETURN
  273 FC=2.0D 0*X(1)*X(5)*X(7)**(4.0D 0/3.0D 0)/(X(3)**1.5D 0*X(6))+
     * 0.1D 0*X(2)*X(5)/(X(3)**0.5D 0*X(6)*X(7)**0.5D 0)+
     * X(2)*X(3)**0.5D 0*X(5)/X(1)+
     * 0.65D 0*X(3)*X(5)*X(7)/(X(2)**2*X(6))
      RETURN
  274 FC=0.2D 0*X(2)*X(5)**0.5D 0*X(7)**(1.0D 0/3.0D 0)/
     * (X(1)**2*X(4))+
     * 0.3D 0*X(1)**0.5D 0*X(2)**2*X(3)*X(4)**(1.0D 0/3.0D 0)*
     * X(7)**0.25D 0/X(5)**(2.0D 0/3.0D 0)+
     * 0.4D 0*X(3)*X(5)*X(7)**0.75D 0/(X(1)**3*X(2)**2)+
     * 0.5D 0*X(4)*X(7)**0.5D 0/X(3)**2
      RETURN
  275 FC=10.0D 0*X(1)*X(4)**2/(X(2)*X(6)**3*X(7)**0.25D 0)+
     * 15.0D 0*X(3)*X(4)/(X(1)*X(2)**2*X(5)*X(7)**0.5D 0)+
     * 20.0D 0*X(2)*X(6)/(X(1)**2*X(4)*X(5)**2)+
     * 25.0D 0*X(1)**2*X(2)**2*X(5)**0.5D 0*X(7)/(X(3)*X(6)**2)
      RETURN
  280 GO TO (281,282,283,284,285),KC
  281 FC=0.0588D 0*X(5)*X(7)+0.1D 0*X(1)
      RETURN
  282 FC=0.0588D 0*X(6)*X(8)+0.1D 0*(X(1)+X(2))
      RETURN
  283 FC=4.0D 0*X(3)/X(5)+2.0D 0/(X(3)**0.71D 0*X(5))+
     * 0.0588D 0*X(7)/X(3)**1.3D 0
      RETURN
  284 FC=4.0D 0*X(4)/X(6)+2.0D 0/(X(4)**0.71D 0*X(6))+
     * 0.0588D 0*X(8)/X(4)**1.3D 0
      RETURN
  285 FC=0.4D 0*(X(1)/X(7))**0.67D 0+0.4D 0*(X(2)/X(8))**0.67D 0+
     * 10.0D 0-X(1)-X(2)
      RETURN
  290 GO TO (291,292,293,294,295,296),KC
  291 FC=0.0025D 0*(X(4)+X(6))
      RETURN
  292 FC=0.0025D 0*(X(5)+X(7)-X(4))
      RETURN
  293 FC=0.01D 0*(X(8)-X(5))
      RETURN
  294 FC=X(1)*X(6)-833.33252D 0*X(4)-100.0D 0*X(1)+83333.333D 0
      RETURN
  295 FC=X(2)*X(7)-1250.0D 0*X(5)-X(2)*X(4)+1250.0D 0*X(4)
      RETURN
  296 FC=X(3)*X(8)-1250000.0D 0-X(3)*X(5)+2500.0D 0*X(5)
      RETURN
  300 GO TO (301,302,303,304,305,306,307,308,309,110,111,112,113),KC
  301 FC=X(3)**2+X(4)**2
      RETURN
  302 FC=X(5)**2+X(6)**2
      RETURN
  303 FC=(X(1)-X(5))**2+(X(2)-X(6))**2
      RETURN
  304 FC=(X(1)-X(7))**2+(X(2)-X(8))**2
      RETURN
  305 FC=(X(3)-X(5))**2+(X(4)-X(6))**2
      RETURN
  306 FC=(X(3)-X(7))**2+(X(4)-X(8))**2
      RETURN
  307 FC=X(7)**2+(X(8)-X(9))**2
      RETURN
  308 FC=X(1)**2+(X(2)-X(9))**2
      RETURN
  309 FC=X(9)**2
      RETURN
  110 FC=X(3)*X(9)
      RETURN
  111 FC=X(5)*X(8)-X(6)*X(7)
      RETURN
  112 FC=X(1)*X(4)-X(2)*X(3)
      RETURN
  113 FC=X(5)*X(9)
      RETURN
  310 GO TO (311,312,313,314,315,316,317,318),KC
  311 FC=3.0D 0*(X(1)-2.0D 0)**2+4.0D 0*(X(2)-
     &3.0D 0)**2+2.0D 0*X(3)**2-7.0D 0*X(4)
      RETURN
  312 FC=5.0D 0*X(1)**2+8.0D 0*X(2)+(X(3)-6.0D 0)**2-
     &2.0D 0*X(4)
      RETURN
  313 FC=0.5D 0*(X(1)-8.0D 0)**2+2.0D 0*(X(2)-
     &4.0D 0)**2+3.0D 0*X(5)**2-X(6)
      RETURN
  314 FC=X(1)**2+2.0D 0*(X(2)-2.0D 0)**2-
     &2.0D 0*X(1)*X(2)+1.4D 1*X(5)-6.0D 0*X(6)
      RETURN
  315 FC=4.0D 0*X(1)+5.0D 0*X(2)-3.0D 0*X(7)+
     &9.0D 0*X(8)
      RETURN
  316 FC=1.0D 1*X(1)-8.0D 0*X(2)-1.7D 1*X(7)+
     &2.0D 0*X(8)
      RETURN
  317 FC=6.0D 0*X(2)-3.0D 0*X(1)+1.2D 1*(X(9)-
     &8.0D 0)**2-7.0D 0*X(10)
      RETURN
  318 FC=2.0D 0*X(2)-8.0D 0*X(1)+5.0D 0*X(9)-
     &2.0D 0*X(10)
      RETURN
  320 FC=3.0D 0*Y(KC+85)*X(KC+10)**2
      K=(KC-1)*5
      DO 321 J=1,5
      FC=FC+2.0D 0*Y(K+J+60)*X(J+10)
  321 CONTINUE
      K=(KC-1)*10
      DO 322 I=1,10
      FC=FC-Y(K+I)*X(I)
  322 CONTINUE
      RETURN
  330 IF (KC.LT.6) THEN
      K=11+3*(KC-1)
      FC=(Y(K)+Y(K+1)*X(KC+5)+Y(K+2)*X(KC))*X(KC)/X(KC+5)
      RETURN
      ELSE IF (KC.LT.14) THEN
      GO TO (331,332,333,334,335,336,337,338),KC-5
  331 FC=(Y(26)*X(6)+(Y(27)*X(1)+Y(28)*X(6))*X(12)/X(11))/X(7)
      RETURN
  332 FC=(Y(29)*X(7)+(Y(30)*X(7)+Y(33)*X(1))*X(12)+Y(31)*X(2)*X(13))/
     * X(8)+Y(32)*X(13)
      RETURN
  333 FC=Y(34)*X(8)+Y(37)*X(9)+(Y(35)*X(8)+Y(38)*X(2))*X(13)+
     * (Y(36)*X(3)+Y(39)*X(9))*X(14)
      RETURN
  334 FC=((Y(40)*X(14)+Y(43))*X(9)+(Y(41)*X(4)+Y(44)*X(8))*X(15)+
     * Y(42)*X(10))/(X(3)*X(14))
      RETURN
  335 FC=(Y(45)*X(5)*X(16)+Y(49)*X(10))/(X(4)*X(15))+
     * Y(46)*X(10)/X(4)+(Y(47)+Y(48)*X(16))/X(15)
      RETURN
  336 FC=(Y(50)+Y(52)*X(5)*X(16))/X(4)+Y(51)*X(16)
      RETURN
  337 FC=Y(53)*X(11)+Y(54)*X(12)
      RETURN
  338 FC=Y(55)*X(12)/X(11)
      RETURN
      ELSE IF (KC.LT.18) THEN
      I=KC-13
      K=5-I
      FC=Y(55+I)*X(K)/X(K+1)
      RETURN
      ELSE IF (KC.EQ.18) THEN
      FC=Y(60)*X(9)/X(10)
      RETURN
      ELSE
      FC=Y(61)*X(8)/X(9)
      RETURN
      ENDIF
  340 GO TO (311,312,313,314,315,316,317,318,341,342,343,344,345,346,
     * 347,348,349),KC
  341 FC=X(1)+X(2)+4.0D 0*X(11)-2.1D 1*X(12)
      RETURN
  342 FC=X(1)**2+1.5D 1*X(11)-8.0D 0*X(12)
      RETURN
  343 FC=4.0D 0*X(1)+9.0D 0*X(2)+5.0D 0*X(13)**2-9.0D 0*X(14)
      RETURN
  344 FC=3.0D 0*X(1)+4.0D 0*X(2)+3.0D 0*(X(13)-
     &6.0D 0)**2-1.4D 1*X(14)
      RETURN
  345 FC=1.4D 1*X(1)**2+3.5D 1*X(15)-7.9D 1*X(16)
      RETURN
  346 FC=1.5D 1*X(2)**2+1.1D 1*X(15)-6.1D 1*X(16)
      RETURN
  347 FC=5.0D 0*X(1)**2+2.0D 0*X(2)+9.0D 0*X(17)**4-X(18)
      RETURN
  348 FC=X(1)**2-X(2)+1.9D 1*X(19)-2.0D 1*X(20)
      RETURN
  349 FC=7.0D 0*X(1)**2+5.0D 0*X(2)**2+X(19)**2-3.0D 1*X(20)
      RETURN
      END
* SUBROUTINE TCGU07             ALL SYSTEMS                 90/12/01
C PORTABILITY : ALL SYSTEMS
C 90/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  GRADIENTS OF TEST FUNCTIONS FOR NONLINEAR PROGRAMMING.
*  UNIVERSAL VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  KC  INDEX OF THE CONSTRAINT.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  GC GRADIENT OF THE CONSTRAINT FUNCTION.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TCGU07(N,KC,X,GC,NEXT)
      INTEGER N,KC,NEXT
      DOUBLE PRECISION X(N),GC(N)
      DOUBLE PRECISION X1,X2,X3,X4,X5,X6,X7,X8
      INTEGER I,J,K
      DOUBLE PRECISION Y(128)
      COMMON /EMPR07/ Y
      GO TO(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,
     * 180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,
     * 340),NEXT
   10 GC(1)=2.0D 0*X(1)
      GC(2)=2.0D 0*X(2)
      RETURN
   20 GC(1)=6.0D 0*X(1)-2.0D 0*X(2)
      GC(2)=2.0D 0*(X(2)-X(1))
      RETURN
   30 GC(1)=8.0D 0*X(1)
      GC(2)=2.0D 0*X(2)
      RETURN
   40 GC(1)=-3.0D 0*(1.0D 0-X(1))**2
      GC(2)=-1.0D 0
      RETURN
   50 GC(1)=-X(2)
      GC(2)=0.49D 0-X(1)
      RETURN
   60 GO TO (61,62),KC
   61 GC(1)=-1.0D 0
      GC(2)= 2.0D 0*X(2)
      RETURN
   62 GC(1)=2.0D 0*X(1)
      GC(2)=-1.0D 0
      RETURN
   70 GO TO (71,72),KC
   71 GC(1)=2.0D 0*(X(1)-5.0D 0)
      GC(2)=2.0D 0*(X(2)-5.0D 0)
      RETURN
   72 GC(1)=2.0D 0*(X(1)-6.0D 0)
      GC(2)=2.0D 0*(X(2)-5.0D 0)
      RETURN
   80 GO TO (81,82,10),KC
   81 GC(1)=1.0D 0
      GC(2)=2.0D 0*X(2)
      RETURN
   82 GC(1)=2.0D 0*X(1)
      GC(2)=1.0D 0
      RETURN
   90 GO TO (91,92,93),KC
   91 GC(1)=X(2)
      GC(2)=X(1)
      RETURN
   92 GC(1)=-2.0D 0*X(1)/125.0D 0
      GC(2)= 1.0D 0
      RETURN
   93 GC(1)=-5.0D 0
      GC(2)= 2.0D 0*(X(2)-5.0D 1)
      RETURN
  100 GC(1)= 2.0D 0*X(1)
      GC(2)= 4.0D 0*X(2)
      GC(3)= 8.0D 0*X(3)
      RETURN
  110 GC(1)=2.0D 0*X(1)
      GC(2)=2.0D 0*X(2)
      GC(3)=0.0D 0
      RETURN
  120 GC(1)=X(2)
      GC(2)=X(1)
      GC(3)=0.0D 0
      RETURN
  130 GC(1)=-4.0D 0/X(1)**2
      GC(2)=-3.2D 1/X(2)**2
      GC(3)=-1.2D 2/X(3)**2
      RETURN
  140 GO TO (141,142),KC
  141 GC(1)=EXP(X(1))
      GC(2)=-1.0D 0
      GC(3)= 0.0D 0
      RETURN
  142 GC(1)= 0.0D 0
      GC(2)=EXP(X(2))
      GC(3)=-1.0D 0
      RETURN
  150 X2=1.6D 0*X(1)
      Y(4)=1.6D 0
      Y(5)=0.0D 0
      Y(6)=0.0D 0
  151 X3=1.22D 0*X2-X(1)
      Y(7)=1.22D 0*Y(4)-1.0D 0
      Y(8)=1.22D 0*Y(5)
      Y(9)=1.22D 0*Y(6)
      X6=(X(2)+X3)/X(1)
      Y(16)=Y(7)/X(1)-(X(2)+X3)/X(1)**2
      Y(17)=(1.0D 0+Y(8))/X(1)
      Y(18)=Y(9)/X(1)
      X1=0.01D 0*X(1)*(112.0D 0+13.167D 0*X6-0.6667D 0*X6**2)
      Y(1)=0.01D 0*(112.0D 0+13.167D 0*X6-0.6667D 0*X6**2)+
     * 0.01D 0*X(1)*(13.167*Y(16)-1.3334D 0*X6*Y(16))
      Y(2)=0.01D 0*X(1)*(13.167D 0*Y(17)-1.3334D 0*X6*Y(17))
      Y(3)=0.01D 0*X(1)*(13.167D 0*Y(18)-1.3334D 0*X6*Y(18))
      IF (ABS(X1-X2).GT.0.001D 0) THEN
      X2=X1
      Y(4)=Y(1)
      Y(5)=Y(2)
      Y(6)=Y(3)
      GO TO 151
      ENDIF
      X4=93.0D 0
      Y(10)=0.0D 0
      Y(11)=0.0D 0
      Y(12)=0.0D 0
  152 X5=86.35D 0+1.098D 0*X6-0.038D 0*X6**2+0.325D 0*(X4-89.0D 0)
      Y(13)=1.098D 0*Y(16)-0.076D 0*X6*Y(16)+0.325D 0*Y(10)
      Y(14)=1.098D 0*Y(17)-0.076D 0*X6*Y(17)+0.325D 0*Y(11)
      Y(15)=1.098D 0*Y(18)-0.076D 0*X6*Y(18)+0.325D 0*Y(12)
      X8=3.0D 0*X5-133.0D 0
      Y(22)=3.0D 0*Y(13)
      Y(23)=3.0D 0*Y(14)
      Y(24)=3.0D 0*Y(15)
      X7=35.82D 0-0.222D 0*X8
      Y(19)=-0.222D 0*Y(22)
      Y(20)=-0.222D 0*Y(23)
      Y(21)=-0.222D 0*Y(24)
      Y(89)=X2*X7+1000.0D 0*X(3)
      Y(90)=Y(89)**2
      X1=98000.0D 0*X(3)/Y(89)
      Y(1)=-98000.0D 0*X(3)*(Y(4)*X7+X2*Y(19))/Y(90)
      Y(2)=-98000.0D 0*X(3)*(Y(5)*X7+X2*Y(20))/Y(90)
      Y(3)=-98000.0D 0*X(3)*(Y(6)*X7+X2*Y(21)+1000.0D 0)/Y(90)+
     * 98000.0D 0/Y(89)
      IF (ABS(X1-X4).GT.0.001D 0) THEN
      X4=X1
      Y(10)=Y(1)
      Y(11)=Y(2)
      Y(12)=Y(3)
      GO TO 152
      ENDIF
      GO TO (153,154,155,156,157,158,159),KC
  153 CONTINUE
      GC(1)=Y(4)
      GC(2)=Y(5)
      GC(3)=Y(6)
      RETURN
  154 CONTINUE
      GC(1)=Y(7)
      GC(2)=Y(8)
      GC(3)=Y(9)
      RETURN
  155 CONTINUE
      GC(1)=Y(10)
      GC(2)=Y(11)
      GC(3)=Y(12)
      RETURN
  156 CONTINUE
      GC(1)=Y(13)
      GC(2)=Y(14)
      GC(3)=Y(15)
      RETURN
  157 CONTINUE
      GC(1)=Y(16)
      GC(2)=Y(17)
      GC(3)=Y(18)
      RETURN
  158 CONTINUE
      GC(1)=Y(19)
      GC(2)=Y(20)
      GC(3)=Y(21)
      RETURN
  159 CONTINUE
      GC(1)=Y(22)
      GC(2)=Y(23)
      GC(3)=Y(24)
      RETURN
  160 GC(1)=0.0D 0
      GC(2)=0.0D 0
      GC(3)=1.0D 0-X(4)
      GC(4)=1.0D 0-X(3)
      RETURN
  170 GO TO (171,172),KC
  171 GC(1)=-4.00D 0/X(1)**2
      GC(2)=-2.25D 0/X(2)**2
      GC(3)=-1.00D 0/X(3)**2
      GC(4)=-0.25D 0/X(4)**2
      RETURN
  172 GC(1)=-0.16D 0/X(1)**2
      GC(2)=-0.36D 0/X(2)**2
      GC(3)=-0.64D 0/X(3)**2
      GC(4)=-0.64D 0/X(4)**2
      RETURN
  180 X5=X(1)+X(1)
      X6=X(2)+X(2)
      X7=X(3)+X(3)
      X8=X(4)+X(4)
  181 GO TO (182,183,184),KC
  182 GC(1)=X5+1.0D 0
      GC(2)=X6-1.0D 0
      GC(3)=X7+1.0D 0
      GC(4)=X8-1.0D 0
      RETURN
  183 GC(1)=X5-1.0D 0
      GC(2)=X6+X6
      GC(3)=X7
      GC(4)=X8+X8-1.0D 0
      RETURN
  184 GC(1)=X5+X5+2.0D 0
      GC(2)=X6-1.0D 0
      GC(3)=X7
      GC(4)=-1.0D 0
      RETURN
  190 GO TO (191,192,193),KC
  191 GC(1)=2.0D 0*X(1)
      GC(2)=2.0D 0*X(2)
      GC(3)=2.0D 0*X(3)
      GC(4)=2.0D 0*X(4)
      GC(5)=2.0D 0*X(5)
      RETURN
  192 GC(1)=2.0D 0*X(1)*X(3)
      GC(2)=0.0D 0
      GC(3)=X(1)**2
      GC(4)=X(5)
      GC(5)=X(4)
      RETURN
  193 GC(1)=10.0D 0*X(5)
      GC(2)=2.0D 0*X(2)*X(4)
      GC(3)=0.0D 0
      GC(4)=X(2)**2
      GC(5)=10.0D 0*X(1)
      RETURN
  200 GO TO (201,202,203),KC
  201 GC(1)= 0.0006262D 0*X(4)
      GC(2)= 0.0056868D 0*X(5)
      GC(3)=-0.0022053D 0*X(5)
      GC(4)= 0.0006262D 0*X(1)
      GC(5)= 0.0056868D 0*X(2)-0.0022053D 0*X(3)
      RETURN
  202 GC(1)=0.0029955D 0*X(2)
      GC(2)=0.0071317D 0*X(5)+0.0029955D 0*X(1)
      GC(3)=0.0043626D 0*X(3)
      GC(4)=0.0D 0
      GC(5)=0.0071317D 0*X(2)
      RETURN
  203 GC(1)=0.0012547D 0*X(3)
      GC(2)=0.0D 0
      GC(3)=0.0047026D 0*X(5)+0.0012547D 0*X(1)+0.0019085D 0*X(4)
      GC(4)=0.0019085D 0*X(3)
      GC(5)=0.0047026D 0*X(3)
      RETURN
  210 GO TO (211,212,213),KC
  211 GC(1)=-145421.4020D 0+2931.15060D 0*X(2)-
     * 40.4279320D 0*X(3)+5106.1920D 0*X(4)+15711.36D 0*X(5)
      GC(2)=2931.15060D 0*X(1)
      GC(3)=-40.4279320D 0*X(1)
      GC(4)=5106.1920D 0*X(1)
      GC(5)=15711.36D 0*X(1)
      RETURN
  212 GC(1)=-155011.1084D 0+4360.53352D 0*X(2)+
     * 12.9492344D 0*X(3)+10236.884D 0*X(4)+13176.786D 0*X(5)
      GC(2)=4360.53352D 0*X(1)
      GC(3)=12.9492344D 0*X(1)
      GC(4)=10236.884D 0*X(1)
      GC(5)=13176.786D 0*X(1)
      RETURN
  213 GC(1)=-326669.5104D 0+7390.68412D 0*X(2)-
     * 27.8986976D 0*X(3)+16643.076D 0*X(4)+30988.146D 0*X(5)
      GC(2)=7390.68412D 0*X(1)
      GC(3)=-27.8986976D 0*X(1)
      GC(4)=16643.076D 0*X(1)
      GC(5)=30988.146D 0*X(1)
      RETURN
  220 IF(KC.EQ.2) THEN
      GC(1)=0.0D 0
      GC(2)=1.5D 0
      GC(3)=-1.0D 0
      GC(4)=0.0D 0
      GC(5)=0.0D 0
      RETURN
      ENDIF
      X1=146.312D 3
      Y(1)=X(2)+X(3)+41.6D 0
      Y(18)=0.024D 0*X(4)-4.62D 0
      Y(2)=12.5D 0/Y(18)+12.0D 0
      Y(19)=0.0003535D 0*X(1)**2+0.5311D 0*X(1)+0.08705D 0*Y(2)*X(1)
      Y(20)=0.052D 0*X(1)+78.0D 0+0.002377D 0*Y(2)*X(1)
      Y(3)=Y(19)/Y(20)
      Y(4)=19.0D 0*Y(3)
      Y(21)=0.04782D 0*(X(1)-Y(3))+0.1956D 0*(X(1)-Y(3))**2/X(2)+
     * 0.6376D 0*Y(4)+1.594D 0*Y(3)
      Y(22)=100.0D 0*X(2)
      Y(23)=X(1)-Y(3)-Y(4)
      Y(24)=0.95D 0-Y(21)/Y(22)
      Y(5)=Y(23)*Y(24)
      Y(6)=X(1)-Y(5)-Y(4)-Y(3)
      Y(25)=(Y(5)+Y(4))*0.995D 0
      Y(7)=Y(25)/Y(1)
      Y(8)=Y(25)/3798.0D 0
      Y(26)=Y(7)-0.0663D 0*Y(7)/Y(8)-0.3153D 0
      Y(9)=96.82D 0/Y(26)+0.321D 0*Y(1)
      Y(10)=1.29D 0*Y(5)+1.258D 0*Y(4)+2.29D 0*Y(3)+1.71D 0*Y(6)
      Y(11)=1.71D 0*X(1)-0.452D 0*Y(4)+0.58D 0*Y(3)
      Y(27)=12.3D 0/752.3D 0
      Y(28)=1.75D 0*Y(2)*0.995D 0*X(1)
      Y(29)=0.995D 0*Y(10)+1998.0D 0
      Y(12)=Y(27)*X(1)+Y(28)/Y(29)
      Y(13)=Y(29)-1.75D 0*Y(2)
      Y(14)=3623.0D 0+64.4D 0*X(2)+58.4D 0*X(3)+X1/(Y(9)+
     * X(5))
      Y(30)=0.995D 0*Y(10)+60.8D 0*X(2)+48.0D 0*X(4)-0.1121D 0*Y(14)-
     * 5095.0D 0
      Y(15)=Y(13)/Y(30)
      Y(16)=148.0D 3-331.0D 3*Y(15)+40.0D 0*Y(13)-61.0D 0*Y(15)*Y(13)
      Y(31)=2324.0D 0*Y(10)-2.874D 7*Y(2)
      Y(17)=1.413D 7-1328.0D 0*Y(10)-531.0D 0*Y(11)+Y(31)/Y(29)
      Y(32)=Y(13)/Y(15)-Y(13)/0.52D 0
      Y(33)=1.104D 0-0.72D 0*Y(15)
      Y(34)=Y(9)+X(5)
      Y(35)=1.0D 0
      Y(36)=0.024D 0
      Y(37)=-12.5D 0*Y(36)/Y(18)**2
      Y(38)=0.000707D 0*X(1)+0.5311D 0+0.08705D 0*Y(2)
      Y(39)=0.08705D 0*X(1)*Y(37)
      Y(40)=0.052D 0+0.002377D 0*Y(2)
      Y(41)=0.002377D 0*X(1)*Y(37)
      Y(42)=(Y(38)*Y(20)-Y(19)*Y(40))/Y(20)**2
      Y(43)=(Y(39)*Y(20)-Y(19)*Y(41))/Y(20)**2
      Y(44)=19.0D 0*Y(42)
      Y(45)=19.0D 0*Y(43)
      Y(46)=0.04782D 0*(1.0D 0-Y(42))+0.3912D 0*(X(1)-Y(3))*(1.0D 0-
     * Y(42))/X(2)+0.6376D 0*Y(44)+1.594D 0*Y(42)
      Y(47)=-0.1956D 0*(X(1)-Y(3))**2/X(2)**2
      Y(48)=-0.04782D 0*Y(43)-0.3912D 0*(X(1)-Y(3))*Y(43)/X(2)+
     * 0.6376D 0*Y(45)+1.594D 0*Y(43)
      Y(49)=100.0D 0
      Y(50)=1.0D 0-Y(42)-Y(44)
      Y(51)=-Y(43)-Y(45)
      Y(52)=-Y(46)/Y(22)
      Y(53)=(Y(21)*Y(49)-Y(47)*Y(22))/Y(22)**2
      Y(54)=-Y(48)/Y(22)
      Y(55)=Y(50)*Y(24)+Y(23)*Y(52)
      Y(56)=Y(23)*Y(53)
      Y(57)=Y(51)*Y(24)+Y(23)*Y(54)
      Y(58)=1.0D 0-Y(55)-Y(44)-Y(42)
      Y(59)=-Y(56)
      Y(60)=-Y(57)-Y(45)-Y(43)
      Y(61)=0.995D 0*(Y(55)+Y(44))
      Y(62)=0.995D 0*Y(56)
      Y(63)=0.995D 0*(Y(57)+Y(45))
      Y(64)=Y(61)/Y(1)
      Y(65)=(Y(62)*Y(1)-Y(25)*Y(35))/Y(1)**2
      Y(66)=-Y(25)*Y(35)/Y(1)**2
      Y(67)=Y(63)/Y(1)
      Y(68)=Y(61)/3798.0D 0
      Y(69)=Y(62)/3798.0D 0
      Y(70)=Y(63)/3798.0D 0
      Y(71)=Y(64)-0.0663D 0*(Y(64)*Y(8)-Y(7)*Y(68))/Y(8)**2
      Y(72)=Y(65)-0.0663D 0*(Y(65)*Y(8)-Y(7)*Y(69))/Y(8)**2
      Y(73)=Y(66)-0.0663D 0*Y(66)/Y(8)
      Y(74)=Y(67)-0.0663D 0*(Y(67)*Y(8)-Y(7)*Y(70))/Y(8)**2
      Y(75)=-96.82D 0*Y(71)/Y(26)**2
      Y(76)=-96.82D 0*Y(72)/Y(26)**2+0.321D 0*Y(35)
      Y(77)=-96.82D 0*Y(73)/Y(26)**2+0.321D 0*Y(35)
      Y(78)=-96.82D 0*Y(74)/Y(26)**2
      Y(79)=1.29D 0*Y(55)+1.258D 0*Y(44)+2.29D 0*Y(42)+1.71D 0*Y(58)
      Y(80)=1.29D 0*Y(56)+1.71D 0*Y(59)
      Y(81)=1.29D 0*Y(57)+1.258D 0*Y(45)+2.29D 0*Y(43)+1.71D 0*Y(60)
      Y(82)=1.71D 0-0.452D 0*Y(44)+0.58D 0*Y(42)
      Y(83)=-0.452D 0*Y(45)+0.58D 0*Y(43)
      Y(84)=1.75D 0*Y(2)*0.995D 0
      Y(85)=1.75D 0*Y(37)*0.995D 0*X(1)
      Y(86)=0.995D 0*Y(79)
      Y(87)=0.995D 0*Y(80)
      Y(88)=0.995D 0*Y(81)
      Y(89)=Y(27)+(Y(84)*Y(29)-Y(28)*Y(86))/Y(29)**2
      Y(90)=-Y(28)*Y(87)/Y(29)**2
      Y(91)=(Y(85)*Y(29)-Y(28)*Y(88))/Y(29)**2
      Y(92)=Y(88)-1.75D 0*Y(37)
      Y(93)=-X1*Y(75)/(Y(9)+X(5))**2
      Y(94)=64.4D 0-X1*Y(76)/(Y(9)+X(5))**2
      Y(95)=58.4D 0-X1*Y(77)/(Y(9)+X(5))**2
      Y(96)=-X1*Y(78)/(Y(9)+X(5))**2
      Y(97)=-X1/(Y(9)+X(5))**2
      Y(98)=0.995D 0*Y(79)-0.1121D 0*Y(93)
      Y(99)=0.995D 0*Y(80)+60.8D 0-0.1121D 0*Y(94)
      Y(100)=-0.1121D 0*Y(95)
      Y(101)=0.995D 0*Y(81)+48.0D 0-0.1121D 0*Y(96)
      Y(102)=-0.1121D 0*Y(97)
      Y(103)=(Y(86)*Y(30)-Y(13)*Y(98))/Y(30)**2
      Y(104)=(Y(87)*Y(30)-Y(13)*Y(99))/Y(30)**2
      Y(105)=-Y(13)*Y(100)/Y(30)**2
      Y(106)=(Y(92)*Y(30)-Y(13)*Y(101))/Y(30)**2
      Y(107)=-Y(13)*Y(102)/Y(30)**2
      Y(108)=-3.31D 5*Y(103)+40.0D 0*Y(86)-61.0D 0*(Y(103)*Y(13)+
     * Y(15)*Y(86))
      Y(109)=-3.31D 5*Y(104)+40.0D 0*Y(87)-61.0D 0*(Y(104)*Y(13)+
     * Y(15)*Y(87))
      Y(110)=-3.31D 5*Y(105)-61.0D 0*Y(105)*Y(13)
      Y(111)=-3.31D 5*Y(106)+40.0D 0*Y(92)-61.0D 0*(Y(106)*Y(13)+
     * Y(15)*Y(92))
      Y(112)=-3.31D 5*Y(107)-61.0D 0*Y(107)*Y(13)
      Y(113)=2.324D 3*Y(79)
      Y(114)=2.324D 3*Y(80)
      Y(115)=2.324D 3*Y(81)-2.874D 7*Y(37)
      Y(116)=-1.328D 3*Y(79)-531.0D 0*Y(82)+(Y(113)*Y(29)-Y(31)*
     * Y(86))/Y(29)**2
      Y(117)=-1.328D 3*Y(80)+(Y(114)*Y(29)-Y(31)*Y(87))/Y(29)**2
      Y(118)=-1.328D 3*Y(81)-531.0D 0*Y(83)+(Y(115)*Y(29)-Y(31)*
     * Y(88))/Y(29)**2
      Y(119)=(Y(86)*Y(15)-Y(13)*Y(103))/Y(15)**2-Y(86)/0.52D 0
      Y(120)=(Y(87)*Y(15)-Y(13)*Y(104))/Y(15)**2-Y(87)/0.52D 0
      Y(121)=-Y(13)*Y(105)/Y(15)**2
      Y(122)=(Y(92)*Y(15)-Y(13)*Y(106))/Y(15)**2-Y(92)/0.52D 0
      Y(123)=-Y(13)*Y(107)/Y(15)**2
      Y(124)=-0.72D 0*Y(103)
      Y(125)=-0.72D 0*Y(104)
      Y(126)=-0.72D 0*Y(105)
      Y(127)=-0.72D 0*Y(106)
      Y(128)=-0.72D 0*Y(107)
      IF (KC.EQ.1) THEN
      GC(1)=Y(44)-0.28D 0*Y(55)/0.72D 0
      GC(2)=-0.28D 0*Y(56)/0.72D 0
      GC(3)=0.0D 0
      GC(4)=Y(45)-0.28D 0*Y(57)/0.72D 0
      GC(5)=0.0D 0
      RETURN
      ENDIF
      GO TO (503,504,505,506,507,508,509,510,511,512,513,514,515,516,
     * 517,518,519,520,521),KC-2
  503 GC(1)=-3496.0D 0*Y(2)*Y(86)/Y(29)**2
      GC(2)=-3496.0D 0*Y(2)*Y(87)/Y(29)**2
      GC(3)=0.0D 0
      GC(4)=3496.0D 0*(Y(37)*Y(29)-Y(2)*Y(88))/Y(29)**2
      GC(5)=0.0D 0
      RETURN
  504 GC(1)=-62212.0D 0*Y(75)/Y(34)**2
      GC(2)=-62212.0D 0*Y(76)/Y(34)**2-Y(35)
      GC(3)=-62212.0D 0*Y(77)/Y(34)**2-Y(35)
      GC(4)=-62212.0D 0*Y(78)/Y(34)**2
      GC(5)=-62212.0D 0*Y(35)/Y(34)**2
      RETURN
  505 GC(1)=0.0D 0
      GC(2)=Y(35)
      GC(3)=Y(35)
      GC(4)=0.0D 0
      GC(5)=0.0D 0
      RETURN
  506 GC(1)=0.0D 0
      GC(2)=0.0D 0
      GC(3)=0.0D 0
      GC(4)=Y(37)
      GC(5)=0.0D 0
      RETURN
  507 GC(1)=Y(42)
      GC(2)=0.0D 0
      GC(3)=0.0D 0
      GC(4)=Y(43)
      GC(5)=0.0D 0
      RETURN
  508 GC(1)=Y(44)
      GC(2)=0.0D 0
      GC(3)=0.0D 0
      GC(4)=Y(45)
      GC(5)=0.0D 0
      RETURN
  509 GC(1)=Y(55)
      GC(2)=Y(56)
      GC(3)=0.0D 0
      GC(4)=Y(57)
      GC(5)=0.0D 0
      RETURN
  510 GC(1)=Y(58)
      GC(2)=Y(59)
      GC(3)=0.0D 0
      GC(4)=Y(60)
      GC(5)=0.0D 0
      RETURN
  511 GC(1)=Y(64)
      GC(2)=Y(65)
      GC(3)=Y(66)
      GC(4)=Y(67)
      GC(5)=0.0D 0
      RETURN
  512 GC(1)=Y(68)
      GC(2)=Y(69)
      GC(3)=0.0D 0
      GC(4)=Y(70)
      GC(5)=0.0D 0
      RETURN
  513 GC(1)=Y(75)
      GC(2)=Y(76)
      GC(3)=Y(77)
      GC(4)=Y(78)
      GC(5)=0.0D 0
      RETURN
  514 GC(1)=Y(79)
      GC(2)=Y(80)
      GC(3)=0.0D 0
      GC(4)=Y(81)
      GC(5)=0.0D 0
      RETURN
  515 GC(1)=Y(82)
      GC(2)=0.0D 0
      GC(3)=0.0D 0
      GC(4)=Y(83)
      GC(5)=0.0D 0
      RETURN
  516 GC(1)=Y(89)
      GC(2)=Y(90)
      GC(3)=0.0D 0
      GC(4)=Y(91)
      GC(5)=0.0D 0
      RETURN
  517 GC(1)=Y(86)
      GC(2)=Y(87)
      GC(3)=0.0D 0
      GC(4)=Y(92)
      GC(5)=0.0D 0
      RETURN
  518 GC(1)=Y(93)
      GC(2)=Y(94)
      GC(3)=Y(95)
      GC(4)=Y(96)
      GC(5)=Y(97)
      RETURN
  519 GC(1)=Y(103)
      GC(2)=Y(104)
      GC(3)=Y(105)
      GC(4)=Y(106)
      GC(5)=Y(107)
      RETURN
  520 GC(1)=Y(108)
      GC(2)=Y(109)
      GC(3)=Y(110)
      GC(4)=Y(111)
      GC(5)=Y(112)
      RETURN
  521 GC(1)=Y(116)
      GC(2)=Y(117)
      GC(3)=0.0D 0
      GC(4)=Y(118)
      GC(5)=0.0D 0
      RETURN
  230 GO TO (231,232),KC
  231 GC(1)=2.0D 0*X(1)
      GC(2)=2.0D 0*X(2)
      GC(3)=2.0D 0*X(3)
      GC(4)=0.0D 0
      GC(5)=0.0D 0
      GC(6)=0.0D 0
      RETURN
  232 GC(1)=0.0D 0
      GC(2)=0.0D 0
      GC(3)=0.0D 0
      GC(4)=2.0D 0*(X(4)-3.0D 0)
      GC(5)=2.0D 0*X(5)
      GC(6)=0.0D 0
      RETURN
  240 GO TO (241,242),KC
  241 GC(1)=0.001D 0*X(2)*X(3)*X(4)*X(5)*X(6)
      GC(2)=0.001D 0*X(1)*X(3)*X(4)*X(5)*X(6)
      GC(3)=0.001D 0*X(1)*X(2)*X(4)*X(5)*X(6)
      GC(4)=0.001D 0*X(1)*X(2)*X(3)*X(5)*X(6)
      GC(5)=0.001D 0*X(1)*X(2)*X(3)*X(4)*X(6)
      GC(6)=0.001D 0*X(1)*X(2)*X(3)*X(4)*X(5)
      RETURN
  242 GC(1)=0.00124D 0*X(1)*X(4)*X(5)**2+
     * 0.00062D 0*X(4)*X(5)**2*(X(2)+X(3))+
     * 0.00058D 0*X(2)*X(3)*X(6)**2
      GC(2)=0.00062D 0*X(1)*X(4)*X(5)**2+
     * 0.00116D 0*1.57D 0*X(2)*X(3)*X(6)**2+
     * 0.00058D 0*X(3)*X(6)**2*(X(1)+X(4))
      GC(3)=0.00062D 0*X(1)*X(4)*X(5)**2+
     * 0.00058D 0*X(2)*X(6)**2*(X(1)+1.57D 0*X(2)+X(4))
      GC(4)=0.00062D 0*X(1)*X(5)**2*(X(1)+X(2)+X(3))+
     * 0.00058D 0*X(2)*X(3)*X(6)**2
      GC(5)=0.00124D 0*X(1)*X(4)*X(5)*(X(1)+X(2)+X(3))
      GC(6)=0.00116D 0*X(2)*X(3)*X(6)*(X(1)+1.57D 0*X(2)+X(4))
      RETURN
  250 GO TO (251,252,253,254),KC
  251 GC(1)=17.1D 0-169.0D 0*X(3)
      GC(2)=38.2D 0
      GC(3)=204.2D 0-169.0D 0*X(1)-3580.0D 0*X(5)
      GC(4)=212.3D 0-3810.0D 0*X(5)-18500.0D 0*X(6)
      GC(5)=623.4D 0-3580.0D 0*X(3)-3810.0D 0*X(4)-24300.0D 0*X(6)
      GC(6)=1495.5D 0-18500.0D 0*X(4)-24300.0D 0*X(5)
      RETURN
  252 GC(1)=17.9D 0-139.0D 0*X(3)
      GC(2)=36.8D 0
      GC(3)=113.9D 0-139.0D 0*X(1)
      GC(4)=169.7D 0-2450.0D 0*X(5)-16600.0D 0*X(6)
      GC(5)=337.8D 0-2450.0D 0*X(4)-17200.0D 0*X(6)
      GC(6)=1385.2D 0-16600.0D 0*X(4)-17200.0D 0*X(5)
      RETURN
  253 GC(1)=0.0D 0
      GC(2)=-273.0D 0
      GC(3)=0.0D 0
      GC(4)=-70.0D 0+26000.0D 0*X(5)
      GC(5)=-819.0D 0+26000.0D 0*X(4)
      GC(6)=0 0D 0
      RETURN
  254 GC(1)=159.9D 0-14000.0D 0*X(6)
      GC(2)=-311.0D 0
      GC(3)=0.0D 0
      GC(4)=587.0D 0
      GC(5)=391.0D 0
      GC(6)=2198.0D 0-14000.0D 0*X(1)
      RETURN
  260 GC(1)=0.0D 0
      GC(2)=0.0D 0
      GC(3)=0.0D 0
      GC(4)=0.0D 0
      GC(5)=0.0D 0
      GC(6)=0.0D 0
      GC(7)=0.0D 0
      GO TO (261,262,263,264),KC
  261 GC(1)=4.0D 0*X(1)
      GC(2)=1.2D 1*X(2)**3
      GC(3)=1.0D 0
      GC(4)=8.0D 0*X(4)
      GC(5)=5.0D 0
      RETURN
  262 GC(1)=7.0D 0
      GC(2)=3.0D 0
      GC(3)=2.0D 1*X(3)
      GC(4)=1.0D 0
      GC(5)=-1.0D 0
      RETURN
  263 GC(1)=2.3D 1
      GC(2)=2.0D 0*X(2)
      GC(6)=1.2D 1*X(6)
      GC(7)=-8.0D 0
      RETURN
  264 GC(1)=8.0D 0*X(1)-3.0D 0*X(2)
      GC(2)=2.0D 0*X(2)-3.0D 0*X(1)
      GC(3)=4.0D 0*X(3)
      GC(6)=5.0D 0
      GC(7)=-1.1D 1
      RETURN
  270 GO TO (271,272,273,274,275),KC
  271 GC(1)=0.25D 0*X(7)/(X(1)**0.5D 0*X(3)*X(6)**2)+
     * 2.1D 0*X(1)**2*X(2)*X(6)*X(7)**0.5D 0/X(3)**2
      GC(2)=0.7D 0*X(1)**3*X(6)*X(7)**0.5D 0/X(3)**2-
     * 0.2D 0*X(3)*X(6)**(2.0D 0/3.0D 0)*X(7)**(1.0D 0/4.0D 0)/
     * (X(2)**2*X(4)**0.5D 0)
      GC(3)=-0.5D 0*X(1)**0.5D 0*X(7)/(X(3)*X(6))**2-
     * 1.4D 0*X(1)**3*X(2)*X(6)*X(7)**0.5D 0/X(3)**3+
     * 0.2D 0*X(6)**(2.0D 0/3.0D 0)*X(7)**(1.0D 0/4.0D 0)/
     * (X(2)*X(4)**0.5D 0)
      GC(4)=-0.1D 0*X(3)*X(6)**(2.0D 0/3.0D 0)*
     * X(7)**(1.0D 0/4.0D 0)/(X(2)*X(4)**1.5D 0)
      GC(5)=0.0D 0
      GC(6)=-1.0D 0*X(1)**0.5D 0*X(7)/(X(3)*X(6)**3)+
     * 0.7D 0*X(1)**3*X(2)*X(7)**0.5D 0/X(3)**2+
     * (0.4D 0/3.0D 0)*X(3)*X(7)**(1.0D 0/4.0D 0)/
     * (X(2)*X(4)**0.5D 0*X(6)**(1.0D 0/3.0D 0))
      GC(7)=0.5D 0*X(1)**0.5D 0/(X(3)*X(6)**2)+
     * 0.35D 0*X(1)**3*X(2)*X(6)/(X(3)**2*X(7)**0.5D 0)+
     * 0.05D 0*X(3)*X(6)**(2.0D 0/3.0D 0)/
     * (X(2)*X(4)**0.5D 0*X(7)**0.75D 0)
      RETURN
  272 GC(1)=-0.65D 0*X(2)*X(6)/(X(1)**1.5D 0*X(3)*X(5))-
     * 3.1D 0*X(2)**0.5D 0*X(6)**(1.0D 0/3.0D 0)/
     * ((X(1)*X(4))**2*X(5))
      GC(2)=1.3D 0*X(6)/(X(1)**0.5D 0*X(3)*X(5))+
     * 1.55D 0*X(6)**(1.0D 0/3.0D 0)/
     * (X(1)*X(2)**0.5D 0*X(4)**2*X(5))
      GC(3)=-1.3D 0*X(2)*X(6)/(X(1)**0.5D 0*X(3)**2*X(5))+
     * 0.8D 0*X(6)**2/(X(4)*X(5))
      GC(4)=-0.8D 0*X(3)*X(6)**2/(X(4)**2*X(5))-
     * 6.2D 0*X(2)**0.5D 0*X(6)**(1.0D 0/3.0D 0)/
     * (X(1)*X(4)**3*X(5))
      GC(5)=-1.3D 0*X(2)*X(6)/(X(1)**0.5D 0*X(3)*X(5)**2)-
     * 0.8D 0*X(3)*X(6)**2/(X(4)*X(5)**2)-
     * 3.1D 0*X(2)**0.5D 0*X(6)**(1.0D 0/3.0D 0)/
     * (X(1)*(X(4)*X(5))**2)
      GC(6)=1.3D 0*X(2)/(X(1)**0.5D 0*X(3)*X(5))+
     * 1.6D 0*X(3)*X(6)/(X(4)*X(5))+
     * (3.1D 0/3.0D 0)*X(2)**0.5D 0/
     * (X(1)*X(4)**2*X(5)*X(6)**(2.0D 0/3.0D 0))
      GC(7)=0.0D 0
      RETURN
  273 GC(1)=2.0D 0*X(5)*X(7)**(4.0D 0/3.0D 0)/(X(3)**1.5D 0*X(6))-
     * X(2)*X(3)**0.5D 0*X(5)/X(1)**2
      GC(2)=0.1D 0*X(5)/(X(3)**0.5D 0*X(6)*X(7)**0.5D 0)+
     * X(3)**0.5D 0*X(5)/X(1)-
     * 1.3D 0*X(3)*X(5)*X(7)/(X(2)**3*X(6))
      GC(3)=-3.0D 0*X(1)*X(5)*X(7)**(4.0D 0/3.0D 0)/(X(3)**2.5D 0*
     * X(6))-0.05D 0*X(2)*X(5)/(X(3)**1.5D 0*X(6)*X(7)**0.5D 0)+
     * 0.5D 0*X(2)*X(5)/(X(1)*X(3)**0.5D 0)+
     * 0.65D 0*X(5)*X(7)/(X(2)**2*X(6))
      GC(4)=0.0D 0
      GC(5)=2.0D 0*X(1)*X(7)**(4.0D 0/3.0D 0)/(X(3)**1.5D 0*X(6))+
     * 0.1D 0*X(2)/(X(3)**0.5D 0*X(6)*X(7)**0.5D 0)+
     * X(2)*X(3)**0.5D 0/X(1)+
     * 0.65D 0*X(3)*X(7)/(X(2)**2*X(6))
      GC(6)=-2.0D 0*X(1)*X(5)*X(7)**(4.0D 0/3.0D 0)/(X(3)**1.5D 0*
     * X(6)**2)-0.1D 0*X(2)*X(5)/(X(3)**0.5D 0*X(6)**2*X(7)**0.5D 0)-
     * 0.65D 0*X(3)*X(5)*X(7)/(X(2)**2*X(6)**2)
      GC(7)=(8.0D 0/3.0D 0)*X(1)*X(5)*X(7)**(1.0D 0/3.0D 0)/
     * (X(3)**1.5D 0*X(6))-0.05D 0*X(2)*X(5)/(X(3)**0.5D 0*X(6)*
     * X(7)**1.5D 0)+0.65D 0*X(3)*X(5)/(X(2)**2*X(6))
      RETURN
  274 GC(1)=-0.4D 0*X(2)*X(5)**0.5D 0*X(7)**(1.0D 0/3.0D 0)/
     * (X(1)**3*X(4))+
     * 0.15D 0*X(2)**2*X(3)*X(4)**(1.0D 0/3.0D 0)*
     * X(7)**(1.0D 0/4.0D 0)/(X(1)**0.5D 0*X(5)**(2.0D 0/3.0D 0))-
     * 1.2D 0*X(3)*X(5)*X(7)**(3.0D 0/4.0D 0)/(X(1)**4*X(2)**2)
      GC(2)=0.2D 0*X(5)**0.5D 0*X(7)**(1.0D 0/3.0D 0)/
     * (X(1)**2*X(4))+
     * 0.6D 0*X(1)**0.5D 0*X(2)*X(3)*X(4)**(1.0D 0/3.0D 0)*
     * X(7)**0.25D 0/X(5)**(2.0D 0/3.0D 0)-
     * 0.8D 0*X(3)*X(5)*X(7)**0.75D 0/(X(1)*X(2))**3
      GC(3)=0.3D 0*X(1)**0.5D 0*X(2)**2*X(4)**(1.0D 0/3.0D 0)*
     * X(7)**(1.0D 0/4.0D 0)/X(5)**(2.0D 0/3.0D 0)+
     * 0.4D 0*X(5)*X(7)**(3.0D 0/4.0D 0)/(X(1)**3*X(2)**2)-
     * X(4)*X(7)**0.5D 0/X(3)**3
      GC(4)=-0.2D 0*X(2)*X(5)**0.5D 0*X(7)**(1.0D 0/3.0D 0)/
     * (X(1)*X(4))**2+
     * (0.3D 0/3.0D 0)*X(1)**0.5D 0*X(2)**2*X(3)*
     * X(7)**(1.0D 0/4.0D 0)/(X(4)**(2.0D 0/3.0D 0)*X(5)**
     * (2.0D 0/3.0D 0))+0.5D 0*X(7)**0.5D 0/X(3)**2
      GC(5)=0.1D 0*X(2)*X(7)**(1.0D 0/3.0D 0)/
     * (X(1)**2*X(4)*X(5)**0.5D 0)-
     * (0.6D 0/3.0D 0)*X(1)**0.5D 0*X(2)**2*X(3)*X(4)**(1.0D 0/
     * 3.0D 0)*X(7)**(1.0D 0/4.0D 0)/X(5)**(5.0D 0/3.0D 0)+
     * 0.4D 0*X(3)*X(7)**(3.0D 0/4.0D 0)/(X(1)**3*X(2)**2)
      GC(6)=0.0D 0
      GC(7)=(0.2D 0/3.0D 0)*X(2)*X(5)**0.5D 0/
     * (X(1)**2*X(4)*X(7)**(2.0D 0/3.0D 0))+
     * 0.075D 0*X(1)**0.5D 0*X(2)**2*X(3)*X(4)**(1.0D 0/3.0D 0)/
     * (X(5)**(2.0D 0/3.0D 0)*X(7)**(3.0D 0/4.0D 0))+
     * 0.3D 0*X(3)*X(5)/(X(1)**3*X(2)**2*X(7)**(1.0D 0/4.0D 0))+
     * 0.25D 0*X(4)/(X(3)**2*X(7)**0.5D 0)
      RETURN
  275 GC(1)=10.0D 0*X(4)**2/(X(2)*X(6)**3*X(7)**0.25D 0)-
     * 15.0D 0*X(3)*X(4)/((X(1)*X(2))**2*X(5)*X(7)**0.5D 0)-
     * 40.0D 0*X(2)*X(6)/(X(1)**3*X(4)*X(5)**2)+
     * 50.0D 0*X(1)*X(2)**2*X(5)**0.5D 0*X(7)/(X(3)*X(6)**2)
      GC(2)=-10.0D 0*X(1)*X(4)**2/(X(2)**2*X(6)**3*X(7)**0.25D 0)-
     * 30.0D 0*X(3)*X(4)/(X(1)*X(2)**3*X(5)*X(7)**0.5D 0)+
     * 20.0D 0*X(6)/(X(1)**2*X(4)*X(5)**2)+
     * 50.0D 0*X(1)**2*X(2)*X(5)**0.5D 0*X(7)/(X(3)*X(6)**2)
      GC(3)=15.0D 0*X(4)/(X(1)*X(2)**2*X(5)*X(7)**0.5D 0)-
     * 25.0D 0*X(1)**2*X(2)**2*X(5)**0.5D 0*X(7)/(X(3)*X(6))**2
      GC(4)=20.0D 0*X(1)*X(4)/(X(2)*X(6)**3*X(7)**0.25D 0)+
     * 15.0D 0*X(3)/(X(1)*X(2)**2*X(5)*X(7)**0.5D 0)-
     * 20.0D 0*X(2)*X(6)/(X(1)*X(4)*X(5))**2
      GC(5)=-15.0D 0*X(3)*X(4)/(X(1)*(X(2)*X(5))**2*X(7)**0.5D 0)-
     * 40.0D 0*X(2)*X(6)/(X(1)**2*X(4)*X(5)**3)+
     * 12.5D 0*X(1)**2*X(2)**2*X(7)/(X(3)*X(5)**0.5D 0*X(6)**2)
      GC(6)=-30.0D 0*X(1)*X(4)**2/(X(2)*X(6)**4*X(7)**0.25D 0)+
     * 20.0D 0*X(2)/(X(1)**2*X(4)*X(5)**2)-
     * 50.0D 0*X(1)**2*X(2)**2*X(5)**0.5D 0*X(7)/(X(3)*X(6)**3)
      GC(7)=-2.5D 0*X(1)*X(4)**2/(X(2)*X(6)**3*X(7)**1.25D 0)-
     * 7.5D 0*X(3)*X(4)/(X(1)*X(2)**2*X(5)*X(7)**1.5D 0)+
     * 25.0D 0*X(1)**2*X(2)**2*X(5)**0.5D 0/(X(3)*X(6)**2)
      RETURN
  280 GC(1)=0.0D 0
      GC(2)=0.0D 0
      GC(3)=0.0D 0
      GC(4)=0.0D 0
      GC(5)=0.0D 0
      GC(6)=0.0D 0
      GC(7)=0.0D 0
      GC(8)=0.0D 0
      GO TO (281,282,283,284,285),KC
  281 GC(1)=0.1D 0
      GC(5)=0.0588D 0*X(7)
      GC(7)=0.0588D 0*X(5)
      RETURN
  282 GC(1)=0.1D 0
      GC(2)=0.1D 0
      GC(6)=0.0588D 0*X(8)
      GC(8)=0.0588D 0*X(6)
      RETURN
  283 GC(3)=4.0D 0/X(5)-1.42D 0/(X(3)**1.71D 0*X(5))-
     * 0.0588D 0*1.3D 0*X(7)/X(3)**2.3D 0
      GC(5)=-4.0D 0*X(3)/X(5)**2-2.0D 0/(X(3)**0.71D 0*X(5)**2)
      GC(7)=0.0588D 0/X(3)**1.3D 0
      RETURN
  284 GC(4)=4.0D 0/X(6)-1.42D 0/(X(4)**1.71D 0*X(6))-
     * 0.0588D 0*1.3D 0*X(8)/X(4)**2.3D 0
      GC(6)=-4.0D 0*X(4)/X(6)**2-2.0D 0/(X(4)**0.71D 0*X(6)**2)
      GC(8)=0.0588D 0/X(4)**1.3D 0
      RETURN
  285 GC(1)=0.268D 0*(X(7)/X(1))**0.33D 0/X(7)-1.0D 0
      GC(2)=0.268D 0*(X(8)/X(2))**0.33D 0/X(8)-1.0D 0
      GC(7)=-0.268D 0*(X(1)/X(7))**0.67D 0/X(7)
      GC(8)=-0.268D 0*(X(2)/X(8))**0.67D 0/X(8)
      RETURN
  290 GC(1)=0.0D 0
      GC(2)=0.0D 0
      GC(3)=0.0D 0
      GC(4)=0.0D 0
      GC(5)=0.0D 0
      GC(6)=0.0D 0
      GC(7)=0.0D 0
      GC(8)=0.0D 0
      GO TO (291,292,293,294,295,296),KC
  291 GC(4)=0.0025D 0
      GC(6)=0.0025D 0
      RETURN
  292 GC(4)=-0.0025D 0
      GC(5)= 0.0025D 0
      GC(7)= 0.0025D 0
      RETURN
  293 GC(5)=-0.01D 0
      GC(8)= 0.01D 0
      RETURN
  294 GC(1)=X(6)-100.0D 0
      GC(4)=-833.33252D 0
      GC(6)=X(1)
      RETURN
  295 GC(2)=X(7)-X(4)
      GC(4)=-X(2)+1250.0D 0
      GC(5)=-1250.0D 0
      GC(7)=X(2)
      RETURN
  296 GC(3)=X(8)-X(5)
      GC(5)=-X(3)+2500.0D 0
      GC(8)=X(3)
      RETURN
  300 GC(1)=0.0D 0
      GC(2)=0.0D 0
      GC(3)=0.0D 0
      GC(4)=0.0D 0
      GC(5)=0.0D 0
      GC(6)=0.0D 0
      GC(7)=0.0D 0
      GC(8)=0.0D 0
      GC(9)=0.0D 0
      GO TO (301,302,303,304,305,306,307,308,309,910,911,912,913),KC
  301 GC(3)=2.0D 0*X(3)
      GC(4)=2.0D 0*X(4)
      RETURN
  302 GC(5)=2.0D 0*X(5)
      GC(6)=2.0D 0*X(6)
      RETURN
  303 GC(1)= 2.0D 0*(X(1)-X(5))
      GC(2)= 2.0D 0*(X(2)-X(6))
      GC(5)=-2.0D 0*(X(1)-X(5))
      GC(6)=-2.0D 0*(X(2)-X(6))
      RETURN
  304 GC(1)= 2.0D 0*(X(1)-X(7))
      GC(2)= 2.0D 0*(X(2)-X(8))
      GC(7)=-2.0D 0*(X(1)-X(7))
      GC(8)=-2.0D 0*(X(2)-X(8))
      RETURN
  305 GC(3)= 2.0D 0*(X(3)-X(5))
      GC(4)= 2.0D 0*(X(4)-X(6))
      GC(5)=-2.0D 0*(X(3)-X(5))
      GC(6)=-2.0D 0*(X(4)-X(6))
      RETURN
  306 GC(3)= 2.0D 0*(X(3)-X(7))
      GC(4)= 2.0D 0*(X(4)-X(8))
      GC(7)=-2.0D 0*(X(3)-X(7))
      GC(8)=-2.0D 0*(X(4)-X(8))
      RETURN
  307 GC(7)= 2.0D 0*X(7)
      GC(8)= 2.0D 0*(X(8)-X(9))
      GC(9)=-2.0D 0*(X(8)-X(9))
      RETURN
  308 GC(1)= 2.0D 0*X(1)
      GC(2)= 2.0D 0*(X(2)-X(9))
      GC(9)=-2.0D 0*(X(2)-X(9))
      RETURN
  309 GC(9)=2.0D 0*X(9)
      RETURN
  910 GC(3)=X(9)
      GC(9)=X(3)
      RETURN
  911 GC(5)= X(8)
      GC(6)=-X(7)
      GC(7)=-X(6)
      GC(8)= X(5)
      RETURN
  912 GC(1)= X(4)
      GC(2)=-X(3)
      GC(3)=-X(2)
      GC(4)= X(1)
      RETURN
  913 GC(5)=X(9)
      GC(9)=X(5)
      RETURN
  310 GC(1)=0.0D 0
      GC(2)=0.0D 0
      GC(3)=0.0D 0
      GC(4)=0.0D 0
      GC(5)=0.0D 0
      GC(6)=0.0D 0
      GC(7)=0.0D 0
      GC(8)=0.0D 0
      GC(9)=0.0D 0
      GC(10)=0.0D 0
      GO TO (311,312,313,314,315,316,317,318),KC
  311 GC(1)=6.0D 0*(X(1)-2.0D 0)
      GC(2)=8.0D 0*(X(2)-3.0D 0)
      GC(3)=4.0D 0*X(3)
      GC(4)=-7.0D 0
      RETURN
  312 GC(1)=1.0D 1*X(1)
      GC(2)=8.0D 0
      GC(3)=2.0D 0*(X(3)-6.0D 0)
      GC(4)=-2.0D 0
      RETURN
  313 GC(1)=1.0D 0*(X(1)-8.0D 0)
      GC(2)=4.0D 0*(X(2)-4.0D 0)
      GC(5)=6.0D 0*X(5)
      GC(6)=-1.0D 0
      RETURN
  314 GC(1)=2.0D 0*X(1)-2.0D 0*X(2)
      GC(2)=4.0D 0*(X(2)-2.0D 0)-2.0D 0*X(1)
      GC(5)=1.4D 1
      GC(6)=-6.0D 0
      RETURN
  315 GC(1)=4.0D 0
      GC(2)=5.0D 0
      GC(7)=-3.0D 0
      GC(8)=9.0D 0
      RETURN
  316 GC(1)=1.0D 1
      GC(2)=-8.0D 0
      GC(7)=-1.7D 1
      GC(8)=2.0D 0
      RETURN
  317 GC(1)=-3.0D 0
      GC(2)=6.0D 0
      GC(9)=2.4D 1*(X(9)-8.0D 0)
      GC(10)=-7.0D 0
      RETURN
  318 GC(1)=-8.0D 0
      GC(2)=2.0D 0
      GC(9)=5.0D 0
      GC(10)=-2.0D 0
      RETURN
  320 K=(KC-1)*5
      DO 321 J=1,5
      GC(J+10)=2.0D 0*Y(K+J+60)
  321 CONTINUE
      GC(KC+10)=GC(KC+10)+6.0D 0*Y(KC+85)*X(KC+10)
      K=(KC-1)*10
      DO 322 I=1,10
      GC(I)=-Y(K+I)
  322 CONTINUE
      RETURN
  330 DO 331 I=1,16
      GC(I)=0.0D 0
  331 CONTINUE
      IF (KC.LT.6) THEN
      K=11+3*(KC-1)
      GC(KC)=Y(K+1)+(Y(K)+2.0D 0*X(KC)*Y(K+2))/X(KC+5)
      GC(KC+5)=-(Y(K)+Y(K+2)*X(KC))*X(KC)/X(KC+5)**2
      RETURN
      ELSE IF (KC.LT.14) THEN
      GO TO (332,333,334,335,336,337,338,339),KC-5
  332 GC(1)=Y(27)*X(12)/(X(7)*X(11))
      GC(6)=(Y(26)+Y(28)*X(12)/X(11))/X(7)
      GC(7)=-(Y(26)*X(6)+(Y(27)*X(1)+Y(28)*X(6))*X(12)/X(11))/X(7)**2
      GC(11)=-((Y(27)*X(1)+Y(28)*X(6))*X(12)/X(7))/X(11)**2
      GC(12)=(Y(27)*X(1)+Y(28)*X(6))/(X(7)*X(11))
      RETURN
  333 GC(1)=Y(33)*X(12)/X(8)
      GC(2)=Y(31)*X(13)/X(8)
      GC(7)=(Y(29)+Y(30)*X(12))/X(8)
      GC(8)=-(Y(29)*X(7)+(Y(30)*X(7)+Y(33)*X(1))*X(12)+
     * Y(31)*X(2)*X(13))/X(8)**2
      GC(12)=(Y(30)*X(7)+Y(33)*X(1))/X(8)
      GC(13)=Y(32)+Y(31)*X(2)/X(8)
      RETURN
  334 GC(2)=Y(38)*X(13)
      GC(3)=Y(36)*X(14)
      GC(8)=Y(35)*X(13)+Y(34)
      GC(9)=Y(39)*X(14)+Y(37)
      GC(13)=Y(35)*X(8)+Y(38)*X(2)
      GC(14)=Y(36)*X(3)+Y(39)*X(9)
      RETURN
  335 GC(3)=-((Y(40)+Y(43)/X(14))*X(9)+((Y(41)*X(4)+Y(44)*X(8))*X(15)+
     * Y(42)*X(10))/X(14))/X(3)**2
      GC(4)=Y(41)*X(15)/(X(3)*X(14))
      GC(8)=Y(44)*X(15)/(X(3)*X(14))
      GC(9)=(Y(40)+Y(43)/X(14))/X(3)
      GC(10)=Y(42)/(X(3)*X(14))
      GC(14)=-((Y(41)*X(4)+Y(44)*X(8))*X(15)+Y(42)*X(10)+
     * Y(43)*X(9))/(X(3)*X(14)**2)
      GC(15)=(Y(41)*X(4)+Y(44)*X(8))/(X(3)*X(14))
      RETURN
  336 GC(4)=-((Y(45)*X(5)*X(16)+Y(49)*X(10))/X(15)+Y(46)*X(10))/X(4)**2
      GC(5)=Y(45)*X(16)/(X(4)*X(15))
      GC(10)=(Y(46)+Y(49)/X(15))/X(4)
      GC(15)=-((Y(45)*X(5)*X(16)+Y(49)*X(10))/X(4)+Y(47)+
     * Y(48)*X(16))/X(15)**2
      GC(16)=(Y(45)*X(5)/X(4)+Y(48))/X(15)
      RETURN
  337 GC(4)=-(Y(50)+Y(52)*X(5)*X(16))/X(4)**2
      GC(5)=Y(52)*X(16)/X(4)
      GC(16)=Y(51)+Y(52)*X(5)/X(4)
      RETURN
  338 GC(11)=Y(53)
      GC(12)=Y(54)
      RETURN
  339 GC(11)=-Y(55)*X(12)/X(11)**2
      GC(12)=Y(55)/X(11)
      RETURN
      ELSE IF (KC.LT.18) THEN
      I=KC-13
      K=5-I
      GC(K)=Y(I+55)/X(K+1)
      GC(K+1)=-Y(I+55)*X(K)/X(K+1)**2
      RETURN
      ELSE IF (KC.EQ.18) THEN
      GC(9)=Y(60)/X(10)
      GC(10)=-Y(60)*X(9)/X(10)**2
      RETURN
      ELSE
      GC(8)=Y(61)/X(9)
      GC(9)=-Y(61)*X(8)/X(9)**2
      RETURN
      ENDIF
  340 GC(1)=0.0D 0
      GC(2)=0.0D 0
      GC(3)=0.0D 0
      GC(4)=0.0D 0
      GC(5)=0.0D 0
      GC(6)=0.0D 0
      GC(7)=0.0D 0
      GC(8)=0.0D 0
      GC(9)=0.0D 0
      GC(10)=0.0D 0
      GC(11)=0.0D 0
      GC(12)=0.0D 0
      GC(13)=0.0D 0
      GC(14)=0.0D 0
      GC(15)=0.0D 0
      GC(16)=0.0D 0
      GC(17)=0.0D 0
      GC(18)=0.0D 0
      GC(19)=0.0D 0
      GC(20)=0.0D 0
      GO TO (311,312,313,314,315,316,317,318,341,342,343,344,345,346,
     * 347,348,349),KC
  341 GC(1)=1.0D 0
      GC(2)=1.0D 0
      GC(11)=4.0D 0
      GC(12)=-2.1D 1
      RETURN
  342 GC(1)=2.0D 0*X(1)
      GC(11)=1.5D 1
      GC(12)=-8.0D 0
      RETURN
  343 GC(1)=4.0D 0
      GC(2)=9.0D 0
      GC(13)=1.0D 1*X(13)
      GC(14)=-9.0D 0
      RETURN
  344 GC(1)=3.0D 0
      GC(2)=4.0D 0
      GC(13)=6.0D 0*(X(13)-6.0D 0)
      GC(14)=-1.4D 1
      RETURN
  345 GC(1)=2.8D 1*X(1)
      GC(15)=3.5D 1
      GC(16)=-7.9D 1
      RETURN
  346 GC(2)=3.0D 1*X(2)
      GC(15)=1.1D 1
      GC(16)=-6.1D 1
      RETURN
  347 GC(1)=1.0D 1*X(1)
      GC(2)=2.0D 0
      GC(17)=3.6D 1*X(17)**3
      GC(18)=-1.0D 0
      RETURN
  348 GC(1)=2.0D 0*X(1)
      GC(2)=-1.0D 0
      GC(19)=1.9D 1
      GC(20)=-2.0D 1
      RETURN
  349 GC(1)=1.4D 1*X(1)
      GC(2)=1.0D 1*X(2)
      GC(19)=2.0D 0*X(19)
      GC(20)=-3.0D 1
      RETURN
      END
* SUBROUTINE TFFU07             ALL SYSTEMS                 90/12/01
C PORTABILITY : ALL SYSTEMS
C 90/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  VALUES OF TEST FUNCTIONS FOR NONLINEAR APPROXIMATION.
*  UNIVERSAL VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  FF  VALUE OF THE OBJECTIVE FUNCTION.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TFFU07(N,X,FF,NEXT)
      INTEGER N,NEXT
      DOUBLE PRECISION X(N),FF
      DOUBLE PRECISION X1,X2,X3,X4,X5,X6,X7,X8,A,A1,A2,A3,A4,A5,A6,A7,
     * A8,B,B1,B3,B4,B5,B6,C,P,Q
      INTEGER I,J,K
      DOUBLE PRECISION Y(128)
      COMMON /EMPR07/ Y
      GO TO(10,20,30,40,50,60,70,60,90,100,110,120,130,140,150,160,170,
     * 180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,
     * 340),NEXT
   10 FF=-X(1)-X(2)
      RETURN
   20 FF=X(1)-X(2)
      RETURN
   30 FF=0.5D 0*X(1)**2+X(2)**2-X(1)*X(2)-7.0D 0*X(1)-7.0D 0*X(2)
      RETURN
   40 FF=(X(1)-2.0D 0)**2+X(2)**2
      RETURN
   50 FF=0.0D 0
      DO 51 I=1,44
      X1=8.0D 0-Y(I)
      X2=EXP(X1*X(2))
      X3=X(1)-0.49D 0
      X4=Y(I+44)-X(1)+X2*X3
      FF=FF+X4*X4
   51 CONTINUE
      RETURN
   60 FF=100.0D 0*(X(2)-X(1)**2)**2+(1.0D 0-X(1))**2
      RETURN
   70 FF=(X(1)-10.0D 0)**3+(X(2)-20.0D 0)**3
      RETURN
   90 FF=-75.196D 0+3.8112D 0*X(1)+0.0020567D 0*X(1)**3-
     * 1.0345D-5*X(1)**4+6.8306D 0*X(2)-0.030234D 0*X(1)*X(2)+
     * 1.28134D-3*X(2)*X(1)**2+2.266D-7*X(2)*X(1)**4-
     * 0.25645D 0*X(2)**2+0.0034604D 0*X(2)**3-1.3514D-5*X(2)**4+
     * 28.106D 0/(X(2)+1.0D 0)+5.2375D-6*(X(1)*X(2))**2+
     *  6.3D-8*X(1)*(X(1)*X(2))**2-7.0D-10*(X(1)*X(2))**3-
     * 3.405D-4*X(1)*X(2)**2+1.6638D-6*X(1)*X(2)**3+
     *2.8673D 0*EXP(0.0005D 0*X(1)*X(2))-3.5256D-5*X(2)*X(1)**3
      RETURN
  100 FF=-X(1)*X(2)*X(3)
      RETURN
  110 FF=X(1)**2+X(2)**2+X(3)**2
      RETURN
  120 FF=9.0D 0*X(1)**2+X(2)**2+9.0D 0*X(3)**2
      RETURN
  130 FF=5.0D 0*X(1)+5.0D 4/X(1)+20.0D 0*X(2)+72.0D 3/X(2)+
     * 10.0D 0*X(3)+14.4D 4/X(3)
      RETURN
  140 FF=-X(1)
      RETURN
  150 X2=1.6D 0*X(1)
  151 X3=1.22D 0*X2-X(1)
      X6=(X(2)+X3)/X(1)
      X1=0.01D 0*X(1)*(112.0D 0+13.167D 0*X6-0.6667D 0*X6**2)
      IF (ABS(X1-X2).GT.0.001D 0) THEN
      X2=X1
      GO TO 151
      ENDIF
      X4=93.0D 0
  152 X5=86.35D 0+1.098D 0*X6-0.038D 0*X6**2+0.325D 0*(X4-89.0D 0)
      X8=3.0D 0*X5-133.0D 0
      X7=35.82D 0-0.222D 0*X8
      X1=98000.0D 0*X(3)/(X2*X7+1000.0D 0*X(3))
      IF (ABS(X1-X4).GT.0.001D 0) THEN
      X4=X1
      GO TO 152
      ENDIF
      FF=-(0.063D 0*X2*X5-5.04D 0*X(1)-3.36D 0*X3-0.035D 0*X(2)-
     * 10.0D 0*X(3))
      RETURN
  160 FF=0.0D 0
      B3=X(3)
      B4=1.0D 0-X(3)
      C=1.0D 0/X(4)
      B=B3+B4*X(4)
      A=B*C
      A1=LOG(A)+1.0D 0
      B1=LOG(B)+1.0D 0
      X1=1.0D 0/(12.0D 0*X(1)+1.0D 0)
      X2=1.0D 0/(12.0D 0*X(2)+1.0D 0)
      X3=SQRT(0.1591545D 0*X(1))*12.0D 0*X1
      X4=SQRT(0.1591545D 0*X(2))*12.0D 0*X2
      X5=X3*B4
      X6=X4*B3
      X7=X5*X(1)
      X8=X6*X(2)
      DO 161 I=1,19
      A5=X(1)*(Y(I+38)+A1-A*Y(I))
      A6=EXP(A5-Y(I+38))
      B5=X(2)*(Y(I+38)+B1-B*Y(I))
      B6=EXP(B5-Y(I+38))
      P=X7*A6+X8*B6
      Q=P-Y(I+19)
      FF=FF+Q*Q
  161 CONTINUE
      RETURN
  170 FF=1.0D 0+X(1)+X(2)+X(3)+X(4)
      RETURN
  180 X1=X(1)*X(1)
      X2=X(2)*X(2)
      X3=X(3)*X(3)
      X4=X(4)*X(4)
      X5=X(1)+X(1)
      X6=X(2)+X(2)
      X7=X(3)+X(3)
      X8=X(4)+X(4)
      FF=X1+X2+X3+X3+X4-5.0D 0*(X(1)+X(2))-2.1D 1*X(3)+7.0D 0*X(4)
      RETURN
  190 FF=10.0D 0*X(1)*X(4)-6.0D 0*X(3)*X(2)**2+X(2)*X(1)**3+
     * 9.0D 0*SIN(X(5)-X(3))+X(2)**3*X(4)**2*X(5)**4
      RETURN
  200 FF=5.3578547D 0*X(3)**2+0.8356891D 0*X(1)*X(5)+
     * 37.293239D 0*X(1)-40792.141D 0
      RETURN
  210 FF=24345.0D 0+8720288.849D 0*X(1)-150512.5253D 0*X(1)*X(2)+
     * 156.6950325D 0*X(1)*X(3)-476470.3222D 0*X(1)*X(4)-
     * 729482.8271D 0*X(1)*X(5)
      RETURN
  220 A1=-5.843D-7
      A2=1.17D-4
      A3=2.358D-5
      A4=1.502D-6
      A5=0.0321D 0
      A6=0.004324D 0
      A7=1.0D-4
      A8=37.48D 0
      X1=146.312D 3
      Y(1)=X(2)+X(3)+41.6D 0
      Y(18)=0.024D 0*X(4)-4.62D 0
      Y(2)=12.5D 0/Y(18)+12.0D 0
      Y(19)=0.0003535D 0*X(1)**2+0.5311D 0*X(1)+0.08705D 0*Y(2)*X(1)
      Y(20)=0.052D 0*X(1)+78.0D 0+0.002377D 0*Y(2)*X(1)
      Y(3)=Y(19)/Y(20)
      Y(4)=19.0D 0*Y(3)
      Y(21)=0.04782D 0*(X(1)-Y(3))+0.1956D 0*(X(1)-Y(3))**2/X(2)+
     * 0.6376D 0*Y(4)+1.594D 0*Y(3)
      Y(22)=100.0D 0*X(2)
      Y(23)=X(1)-Y(3)-Y(4)
      Y(24)=0.95D 0-Y(21)/Y(22)
      Y(5)=Y(23)*Y(24)
      Y(6)=X(1)-Y(5)-Y(4)-Y(3)
      Y(25)=(Y(5)+Y(4))*0.995D 0
      Y(7)=Y(25)/Y(1)
      Y(8)=Y(25)/3798.0D 0
      Y(26)=Y(7)-0.0663D 0*Y(7)/Y(8)-0.3153D 0
      Y(9)=96.82D 0/Y(26)+0.321D 0*Y(1)
      Y(10)=1.29D 0*Y(5)+1.258D 0*Y(4)+2.29D 0*Y(3)+1.71D 0*Y(6)
      Y(11)=1.71D 0*X(1)-0.452D 0*Y(4)+0.58D 0*Y(3)
      Y(27)=12.3D 0/752.3D 0
      Y(28)=1.75D 0*Y(2)*0.995D 0*X(1)
      Y(29)=0.995D 0*Y(10)+1998.0D 0
      Y(12)=Y(27)*X(1)+Y(28)/Y(29)
      Y(13)=Y(29)-1.75D 0*Y(2)
      Y(14)=3623.0D 0+64.4D 0*X(2)+58.4D 0*X(3)+X1/(Y(9)+
     * X(5))
      Y(30)=0.995D 0*Y(10)+60.8D 0*X(2)+48.0D 0*X(4)-0.1121D 0*Y(14)-
     * 5095.0D 0
      Y(15)=Y(13)/Y(30)
      Y(16)=148.0D 3-331.0D 3*Y(15)+40.0D 0*Y(13)-61.0D 0*Y(15)*Y(13)
      Y(31)=2324.0D 0*Y(10)-2.874D 7*Y(2)
      Y(17)=1.413D 7-1328.0D 0*Y(10)-531.0D 0*Y(11)+Y(31)/Y(29)
      Y(32)=Y(13)/Y(15)-Y(13)/0.52D 0
      Y(33)=1.104D 0-0.72D 0*Y(15)
      Y(34)=Y(9)+X(5)
      FF=A1*Y(17)+A2*Y(14)+A3*Y(13)+A4*Y(16)+A5*Y(12)+A6*Y(5)+
     * A7*Y(32)/Y(33)+A8*Y(2)/Y(29)+0.1365D 0
      RETURN
  230 FF=(X(1)-X(4))**2+(X(2)-X(5))**2+(X(3)-X(6))**2
      RETURN
  240 FF=0.0204D 0*X(1)*X(4)*(X(1)+X(2)+X(3))+0.0187D 0*X(2)*X(3)*
     * (X(1)+1.57D 0*X(2)+X(4))+0.0607D 0*X(1)*X(4)*X(5)**2*(X(1)+
     * X(2)+X(3))+0.0437D 0*X(2)*X(3)*X(6)**2*(X(1)+1.57D 0*X(2)+X(4))
      RETURN
  250 FF=4.3D 0*X(1)+31.8D 0*X(2)+63.3D 0*X(3)+15.8D 0*X(4)+
     * 68.5D 0*X(5)+4.7D 0*X(6)
      RETURN
  260 FF=(X(1)-1.0D 1)**2+5.0D 0*(X(2)-1.2D 1)**2+X(3)**4+3.0D 0*
     &(X(4)-1.1D 1)**2+1.0D 1*X(5)**6+7.0D 0*X(6)**2+X(7)**4-4.0D 0*
     &X(6)*X(7)-1.0D 1*X(6)-8.0D 0*X(7)
      RETURN
  270 FF=10.0D 0*X(1)*X(4)**2/(X(2)*X(6)**3*X(7)**0.25D 0)+
     * 15.0D 0*X(3)*X(4)/(X(1)*X(2)**2*X(5)*X(7)**0.5D 0)+
     * 20.0D 0*X(2)*X(6)/(X(1)**2*X(4)*X(5)**2)+
     * 25.0D 0*X(1)**2*X(2)**2*X(5)**0.5D 0*X(7)/(X(3)*X(6)**2)
      RETURN
  280 FF=0.4D 0*(X(1)/X(7))**0.67D 0+0.4D 0*(X(2)/X(8))**0.67D 0+
     * 10.0D 0-X(1)-X(2)
      RETURN
  290 FF=X(1)+X(2)+X(3)
      RETURN
  300 FF=-0.5D 0*(X(1)*X(4)-X(2)*X(3)+X(3)*X(9)-X(5)*X(9)+X(5)*X(8)-
     * X(6)*X(7))
      RETURN
  310 FF=X(1)**2+X(2)**2+X(1)*X(2)-1.4D 1*X(1)-1.6D 1*X(2)+
     &(X(3)-1.0D 1)**2+4.0D 0*(X(4)-5.0D 0)**2+(X(5)-3.0D 0)**2+
     &2.0D 0*(X(6)-1.0D 0)**2+5.0D 0*X(7)**2+7.0D 0*(X(8)-
     &1.1D 1)**2+2.0D 0*(X(9)-1.0D 1)**2+(X(10)-7.0D 0)**2+4.5D 1
      RETURN
  320 FF=0.0D 0
      DO 321 I=1,10
      FF=FF-Y(I+50)*X(I)
  321 CONTINUE
      K=0
      DO 323 I=1,5
      X1=2.0D 0*Y(I+85)*X(I+10)**2
      X2=0.0D 0
      DO 322 J=1,5
      X2=X2+Y(K+J+60)*X(J+10)
  322 CONTINUE
      FF=FF+(X1+X2)*X(I+10)
      K=K+5
  323 CONTINUE
      RETURN
  330 FF=0.0D 0
      DO 331 I=1,5
      FF=FF+Y(I)*X(I+11)+Y(I+5)*X(I)*X(I+11)
  331 CONTINUE
      RETURN
  340 FF=X(1)**2+X(2)**2+X(1)*X(2)-1.4D 1*X(1)-1.6D 1*X(2)+(X(3)-
     &1.0D 1)**2+4.0D 0*(X(4)-5.0D 0)**2+(X(5)-3.0D 0)**2+2.0D 0*
     &(X(6)-1.0D 0)**2+5.0D 0*X(7)**2+7.0D 0*(X(8)-1.1D 1)**2+
     &2.0D 0*(X(9)-1.0D 1)**2+(X(10)-7.0D 0)**2+(X(11)-9.0D 0)**2+
     &1.0D 1*(X(12)-1.0D 0)**2+5.0D 0*(X(13)-7.0D 0)**2+4.0D 0*
     &(X(14)-1.4D 1)**2+2.7D 1*(X(15)-1.0D 0)**2+X(16)**4+(X(17)-
     &2.0D 0)**2+1.3D 1*(X(18)-2.0D 0)**2+(X(19)-3.D 0)**2+X(20)**2+
     &9.5D 1
      RETURN
      END
* SUBROUTINE TFGU07             ALL SYSTEMS                 90/12/01
C PORTABILITY : ALL SYSTEMS
C 90/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  GRADIENTS OF TEST FUNCTIONS FOR NONLINEAR APPROXIMATION.
*  UNIVERSAL VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  GF(N)  GRADIENT OF THE OBJECTIVE FUNCTION.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TFGU07(N,X,GF,NEXT)
      INTEGER N,NEXT
      DOUBLE PRECISION X(N),GF(N)
      DOUBLE PRECISION X1,X2,X3,X4,X5,X6,X7,X8,A,A1,A2,A3,A4,A5,A6,A7,
     * A8,B,B1,B2,B3,B4,B5,B6,B7,C,P,Q
      INTEGER I,J,K
      DOUBLE PRECISION Y(128)
      COMMON /EMPR07/ Y
      GO TO(10,20,30,40,50,60,70,60,90,100,110,120,130,140,150,160,170,
     * 180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,
     * 340),NEXT
   10 GF(1)=-1.0D 0
      GF(2)=-1.0D 0
      RETURN
   20 GF(1)= 1.0D 0
      GF(2)=-1.0D 0
      RETURN
   30 GF(1)=X(1)-X(2)-7.0D 0
      GF(2)=2.0D 0*X(2)-X(1)-7.0D 0
      RETURN
   40 GF(1)=2.0D 0*(X(1)-2.0D 0)
      GF(2)=2.0D 0*X(2)
      RETURN
   50 GF(1)=0.0D 0
      GF(2)=0.0D 0
      DO 51 I=1,44
      X1=8.0D 0-Y(I)
      X2=EXP(X1*X(2))
      X3=X(1)-0.49D 0
      X4=Y(I+44)-X(1)+X2*X3
      GF(1)=GF(1)+X4*(X2-1.0D 0)
      GF(2)=GF(2)+X1*X2*X3*X4
   51 CONTINUE
      GF(1)=2.0D 0*GF(1)
      GF(2)=2.0D 0*GF(2)
      RETURN
   60 GF(1)=-400.0D 0*(X(2)-X(1)**2)*X(1)-2.0D 0*(1.0D 0-X(1))
      GF(2)=200.0D 0*(X(2)-X(1)**2)
      RETURN
   70 GF(1)=3.0D 0*(X(1)-10.0D 0)**2
      GF(2)=3.0D 0*(X(2)-20.0D 0)**2
      RETURN
   90 GF(1)=3.8112D 0+0.0061701D 0*X(1)**2-4.1380D-5*X(1)**3-
     * 0.030234D 0*X(2)+2.56268D-3*X(1)*X(2)+9.064D-7*X(2)*X(1)**3+
     * 10.4750D-6*X(1)*X(2)**2+18.9D-8*(X(1)*X(2))**2-
     * 21.0D-10*(X(1)*X(2))**2*X(2)-3.405D-4*X(2)**2+
     * 1.6638D-6*X(2)**3-10.5768D-5*X(2)*X(1)**2+
     * 2.8673D 0*0.0005D 0*EXP(0.0005D 0*X(1)*X(2))*X(2)
      GF(2)=6.8306D 0-0.030234D 0*X(1)+1.28134D-3*X(1)**2+
     * 2.266D-7*X(1)**4-0.51290D 0*X(2)+0.0103812D 0*X(2)**2-
     * 5.4056D-5*X(2)**3-28.106D 0/(X(2)+1.0D 0)**2+
     * 10.4750D-6*X(2)*X(1)**2+12.6D-8*X(2)*X(1)**3-
     * 21.0D-10*X(1)*(X(1)*X(2))**2-6.810D-4*X(1)*X(2)+
     * 4.9914D-6*X(1)*X(2)**2-3.5256D-5*X(1)**3+
     * 2.8673D 0*0.0005D 0*EXP(0.0005D 0*X(1)*X(2))*X(1)
      RETURN
  100 GF(1)=-X(2)*X(3)
      GF(2)=-X(1)*X(3)
      GF(3)=-X(1)*X(2)
      RETURN
  110 GF(1)=2.0D 0*X(1)
      GF(2)=2.0D 0*X(2)
      GF(3)=2.0D 0*X(3)
      RETURN
  120 GF(1)=18.0D 0*X(1)
      GF(2)= 2.0D 0*X(2)
      GF(3)=18.0D 0*X(3)
      RETURN
  130 GF(1)=5.0D 0-5.0D 4/X(1)**2
      GF(2)=2.0D 1-72.0D 3/X(2)**2
      GF(3)=1.0D 1-14.4D 4/X(3)**2
      RETURN
  140 GF(1)=-1.0D 0
      GF(2)= 0.0D 0
      GF(3)= 0.0D 0
      RETURN
  150 X2=1.6D 0*X(1)
      Y(4)=1.6D 0
      Y(5)=0.0D 0
      Y(6)=0.0D 0
  151 X3=1.22D 0*X2-X(1)
      Y(7)=1.22D 0*Y(4)-1.0D 0
      Y(8)=1.22D 0*Y(5)
      Y(9)=1.22D 0*Y(6)
      X6=(X(2)+X3)/X(1)
      Y(16)=Y(7)/X(1)-(X(2)+X3)/X(1)**2
      Y(17)=(1.0D 0+Y(8))/X(1)
      Y(18)=Y(9)/X(1)
      X1=0.01D 0*X(1)*(112.0D 0+13.167D 0*X6-0.6667D 0*X6**2)
      Y(1)=0.01D 0*(112.0D 0+13.167D 0*X6-0.6667D 0*X6**2)+
     * 0.01D 0*X(1)*(13.167*Y(16)-1.3334D 0*X6*Y(16))
      Y(2)=0.01D 0*X(1)*(13.167D 0*Y(17)-1.3334D 0*X6*Y(17))
      Y(3)=0.01D 0*X(1)*(13.167D 0*Y(18)-1.3334D 0*X6*Y(18))
      IF (ABS(X1-X2).GT.0.001D 0) THEN
      X2=X1
      Y(4)=Y(1)
      Y(5)=Y(2)
      Y(6)=Y(3)
      GO TO 151
      ENDIF
      X4=93.0D 0
      Y(10)=0.0D 0
      Y(11)=0.0D 0
      Y(12)=0.0D 0
  152 X5=86.35D 0+1.098D 0*X6-0.038D 0*X6**2+0.325D 0*(X4-89.0D 0)
      Y(13)=1.098D 0*Y(16)-0.076D 0*X6*Y(16)+0.325D 0*Y(10)
      Y(14)=1.098D 0*Y(17)-0.076D 0*X6*Y(17)+0.325D 0*Y(11)
      Y(15)=1.098D 0*Y(18)-0.076D 0*X6*Y(18)+0.325D 0*Y(12)
      X8=3.0D 0*X5-133.0D 0
      Y(22)=3.0D 0*Y(13)
      Y(23)=3.0D 0*Y(14)
      Y(24)=3.0D 0*Y(15)
      X7=35.82D 0-0.222D 0*X8
      Y(19)=-0.222D 0*Y(22)
      Y(20)=-0.222D 0*Y(23)
      Y(21)=-0.222D 0*Y(24)
      Y(89)=X2*X7+1000.0D 0*X(3)
      Y(90)=Y(89)**2
      X1=98000.0D 0*X(3)/Y(89)
      Y(1)=-98000.0D 0*X(3)*(Y(4)*X7+X2*Y(19))/Y(90)
      Y(2)=-98000.0D 0*X(3)*(Y(5)*X7+X2*Y(20))/Y(90)
      Y(3)=-98000.0D 0*X(3)*(Y(6)*X7+X2*Y(21)+1000.0D 0)/Y(90)+
     * 98000.0D 0/Y(89)
      IF (ABS(X1-X4).GT.0.001D 0) THEN
      X4=X1
      Y(10)=Y(1)
      Y(11)=Y(2)
      Y(12)=Y(3)
      GO TO 152
      ENDIF
      GF(1)=-(0.063D 0*(Y(4)*X5+X2*Y(13))-3.36D 0*Y(7)-5.04D 0)
      GF(2)=-(0.063D 0*(Y(5)*X5+X2*Y(14))-3.36D 0*Y(8)-0.35D-1)
      GF(3)=-(0.063D 0*(Y(6)*X5+X2*Y(15))-3.36D 0*Y(9)-10.0D 0)
      RETURN
  160 GF(1)=0.0D 0
      GF(2)=0.0D 0
      GF(3)=0.0D 0
      GF(4)=0.0D 0
      B3=X(3)
      B4=1.0D 0-X(3)
      C=1.0D 0/X(4)
      B=B3+B4*X(4)
      A=B*C
      A1=LOG(A)+1.0D 0
      B1=LOG(B)+1.0D 0
      X1=1.0D 0/(12.0D 0*X(1)+1.0D 0)
      X2=1.0D 0/(12.0D 0*X(2)+1.0D 0)
      X3=SQRT(0.1591545D 0*X(1))*12.0D 0*X1
      X4=SQRT(0.1591545D 0*X(2))*12.0D 0*X2
      X5=X3*B4
      X6=X4*B3
      X7=X5*X(1)
      X8=X6*X(2)
      A2=1.0D 0/A
      B2=1.0D 0/B
      B3=1.0D 0-X(4)
      A3=C*B3
      A4=-C*C*X(3)
      X1=X1+0.5D 0
      X2=X2+0.5D 0
      X3=X3*X(1)
      X4=X4*X(2)
      A5=X7*X(1)
      B5=X8*X(2)
      A3=A5*A3
      B3=B5*B3
      A4=A5*A4
      B4=B5*B4
      DO 161 I=1,19
      A5=X(1)*(Y(I+38)+A1-A*Y(I))
      A6=EXP(A5-Y(I+38))
      B5=X(2)*(Y(I+38)+B1-B*Y(I))
      B6=EXP(B5-Y(I+38))
      P=X7*A6+X8*B6
      Q=P-Y(I+19)
      A7=A6*(A2-Y(I))
      B7=B6*(B2-Y(I))
      GF(1)=GF(1)+Q*A6*X5*(X1+A5)
      GF(2)=GF(2)+Q*B6*X6*(X2+B5)
      GF(3)=GF(3)+Q*(A7*A3+B7*B3-A6*X3+B6*X4)
      GF(4)=GF(4)+Q*(A7*A4+B7*B4)
  161 CONTINUE
      GF(1)=2.0D 0*GF(1)
      GF(2)=2.0D 0*GF(2)
      GF(3)=2.0D 0*GF(3)
      GF(4)=2.0D 0*GF(4)
      RETURN
  170 GF(1)=1.0D 0
      GF(2)=1.0D 0
      GF(3)=1.0D 0
      GF(4)=1.0D 0
      RETURN
  180 X5=X(1)+X(1)
      X6=X(2)+X(2)
      X7=X(3)+X(3)
      X8=X(4)+X(4)
      GF(1)=X5-5.0D 0
      GF(2)=X6-5.0D 0
      GF(3)=X7+X7-2.1D 1
      GF(4)=X8+7.0D 0
      RETURN
  190 GF(1)=10.0D 0*X(4)+3.0D 0*X(2)*X(1)**2
      GF(2)=-12.0D 0*X(2)*X(3)+X(1)**3+
     * 3.0D 0*(X(2)*X(4))**2*X(5)**4
      GF(3)=-6.0D 0*X(2)**2-9.0D 0*COS(X(5)-X(3))
      GF(4)=10.0D 0*X(1)+2.0D 0*X(4)*X(2)**3*X(5)**4
      GF(5)=9.0D 0*COS(X(5)-X(3))+4.0D 0*X(2)**3*X(4)**2*X(5)**3
      RETURN
  200 GF(1)=0.8356891D 0*X(5)+37.293239D 0
      GF(2)=0.0D 0
      GF(3)=10.7157094D 0*X(3)
      GF(4)=0.0D 0
      GF(5)=0.8356891D 0*X(1)
      RETURN
  210 GF(1)= 8720288.849D 0-150512.5253D 0*X(2)+156.6950325D 0*X(3)-
     * 476470.3222D 0*X(4)-729482.8271D 0*X(5)
      GF(2)=-150512.5253D 0*X(1)
      GF(3)= 156.6950325D 0*X(1)
      GF(4)=-476470.3222D 0*X(1)
      GF(5)=-729482.8271D 0*X(1)
      RETURN
  220 A1=-5.843D-7
      A2=1.17D-4
      A3=2.358D-5
      A4=1.502D-6
      A5=0.0321D 0
      A6=0.004324D 0
      A7=1.0D-4
      A8=37.48D 0
      X1=146.312D 3
      Y(1)=X(2)+X(3)+41.6D 0
      Y(18)=0.024D 0*X(4)-4.62D 0
      Y(2)=12.5D 0/Y(18)+12.0D 0
      Y(19)=0.0003535D 0*X(1)**2+0.5311D 0*X(1)+0.08705D 0*Y(2)*X(1)
      Y(20)=0.052D 0*X(1)+78.0D 0+0.002377D 0*Y(2)*X(1)
      Y(3)=Y(19)/Y(20)
      Y(4)=19.0D 0*Y(3)
      Y(21)=0.04782D 0*(X(1)-Y(3))+0.1956D 0*(X(1)-Y(3))**2/X(2)+
     * 0.6376D 0*Y(4)+1.594D 0*Y(3)
      Y(22)=100.0D 0*X(2)
      Y(23)=X(1)-Y(3)-Y(4)
      Y(24)=0.95D 0-Y(21)/Y(22)
      Y(5)=Y(23)*Y(24)
      Y(6)=X(1)-Y(5)-Y(4)-Y(3)
      Y(25)=(Y(5)+Y(4))*0.995D 0
      Y(7)=Y(25)/Y(1)
      Y(8)=Y(25)/3798.0D 0
      Y(26)=Y(7)-0.0663D 0*Y(7)/Y(8)-0.3153D 0
      Y(9)=96.82D 0/Y(26)+0.321D 0*Y(1)
      Y(10)=1.29D 0*Y(5)+1.258D 0*Y(4)+2.29D 0*Y(3)+1.71D 0*Y(6)
      Y(11)=1.71D 0*X(1)-0.452D 0*Y(4)+0.58D 0*Y(3)
      Y(27)=12.3D 0/752.3D 0
      Y(28)=1.75D 0*Y(2)*0.995D 0*X(1)
      Y(29)=0.995D 0*Y(10)+1998.0D 0
      Y(12)=Y(27)*X(1)+Y(28)/Y(29)
      Y(13)=Y(29)-1.75D 0*Y(2)
      Y(14)=3623.0D 0+64.4D 0*X(2)+58.4D 0*X(3)+X1/(Y(9)+
     * X(5))
      Y(30)=0.995D 0*Y(10)+60.8D 0*X(2)+48.0D 0*X(4)-0.1121D 0*Y(14)-
     * 5095.0D 0
      Y(15)=Y(13)/Y(30)
      Y(16)=148.0D 3-331.0D 3*Y(15)+40.0D 0*Y(13)-61.0D 0*Y(15)*Y(13)
      Y(31)=2324.0D 0*Y(10)-2.874D 7*Y(2)
      Y(17)=1.413D 7-1328.0D 0*Y(10)-531.0D 0*Y(11)+Y(31)/Y(29)
      Y(32)=Y(13)/Y(15)-Y(13)/0.52D 0
      Y(33)=1.104D 0-0.72D 0*Y(15)
      Y(34)=Y(9)+X(5)
      Y(35)=1.0D 0
      Y(36)=0.024D 0
      Y(37)=-12.5D 0*Y(36)/Y(18)**2
      Y(38)=0.000707D 0*X(1)+0.5311D 0+0.08705D 0*Y(2)
      Y(39)=0.08705D 0*X(1)*Y(37)
      Y(40)=0.052D 0+0.002377D 0*Y(2)
      Y(41)=0.002377D 0*X(1)*Y(37)
      Y(42)=(Y(38)*Y(20)-Y(19)*Y(40))/Y(20)**2
      Y(43)=(Y(39)*Y(20)-Y(19)*Y(41))/Y(20)**2
      Y(44)=19.0D 0*Y(42)
      Y(45)=19.0D 0*Y(43)
      Y(46)=0.04782D 0*(1.0D 0-Y(42))+0.3912D 0*(X(1)-Y(3))*(1.0D 0-
     * Y(42))/X(2)+0.6376D 0*Y(44)+1.594D 0*Y(42)
      Y(47)=-0.1956D 0*(X(1)-Y(3))**2/X(2)**2
      Y(48)=-0.04782D 0*Y(43)-0.3912D 0*(X(1)-Y(3))*Y(43)/X(2)+
     * 0.6376D 0*Y(45)+1.594D 0*Y(43)
      Y(49)=100.0D 0
      Y(50)=1.0D 0-Y(42)-Y(44)
      Y(51)=-Y(43)-Y(45)
      Y(52)=-Y(46)/Y(22)
      Y(53)=(Y(21)*Y(49)-Y(47)*Y(22))/Y(22)**2
      Y(54)=-Y(48)/Y(22)
      Y(55)=Y(50)*Y(24)+Y(23)*Y(52)
      Y(56)=Y(23)*Y(53)
      Y(57)=Y(51)*Y(24)+Y(23)*Y(54)
      Y(58)=1.0D 0-Y(55)-Y(44)-Y(42)
      Y(59)=-Y(56)
      Y(60)=-Y(57)-Y(45)-Y(43)
      Y(61)=0.995D 0*(Y(55)+Y(44))
      Y(62)=0.995D 0*Y(56)
      Y(63)=0.995D 0*(Y(57)+Y(45))
      Y(64)=Y(61)/Y(1)
      Y(65)=(Y(62)*Y(1)-Y(25)*Y(35))/Y(1)**2
      Y(66)=-Y(25)*Y(35)/Y(1)**2
      Y(67)=Y(63)/Y(1)
      Y(68)=Y(61)/3798.0D 0
      Y(69)=Y(62)/3798.0D 0
      Y(70)=Y(63)/3798.0D 0
      Y(71)=Y(64)-0.0663D 0*(Y(64)*Y(8)-Y(7)*Y(68))/Y(8)**2
      Y(72)=Y(65)-0.0663D 0*(Y(65)*Y(8)-Y(7)*Y(69))/Y(8)**2
      Y(73)=Y(66)-0.0663D 0*Y(66)/Y(8)
      Y(74)=Y(67)-0.0663D 0*(Y(67)*Y(8)-Y(7)*Y(70))/Y(8)**2
      Y(75)=-96.82D 0*Y(71)/Y(26)**2
      Y(76)=-96.82D 0*Y(72)/Y(26)**2+0.321D 0*Y(35)
      Y(77)=-96.82D 0*Y(73)/Y(26)**2+0.321D 0*Y(35)
      Y(78)=-96.82D 0*Y(74)/Y(26)**2
      Y(79)=1.29D 0*Y(55)+1.258D 0*Y(44)+2.29D 0*Y(42)+1.71D 0*Y(58)
      Y(80)=1.29D 0*Y(56)+1.71D 0*Y(59)
      Y(81)=1.29D 0*Y(57)+1.258D 0*Y(45)+2.29D 0*Y(43)+1.71D 0*Y(60)
      Y(82)=1.71D 0-0.452D 0*Y(44)+0.58D 0*Y(42)
      Y(83)=-0.452D 0*Y(45)+0.58D 0*Y(43)
      Y(84)=1.75D 0*Y(2)*0.995D 0
      Y(85)=1.75D 0*Y(37)*0.995D 0*X(1)
      Y(86)=0.995D 0*Y(79)
      Y(87)=0.995D 0*Y(80)
      Y(88)=0.995D 0*Y(81)
      Y(89)=Y(27)+(Y(84)*Y(29)-Y(28)*Y(86))/Y(29)**2
      Y(90)=-Y(28)*Y(87)/Y(29)**2
      Y(91)=(Y(85)*Y(29)-Y(28)*Y(88))/Y(29)**2
      Y(92)=Y(88)-1.75D 0*Y(37)
      Y(93)=-X1*Y(75)/(Y(9)+X(5))**2
      Y(94)=64.4D 0-X1*Y(76)/(Y(9)+X(5))**2
      Y(95)=58.4D 0-X1*Y(77)/(Y(9)+X(5))**2
      Y(96)=-X1*Y(78)/(Y(9)+X(5))**2
      Y(97)=-X1/(Y(9)+X(5))**2
      Y(98)=0.995D 0*Y(79)-0.1121D 0*Y(93)
      Y(99)=0.995D 0*Y(80)+60.8D 0-0.1121D 0*Y(94)
      Y(100)=-0.1121D 0*Y(95)
      Y(101)=0.995D 0*Y(81)+48.0D 0-0.1121D 0*Y(96)
      Y(102)=-0.1121D 0*Y(97)
      Y(103)=(Y(86)*Y(30)-Y(13)*Y(98))/Y(30)**2
      Y(104)=(Y(87)*Y(30)-Y(13)*Y(99))/Y(30)**2
      Y(105)=-Y(13)*Y(100)/Y(30)**2
      Y(106)=(Y(92)*Y(30)-Y(13)*Y(101))/Y(30)**2
      Y(107)=-Y(13)*Y(102)/Y(30)**2
      Y(108)=-3.31D 5*Y(103)+40.0D 0*Y(86)-61.0D 0*(Y(103)*Y(13)+
     * Y(15)*Y(86))
      Y(109)=-3.31D 5*Y(104)+40.0D 0*Y(87)-61.0D 0*(Y(104)*Y(13)+
     * Y(15)*Y(87))
      Y(110)=-3.31D 5*Y(105)-61.0D 0*Y(105)*Y(13)
      Y(111)=-3.31D 5*Y(106)+40.0D 0*Y(92)-61.0D 0*(Y(106)*Y(13)+
     * Y(15)*Y(92))
      Y(112)=-3.31D 5*Y(107)-61.0D 0*Y(107)*Y(13)
      Y(113)=2.324D 3*Y(79)
      Y(114)=2.324D 3*Y(80)
      Y(115)=2.324D 3*Y(81)-2.874D 7*Y(37)
      Y(116)=-1.328D 3*Y(79)-531.0D 0*Y(82)+(Y(113)*Y(29)-Y(31)*
     * Y(86))/Y(29)**2
      Y(117)=-1.328D 3*Y(80)+(Y(114)*Y(29)-Y(31)*Y(87))/Y(29)**2
      Y(118)=-1.328D 3*Y(81)-531.0D 0*Y(83)+(Y(115)*Y(29)-Y(31)*
     * Y(88))/Y(29)**2
      Y(119)=(Y(86)*Y(15)-Y(13)*Y(103))/Y(15)**2-Y(86)/0.52D 0
      Y(120)=(Y(87)*Y(15)-Y(13)*Y(104))/Y(15)**2-Y(87)/0.52D 0
      Y(121)=-Y(13)*Y(105)/Y(15)**2
      Y(122)=(Y(92)*Y(15)-Y(13)*Y(106))/Y(15)**2-Y(92)/0.52D 0
      Y(123)=-Y(13)*Y(107)/Y(15)**2
      Y(124)=-0.72D 0*Y(103)
      Y(125)=-0.72D 0*Y(104)
      Y(126)=-0.72D 0*Y(105)
      Y(127)=-0.72D 0*Y(106)
      Y(128)=-0.72D 0*Y(107)
      GF(1)=A1*Y(116)+A2*Y(93)+A3*Y(86)+A4*Y(108)+A5*Y(89)+A6*Y(55)+
     * A7*(Y(119)*Y(33)-Y(32)*Y(124))/Y(33)**2-A8*Y(2)*Y(86)/Y(29)**2
      GF(2)=A1*Y(117)+A2*Y(94)+A3*Y(87)+A4*Y(109)+A5*Y(90)+A6*Y(56)+
     * A7*(Y(120)*Y(33)-Y(32)*Y(125))/Y(33)**2-A8*Y(2)*Y(87)/Y(29)**2
      GF(3)=A2*Y(95)+A4*Y(110)+A7*(Y(121)*Y(33)-Y(32)*Y(126))/Y(33)**2
      GF(4)=A1*Y(118)+A2*Y(96)+A3*Y(92)+A4*Y(111)+A5*Y(91)+A6*Y(57)+
     * A7*(Y(122)*Y(33)-Y(32)*Y(127))/Y(33)**2+A8*(Y(37)*Y(29)-Y(2)*
     * Y(88))/Y(29)**2
      GF(5)=A2*Y(97)+A4*Y(112)+A7*(Y(123)*Y(33)-Y(32)*Y(128))/Y(33)**2
      RETURN
  230 GF(1)= 2.0D 0*(X(1)-X(4))
      GF(2)= 2.0D 0*(X(2)-X(5))
      GF(3)= 2.0D 0*(X(3)-X(6))
      GF(4)=-2.0D 0*(X(1)-X(4))
      GF(5)=-2.0D 0*(X(2)-X(5))
      GF(6)=-2.0D 0*(X(3)-X(6))
      RETURN
  240 GF(1)=0.0408D 0*X(1)*X(4)+0.0204D 0*X(4)*(X(2)+X(3))+
     * 0.0187D 0*X(2)*X(3)+0.1214D 0*X(1)*X(4)*X(5)**2+
     * 0.0607D 0*X(4)*X(5)**2*(X(2)+X(3))+0.0437D 0*X(2)*X(3)*X(6)**2
      GF(2)=0.0204D 0*X(1)*X(4)+0.0374D 0*1.57D 0*X(2)*X(3)+
     * 0.0187D 0*X(3)*(X(1)+X(4))+0.0607D 0*X(1)*X(4)*X(5)**2+
     * 0.0874D 0*1.57D 0*X(2)*X(3)*X(6)**2+
     * 0.0437D 0*X(3)*X(6)**2*(X(1)+X(4))
      GF(3)=0.0204D 0*X(1)*X(4)+
     * 0.0187D 0*X(2)*(X(1)+1.57D 0*X(2)+X(4))+
     * 0.0607D 0*X(1)*X(4)*X(5)**2+
     * 0.0437D 0*X(2)*X(6)**2*(X(1)+1.57D 0*X(2)+X(4))
      GF(4)=0.0204D 0*X(1)*(X(1)+X(2)+X(3))+
     * 0.0187D 0*X(2)*X(3)+
     * 0.0607D 0*X(1)*X(5)**2*(X(1)+X(2)+X(3))+
     * 0.0437D 0*X(2)*X(3)*X(6)**2
      GF(5)=0.1214D 0*X(1)*X(4)*X(5)*(X(1)+X(2)+X(3))
      GF(6)=0.0874D 0*X(2)*X(3)*X(6)*(X(1)+1.57D 0*X(2)+X(4))
      RETURN
  250 GF(1)= 4.3D 0
      GF(2)=31.8D 0
      GF(3)=63.3D 0
      GF(4)=15.8D 0
      GF(5)=68.5D 0
      GF(6)= 4.7D 0
      RETURN
  260 GF(1)=2.0D 0*(X(1)-1.0D 1)
      GF(2)=1.0D 1*(X(2)-1.2D 1)
      GF(3)=4.0D 0*X(3)**3
      GF(4)=6.0D 0*(X(4)-1.1D 1)
      GF(5)=6.0D 1*X(5)**5
      GF(6)=1.4D 1*X(6)-4.0D 0*X(7)-1.0D 1
      GF(7)=4.0D 0*X(7)**3-4.0D 0*X(6)-8.0D 0
      RETURN
  270 GF(1)=10.0D 0*X(4)**2/(X(2)*X(6)**3*X(7)**0.25D 0)-
     * 15.0D 0*X(3)*X(4)/((X(1)*X(2))**2*X(5)*X(7)**0.5D 0)-
     * 40.0D 0*X(2)*X(6)/(X(1)**3*X(4)*X(5)**2)+
     * 50.0D 0*X(1)*X(2)**2*X(5)**0.5D 0*X(7)/(X(3)*X(6)**2)
      GF(2)=-10.0D 0*X(1)*X(4)**2/(X(2)**2*X(6)**3*X(7)**0.25D 0)-
     * 30.0D 0*X(3)*X(4)/(X(1)*X(2)**3*X(5)*X(7)**0.5D 0)+
     * 20.0D 0*X(6)/(X(1)**2*X(4)*X(5)**2)+
     * 50.0D 0*X(1)**2*X(2)*X(5)**0.5D 0*X(7)/(X(3)*X(6)**2)
      GF(3)=15.0D 0*X(4)/(X(1)*X(2)**2*X(5)*X(7)**0.5D 0)-
     * 25.0D 0*X(1)**2*X(2)**2*X(5)**0.5D 0*X(7)/(X(3)*X(6))**2
      GF(4)=20.0D 0*X(1)*X(4)/(X(2)*X(6)**3*X(7)**0.25D 0)+
     * 15.0D 0*X(3)/(X(1)*X(2)**2*X(5)*X(7)**0.5D 0)-
     * 20.0D 0*X(2)*X(6)/(X(1)*X(4)*X(5))**2
      GF(5)=-15.0D 0*X(3)*X(4)/(X(1)*(X(2)*X(5))**2*X(7)**0.5D 0)-
     * 40.0D 0*X(2)*X(6)/(X(1)**2*X(4)*X(5)**3)+
     * 12.5D 0*X(1)**2*X(2)**2*X(7)/(X(3)*X(5)**0.5D 0*X(6)**2)
      GF(6)=-30.0D 0*X(1)*X(4)**2/(X(2)*X(6)**4*X(7)**0.25D 0)+
     * 20.0D 0*X(2)/(X(1)**2*X(4)*X(5)**2)-
     * 50.0D 0*X(1)**2*X(2)**2*X(5)**0.5D 0*X(7)/(X(3)*X(6)**3)
      GF(7)=-2.5D 0*X(1)*X(4)**2/(X(2)*X(6)**3*X(7)**1.25D 0)-
     * 7.5D 0*X(3)*X(4)/(X(1)*X(2)**2*X(5)*X(7)**1.5D 0)+
     * 25.0D 0*X(1)**2*X(2)**2*X(5)**0.5D 0/(X(3)*X(6)**2)
      RETURN
  280 GF(1)=0.268D 0*(X(7)/X(1))**0.33D 0/X(7)-1.0D 0
      GF(2)=0.268D 0*(X(8)/X(2))**0.33D 0/X(8)-1.0D 0
      GF(3)=0.0D 0
      GF(4)=0.0D 0
      GF(5)=0.0D 0
      GF(6)=0.0D 0
      GF(7)=-0.268D 0*(X(1)/X(7))**0.67D 0/X(7)
      GF(8)=-0.268D 0*(X(2)/X(8))**0.67D 0/X(8)
      RETURN
  290 GF(1)=1.0D 0
      GF(2)=1.0D 0
      GF(3)=1.0D 0
      GF(4)=0.0D 0
      GF(5)=0.0D 0
      GF(6)=0.0D 0
      GF(7)=0.0D 0
      GF(8)=0.0D 0
      RETURN
  300 GF(1)=-0.5D 0*X(4)
      GF(2)= 0.5D 0*X(3)
      GF(3)=-0.5D 0*(X(9)-X(2))
      GF(4)=-0.5D 0*X(1)
      GF(5)=-0.5D 0*(X(8)-X(9))
      GF(6)= 0.5D 0*X(7)
      GF(7)= 0.5D 0*X(6)
      GF(8)=-0.5D 0*X(5)
      GF(9)=-0.5D 0*(X(3)-X(5))
      RETURN
  310 GF(1)=2.0D 0*X(1)+X(2)-1.4D 1
      GF(2)=2.0D 0*X(2)+X(1)-1.6D 1
      GF(3)=2.0D 0*(X(3)-1.0D 1)
      GF(4)=8.0D 0*(X(4)-5.0D 0)
      GF(5)=2.0D 0*(X(5)-3.0D 0)
      GF(6)=4.0D 0*(X(6)-1.0D 0)
      GF(7)=1.0D 1*X(7)
      GF(8)=1.4D 1*(X(8)-1.1D 1)
      GF(9)=4.0D 0*(X(9)-1.0D 1)
      GF(10)=2.0D 0*(X(10)-7.0D 0)
      RETURN
  320 DO 321 I=1,10
      GF(I)=-Y(I+50)
  321 CONTINUE
      K=0
      DO 323 I=1,5
      X1=2.0D 0*Y(I+85)*X(I+10)**2
      X2=0.0D 0
      DO 322 J=1,5
      X2=X2+Y(K+J+60)*X(J+10)
  322 CONTINUE
      GF(I+10)=3.0D 0*X1+2.0D 0*X2
      K=K+5
  323 CONTINUE
      RETURN
  330 DO 331 I=1,5
      GF(I)=Y(I+5)*X(I+11)
      GF(I+5)=0.0D 0
      GF(I+11)=Y(I)+Y(I+5)*X(I)
  331 CONTINUE
      GF(11)=0.0D 0
      RETURN
  340 GF(1)=2.0D 0*X(1)+X(2)-1.4D 1
      GF(2)=2.0D 0*X(2)+X(1)-1.6D 1
      GF(3)=2.0D 0*(X(3)-1.0D 1)
      GF(4)=8.0D 0*(X(4)-5.0D 0)
      GF(5)=2.0D 0*(X(5)-3.0D 0)
      GF(6)=4.0D0*(X(6)-1.0D 0)
      GF(7)=1.0D 1*X(7)
      GF(8)=1.4D 1*(X(8)-1.1D 1)
      GF(9)=4.0D 0*(X(9)-1.0D 1)
      GF(10)=2.0D 0*(X(10)-7.0D0)
      GF(11)=2.0D 0*(X(11)-9.0D0)
      GF(12)=2.0D 1*(X(12)-1.0D 0)
      GF(13)=1.0D 1*(X(13)-7.0D 0)
      GF(14)=8.0D 0*(X(14)-1.4D 1)
      GF(15)=5.4D 1*(X(15)-1.0D 0)
      GF(16)=4.0D 0*X(16)**3
      GF(17)=2.0D 0*(X(17)-2.0D 0)
      GF(18)=2.6D 1*(X(18)-2.0D 0)
      GF(19)=2.0D 0*(X(19)-3.0D 0)
      GF(20)=2.0D 0*X(20)
      RETURN
      END
* SUBROUTINE TIUD15                ALL SYSTEMS                00/12/01
C PORTABILITY : ALL SYSTEMS
C 92/12/01 RA : ORIGINAL VERSION
*
* PURPOSE :
*  INITIAL VALUES OF THE VARIABLES AND STRUCTURE OF THE SPARSE HESSIAN
*  MATRIX FOR UNCONSTRAINED MINIMIZATION.
*  SPARSE VERSION WITH CHANGED TESTS 7-10, 12.
*  CHANGED FOR THE TESTS.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  NB  NUMBER OF ELEMENTS OF THE SPARSE MATRIX.
*  RO  X(N)  VECTOR OF VARIABLES.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*  IO  IERR  ERROR INDICATOR.
*
      SUBROUTINE TIUD15(N,NB,X,XMAX,NEXT,IERR)
      INTEGER N,NB,NEXT,IERR
      DOUBLE PRECISION X(N),XMAX
      INTEGER I
      DOUBLE PRECISION Y(20)
      COMMON /EMPR15/ Y
      REAL*8 ETA9
      PARAMETER (ETA9=1.0D 60)
      XMAX=1.0D 3
      IERR=0
      GOTO(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,
     & 170,180,190,200,210,220,260,270),NEXT
   10 IF (N.LT.2) GO TO 999
      N=N-MOD(N,2)
      DO 11 I=1,N
        IF(MOD(I,2).EQ.1) THEN
          X(I)=-1.2D 0
        ELSE
          X(I)=1.0D 0
        ENDIF
   11 CONTINUE
      NB=2*(N-1)
      RETURN
   20 IF (N.LT.4) GO TO 999
      N=N-MOD(N,2)
      DO 21 I=1,N
        IF(MOD(I,2).EQ.1) THEN
          X(I)=-2.0D 0
          IF(I.LE.4) X(I)=-3.0D 0
        ELSE
          X(I)=0.0D 0
          IF(I.LE.4) X(I)=-1.0D 0
        ENDIF
   21 CONTINUE
      NB=3*(N-2)
      RETURN
   30 IF (N.LT.4) GO TO 999
      N=N-MOD(N,2)
      DO 31 I=1,N
        IF(MOD(I,4).EQ.1) THEN
          X(I)=3.0D 0
        ELSEIF(MOD(I,4).EQ.2) THEN
          X(I)=-1.0D 0
        ELSEIF(MOD(I,4).EQ.3) THEN
          X(I)=0.0D 0
        ELSE
          X(I)=1.0D 0
        ENDIF
   31 CONTINUE
      NB=2*(N-2)
      RETURN
   40 IF (N.LT.4) GO TO 999
      N=N-MOD(N,2)
      DO 41 I=1,N
        X(I)=2.0D 0
   41 CONTINUE
      X(1)=1.0D 0
      NB=5*(N-2)/2
      XMAX=1.0D 1
      RETURN
   50 IF (N.LT.3) GO TO 999
      DO 51 I=1,N
        X(I)=-1.0D 0
   51 CONTINUE
      NB=N
      RETURN
   60 IF (N.LT.6) GO TO 999
      DO 61 I=1,N
        X(I)=-1.0D 0
   61 CONTINUE
      NB=N
      RETURN
   70 IF (N.LT.2) GO TO 999
      DO 71 I=1,N-1
        X(I)=0.5D 0
   71 CONTINUE
      X(N)=-2.0D 0
      NB=2*(N-1)
      RETURN
   80 IF (N.LT.4) GO TO 999
      N=N-MOD(N,4)
      DO 81 I=1,N
        X(I)=SIN(DBLE(I))**2
   81 CONTINUE
      NB=5*N
      RETURN
   90 IF (N.LT.4) GO TO 999
      N=N-MOD(N,2)
      DO 91 I=1,N
        X(I)=5.0D 0
   91 CONTINUE
      NB=3*(N-2)
      RETURN
  100 IF (N.LT.2) GO TO 999
      DO 101 I=1,N
        X(I)=0.2D 0
  101 CONTINUE
      NB=2*(N-1)
      RETURN
  110 CONTINUE
      IF (N.LT.2) GO TO 999
      N=N-MOD(N,2)
      DO 111 I=1,N
        IF(MOD(I,2).EQ.1) THEN
          X(I)=-0.8D 0
        ELSE
          X(I)=-0.8D 0
        ENDIF
  111 CONTINUE
      NB=2*(N-1)
      RETURN
  120 CONTINUE
      IF (N.LT.5) GO TO 999
      IF (MOD(N-5,3).NE.0) N=N-MOD(N-5,3)
      DO 121 I=1,N
        X(I)=-1.0D 0
  121 CONTINUE
      NB=6*((N-5)/3+1)
      RETURN
  130 CONTINUE
      IF (N.LT.5) GO TO 999
      IF (MOD(N-5,3).NE.0) N=N-MOD(N-5,3)
      DO 131 I=1,N
        X(I)=-1.0D 0
  131 CONTINUE
      NB=7*((N-5)/3+1)
      RETURN
  140 CONTINUE
      IF (N.LT.4) GO TO 999
      DO 141 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=-0.8D 0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)= 1.2D 0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=-1.2D 0
        ELSE
          X(I)= 0.8D 0
        ENDIF
  141 CONTINUE
      Y(1)=14.4D 0
      Y(2)=6.8D 0
      Y(3)=4.2D 0
      Y(4)=3.2D 0
  142 IF (MOD(N-4,2).NE.0) N=N-MOD(N-4,2)
      NB=4*((N-4)/2+1)
      RETURN
  150 CONTINUE
      IF (N.LT.4) GO TO 999
      DO 151 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=-0.8D 0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)= 1.2D 0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=-1.2D 0
        ELSE
          X(I)= 0.8D 0
        ENDIF
  151 CONTINUE
      Y(1)=35.8D 0
      Y(2)=11.2D 0
      Y(3)=6.2D 0
      Y(4)=4.4D 0
      GO TO 142
  160 CONTINUE
      IF (N.LT.4) GO TO 999
      DO 161 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=-0.8D 0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)= 1.2D 0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=-1.2D 0
        ELSE
          X(I)= 0.8D 0
        ENDIF
  161 CONTINUE
      Y(1)=30.6D 0
      Y(2)=72.2D 0
      Y(3)=124.4D 0
      Y(4)=187.4D 0
      GO TO 142
  170 IF (N.LT.4) GO TO 999
      N=N-MOD(N,2)
      NB=N
      DO 171 I=1,N
        IF (MOD(I,8).EQ.1) X(I)=1.0D-1
        IF (MOD(I,8).EQ.2.OR.MOD(I,8).EQ.0) X(I)=2.0D-1
        IF (MOD(I,8).EQ.3.OR.MOD(I,8).EQ.7) X(I)=3.0D-1
        IF (MOD(I,8).EQ.4.OR.MOD(I,8).EQ.6) X(I)=4.0D-1
        IF (MOD(I,8).EQ.5) X(I)=5.0D-1
  171 CONTINUE
      RETURN
  180 IF (N.LT.3) GO TO 999
      NB=N
      DO 181 I=1,N
        X(I)=1.2D 1
  181 CONTINUE
      RETURN
  190 IF (N.LT.7) GO TO 999
      NB=N
      DO 191 I=1,N
        X(I)=-1.0D 0
  191 CONTINUE
      RETURN
  200 IF (N.LT.3) GO TO 999
      NB=N
      DO 201 I=1,N
        X(I)=DBLE(I)/DBLE(N+1)
        X(I)=X(I)*(X(I)-1.0D 0)
  201 CONTINUE
      RETURN
  210 CONTINUE
      IF (N.LT.5) GO TO 999
      IF (MOD(N-5,3).NE.0) N=N-MOD(N-5,3)
      DO 211 I=1,N
        X(I)=-1.0D 0
  211 CONTINUE
      NB=7*((N-5)/3+1)
      RETURN
  220 IF (N.LT.3) GO TO 999
      DO 221 I=1,N
        IF(MOD(I,2).EQ.1) THEN
          X(I)=-1.2D 0
        ELSE
          X(I)= 1.0D 0
        ENDIF
  221 CONTINUE
      NB=2*(N-1)
      RETURN
  260 IF (N.LT.4) GO TO 999
      N=N-MOD(N,4)
      DO 261 I=1,N
        X(I)=1.0D 0
  261 CONTINUE
      NB=5*N
      RETURN
  270 IF (N.LT.7) GO TO 999
      DO 271 I=1,N
        X(I)=5.0D 0
  271 CONTINUE
      Y(1)=SIN(1.0D 0)
      NB=13*(N-6)
      RETURN
  999 IERR=1
      RETURN
      END
* SUBROUTINE TAFU15             ALL SYSTEMS                92/12/01
C PORTABILITY : ALL SYSTEMS
C 92/12/01 RA : ORIGINAL VERSION
*
* PURPOSE :
*  VALUES OF PARTIAL FUNCTIONS IN THE SUM OF SQUARES.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  KA  INDEX OF THE GIVEN PARTIAL FUNCTION.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  FA  VALUE OF THE KA-TH PARTIAL FUNCTION AT THE POINT X.
*  II  NEXT  NUMBER OF THE SELECTED TEST PROBLEM.
*
      SUBROUTINE TAFU15(N,KA,X,FA,NEXT)
      INTEGER N,KA,NEXT
      DOUBLE PRECISION X(*),FA
      DOUBLE PRECISION A,C,D,P,ALFA,U,V
      INTEGER I,J,K,L,M,IA,IB,IC
      DOUBLE PRECISION Y(20)
      COMMON /EMPR15/ Y
      GO TO(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,
     & 170,180,190,200,210,220),NEXT
   10 I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
      FA=1.0D 1*(X(I)**2-X(I+1))
      ELSE
      FA=X(I)-1.0D 0
      END IF
      RETURN
   20 A=SQRT(10.0D 0)
      D=SQRT(90.0D 0)
      I=2*((KA+5)/6)
      IF (MOD(KA,6).EQ.1) THEN
      FA=1.0D 1*(X(I-1)**2-X(I))
      ELSE IF (MOD(KA,6).EQ.2) THEN
      FA=X(I-1)-1.0D 0
      ELSE IF (MOD(KA,6).EQ.3) THEN
      FA=D*(X(I+1)**2-X(I+2))
      ELSE IF (MOD(KA,6).EQ.4) THEN
      FA=X(I+1)-1.0D 0
      ELSE IF (MOD(KA,6).EQ.5) THEN
      FA=A*(X(I)+X(I+2)-2.0D 0)
      ELSE
      FA=(X(I)-X(I+2))/A
      END IF
      RETURN
   30 A=SQRT(1.0D 1)
      C=SQRT(5.0D 0)
      I=2*((KA+3)/4)
      IF (MOD(KA,4).EQ.1) THEN
      FA=X(I-1)+1.0D 1*X(I)
      ELSE IF (MOD(KA,4).EQ.2) THEN
      FA=C*(X(I+1)-X(I+2))
      ELSE IF (MOD(KA,4).EQ.3) THEN
      FA=(X(I)-2.0D 0*X(I+1))**2
      ELSE
      FA=A*(X(I-1)-X(I+2))**2
      END IF
      RETURN
   40 I=2*((KA+4)/5)
      IF (MOD(KA,5).EQ.1) THEN
      FA=(EXP(X(I-1))-X(I))**2
      ELSE IF (MOD(KA,5).EQ.2) THEN
      FA=1.0D 1*(X(I)-X(I+1))**3
      ELSE IF (MOD(KA,5).EQ.3) THEN
      P=X(I+1)-X(I+2)
      FA=(SIN(P)/COS(P))**2
      ELSE IF (MOD(KA,5).EQ.4) THEN
      FA=X(I-1)**4
      ELSE
      FA=X(I+2)-1.0D 0
      END IF
      RETURN
   50 I=KA
      FA=(3.0D 0-2.0D 0*X(I))*X(I)+1.0D 0
      IF (I.GT.1) FA=FA-X(I-1)
      IF (I.LT.N) FA=FA-X(I+1)
      RETURN
   60 I=KA
      FA=(2.0D 0+5.0D 0*X(I)**2)*X(I)+1.0D 0
      DO 61 J=MAX(1,I-5),MIN(N,I+1)
      FA=FA+X(J)*(1.0D 0+X(J))
   61 CONTINUE
      RETURN
   70 I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
      FA=X(I)+X(I+1)*((5.0D 0-X(I+1))*X(I+1)-2.0D 0)-1.3D 1
      ELSE
      FA=X(I)+X(I+1)*((1.0D 0+X(I+1))*X(I+1)-1.4D 1)-2.9D 1
      END IF
      RETURN
   80 I=MOD(KA,N/2)+1
      J=I+N/2
      M=5*N
      IF (KA.LE.M/2) THEN
      IA=1
      ELSE
      IA=2
      END IF
      IB=5-KA/(M/4)
      IC=MOD(KA,5)+1
      FA=(X(I)**IA-X(J)**IB)**IC
      RETURN
   90 I=2*((KA+5)/6)-1
      IF (MOD(KA,6).EQ.1) THEN
      FA=X(I)+3.0D 0*X(I+1)*(X(I+2)-1.0D 0)+X(I+3)**2-1.0D 0
      ELSE IF (MOD(KA,6).EQ.2) THEN
      FA=(X(I)+X(I+1))**2+(X(I+2)-1.0D 0)**2-X(I+3)-3.0D 0
      ELSE IF (MOD(KA,6).EQ.3) THEN
      FA=X(I)*X(I+1)-X(I+2)*X(I+3)
      ELSE IF (MOD(KA,6).EQ.4) THEN
      FA=2.0D 0*X(I)*X(I+2)+X(I+1)*X(I+3)-3.0D 0
      ELSE IF (MOD(KA,6).EQ.5) THEN
      FA=(X(I)+X(I+1)+X(I+2)+X(I+3))**2+(X(I)-1.0D 0)**2
      ELSE
      FA=X(I)*X(I+1)*X(I+2)*X(I+3)+(X(I+3)-1.0D 0)**2-1.0D 0
      END IF
      RETURN
  100 I=(KA+1)/2
      J=MOD(KA,2)
      IF (J.EQ.0) THEN
      FA=6.0D 0-EXP(2.0D 0*X(I))-EXP(2.0D 0*X(I+1))
      ELSE IF (I.EQ.1) THEN
      FA=4.0D 0-EXP(X(I))-EXP(X(I+1))
      ELSE IF (I.EQ.N) THEN
      FA=8.0D 0-EXP(3.0D 0*X(I-1))-EXP(3.0D 0*X(I))
      ELSE
      FA=8.0D 0-EXP(3.0D 0*X(I-1))-EXP(3.0D 0*X(I))+
     & 4.0D 0-EXP(X(I))-EXP(X(I+1))
      END IF
      RETURN
  110 I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
      FA=1.0D 1*(2.0D 0*X(I)/(1.0D 0+X(I)**2)-X(I+1))
      ELSE
      FA=X(I)-1.0D 0
      END IF
      RETURN
  120 I=3*((KA+5)/6)-2
      IF (MOD(KA,6).EQ.1) THEN
      FA=1.0D 1*(X(I)**2-X(I+1))
      ELSE IF (MOD(KA,6).EQ.2) THEN
      FA=X(I+2)-1.0D 0
      ELSE IF (MOD(KA,6).EQ.3) THEN
      FA=(X(I+3)-1.0D 0)**2
      ELSE IF (MOD(KA,6).EQ.4) THEN
      FA=(X(I+4)-1.0D 0)**3
      ELSE IF (MOD(KA,6).EQ.5) THEN
      FA=X(I)**2*X(I+3)+SIN(X(I+3)-X(I+4))-1.0D 1
      ELSE
      FA=X(I+1)+(X(I+2)**2*X(I+3))**2-2.0D 1
      END IF
      RETURN
  130 I=3*((KA+6)/7)-2
      IF (MOD(KA,7).EQ.1) THEN
      FA=1.0D 1*(X(I)**2-X(I+1))
      ELSE IF (MOD(KA,7).EQ.2) THEN
      FA=1.0D 1*(X(I+1)**2-X(I+2))
      ELSE IF (MOD(KA,7).EQ.3) THEN
      FA=(X(I+2)-X(I+3))**2
      ELSE IF (MOD(KA,7).EQ.4) THEN
      FA=(X(I+3)-X(I+4))**2
      ELSE IF (MOD(KA,7).EQ.5) THEN
      FA=X(I)+X(I+1)**2+X(I+2)-3.0D 1
      ELSE IF (MOD(KA,7).EQ.6) THEN
      FA=X(I+1)-X(I+2)**2+X(I+3)-1.0D 1
      ELSE
      FA=X(I)*X(I+4)-1.0D 1
      END IF
      RETURN
  140 I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      FA=-Y(L)
      DO 142 K=1,3
      A=DBLE(K*K)/DBLE(L)
      DO 141 J=1,4
      IF (X(I+J).EQ.0) X(I+J)=1.0D-16
      A=A*SIGN(1.0D 0,X(I+J))*ABS(X(I+J))**(DBLE(J)/DBLE(K*L))
  141 CONTINUE
      FA=FA+A
  142 CONTINUE
      RETURN
  150 I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      FA=-Y(L)
      DO 152 K=1,3
      A=0.0D 0
      DO 151 J=1,4
      A=A+X(I+J)*(DBLE(J)/DBLE(K*L))
  151 CONTINUE
      FA=FA+EXP(A)*DBLE(K*K)/DBLE(L)
  152 CONTINUE
      RETURN
  160 I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      FA=-Y(L)
      DO 161 J=1,4
      FA=FA+DBLE((1-2*MOD(J,2))*L*J*J)*SIN(X(I+J))+
     & DBLE(L*L*J)*COS(X(I+J))
  161 CONTINUE
      RETURN
  170 ALFA=0.5D 0
      IF (KA.EQ.1) THEN
      FA=ALFA-(1.0D 0-ALFA)*X(3)-X(1)*(1.0D 0+4.0D 0*X(2))
      ELSE IF (KA.EQ.2) THEN
      FA=-(2.0D 0-ALFA)*X(4)-X(2)*(1.0D 0+4.0D 0*X(1))
      ELSE IF (KA.EQ.N-1) THEN
      FA=ALFA*X(N-3)-X(N-1)*(1.0D 0+4.0D 0*X(N))
      ELSE IF (KA.EQ.N) THEN
      FA=ALFA*X(N-2)-(2.0D 0-ALFA)-X(N)*(1.0D 0+4.0D 0*X(N-1))
      ELSE IF (MOD(KA,2).EQ.1) THEN
      FA=ALFA*X(KA-2)-(1.0D 0-ALFA)*X(KA+2)-
     & X(KA)*(1.0D 0+4.0D 0*X(KA+1))
      ELSE
      FA=ALFA*X(KA-2)-(2.0D 0-ALFA)*X(KA+2)-
     & X(KA)*(1.0D 0+4.0D 0*X(KA-1))
      END IF
      RETURN
  180 IF (KA.LT.2) THEN
      FA=4.0D 0*(X(KA)-X(KA+1)**2)
      ELSE IF (KA.LT.N) THEN
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))+
     & 4.0D 0*(X(KA)-X(KA+1)**2)
      ELSE
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))
      END IF
      RETURN
  190 IF (KA.EQ.1) THEN
      FA=-2.0D 0*X(KA)**2+3.0D 0*X(KA)-2.0D 0*X(KA+1)+
     & 3.0D 0*X(N-4)-X(N-3)-X(N-2)+0.5D 0*X(N-1)-X(N)+1.0D 0
      ELSE IF (KA.LE.N-1) THEN
      FA=-2.0D 0*X(KA)**2+3.0D 0*X(KA)-X(KA-1)-2.0D 0*X(KA+1)+
     & 3.0D 0*X(N-4)-X(N-3)-X(N-2)+0.5D 0*X(N-1)-X(N)+1.0D 0
      ELSE
      FA=-2.0D 0*X(N)**2+3.0D 0*X(N)-X(N-1)+
     & 3.0D 0*X(N-4)-X(N-3)-X(N-2)+0.5D 0*X(N-1)-X(N)+1.0D 0
      END IF
      RETURN
  200 U=1.0D 0/DBLE(N+1)
      V=DBLE(KA)*U
      FA=2.0D 0*X(KA)+0.5D 0*U*U*(X(KA)+V+1.0D 0)**3+1.0D 0
      IF (KA.GT.1) FA=FA-X(KA-1)
      IF (KA.LT.N) FA=FA-X(KA+1)
      RETURN
  210 I=3*((KA+6)/7)-2
      IF (MOD(KA,7).EQ.1) THEN
      FA=1.0D 1*(X(I)**2-X(I+1))
      ELSE IF (MOD(KA,7).EQ.2) THEN
      FA=X(I+1)+X(I+2)-2.0D 0
      ELSE IF (MOD(KA,7).EQ.3) THEN
      FA=X(I+3)-1.0D 0
      ELSE IF (MOD(KA,7).EQ.4) THEN
      FA=X(I+4)-1.0D 0
      ELSE IF (MOD(KA,7).EQ.5) THEN
      FA=X(I)+3.0D 0*X(I+1)
      ELSE IF (MOD(KA,7).EQ.6) THEN
      FA=X(I+2)+X(I+3)-2.0D 0*X(I+4)
      ELSE
      FA=1.0D 1*(X(I+1)**2-X(I+4))
      END IF
      RETURN
  220 I=KA/2
      IF (KA.EQ.1) THEN
      FA=X(KA)-1.0D 0
      ELSE IF (MOD(KA,2).EQ.0) THEN
      FA=1.0D 1*(X(I)**2-X(I+1))
      ELSE
      FA=2.0D 0*EXP(-(X(I)-X(I+1))**2)+
     & EXP(-2.0D 0*(X(I+1)-X(I+2))**2)
      END IF
      RETURN
      END
* SUBROUTINE TAGU15                ALL SYSTEMS                92/12/01
C PORTABILITY : ALL SYSTEMS
C 92/12/01 RA : ORIGINAL VERSION
*
* PURPOSE :
*  GRADIENTS OF PARTIAL FUNCTIONS IN THE SUM OF SQUARES.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  KA  INDEX OF THE GIVEN PARTIAL FUNCTION.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  GA(N)  GRADIENT OF THE KA-TH PARTIAL FUNCTION AT THE POINT X.
*  II  NEXT  NUMBER OF THE SELECTED TEST PROBLEM.
*
      SUBROUTINE TAGU15(N,KA,X,GA,NEXT)
      INTEGER N,KA,NEXT
      DOUBLE PRECISION X(*),GA(*)
      DOUBLE PRECISION A,B,C,D,Q,R,E,ALFA,U,V
      INTEGER I,J,K,L,M,IA,IB,IC
      DOUBLE PRECISION Y(20)
      COMMON /EMPR15/ Y
      GO TO(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,
     & 170,180,190,200,210,220,260,270),NEXT
   10 I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
      GA(I)=2.0D 1*X(I)
      GA(I+1)=-1.0D 1
      ELSE
      GA(I)=1.0D 0
      END IF
      RETURN
   20 I=2*((KA+5)/6)
      A=SQRT(9.0D 1)
      B=SQRT(1.0D 1)
      IF (MOD(KA,6).EQ.1) THEN
      GA(I-1)=2.0D 1 *X(I-1)
      GA(I)=-1.0D 1
      ELSE IF (MOD(KA,6).EQ.2) THEN
      GA(I-1)=1.0D 0
      ELSE IF (MOD(KA,6).EQ.3) THEN
      GA(I+1)=2.0D 0*A*X(I+1)
      GA(I+2)=-A
      ELSE IF (MOD(KA,6).EQ.4) THEN
      GA(I+1)=1.0D 0
      ELSE IF (MOD(KA,6).EQ.5) THEN
      GA(I)=B
      GA(I+2)=B
      ELSE
      GA(I)=1.0D 0/B
      GA(I+2)=-1.0D 0/B
      END IF
      RETURN
   30 I=2*((KA+3)/4)
      A=SQRT(5.0D 0)
      B=SQRT(1.0D 1)
      IF (MOD(KA,4).EQ.1) THEN
      GA(I-1)=1.0D 0
      GA(I)=1.0D 1
      ELSE IF (MOD(KA,4).EQ.2) THEN
      GA(I+1)=A
      GA(I+2)=-A
      ELSE IF (MOD(KA,4).EQ.3) THEN
      C=X(I)-2.0D 0*X(I+1)
      GA(I)=2.0D 0*C
      GA(I+1)=-4.0D 0*C
      ELSE
      C=X(I-1)-X(I+2)
      D=2.0D 0*B*C
      GA(I-1)=D
      GA(I+2)=-D
      END IF
      RETURN
   40 I=2*((KA+4)/5)
      IF (MOD(KA,5).EQ.1) THEN
      A=EXP(X(I-1))
      B=A-X(I)
      C=2.0D 0*B
      GA(I-1)=C*A
      GA(I)=-C
      ELSE IF (MOD(KA,5).EQ.2) THEN
      A=X(I)-X(I+1)
      B=3.0D 1*A**2
      GA(I)=B
      GA(I+1)=-B
      ELSE IF (MOD(KA,5).EQ.3) THEN
      C=X(I+1)-X(I+2)
      Q=SIN(C)/COS(C)
      R=COS(C)
      D=2.0D 0*Q/R**2
      GA(I+1)=D
      GA(I+2)=-D
      ELSE IF (MOD(KA,5).EQ.4) THEN
      GA(I-1)=4.0D 0*X(I-1)**3
      ELSE
      GA(I+2)=1.0D 0
      END IF
      RETURN
   50 I=KA
      GA(I)=3.0D 0-4.0D 0*X(I)
      IF (I.GT.1) GA(I-1)=-1.0D 0
      IF (I.LT.N) GA(I+1)=-1.0D 0
      RETURN
   60 I=KA
      DO 61 J=MAX(1,I-5),MIN(N,I+1)
      GA(J)=1.0D 0+2.0D 0*X(J)
   61 CONTINUE
      GA(I)=GA(I)+2.0D 0+1.5D 1*X(I)**2
      RETURN
   70 I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
      GA(I)=1.0D 0
      GA(I+1)=1.0D 1*X(I+1)-3.0D 0*X(I+1)**2-2.0D 0
      ELSE
      GA(I)=1.0D 0
      GA(I+1)=2.0D 0*X(I+1)+3.0D 0*X(I+1)**2-1.4D 1
      END IF
      RETURN
   80 I=MOD(KA,N/2)+1
      J=I+N/2
      M=5*N
      IF (KA.LE.M/2) THEN
      IA=1
      ELSE
      IA=2
      END IF
      IB=5-KA/(M/4)
      IC=MOD(KA,5)+1
      A=DBLE(IA)
      B=DBLE(IB)
      C=DBLE(IC)
      D=X(I)**IA-X(J)**IB
      IF (D.EQ.0.0D 0) THEN
      GA(I)=0.0D 0
      GA(J)=0.0D 0
      ELSE
      E=C*D**(IC-1)
      IF (X(I).EQ.0.0D 0.AND.IA.LE.1) THEN
      GA(I)=0.0D 0
      ELSE
      GA(I)= E*A*X(I)**(IA-1)
      END IF
      IF (X(J).EQ.0.0D 0.AND.IB.LE.1) THEN
      GA(J)=0.0D 0
      ELSE
      GA(J)=-E*B*X(J)**(IB-1)
      END IF
      END IF
      RETURN
   90 I=2*((KA+5)/6)-1
      IF (MOD(KA,6).EQ.1) THEN
      GA(I)=1.0D 0
      GA(I+1)=3.0D 0*(X(I+2)-1.0D 0)
      GA(I+2)=3.0D 0*X(I+1)
      GA(I+3)=2.0D 0*X(I+3)
      ELSE IF (MOD(KA,6).EQ.2) THEN
      GA(I)=2.0D 0*(X(I)+X(I+1))
      GA(I+1)=2.0D 0*(X(I)+X(I+1))
      GA(I+2)=2.0D 0*(X(I+2)-1.0D 0)
      GA(I+3)=-1.0D 0
      ELSE IF (MOD(KA,6).EQ.3) THEN
      GA(I)=X(I+1)
      GA(I+1)=X(I)
      GA(I+2)=-X(I+3)
      GA(I+3)=-X(I+2)
      ELSE IF (MOD(KA,6).EQ.4) THEN
      GA(I)=2.0D 0*X(I+2)
      GA(I+1)=X(I+3)
      GA(I+2)=2.0D 0*X(I)
      GA(I+3)=X(I+1)
      ELSE IF (MOD(KA,6).EQ.5) THEN
      GA(I)=2.0D 0*(X(I)+X(I+1)+X(I+2)+X(I+3))+2.0D 0*(X(I)-1.0D 0)
      GA(I+1)=2.0D 0*(X(I)+X(I+1)+X(I+2)+X(I+3))
      GA(I+2)=2.0D 0*(X(I)+X(I+1)+X(I+2)+X(I+3))
      GA(I+3)=2.0D 0*(X(I)+X(I+1)+X(I+2)+X(I+3))
      ELSE
      GA(I)=X(I+1)*X(I+2)*X(I+3)
      GA(I+1)=X(I)*X(I+2)*X(I+3)
      GA(I+2)=X(I)*X(I+1)*X(I+3)
      GA(I+3)=X(I)*X(I+1)*X(I+2)+2.0D 0*(X(I+3)-1.0D 0)
      END IF
      RETURN
  100 IF (N.GE.2) THEN
      I=(KA+1)/2
      J=MOD(KA,2)
      IF (J.EQ.0) THEN
      GA(I)=-2.0D 0*EXP(2.0D 0*X(I))
      GA(I+1)=-2.0D 0*EXP(2.0D 0*X(I+1))
      ELSE IF (I.EQ.1) THEN
      GA(I)=-EXP(X(I))
      GA(I+1)=-EXP(X(I+1))
      ELSE IF (I.EQ.N) THEN
      GA(I-1)=-3.0D 0*EXP(3.0D 0*X(I-1))
      GA(I)=-3.0D 0*EXP(3.0D 0*X(I))
      ELSE
      GA(I-1)=-3.0D 0*EXP(3.0D 0*X(I-1))
      GA(I)=-3.0D 0*EXP(3.0D 0*X(I))-EXP(X(I))
      GA(I+1)=-EXP(X(I+1))
      END IF
      END IF
      RETURN
  110 I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
      GA(I)=2.0D 1*(1.0D 0-X(I)**2)/(1.0D 0+X(I)**2)**2
      GA(I+1)=-1.0D 1
      ELSE
      GA(I)=1.0D 0
      END IF
      RETURN
  120 I=3*((KA+5)/6)-2
      IF (MOD(KA,6).EQ.1) THEN
      GA(I)=2.0D 1*X(I)
      GA(I+1)=-1.0D 1
      ELSE IF (MOD(KA,6).EQ.2) THEN
      GA(I+2)=1.0D 0
      ELSE IF (MOD(KA,6).EQ.3) THEN
      GA(I+3)=2.0D 0*(X(I+3)-1.0D 0)
      ELSE IF (MOD(KA,6).EQ.4) THEN
      GA(I+4)=3.0D 0*(X(I+4)-1.0D 0)**2
      ELSE IF (MOD(KA,6).EQ.5) THEN
      GA(I)=2.0D 0*X(I)*X(I+3)
      GA(I+3)=X(I)**2+COS(X(I+3)-X(I+4))
      GA(I+4)=-COS(X(I+3)-X(I+4))
      ELSE
      GA(I+1)=1.0D 0
      GA(I+2)=4.0D 0*X(I+2)*(X(I+2)*X(I+3))**2
      GA(I+3)=2.0D 0*X(I+2)**4*X(I+3)
      END IF
      RETURN
  130 I=3*((KA+6)/7)-2
      IF (MOD(KA,7).EQ.1) THEN
      GA(I)=2.0D 1*X(I)
      GA(I+1)=-1.0D 1
      ELSE IF (MOD(KA,7).EQ.2) THEN
      GA(I+1)=2.0D 1*X(I+1)
      GA(I+2)=-1.0D 1
      ELSE IF (MOD(KA,7).EQ.3) THEN
      GA(I+2)= 2.0D 0*(X(I+2)-X(I+3))
      GA(I+3)=-2.0D 0*(X(I+2)-X(I+3))
      ELSE IF (MOD(KA,7).EQ.4) THEN
      GA(I+3)= 2.0D 0*(X(I+3)-X(I+4))
      GA(I+4)=-2.0D 0*(X(I+3)-X(I+4))
      ELSE IF (MOD(KA,7).EQ.5) THEN
      GA(I)=1.0D 0
      GA(I+1)=2.0D 0*X(I+1)
      GA(I+2)=1.0D 0
      ELSE IF (MOD(KA,7).EQ.6) THEN
      GA(I+1)=1.0D 0
      GA(I+2)=-2.0D 0*X(I+2)
      GA(I+3)=1.0D 0
      ELSE
      GA(I)=X(I+4)
      GA(I+4)=X(I)
      END IF
      RETURN
  140 I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      DO 141 J=1,4
      GA(I+J)=0.0D 0
  141 CONTINUE
      DO 144 K=1,3
      A=DBLE(K*K)/DBLE(L)
      DO 142 J=1,4
      A=A*SIGN(1.0D 0,X(I+J))*ABS(X(I+J))**(DBLE(J)/DBLE(K*L))
  142 CONTINUE
      DO 143 J=1,4
      GA(I+J)=GA(I+J)+(DBLE(J)/DBLE(K*L))*A/X(I+J)
  143 CONTINUE
  144 CONTINUE
      RETURN
  150 I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      DO 151 J=1,4
      GA(I+J)=0.0D 0
  151 CONTINUE
      DO 154 K=1,3
      A=0.0D 0
      DO 152 J=1,4
      A=A+X(I+J)*(DBLE(J)/DBLE(K*L))
  152 CONTINUE
      A=EXP(A)*DBLE(K*K)/DBLE(L)
      DO 153 J=1,4
      GA(I+J)=GA(I+J)+A*(DBLE(J)/DBLE(K*L))
  153 CONTINUE
  154 CONTINUE
      RETURN
  160 I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      DO 161 J=1,4
      GA(I+J)=DBLE((1-2*MOD(J,2))*L*J*J)*COS(X(I+J))-
     & DBLE(L*L*J)*SIN(X(I+J))
  161 CONTINUE
      RETURN
  170 ALFA=0.5D 0
      IF (KA.EQ.1) THEN
      GA(1)=-1.0D 0-4.0D 0*X(2)
      GA(2)=-4.0D 0*X(1)
      GA(3)=ALFA-1.0D 0
      ELSE IF (KA.EQ.2) THEN
      GA(1)=-4.0D 0*X(2)
      GA(2)=-1.0D 0-4.0D 0*X(1)
      GA(4)=-2.0D 0+ALFA
      ELSE IF (KA.EQ.N-1) THEN
      GA(N-3)=ALFA
      GA(N-1)=-1.0D 0-4.0D 0*X(N)
      GA(N)=-4.0D 0*X(N-1)
      ELSE IF (KA.EQ.N) THEN
      GA(N-2)=ALFA
      GA(N-1)=-4.0D 0*X(N)
      GA(N)=-1.0D 0-4.0D 0*X(N-1)
      ELSE IF (MOD(KA,2).EQ.1) THEN
      GA(KA-2)=ALFA
      GA(KA)=-1.0D 0-4.0D 0*X(KA+1)
      GA(KA+1)=-4.0D 0*X(KA)
      GA(KA+2)=-1.0D 0+ALFA
      ELSE
      GA(KA-2)=ALFA
      GA(KA-1)=-4.0D 0*X(KA)
      GA(KA)=-1.0D 0-4.0D 0*X(KA-1)
      GA(KA+2)=-2.0D 0+ALFA
      END IF
      RETURN
  180 IF (KA.LT.2) THEN
      GA(KA)=4.0D 0
      GA(KA+1)=-8.0D 0*X(KA+1)
      ELSE IF (KA.LT.N) THEN
      GA(KA-1)=-8.0D 0*X(KA)
      GA(KA)=24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+6.0D 0
      GA(KA+1)=-8.0D 0*X(KA+1)
      ELSE
      GA(KA-1)=-8.0D 0*X(KA)
      GA(KA)=24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+2.0D 0
      END IF
      RETURN
  190 IF (KA.EQ.1) THEN
      GA(N-4)=3.0D 0
      GA(N-3)=-1.0D 0
      GA(N-2)=-1.0D 0
      GA(N-1)=0.50D 0
      GA(N)=-1.0D 0
      GA(KA)=-4.0D 0*X(KA)+3.0D 0
      GA(KA+1)=-2.0D 0
      ELSE IF (KA.LE.N-1) THEN
      GA(KA-1)=0.0D 0
      GA(KA)=0.0D 0
      GA(KA+1)=0.0D 0
      GA(N-4)=3.0D 0
      GA(N-3)=-1.0D 0
      GA(N-2)=-1.0D 0
      GA(N-1)=0.50D 0
      GA(N)=-1.0D 0
      GA(KA-1)=GA(KA-1)-1.0D 0
      GA(KA)=GA(KA)-4.0D 0*X(KA)+3.0D 0
      GA(KA+1)=GA(KA+1)-2.0D 0
      ELSE
      GA(N-4)=3.0D 0
      GA(N-3)=-1.0D 0
      GA(N-2)=-1.0D 0
      GA(N-1)=0.50D 0
      GA(N)=-4.0D 0*X(N)+2.0D 0
      END IF
      RETURN
  200 U=1.0D 0/DBLE(N+1)
      V=DBLE(KA)*U
      GA(KA)=2.0D 0+1.5D 0*U**2*(X(KA)+V+1.0D 0)**2
      IF (KA.GT.1) GA(KA-1)=-1.0D 0
      IF (KA.LT.N) GA(KA+1)=-1.0D 0
      RETURN
  210 I=3*((KA+6)/7)-2
      IF (MOD(KA,7).EQ.1) THEN
      GA(I)=2.0D 1*X(I)
      GA(I+1)=-1.0D 1
      ELSE IF (MOD(KA,7).EQ.2) THEN
      GA(I+1)=1.0D 0
      GA(I+2)=1.0D 0
      ELSE IF (MOD(KA,7).EQ.3) THEN
      GA(I+3)=1.0D 0
      ELSE IF (MOD(KA,7).EQ.4) THEN
      GA(I+4)=1.0D 0
      ELSE IF (MOD(KA,7).EQ.5) THEN
      GA(I)=1.0D 0
      GA(I+1)=3.0D 0
      ELSE IF (MOD(KA,7).EQ.6) THEN
      GA(I+2)=1.0D 0
      GA(I+3)=1.0D 0
      GA(I+4)=-2.0D 0
      ELSE
      GA(I+1)=2.0D 1*X(I+1)
      GA(I+4)=-1.0D 1
      END IF
      RETURN
  220 I=KA/2
      IF (KA.EQ.1) THEN
      GA(KA)=1.0D 0
      ELSE IF (MOD(KA,2).EQ.0) THEN
      GA(I)=2.0D 1*X(I)
      GA(I+1)=-1.0D 1
      ELSE
      A=2.0D 0*EXP(-(X(I)-X(I+1))**2)
      B=EXP(-2.0D 0*(X(I+1)-X(I+2))**2)
      GA(I)= -2.0D 0*(X(I)-X(I+1))*A
      GA(I+1)=2.0D 0*(X(I)-X(I+1))*A-4.0D 0*(X(I+1)-X(I+2))*B
      GA(I+2)=4.0D 0*(X(I+1)-X(I+2))*B
      END IF
      RETURN
  260 I=MOD(KA,N/2)+1
      J=I+N/2
      M=5*N
      IA=KA/(M/4)+1
      A=DBLE(IA)
      B=DBLE(KA/(M/5)+1)
      GA(I)=A*X(I)**(IA-1)*EXP(B*X(J))
      GA(J)=B*X(I)**IA*EXP(B*X(J))+1.0D 0
      RETURN
  270 IA=MIN(MAX(MOD(KA,13)-2,1),7)
      IB=(KA+12)/13
      I=IA+IB-1
      IF (IA.EQ.7) THEN
      J=IB
      ELSE
      J=IA+IB
      END IF
      C=3.0D 0*DBLE(IA)/1.0D 1
      A=COS(C)
      B=EXP(SIN(C*X(J)))
      GA(I)=-COS(X(I))*(1.0D 0+A)+5.0D 0*B
      GA(J)=(1.0D 0+A)+5.0D 0*(X(I)-2.0D 0)*C*COS(C*X(J))*B
      DO 271 L=0,6
      IF (IB+L.NE.I.AND.IB+L.NE.J) GA(IB+L)=0.5D 0*COS(X(IB+L))
  271 CONTINUE
      RETURN
      END
* SUBROUTINE TIUD16             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 99/12/01 GA : ORIGINAL VERSION
*
* PURPOSE :
*  INITIAL VALUES OF THE VARIABLES FOR NONLINEAR EQUATIONS.
*  UNCONSTRAINED AND DENSE VERSION.
*
* PARAMETERS :
*  IO  N  NUMBER OF VARIABLES.
*  IO  NA  NUMBER OF EQUATIONS.
*  RO  X(N)  VECTOR OF VARIABLES.
*  RO  FMIN  LOWER BOUND FOR VALUE OF THE OBJECTIVE FUNCTION.
*  RO  XMAX  MAXIMUM STEPSIZE.
*  IO  NEXT  NUMBER OF THE TEST PROBLEM.
*  IO  IERR  ERROR INDICATOR.
*
      SUBROUTINE TIUD16(N,NA,X,FMIN,XMAX,NEXT,IERR)
      INTEGER N,NA,NEXT,IERR
      DOUBLE PRECISION X(N),FMIN,XMAX
      DOUBLE PRECISION Y(1000),F,S,T
      INTEGER I,N1,ALF,BET
      NA=N
      FMIN=0.0D 0
      XMAX=1.0D 3
      IERR=0
      GO TO (201,202,203,204,205,206,207,208,209,210,211,212,213,214,
     & 217,218,219,220,221,222,223,224,225,226,227,228,
     & 229,230,231,231),NEXT
201   N1=N-1
      DO 3030 I=1,N1
        X(I)=-1.2D 0
3030  CONTINUE
      X(N)=-1.0D 0
      RETURN
202   DO 3040 I=1,N
        X(I)=2.0D 0
3040  CONTINUE
      RETURN
203   DO 3050 I=1,N
        X(I)=0.5D 0
3050  CONTINUE
      RETURN
204   N1=N-1
      DO 3060 I=1,N1,2
        X(I)=-1.2D 0
        X(I+1)=1.0D 0
3060  CONTINUE
      RETURN
205   DO 3070 I=1,N
        X(I)=1.5D 0
3070  CONTINUE
      RETURN
206   DO 3080 I=1,N
        X(I)=0.0D 0
3080  CONTINUE
      RETURN
207   DO 3090 I=1,N
        X(I)=-1.0D 0
3090  CONTINUE
      RETURN
208   DO 4000 I=1,N
        X(I)=-1.0D 0
4000  CONTINUE
      RETURN
209   DO 4010 I=1,N
        X(I)=0.5D 0
4010  CONTINUE
      RETURN
210   DO 4020 I=1,N
        X(I)=0.0D 0
4020  CONTINUE
      RETURN
211   DO 4030 I=1,N
        X(I)=0.0D 0
4030  CONTINUE
      RETURN
212   DO 4040 I=1,N
        X(I)=0.5D 0
4040  CONTINUE
      RETURN
213   DO 4050 I=1,N
        X(I)=1.0D 0
4050  CONTINUE
      RETURN
214   DO 4060 I=1,N
        X(I)=-1.0D 0
4060  CONTINUE
      RETURN
215   DO 4070 I=1,N
        X(I)=0.0D 0
4070  CONTINUE
      RETURN
216   DO 4080 I=1,N
         X(I)=DBLE(I)/DBLE(N+1)
4080  CONTINUE
      RETURN
217   DO 4090 I=1,N
        X(I)=10.0D 0
4090  CONTINUE
      RETURN
218   DO 5000 I=1,N
        T=DBLE(I)/DBLE(N+1)
        X(I)=T*(T-1.0D 0)
5000  CONTINUE
      RETURN
219   DO 5010 I=1,N
        T=DBLE(I)/DBLE(N+1)
        X(I)=T*(T-1.D 0)
5010  CONTINUE
      RETURN
220   DO 5020 I=1,N
        X(I)=1.0D 0/DBLE(N)
5020  CONTINUE
      RETURN
221   DO 5030 I=1,N
        X(I)=1.0D 0-DBLE(I)/DBLE(N)
5030  CONTINUE
      RETURN
222   DO 5040 I=1,N
        X(I)=-1.0D 0
5040  CONTINUE
      RETURN
223   ALF=5
      BET=14
      DO 5051 I=1,N
        X(I)=0.0D 0
5051  CONTINUE
      DO 5052 I=1,N
        CALL EAFU16(N,I,X,F,NEXT)
        Y(I)=-F
5052  CONTINUE
      S=DBLE(BET*N)/DBLE(BET**2*N**2-(ALF+1)**2*(N-1)**2)
      DO 5050 I=1,N
        X(I)=S*Y(I)
5050  CONTINUE
      RETURN
224   DO 5060 I=1,N
        X(I)=1.0D 0
5060  CONTINUE
      RETURN
225   DO 5070 I=1,N
        X(I)=0.0D 0
5070  CONTINUE
      RETURN
226   DO 5080 I=1,N
        X(I)=1.0D 0
5080   CONTINUE
      RETURN
227   DO 5090 I=1,N
        X(I)=1.0D 0
5090   CONTINUE
      RETURN
228   DO 6000 I=1,N
        X(I)=1.0D 0
6000   CONTINUE
      RETURN
229   DO 6010 I=1,N
        X(I)=1.0D 0
6010  CONTINUE
      RETURN
230   T=DBLE(2)/DBLE(N+2)
      N1=N/2
      DO 6020 I=1,N1
        S=DBLE(I)*T
        X(I)=S*(1.0D 0-S)
        X(N1+I)=X(I)
6020  CONTINUE
      RETURN
231   N1=INT(SQRT(DBLE(N)))
      N=N1*N1
      DO 6030 I=1,N
        X(I)=1.0D 0
6030  CONTINUE
      RETURN
      END
* SUBROUTINE EAFU16             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 99/12/01 GA : ORIGINAL VERSION
*
* PURPOSE :
*  VALUES OF TEST FUNCTIONS FOR NONLINEAR EQUATIONS.
*  UNIVERSAL VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  I  INDEX OF THE APPROXIMATED FUNCTION.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  F  VALUE OF THE APPROXIMATED FUNCTION AT THE
*          SELECTED POINT.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE EAFU16(N,I,X,F,NEXT)
      INTEGER N,I,NEXT
      DOUBLE PRECISION X(N),F
      DOUBLE PRECISION S,T,S1,S2,S3,T1,C,H,AL1,AL2,BE1,BE2
      DOUBLE PRECISION AL,BE,A,B,GA,CA,CB,FF,U,LA,H2
      INTEGER J,N1,I1,I2,J1,J2,K,ALF,BET,GAM,L,ND
      GO TO (201,202,203,204,205,206,207,208,209,210,211,212,213,214,
     & 217,218,219,220,221,222,223,224,225,226,227,228,
     & 229,230,231,232),NEXT
201   IF(I.EQ.1) THEN
        F=1.0D 0-X(1)
      ELSE
        F=10.0D 0*DBLE(I-1)*(X(I)-X(I-1))**2
      ENDIF
      RETURN
202   IF(I.EQ.N) THEN
        F=X(I)-0.1D 0*X(1)**2
      ELSE
        F=X(I)-0.1D 0*X(I+1)**2
      ENDIF
      RETURN
203   T=-1.0D 0
      IF(I.LT.N) THEN
        T=T+X(I)
        DO 3010 J=1,N
           T=T+X(J)
3010    CONTINUE
        F=T-DBLE(N)
      ELSE
        S=1.0D 0
        DO 3020 J=1,N
           S=S*X(J)
3020    CONTINUE
        F=T+S
      ENDIF
      RETURN
204   IF(I/2*2.LT.I) THEN
        F=1.D 0-X(I)
      ELSE
        F=10.0D 0*(X(I)-X(I-1)**2)
      ENDIF
      RETURN
205   S=0.0D 0
      DO 3030 J=1,N
        S=S+X(J)**3
3030  CONTINUE
      F=X(I)-1.0D 0/DBLE(2*N)*(S+DBLE(I))
      RETURN
206   S=1.0D 0/DBLE(N+1)
      IF(N.EQ.1) THEN
        F=-2.0D 0*X(I)-S**2*EXP(X(I))
      ELSE IF(I.EQ.1) THEN
        F=-2.0D 0*X(I)+X(I+1)-S**2*EXP(X(I))
      ELSE IF(I.EQ.N) THEN
        F=X(I-1)-2.0D 0*X(I)-S**2*EXP(X(I))
      ELSE
        F=X(I-1)-2.0D 0*X(I)+X(I+1)-S**2*EXP(X(I))
      ENDIF
      RETURN
207   S=0.1D 0
      IF(N.EQ.1) THEN
        F=(3.0D 0-S*X(I))*X(I)+1.0D 0
      ELSE IF(I.EQ.1) THEN
        F=(3.0D 0-S*X(I))*X(I)+1.0D 0-2.0D 0*X(I+1)
      ELSE IF(I.EQ.N) THEN
        F=(3.0D 0-S*X(I))*X(I)+1.0D 0-X(I-1)
      ELSE
        F=(3.0D 0-S*X(I))*X(I)+1.0D 0-X(I-1)-2.0D 0*X(I+1)
      ENDIF
      RETURN
208   S1=1.0D 0
      S2=1.0D 0
      S3=1.0D 0
      J1=3
      J2=3
      IF(I-J1.GT.1) THEN
        I1=I-J1
      ELSE
        I1=1
      ENDIF
      IF(I+J2.LT.N) THEN
        I2=I+J2
      ELSE
        I2=N
      ENDIF
      S=0.0D 0
      DO 3040 J=I1,I2
        IF(J.NE.I) S=S+X(J)+X(J)**2
3040  CONTINUE
      F=(S1+S2*X(I)**2)*X(I)+1.D 0-S3*S
      RETURN
209   IF(I.EQ.1) THEN
        F=X(1)**2-1.0D 0
      ELSE
        F=X(I-1)**2+LOG(X(I))-1.0D 0
      ENDIF
      RETURN
210   S=0.0D 0
      DO 3050 J=1,N
        S=S+X(J)
3050  CONTINUE
      F=EXP(COS(DBLE(I)*S))
      RETURN
211   S=0.0D 0
      DO 3060 J=1,N
        S=S+X(J)**3
3060  CONTINUE
      F=1.0D 0/DBLE(2*N)*(S+DBLE(I))
      RETURN
212   IF(I.EQ.1) THEN
        F=X(1)
      ELSE
        F=COS(X(I-1))+X(I)-1.0D 0
      ENDIF
      RETURN
213   S=(1.0D 0/DBLE(N+1))**2
      IF(N.EQ.1) THEN
        F=2.0D 0*X(1)-1.D 0+S*(X(1)+SIN(X(1)))
      ELSE IF(I.EQ.1) THEN
        F=2.0D 0*X(1)-X(I+1)+S*(X(1)+SIN(X(1)))
      ELSE IF(I.EQ.N) THEN
        F=-X(I-1)+2.0D 0*X(I)-1.D 0+S*(X(I)+SIN(X(I)))
      ELSE
        F=-X(I-1)+2.0D 0*X(I)-X(I+1)+S*(X(I)+SIN(X(I)))
      ENDIF
      RETURN
214   IF(I-5.GT.1) THEN
        I1=I-5
      ELSE
        I1=1
      ENDIF
      IF(I+1.LT.N) THEN
        I2=I+1
      ELSE
        I2=N
      ENDIF
      S=0.D 0
      DO 3070 J=I1,I2
        IF(J.NE.I) S=S+X(J)*(1.0D 0+X(J))
3070  CONTINUE
      F=X(I)*(2.0D 0+5.0D 0*X(I)**2)+1.0D 0-S
      RETURN
215   S3=0.0D 0
      DO 4000 K=1,29
        S1=0.0D 0
        S=DBLE(K)/DBLE(29)
        DO 3080 J=1,N
           S1=S1+S**(J-1)*X(J)
3080    CONTINUE
        T1=0.0D 0
        DO 3090 J=2,N
           T1=T1+(J-1)*S**(J-2)*X(J)
3090    CONTINUE
        S3=S3+S**(I-2)*((I-1.0D 0-2.0D 0*S*S1)*(T1-S1**2-1.0D 0))
4000  CONTINUE
      IF (I.EQ.1) THEN
         F=S3+X(1)*(1.0D 0-2*(X(2)-X(1)*X(1)-1.0D 0))
      ELSE IF(I.EQ.2) THEN
        F=S3+X(2)-X(1)*X(1)-1.0D 0
      ELSE
        F=S3
      ENDIF
      RETURN
216   S=0.D 0
      DO 4020 J=1,N
        S1=2.D 0*X(J)-1.D 0
        T1=S1
        IF(I.GT.1) S2=2.D 0*S1**2-1.D 0
        IF(I.EQ.1) THEN
           S=S+S1
        ELSE IF(I.EQ.2) THEN
           S=S+S2
        ELSE
           DO 4010 K=3,I
              S3=2*T1*S2-S1
              S1=S2
              S2=S3
4010       CONTINUE
           S=S+S3
        ENDIF
4020  CONTINUE
      IF( I/2*2.EQ.I) THEN
         F=S/DBLE(N)+1.0D 0/DBLE(I**2-1)
      ELSE
        F=S/DBLE(N)
      ENDIF
      RETURN
217   IF(N.EQ.1) THEN
        F=3.D 0*X(1)*(20.D 0-2.D 0*X(1))+100.D 0
      ELSE IF(I.EQ.1) THEN
        F=3.D 0*X(I)*(X(I+1)-2.D 0*X(I))+(X(I+1))**2/4.0D 0
      ELSE IF(I.EQ.N) THEN
        F=3.D 0*X(I)*(20.0D 0-2.D 0*X(I)+X(I-1))+
     &   (20.0D 0-X(I-1))**2/4.0D 0
      ELSE
        F=3.D 0*X(I)*(X(I+1)-2.D 0*X(I)+X(I-1))+
     &   (X(I+1)-X(I-1))**2/4.0D 0
      ENDIF
      RETURN
218   S=1.0D 0/DBLE(N+1)
      IF(N.EQ.1) THEN
        F=2.D 0*X(1)+S**2/2.0D 0*(X(1)+S+1.0D 0)**3
      ELSE IF(I.EQ.1) THEN
       F=2.D 0*X(1)-X(I+1)+S**2/2.0D 0*(X(1)+DBLE(I)*S+1.0D 0)**3
      ELSE IF(I.EQ.N) THEN
        F=2.D 0*X(I)-X(I-1)+S**2/2.0D 0*(X(I)+DBLE(I)*S+1.0D 0)**3
      ELSE
        F=2.D 0*X(I)-X(I-1)-X(I+1)+S**2/2.0D 0*(X(I)+
     +  DBLE(I)*S+1.0D 0)**3
      ENDIF
      RETURN
219   S1=1.0D 0/(DBLE(N)+1.0D 0)
      S2=0.0D 0
      S3=0.0D 0
      DO 4030 J=1,I
        S2=S2+DBLE(J)*S1*(X(J)+DBLE(J)*S1+1.0D 0)**3
4030  CONTINUE
      DO 4040 J=I+1,N
        S3=S3+(1.0D 0-DBLE(J)*S1)*(X(J)+DBLE(J)*S1+1.0D 0)**3
4040  CONTINUE
      F=X(I)+S1/2.0D 0*((1.0D 0-DBLE(I)*S1)*S2+DBLE(I)*S1*S3)
      RETURN
220   S=0.0D 0
      DO 4050 J=1,N
        S=S+COS(X(J))
4050  CONTINUE
      F=DBLE(N)-S+DBLE(I)*(1.0D 0-COS(X(I)))-SIN(X(I))
      RETURN
221   S=0.0D 0
      DO 4060 J=1,N
        S=S+DBLE(J)*(X(J)-1.0D 0)
4060  CONTINUE
      F=X(I)-1.0D 0+DBLE(I)*S*(1.0D 0+2.0D 0*S*S)
      RETURN
222   IF(I.EQ.1.AND.N.GT.1) THEN
        F=(3.D 0-2.D 0*X(I))*X(I)-2.D 0*X(I+1)+1.D 0
      ELSE IF(I.EQ.N.AND.N.GT.1) THEN
        F=(3.D 0-2.D 0*X(I))*X(I)-X(I-1)+1.D 0
      ELSE IF(N.GT.1) THEN
        F=(3.D 0-2.D 0*X(I))*X(I)-X(I-1)-2.D 0*X(I+1)+1.D 0
      ELSE
        F=(3.D 0-2.D 0*X(I))*X(I)+1.D 0
      ENDIF
      RETURN
223   ALF=5
      BET=14
      GAM=3
      F=DBLE(BET*N)*X(I)+(DBLE(I)-DBLE(N)/2.0D 0)**GAM
      DO 4070 J=1,N
        IF(J.NE.I) THEN
        T=SQRT(X(J)**2+DBLE(I)/DBLE(J))
        S1=LOG(T)
        F=F+T*(SIN(S1)**ALF+COS(S1)**ALF)
      ENDIF
4070  CONTINUE
      RETURN
224   C=0.5D 0
      H=1.0D 0/DBLE(N)
      F=(1.0D 0-C*H/4.0D 0)
      DO 4080 J=1,N
        S=C*H*DBLE(I)/DBLE(2*(I+J))
        IF(J.EQ.N) S=S/2.0D 0
        F=F-S*X(J)
4080  CONTINUE
      F=-1.0D 0+X(I)*F
      RETURN
225   IF(I.EQ.1) THEN
        F=3.0D 0*X(I)**3+2.0D 0*X(I+1)-5.0D 0+
     +  SIN(X(I)-X(I+1))*SIN(X(I)+X(I+1))
      ELSE IF(I.NE.N) THEN
        F=3.0D 0*X(I)**3+2.0D 0*X(I+1)-5.0D 0+
     +  SIN(X(I)-X(I+1))*SIN(X(I)+X(I+1))+
     +  4.0D 0*X(I)-3.0D 0-X(I-1)*EXP(X(I-1)-X(I))
      ELSE
        F= 4.0D 0*X(I)-3.0D 0-X(I-1)*EXP(X(I-1)-X(I))
      ENDIF
      RETURN
226   T=10.0D 0
      H=T/DBLE(N+1)**2
      S=T*X(I)
      IF(I.EQ.1) THEN
        F=2.0D 0*X(I)-X(I+1)+H*(EXP(S)-EXP(-S))/2.D 0
      endif
      IF(1.LT.I.AND.I.LT.N) THEN
        F=-X(I-1)+2.0D 0*X(I)-X(I+1)+H*(EXP(S)-EXP(-S))/2.D 0
      endif
      IF(I.EQ.N) THEN
        F=-X(I-1)+2.0D 0*X(I)-1.0D 0+H*(EXP(S)-EXP(-S))/2.D 0
      ENDIF
      RETURN
227   S=0.5D 0
      H=1.0D 0/DBLE(N+1)
      T=H**2/S
      H=2.0D 0*H
      IF(I.EQ.1) THEN
        F=2.0D 0*X(I)-X(I+1)-T*(X(I)**2+X(I+1)/H)
      ENDIF
      IF(1.LT.I.AND.I.LT.N) THEN
        F=-X(I-1)+2.0D 0*X(I)-X(I+1)-T*(X(I)**2+(X(I+1)-X(I-1))/H)
      ENDIF
      IF(I.EQ.N) THEN
        F=-X(I-1)+2.0D 0*X(I)-0.5D 0-T*(X(I)**2+(0.5D 0-X(I-1))/H)
      ENDIF
      RETURN
228   S=0.5D 0
      H=1.0D 0/DBLE(N+1)
      T=H**2/S
      T1=2.0D 0*H
      AL=0.0D 0
      BE=0.5D 0
      S1=0.0D 0
      DO 4090 J=1,I
        IF(J.EQ.1) THEN
           S1=S1+DBLE(J)*(X(J)**2+(X(J+1)-AL)/T1)
        ENDIF
        IF(1.LT.J.AND.J.LT.N) THEN
           S1=S1+DBLE(J)*(X(J)**2+(X(J+1)-X(J-1))/T1)
        ENDIF
        IF(J.EQ.N) THEN
           S1=S1+DBLE(J)*(X(J)**2+(BE-X(J-1))/T1)
        ENDIF
4090  CONTINUE
      S1=(1.0D 0-DBLE(I)*H)*S1
      IF(I.EQ.N) GO TO 5010
      S2=0.0D 0
      DO 5000 J=I+1,N
        IF(J.LT.N) THEN
           S2=S2+(1.0D 0-DBLE(J)*H)*(X(J)**2+(X(J+1)-X(J-1))/T1)
        ELSE
           S2=S2+(1.0D 0-DBLE(J)*H)*(X(J)**2+(BE-X(J-1))/T1)
        ENDIF
5000  CONTINUE
      S1=S1+DBLE(I)*S2
5010  F=X(I)-0.5D 0*DBLE(I)*H-T*S1
      RETURN
229   A=-9.0D-3
      B=1.0D-3
      AL=0.0D 0
      BE=25.0D 0
      GA=20.0D 0
      CA=0.3D 0
      CB=0.3D 0
      H=(B-A)/DBLE(N+1)
      T=A+DBLE(I)*H
      H=H**2
      S=DBLE(I)/DBLE(N+1)
      U=AL*(1.0D 0-S)+BE*S+X(I)
      FF=CB*EXP(GA*(U-BE))-CA*EXP(GA*(AL-U))
      IF(T.LE.0) THEN
        FF=FF+CA
      ELSE
        FF=FF-CB
      ENDIF
      IF(N.EQ.1) THEN
        F=-AL+2.0D 0*X(I)-BE+H*FF
      ELSEIF(I.EQ.1) THEN
        F=-AL+2.0D 0*X(I)-X(I+1)+H*FF
      ELSEIF(I.LT.N) THEN
        F=-X(I-1)+2.0D 0*X(I)-X(I+1)+H*FF
      ELSE
        F=-X(I-1)+2.0D 0*X(I)-BE+H*FF
      ENDIF
      RETURN
230   AL1=0.0D 0
      AL2=0.0D 0
      BE1=0.0D 0
      BE2=0.0D 0
      N1=N/2
      H=1.0D 0/DBLE(N1+1)
      T=DBLE(I)*H
      IF(I.EQ.1) THEN
        S1=2.0D 0*X(I)-X(I+1)
        B=AL1
      ELSE IF(I.EQ.N1+1) THEN
        S1=2.0D 0*X(I)-X(I+1)
        B=AL2
      ELSE IF(I.EQ.N1) THEN
        S1=-X(I-1)+2.0D 0*X(I)
        B=BE1
      ELSE IF(I.EQ.N) THEN
        S1=-X(I-1)+2.0D 0*X(I)
        B=BE2
      ELSE
        S1=-X(I-1)+2.0D 0*X(I)-X(I+1)
        B=0.0D 0
      ENDIF
      IF(I.LE.N1) THEN
        S2=X(I)**2+X(I)+0.1D 0*X(N1+I)**2-1.2D 0
      ELSE
        S2=0.2D 0*X(I-N1)**2+X(I)**2+2.0D 0*X(I)-0.6D 0
      ENDIF
      F=S1+H**2*S2-B
      RETURN
231   ND=INT(SQRT(DBLE(N)))
      L=MOD(I,ND)
      IF(L.EQ.0) THEN
         K=I/ND
         L=ND
      ELSE
         K=INT(I/ND)+1
      ENDIF
      LA=1.0D 0
      H=1.0D 0/DBLE(ND+1)
      H2=LA*H*H
      IF(L.EQ.1.AND.K.EQ.1) THEN
         F=4.0D 0*X(1)-X(2)-X(ND+1)+H2*EXP(X(1))
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.1) THEN
         F=4.0D 0*X(L)-X(L-1)-X(L+1)-X(L+ND)+H2*EXP(X(L))
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.1) THEN
         F=4.0D 0*X(ND)-X(ND-1)-X(ND+ND)+H2*EXP(X(ND))
      ENDIF
      IF(L.EQ.1.AND.1.LT.K.AND.K.LT.ND) THEN
         F=4.0D 0*X(I)-X(I+1)-X(I-ND)-X(I+ND)+H2*EXP(X(I))
      ENDIF
      IF(L.EQ.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         F=4.0D 0*X(I)-X(I-ND)-X(I-1)-X(I+ND)+H2*EXP(X(I))
      ENDIF
      IF(L.EQ.1.AND.K.EQ.ND) THEN
         F=4.0D 0*X(I)-X(I+1)-X(I-ND)+H2*EXP(X(I))
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.ND) THEN
         F=4.0D 0*X(I)-X(I-1)-X(I+1)-X(I-ND)+H2*EXP(X(I))
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.ND) THEN
        F=4.0D 0*X(I)-X(I-1)-X(I-ND)+H2*EXP(X(I))
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         F=4.0D 0*X(I)-X(I-1)-X(I+1)-X(I-ND)-X(I+ND)+H2*EXP(X(I))
      ENDIF
      RETURN
232   ND=INT(SQRT(DBLE(N)))
      L=MOD(I,ND)
      IF(L.EQ.0) THEN
         K=I/ND
         L=ND
      ELSE
         K=INT(I/ND)+1
      ENDIF
      H=1.0D 0/DBLE(ND+1)
      H2=H*H
      IF(L.EQ.1.AND.K.EQ.1) THEN
         F=4.0D 0*X(1)-X(2)-X(ND+1)+H2*X(1)**2-24.0D 0/(H+1.0D 0)**2
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.1) THEN
         F=4.0D 0*X(L)-X(L-1)-X(L+1)-X(L+ND)+H2*X(L)**2
     *-12.0D 0/(DBLE(L)*H+1.0D 0)**2
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.1) THEN
         F=4.0D 0*X(ND)-X(ND-1)-X(ND+ND)+H2*X(ND)**2
     *-12.0D 0/(DBLE(ND)*H+1.0D 0)**2-12.0D 0/(H+2.0D 0)**2
      ENDIF
      IF(L.EQ.1.AND.1.LT.K.AND.K.LT.ND) THEN
         F=4.0D 0*X(I)-X(I+1)-X(I-ND)-X(I+ND)+H2*X(I)**2
     *-12.0D 0/(DBLE(K)*H+1.0D 0)**2
      ENDIF
      IF(L.EQ.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         F=4.0D 0*X(I)-X(I-ND)-X(I-1)-X(I+ND)+H2*X(I)**2
     *-12.0D 0/(DBLE(K)*H+2.0D 0)**2
      ENDIF
      IF(L.EQ.1.AND.K.EQ.ND) THEN
         F=4.0D 0*X(I)-X(I+1)-X(I-ND)+H2*X(I)**2
     *-12.0D 0/(DBLE(ND)*H+1.0D 0)**2-12.0D 0/(H+2.0D 0)**2
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.ND) THEN
         F=4.0D 0*X(I)-X(I-1)-X(I+1)-X(I-ND)+H2*X(I)**2
     *-12.0D 0/(DBLE(L)*H+2.0D 0)**2
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.ND) THEN
        F=4.0D 0*X(I)-X(I-1)-X(I-ND)+H2*X(I)**2
     *-24.0D 0/(DBLE(ND)*H+2.0D 0)**2
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         F=4.0D 0*X(I)-X(I-1)-X(I+1)-X(I-ND)-X(I+ND)+H2*X(I)**2
      ENDIF
      RETURN
      END
* SUBROUTINE EAGU16             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 99/12/01 GA : ORIGINAL VERSION
*
* PURPOSE :
*  GRADIENTS OF TEST FUNCTIONS FOR NONLINEAR EQUATIONS.
*  UNIVERSAL VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  I  INDEX OF THE APPROXIMATED FUNCTION.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  G(N)  GRADIENT OF THE APPROXIMATED FUNCTION AT THE
*          SELECTED POINT.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE EAGU16(N,I,X,G,NEXT)
      INTEGER N,I,NEXT
      DOUBLE PRECISION X(N),G(N)
      DOUBLE PRECISION BROW(1000),T,S,S1,S2,S3,T1,T2,C,H,H2,LA,A,B,AL,BE
      DOUBLE PRECISION GA,CA,CB,U,FF,TI,TIM1,TIM2,SI,SIM1,SIM2,V
      INTEGER J,K,I1,I2,J1,J2,L,ALF,BET,N1,ND
      GO TO (201,202,203,204,205,206,207,208,209,210,211,212,213,214,
     & 217,218,219,220,221,222,223,224,225,226,227,228,
     & 229,230,231,232),NEXT
201   IF(I.EQ.1) THEN
      G(1)=-1.0D 0
      DO 3010 J=2,N
        G(J)=0.0D 0
3010  CONTINUE
      ELSE
      S=20.0D 0*DBLE(I-1)*(X(I)-X(I-1))
      DO 3020 J=1,N
        IF(J.EQ.I) THEN
           G(J)=S
        ELSE IF(J.EQ.I-1) THEN
           G(J)=-S
        ELSE
          G(J)=0.0D 0
        ENDIF
3020  CONTINUE
      ENDIF
      RETURN
202   IF(I.EQ.N) THEN
        G(1)=-0.20D 0*X(1)
        N1=N-1
        DO 3030 J=2,N1
           G(J)=0.0D 0
3030    CONTINUE
        G(N)=1.0D 0
      ELSE
        DO 3050 J=1,N
           IF(J.EQ.I) THEN
              G(J)=1.0D 0
           ELSE IF(J.EQ.I+1) THEN
              G(J)=-0.2D 0*X(I+1)
           ELSE
              G(J)=0.0D 0
           ENDIF
3050    CONTINUE
      ENDIF
      RETURN
203   IF(I.EQ.N) THEN
        DO 3070 J=1,N
           T=1.0D 0
           DO 3060 K=1,N
              IF(K.NE.J) T=T*X(K)
3060       CONTINUE
        G(J)=T
3070    CONTINUE
      ELSE
        DO 3080 J=1,N
           IF(J.EQ.I) THEN
              G(J)=2.0D 0
           ELSE
              G(J)=1.0D 0
           ENDIF
3080    CONTINUE
      ENDIF
      RETURN
204   DO 3090 J=1,N
        IF(J.EQ.I.AND.I/2*2.EQ.I) THEN
           G(J)=10.0D 0
        ELSE IF(J.EQ.I) THEN
          G(J)=-1.0D 0
        ELSE IF(J.EQ.I-1.AND.I/2*2.EQ.I) THEN
           G(J)=-20.0D 0*X(I-1)
        ELSE
           G(J)=0.0D 0
        ENDIF
3090  CONTINUE
      RETURN
205   DO 4000 J=1,N
        IF(J.EQ.I) THEN
           G(J)=1.0D 0-3.0D 0*X(I)**2/(2.0D 0*DBLE(N))
        ELSE
           G(J)=-3.0D 0*X(J)**2/(2.0D 0*DBLE(N))
        ENDIF
4000  CONTINUE
      RETURN
206   DO 4010 J=1,N
        IF(J.EQ.I) THEN
           G(J)=-2.0D 0-1.0D 0/DBLE(N+1)**2*EXP(X(I))
        ELSE IF(J.EQ.I-1.OR.J.EQ.I+1) THEN
           G(J)=1.0D 0
        ELSE
           G(J)=0.0D 0
        ENDIF
4010  CONTINUE
      RETURN
207   S=0.1D 0
      DO 4020 J=1,N
        IF(J.EQ.I-1) THEN
           G(J)=-1.0D 0
        ELSE IF(J.EQ.I) THEN
           G(J)=3.0D 0-2.0D 0*S*X(I)
        ELSE IF(J.EQ.I+1) THEN
           G(J)=-2.0D 0
        ELSE
           G(J)=0.0D 0
        ENDIF
4020  CONTINUE
      RETURN
208   S1=1.0D 0
      S2=1.0D 0
      S3=1.0D 0
      J1=3
      J2=3
      IF(I-J1.GT.1) THEN
        I1=I-J1
      ELSE
        I1=1
      ENDIF
      IF(I+J2.LT.N) THEN
        I2=I+J2
      ELSE
        I2=N
      ENDIF
      DO 4030 J=1,N
        IF(J.EQ.I) THEN
           G(J)=S1+3.0D 0*S2*X(I)**2
        ELSE IF(J.GE.I1.AND.J.LE.I2) THEN
           G(J)=-S3*(1.0D 0+2.0D 0*X(J))
        ELSE
           G(J)=0.0D 0
        ENDIF
4030  CONTINUE
      RETURN
209   DO 4040 J=1,N
        IF(I.EQ.1.AND.J.EQ.1) THEN
           G(J)=2.0D 0*X(1)
        ELSE IF(I.EQ.1) THEN
           G(J)=0.0D 0
        ELSE IF(J.EQ.I-1) THEN
           G(J)=2.0D 0*X(I-1)
        ELSE IF(J.EQ.I) THEN
           G(J)=1.0D 0/X(I)
        ELSE
           G(J)=0.0D 0
        ENDIF
4040  CONTINUE
      RETURN
210   S1=0.0D 0
      DO 4050 J=1,N
        S1=S1+X(J)
4050  CONTINUE
      DO 4060 J=1,N
        S2=DBLE(I)
        S3=S2*S1
        G(J)=-S2*SIN(S3)*EXP(COS(S3))
4060  CONTINUE
      RETURN
211   DO 4070 J=1,N
        G(J)=3.0D 0*X(J)**2/(DBLE(N)*2.0D 0)
4070  CONTINUE
      RETURN
212   DO 4080 J=1,N
        IF(J.EQ.I) THEN
           G(J)=1.0D 0
        ELSE IF(J.EQ.I-1) THEN
           G(J)=-SIN(X(I-1))
        ELSE
           G(J)=0.0D 0
        ENDIF
4080  CONTINUE
      RETURN
213   S=(1.0D 0/DBLE(N+1))**2
      DO 4090 J=1,N
        IF(J.EQ.I-1.OR.J.EQ.I+1) THEN
           G(J)=-1.0D 0
        ELSE IF(J.EQ.I) THEN
           G(J)=2.0D 0+S*(1.0D 0+COS(X(J)))
        ELSE
           G(J)=0.0D 0
        ENDIF
4090  CONTINUE
      RETURN
214   IF(I-5.GT.1) THEN
        I1=I-5
      ELSE
        I1=1
      ENDIF
      IF(I+1.LT.N) THEN
        I2=I+1
      ELSE
        I2=N
      ENDIF
      DO 5000 J=1,N
        IF(J.EQ.I) THEN
           G(J)=2.0D 0+15.0D 0*X(J)**2
        ELSE IF(J.GE.I1.AND.J.LE.I2) THEN
           G(J)=-1.0D 0-2.0D 0*X(J)
        ELSE
           G(J)=0.0D 0
        ENDIF
5000  CONTINUE
      RETURN
215   DO 5010 J=1,N
        S3=0.0D 0
        DO 5011 K=1,29
           S1=0.0D 0
           T1=0.D 0
           S=DBLE(K)/29.0D 0
           DO 5012 L=1,N
              S1=S1+S**(L-1)*X(L)
5012       CONTINUE
           DO 5013 L=2,N
              T1=T1+S**(L-2)*(L-1)*X(L)
5013       CONTINUE
         S3=S3+S**(I+J-2)*((DBLE(29*(I-1))/DBLE(K)-2.0D 0*S1)*
     +   (DBLE(29*(J - 1))/DBLE(K)-2.0D 0*S1)-
     +   2.0D 0*(T1-S1*S1-1.0D 0))
         IF ((I.EQ.1).AND.(J.EQ.1)) THEN
             G(J)=S3+6.0D 0*X(1)*X(1)-2.0D 0*X(2)+3.0D 0
         ELSE IF (((I.EQ.1).AND.(J.EQ.2)).OR.((I.EQ.2)
     +    .AND.(J.EQ.1))) THEN
            G(J)=S3-2.0D 0*X(1)
         ELSE IF ((I.EQ.2).AND.(J.EQ.2)) THEN
           G(J)=S3+1.0D 0
         ELSE
           G(J)=S3
         ENDIF
5011   CONTINUE
5010  CONTINUE
      RETURN
216   DO 5020 J=1,N
            SI=2.D 0*X(J)-1.D 0
            V=2.D 0*SI
            TIM1=2.D 0
            TI=2.D 0
            TIM2=0.D 0
            SIM1=1.D 0
            DO 5021 L=2,I
              TI=4.D 0*SI+V*TIM1-TIM2
              TIM2=TIM1
              TIM1=TI
              SIM2=SIM1
              SIM1=SI
              SI=V*SIM1-SIM2
5021    CONTINUE
5022  G(J)=TI/DBLE(N)
5020  CONTINUE
      RETURN
217   IF(I.EQ.1) THEN
        S1=0.0D 0
      ELSE
        S1=X(I-1)
      ENDIF
      IF(I.EQ.N) THEN
        S2=20.0D 0
      ELSE
        S2=X(I+1)
      ENDIF
      DO 5030 J=1,N
        IF(J.EQ.I) THEN
           G(J)=3.0D 0*(S2-4.D 0*X(I)+S1)
        ELSE IF(J.EQ.I-1) THEN
           G(J)=3.0D 0*X(I)-(S2-S1)/2.0D 0
        ELSE IF(J.EQ.I+1) THEN
           G(J)=3.0D 0*X(I)+(S2-S1)/2.0D 0
        ELSE
           G(J)=0.0D 0
        ENDIF
5030  CONTINUE
      RETURN
218   S=1.0D 0/DBLE(N+1)
      DO 5040 J=1,N
        IF(J.EQ.I-1.OR.J.EQ.I+1) THEN
           G(J)=-1.0D 0
        ELSE IF(J.EQ.I) THEN
          G(J)=2.0D 0+3.0D 0*S*S*(1.0D 0+DBLE(I)*S+X(I))**2/2.0D 0
       ELSE
          G(J)=0.0D 0
       ENDIF
5040  CONTINUE
      RETURN
219   S=1.0D 0/DBLE(N+1)
      S1=3.0D 0*S/2.0D 0
      T1=DBLE(I)*S
      DO 5050 J=1,N
        T2=DBLE(J)*S
        IF(J.LT.I) THEN
           G(J)=S1*(1.0D 0-T1)*T2*(X(J)+T2+1.0D 0)**2
        ELSE IF(J.EQ.I) THEN
           G(J)=1.0D 0+S1*(1.0D 0-T1)*T2*(X(J)+T2+1.0D 0)**2
        ELSE
           G(J)=S1*(1.0D 0-T2)*T1*(X(J)+T2+1.0D 0)**2
        ENDIF
5050  CONTINUE
      RETURN
220   DO 5060 J=1,N
        S=SIN(X(J))
        IF(J.EQ.I) THEN
           G(J)=S+DBLE(I)*S-COS(X(J))
        ELSE
           G(J)=S
        ENDIF
5060  CONTINUE
      RETURN
221   S=0.0D 0
      DO 5070 J=1,N
        S=S+DBLE(J)*(X(J)-1.0D 0)
5070  CONTINUE
      DO 5080 J=1,N
        IF(J.NE.I) THEN
           G(J)=DBLE(I*J)*(1.0D 0+6.0D 0*S*S)
        ELSE
           G(J)=1.0D 0+DBLE(I*I)*(1.0D 0+6.0D 0*S*S)
        ENDIF
5080  CONTINUE
      RETURN
222   DO 5090 J=1,N
        IF(J.EQ.I-1) THEN
           G(J)=-1.0D 0
        ELSE IF(J.EQ.I)THEN
           G(J)=3.D 0-4.D 0*X(I)
        ELSE IF(J.EQ.I+1) THEN
           G(J)=-2.0D 0
        ELSE
           G(J)=0.0D 0
        ENDIF
5090  CONTINUE
      RETURN
223   ALF=5
      BET=14
      DO 6000 J=1,N
        IF(J.NE.I) THEN
        T=SQRT(X(J)**2+DBLE(I)/DBLE(J))
        S1=LOG(T)
        G(J)=X(J)*(SIN(S1)**ALF+COS(S1)**ALF+
     +     ALF*SIN(S1)**(ALF-1)*COS(S1)-
     +     ALF*SIN(S1)*COS(S1)**(ALF-1))/T
        ELSE
          G(J)=DBLE(BET*N)
        ENDIF
6000  CONTINUE
      RETURN
224   C=0.5D 0
      H=1.0D 0/DBLE(N)
      DO 6008 J=1,N
        BROW(J)=C*H*DBLE(I)/DBLE(2*(I+J))
6008  CONTINUE
      BROW(N)=0.5D 0*BROW(N)
      DO 6010 J=1,N
        IF(I.NE.J) THEN
           G(J)=-BROW(J)*X(I)
        ELSE
           T=1.0D 0-C*H/4.0D 0
           DO 6011 L=1,N
             IF(L.EQ.I) THEN
               T=T-2.0D 0*BROW(I)*X(I)
             ELSE
               T=T-BROW(L)*X(L)
             ENDIF
6011       CONTINUE
           G(J)=T
        ENDIF
6010  CONTINUE
      RETURN
225   DO 6020 J=1,N
        IF(I.EQ.1.AND.J.EQ.1) THEN
        G(J)=9.0D 0*X(I)**2+COS(X(I)-X(I+1))*SIN(X(I)+X(I+1))+
     +     SIN(X(I)-X(I+1))*COS(X(I)+X(I+1))
        ELSE IF(I.EQ.1.AND.J.EQ.2) THEN
           G(J)=2.0D 0-COS(X(I)-X(I+1))*SIN(X(I)+X(I+1))+
     +     SIN(X(I)-X(I+1))*COS(X(I)+X(I+1))
        ELSE IF(I.LT.N.AND.J.EQ.I) THEN
           G(J)=9.0D 0*X(I)**2+COS(X(I)-X(I+1))*SIN(X(I)+X(I+1))+
     +     SIN(X(I)-X(I+1))*COS(X(I)+X(I+1))+
     +     4.0D 0+X(I-1)*EXP(X(I-1)-X(I))
        ELSE IF(I.LT.N.AND.J.EQ.I+1) THEN
           G(J)=2.0D 0-COS(X(I)-X(I+1))*SIN(X(I)+X(I+1))+
     +     SIN(X(I)-X(I+1))*COS(X(I)+X(I+1))
        ELSE IF(J.EQ.I-1) THEN
           G(J)=-EXP(X(I-1)-X(I))-X(I-1)*EXP(X(I-1)-X(I))
        ELSE IF(I.EQ.N.AND.J.EQ.N) THEN
           G(J)=4.0D 0+X(I-1)*EXP(X(I-1)-X(I))
        ELSE
           G(J)=0.0D 0
        ENDIF
6020  CONTINUE
      RETURN
226   T=10.0D 0
      H=(T/DBLE(N+1))**2
      S=T*X(I)
      DO 6030 J=1,N
        IF(J.EQ.I-1.OR.J.EQ.I+1) THEN
           G(J)=-1.0D 0
        ELSEIF(J.EQ.I) THEN
           G(J)=2.0D 0+H*(EXP(S)+EXP(-S))/2.D 0
        ELSE
           G(J)=0.0D 0
        ENDIF
6030  CONTINUE
      RETURN
227   S=0.5D 0
      H=1.0D 0/DBLE(N+1)
      DO 6040 J=1,N
        IF(J.EQ.I-1) THEN
           G(J)=-1.0D 0+H/(2.0D 0*S)
        ELSEIF(J.EQ.I+1) THEN
           G(J)=-1.0D 0-H/(2.0D 0*S)
        ELSEIF(J.EQ.I) THEN
           G(J)=2.0D 0*(1.0D 0-H**2*X(I)/S)
        ELSE
           G(J)=0.0D 0
        ENDIF
6040  CONTINUE
      RETURN
228   S=0.5D 0
      H=1.0D 0/DBLE(N+1)
      T1=2.0D 0*H
      S1=H**2/S
      S2=1.0D 0-DBLE(I)*H
      DO 6050 J=1,I
        BROW(J)=DBLE(J)*S2
6050  CONTINUE
      IF(I.EQ.N) GO TO 6055
      DO 6051 J=I+1,N
        BROW(J)=DBLE(I)*(1.0D 0-DBLE(J)*H)
6051  CONTINUE
6055  G(1)=-S1*(BROW(1)*2.0D 0*X(1)-BROW(2)/T1)
      G(N)=-S1*(BROW(N-1)/T1+BROW(N)*2.0D 0*X(N))
      DO 6052 J=2,N-1
        G(J)=-S1*((BROW(J-1)-BROW(J+1))/T1+BROW(J)*2.0D 0*X(J))
6052  CONTINUE
      G(I)=1.0D 0+G(I)
      RETURN
229   A=-9.0D-3
      B=1.0D-3
      AL=0.0D 0
      BE=25.0D 0
      GA=20.0D 0
      CA=0.3D 0
      CB=0.3D 0
      H=(B-A)/DBLE(N+1)
      S=DBLE(I)/DBLE(N+1)
      H=H**2
      U=AL*(1.0D 0-S)+BE*S+X(I)
      FF=GA*(CB*EXP(GA*(U-BE))+CA*EXP(GA*(AL-U)))
      DO 6060 J=1,N
        IF(J.EQ.I-1.OR.J.EQ.I+1) THEN
           G(J)=-1.0D 0
        ELSEIF(J.EQ.I) THEN
           G(J)=2.0D 0+H*FF
        ELSE
           G(J)=0.0D 0
        ENDIF
6060  CONTINUE
      RETURN
230   N1=N/2
      DO 6070 J=1,N
        G(J)=0.0D 0
6070  CONTINUE
      H=1.0D 0/DBLE(N1+1)**2
      IF(I.EQ.1) THEN
        G(I)=2.0D 0+H*(2.0D 0*X(I)+1.0D 0)
        G(I+1)=-1.0D 0
        G(N1+1)=H*0.2D 0*X(N1+I)
      ELSE IF(I.EQ.N1+1) THEN
        G(1)=H*0.4D 0*X(1)
        G(I)=2.0D 0+H*(2.0D 0*X(I)+2.0D 0)
        G(I+1)=-1.0D 0
      ELSE IF(I.EQ.N1) THEN
        G(I-1)=-1.0D 0
        G(I)=2.0D 0+H*(2.0D 0*X(I)+1.0D 0)
        G(N1+I)=H*0.2D 0*X(N1+I)
      ELSE IF(I.EQ.N) THEN
        G(N1)=H*0.4D 0*X(N1)
        G(I-1)=-1.0D 0
        G(I)=2.0D 0+H*(2.0D 0*X(I)+2.0D 0)
      ELSE IF(I.LT.N1) THEN
        G(I-1)=-1.0D 0
        G(I)=2.0D 0+H*(2.0D 0*X(I)+1.0D 0)
        G(I+1)=-1.0D 0
        G(N1+I)=H*0.2D 0*X(N1+I)
      ELSE
         G(I-N1)=H*0.4D 0*X(I-N1)
         G(I-1)=-1.0D 0
         G(I)=2.0D 0+H*(2.0D 0*X(I)+2.0D 0)
         G(I+1)=-1.0D 0
      ENDIF
      RETURN
231   ND=INT(SQRT(DBLE(N)))
      L=MOD(I,ND)
      IF(L.EQ.0) THEN
         K=I/ND
         L=ND
      ELSE
         K=INT(I/ND)+1
      ENDIF
      LA=1.0D 0
      H=1.0D 0/DBLE(ND+1)
      H2=LA*H*H
      DO 2010 J=1,N
        G(J)=0.0D 0
2010  CONTINUE
      IF(L.EQ.1.AND.K.EQ.1) THEN
         G(1)=4.0D 0+H2*EXP(X(1))
         G(2)=-1.0D 0
         G(ND+1)=-1.0D 0
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.1) THEN
         G(L)=4.0D 0+H2*EXP(X(L))
         G(L-1)=-1.0D 0
         G(L+1)=-1.0D 0
         G(L+ND)=-1.0D 0
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.1) THEN
         G(ND)=4.0D 0+H2*EXP(X(ND))
         G(ND-1)=-1.0D 0
         G(ND+ND)=-1.0D 0
      ENDIF
      IF(L.EQ.1.AND.1.LT.K.AND.K.LT.ND) THEN
         G(I)=4.0D 0+H2*EXP(X(I))
         G(I-ND)=-1.0D 0
         G(I+1)=-1.0D 0
         G(I+ND)=-1.0D 0
      ENDIF
      IF(L.EQ.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         G(I)=4.0D 0+H2*EXP(X(I))
         G(I-ND)=-1.0D 0
         G(I-1)=-1.0D 0
         G(I+ND)=-1.0D 0
      ENDIF
      IF(L.EQ.1.AND.K.EQ.ND) THEN
         G(I)=4.0D 0+H2*EXP(X(I))
         G(I-ND)=-1.0D 0
         G(I+1)=-1.0D 0
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.ND) THEN
         G(I)=4.0D 0+H2*EXP(X(I))
         G(I-ND)=-1.0D 0
         G(I-1)=-1.0D 0
         G(I+1)=-1.0D 0
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.ND) THEN
         G(I)=4.0D 0+H2*EXP(X(I))
         G(I-ND)=-1.0D 0
         G(I-1)=-1.0D 0
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         G(I)=4.0D 0+H2*EXP(X(I))
         G(I-ND)=-1.0D 0
         G(I-1)=-1.0D 0
         G(I+1)=-1.0D 0
         G(I+ND)=-1.0D 0
      ENDIF
      RETURN
232   ND=INT(SQRT(DBLE(N)))
      L=MOD(I,ND)
      IF(L.EQ.0) THEN
         K=I/ND
         L=ND
      ELSE
         K=INT(I/ND)+1
      ENDIF
      H=1.0D 0/DBLE(ND+1)
      H2=H*H
      DO 2020 J=1,N
        G(J)=0.0D 0
2020  CONTINUE
      IF(L.EQ.1.AND.K.EQ.1) THEN
         G(1)=4.0D 0+H2*X(1)*2.0D 0
         G(2)=-1.0D 0
         G(ND+1)=-1.0D 0
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.1) THEN
         G(L)=4.0D 0+H2*X(L)*2.0D 0
         G(L-1)=-1.0D 0
         G(L+1)=-1.0D 0
         G(L+ND)=-1.0D 0
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.1) THEN
         G(ND)=4.0D 0+H2*X(ND)*2.0D 0
         G(ND-1)=-1.0D 0
         G(ND+ND)=-1.0D 0
      ENDIF
      IF(L.EQ.1.AND.1.LT.K.AND.K.LT.ND) THEN
         G(I)=4.0D 0+H2*X(I)*2.0D 0
         G(I-ND)=-1.0D 0
         G(I+1)=-1.0D 0
         G(I+ND)=-1.0D 0
      ENDIF
      IF(L.EQ.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         G(I)=4.0D 0+H2*X(I)*2.0D 0
         G(I-ND)=-1.0D 0
         G(I-1)=-1.0D 0
         G(I+ND)=-1.0D 0
      ENDIF
      IF(L.EQ.1.AND.K.EQ.ND) THEN
         G(I)=4.0D 0+H2*X(I)*2.0D 0
         G(I-ND)=-1.0D 0
         G(I+1)=-1.0D 0
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.ND) THEN
         G(I)=4.0D 0+H2*X(I)*2.0D 0
         G(I-ND)=-1.0D 0
         G(I-1)=-1.0D 0
         G(I+1)=-1.0D 0
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.ND) THEN
         G(I)=4.0D 0+H2*X(I)*2.0D 0
         G(I-ND)=-1.0D 0
         G(I-1)=-1.0D 0
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         G(I)=4.0D 0+H2*X(I)*2.0D 0
         G(I-ND)=-1.0D 0
         G(I-1)=-1.0D 0
         G(I+1)=-1.0D 0
         G(I+ND)=-1.0D 0
      ENDIF
      RETURN
      END
* SUBROUTINE TIUD28                ALL SYSTEMS                92/12/01
C PORTABILITY : ALL SYSTEMS
C 92/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  INITIAL VALUES OF THE VARIABLES FOR UNCONSTRAINED MINIMIZATION.
*  UNCONSTRAINED AND DENSE VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  NA  NUMBER OF APPROXIMATED FUNCTIONS.
*  RO  X(N)  VECTOR OF VARIABLES.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*  IO  IERR  ERROR INDICATOR.
*
      SUBROUTINE TIUD28(N,X,FMIN,XMAX,NEXT,IERR)
      INTEGER N,NA,NEXT,IERR
      DOUBLE PRECISION X(N),FMIN,XMAX
      DOUBLE PRECISION P,Q,ALF,BET,GAM,F,S,S1,T,Z(1000)
      INTEGER I,J,K,M,N1
      DOUBLE PRECISION Y(20),PAR
      COMMON /EMPR28/ Y,PAR,NA,M
      DOUBLE PRECISION ETA9
      PARAMETER (ETA9=1.0D 60)
      FMIN=0.0D 0
      XMAX=1.0D 3
      IERR=0
      GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,
     &  170,180,190,200,210,220,230,250,310,320,330,350,370,390,400,
     &  450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,
     &  600,610,620,630,720,740,750,760,780,790,810,830,840,860,870,
     &  880,900,910,920,930,940,950,960,970,980,990,800,240,410,420,
     &  650,660,670,680,690,340,360,380,430,440,270,280,290,300,710,
     &  820),NEXT
   10 IF (N.LT.2) GO TO 999
      N=N-MOD(N,2)
      DO 11 I=1,N
        IF(MOD(I,2).EQ.1) THEN
          X(I)=-1.2D 0
        ELSE
          X(I)=1.0D 0
        ENDIF
   11 CONTINUE
      RETURN
   20 IF (N.LT.4) GO TO 999
      N=N-MOD(N,2)
      DO 21 I=1,N
        IF(MOD(I,2).EQ.1) THEN
          X(I)=-2.0D 0
          IF(I.LE.4) X(I)=-3.0D 0
        ELSE
          X(I)=0.0D 0
          IF(I.LE.4) X(I)=-1.0D 0
        ENDIF
   21 CONTINUE
      RETURN
   30 IF (N.LT.4) GO TO 999
      N=N-MOD(N,2)
      DO 31 I=1,N
        IF(MOD(I,4).EQ.1) THEN
          X(I)=3.0D 0
        ELSEIF(MOD(I,4).EQ.2) THEN
          X(I)=-1.0D 0
        ELSEIF(MOD(I,4).EQ.3) THEN
          X(I)=0.0D 0
        ELSE
          X(I)=1.0D 0
        ENDIF
   31 CONTINUE
      RETURN
   40 IF (N.LT.4) GO TO 999
      N=N-MOD(N,2)
      DO 41 I=1,N
        X(I)=2.0D 0
   41 CONTINUE
      X(1)=1.0D 0
      RETURN
   50 IF (N.LT.3) GO TO 999
      DO 51 I=1,N
        X(I)=-1.0D 0
   51 CONTINUE
      RETURN
   60 IF (N.LT.7) GO TO 999
      DO 61 I=1,N
        X(I)=-1.0D 0
   61 CONTINUE
      RETURN
   70 IF (N.LT.4) GO TO 999
      N=N-MOD(N,2)
      DO 71 I=1,N
        X(I)=-1.0D 0
   71 CONTINUE
      RETURN
   80 IF (N.LT.6) GO TO 999
      DO 81 I=1,N
        X(I)=1.0D 0/DBLE(N)
   81 CONTINUE
      RETURN
   90 IF (N.LT.6) GO TO 999
      DO 91 I=1,N
        X(I)=1.0D 0/DBLE(N)
   91 CONTINUE
      FMIN=-ETA9
      RETURN
  100 IF (N.LT.6) GO TO 999
      DO 101 I=1,N
        X(I)=1.0D 0
  101 CONTINUE
      FMIN=-ETA9
      RETURN
  110 IF (N.LT.5) GO TO 999
      N=N-MOD(N,5)
      DO 111 I=0,N-5,5
        X(I+1)=-1.0D 0
        X(I+2)=-1.0D 0
        X(I+3)=2.0D 0
        X(I+4)=-1.0D 0
        X(I+5)=-1.0D 0
  111 CONTINUE
      X(1)=-2.0D 0
      X(2)=2.0D 0
      XMAX=1.0D 0
      RETURN
  120 IF (N.LT.2) GO TO 999
      N=N-MOD(N,2)
      DO 121 I=2,N,2
        X(I-1)=0.0D 0
        X(I)=-1.0D 0
  121 CONTINUE
      XMAX=1.0D 1
      RETURN
  130 IF (N.LT.2) GO TO 999
      N=N-MOD(N,2)
      DO 131 I=2,N,2
        X(I-1)=-1.0D 0
        X(I)=1.0D 0
  131 CONTINUE
      XMAX=1.0D 1
      RETURN
  140 IF (N.LT.3) GO TO 999
      P=1.0D 0/DBLE(N+1)
      DO 141 I=1,N
        Q=P*DBLE(I)
        X(I)=Q*(Q-1.0D 0)
  141 CONTINUE
      RETURN
  150 IF (N.LT.3) GO TO 999
      P=1.0D 0/DBLE(N+1)
      DO 151 I=1,N
        Q=P*DBLE(I)
        X(I)=0.1D 0*Q*DBLE(N+1-I)
  151 CONTINUE
      FMIN=-ETA9
      RETURN
  160 IF (N.LT.3) GO TO 999
      DO 161 I=1,N
        X(I)=1.0D 0
  161 CONTINUE
      FMIN=-ETA9
      RETURN
  170 IF (N.LT.3) GO TO 999
      P=1.0D 0/DBLE(N+1)
      DO 171 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  171 CONTINUE
      FMIN=-ETA9
      RETURN
  180 IF (N.LT.3) GO TO 999
      P=1.0D 0/DBLE(N+1)
      DO 181 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  181 CONTINUE
      FMIN=-ETA9
      RETURN
  190 IF (N.LT.3) GO TO 999
      P=EXP(2.0D 0)/DBLE(N+1)
      DO 191 I=1,N
        X(I)=(P*DBLE(I)+1.0D 0)/3.0D 0
  191 CONTINUE
      FMIN=-ETA9
      RETURN
  200 IF (N.LT.3) GO TO 999
      P=1.0D 0/DBLE(N+1)
      DO 201 I=1,N
        X(I)=P*DBLE(I)
  201 CONTINUE
      FMIN=-ETA9
      RETURN
  210 IF (N.LT.3) GO TO 999
      P=1.0D 0/DBLE(N+1)
      DO 211 I=1,N
        X(I)=P*DBLE(I)+1.0D 0
  211 CONTINUE
      FMIN=-ETA9
      RETURN
  220 IF (N.LT.3) GO TO 999
      P=1.0D 0/DBLE(N+1)
      DO 221 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  221 CONTINUE
      FMIN=-ETA9
      RETURN
  230 IF (N.LT.3) GO TO 999
      P=1.0D 0/DBLE(N+1)
      DO 231 I=1,N
        X(I)=DBLE(I)*P
  231 CONTINUE
      RETURN
  250 IF (N.LT.3) GO TO 999
      P=1.0D 0/DBLE(N+1)
      DO 251 I=1,N
        X(I)=DBLE(I)*P
  251 CONTINUE
      XMAX=1.0D 0
      RETURN
  310 IF(N.GE.2) THEN
      N=N-MOD(N,2)
      NA=N
      DO 311 I=1,N
      IF(MOD(I,2).EQ.1) THEN
      X(I)=-1.2D 0
      ELSE
      X(I)=1.0D 0
      ENDIF
  311 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  320 IF(N.GE.4) THEN
      N=N-MOD(N,4)
      NA=N
      DO 321 I=1,N
      IF (MOD(I,4).EQ.1) THEN
      X(I)=3.0D 0
      ELSE IF(MOD(I,4).EQ.2) THEN
      X(I)=-1.0D 0
      ELSEIF(MOD(I,4).EQ.3) THEN
      X(I)=0.0D 0
      ELSE
      X(I)=1.0D 0
      ENDIF
  321 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  330 IF(N.GE.2) THEN
      NA=N+1
      DO 331 I=1,N
      X(I)=DBLE(I)
  331 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  350 IF(N.GE.2) THEN
      NA=N+2
      DO 351 I=1,N
      X(I)=1.0D 0-DBLE(I)/DBLE(N)
  351 CONTINUE
      XMAX=1.0D 2
      ELSE
      IERR=1
      ENDIF
      RETURN
  370 IF(N.GE.2) THEN
      NA=N
      DO 371 I=1,N
      X(I)=0.5D 0
  371 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  390 IF(N.GE.2) THEN
      NA=N
      DO 391 I=1,N
      X(I)=DBLE(I)/DBLE(N+1)
      X(I)=X(I)*(X(I)-1.0D 0)
  391 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  400 IF(N.GE.2) THEN
      NA=N
      DO 401 I=1,N
      X(I)=-1.0D 0
  401 CONTINUE
      ELSE
      IERR=1
      ENDIF
  450 IF (N.LT.3) GO TO 999
      DO 451 I=1,N
        X(I)=-1.0D 0
  451 CONTINUE
      NA=N
      RETURN
  460 IF (N.LT.6) GO TO 999
      DO 461 I=1,N
        X(I)=-1.0D 0
  461 CONTINUE
      NA=N
      RETURN
  470 IF (N.LT.2) GO TO 999
      DO 471 I=1,N-1
        X(I)=0.5D 0
  471 CONTINUE
      X(N)=-2.0D 0
      NA=2*(N-1)
      RETURN
  480 IF (N.LT.4) GO TO 999
      N=N-MOD(N,4)
      DO 481 I=1,N
        X(I)=SIN(DBLE(I))**2
  481 CONTINUE
      NA=5*N
      RETURN
  490 IF (N.LT.4) GO TO 999
      N=N-MOD(N,2)
      DO 491 I=1,N
        X(I)=5.0D 0
  491 CONTINUE
      NA=3*(N-2)
      RETURN
  500 IF (N.LT.2) GO TO 999
      DO 501 I=1,N
        X(I)=0.2D 0
  501 CONTINUE
      NA=2*(N-1)
      RETURN
  510 CONTINUE
      IF (N.LT.2) GO TO 999
      N=N-MOD(N,2)
      DO 511 I=1,N
        IF(MOD(I,2).EQ.1) THEN
          X(I)=-0.8D 0
        ELSE
          X(I)=-0.8D 0
        ENDIF
  511 CONTINUE
      NA=2*(N-1)
      RETURN
  520 CONTINUE
      IF (N.LT.5) GO TO 999
      IF (MOD(N-5,3).NE.0) N=N-MOD(N-5,3)
      DO 521 I=1,N
        X(I)=-1.0D 0
  521 CONTINUE
      NA=6*((N-5)/3+1)
      RETURN
  530 CONTINUE
      IF (N.LT.5) GO TO 999
      IF (MOD(N-5,3).NE.0) N=N-MOD(N-5,3)
      DO 531 I=1,N
        X(I)=-1.0D 0
  531 CONTINUE
      NA=7*((N-5)/3+1)
      RETURN
  540 CONTINUE
      IF (N.LT.4) GO TO 999
      DO 541 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=-0.8D 0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)= 1.2D 0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=-1.2D 0
        ELSE
          X(I)= 0.8D 0
        ENDIF
  541 CONTINUE
      Y(1)=14.4D 0
      Y(2)=6.8D 0
      Y(3)=4.2D 0
      Y(4)=3.2D 0
  542 IF (MOD(N-4,2).NE.0) N=N-MOD(N-4,2)
      NA=4*((N-4)/2+1)
      RETURN
  550 CONTINUE
      IF (N.LT.4) GO TO 999
      DO 551 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=-0.8D 0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)= 1.2D 0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=-1.2D 0
        ELSE
          X(I)= 0.8D 0
        ENDIF
  551 CONTINUE
      Y(1)=35.8D 0
      Y(2)=11.2D 0
      Y(3)=6.2D 0
      Y(4)=4.4D 0
      GO TO 542
  560 CONTINUE
      IF (N.LT.4) GO TO 999
      DO 561 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=-0.8D 0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)= 1.2D 0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=-1.2D 0
        ELSE
          X(I)= 0.8D 0
        ENDIF
  561 CONTINUE
      Y(1)=30.6D 0
      Y(2)=72.2D 0
      Y(3)=124.4D 0
      Y(4)=187.4D 0
      GO TO 542
  570 IF (N.LT.4) GO TO 999
      N=N-MOD(N,2)
      NA=N
      DO 571 I=1,N
        IF (MOD(I,8).EQ.1) X(I)=1.0D-1
        IF (MOD(I,8).EQ.2.OR.MOD(I,8).EQ.0) X(I)=2.0D-1
        IF (MOD(I,8).EQ.3.OR.MOD(I,8).EQ.7) X(I)=3.0D-1
        IF (MOD(I,8).EQ.4.OR.MOD(I,8).EQ.6) X(I)=4.0D-1
        IF (MOD(I,8).EQ.5) X(I)=5.0D-1
  571 CONTINUE
      RETURN
  580 IF (N.LT.3) GO TO 999
      NA=N
      DO 581 I=1,N
        X(I)=1.2D 1
  581 CONTINUE
      XMAX=1.0D 1
      RETURN
  590 IF (N.LT.7) GO TO 999
      NA=N
      DO 591 I=1,N
        X(I)=-1.0D 0
  591 CONTINUE
      RETURN
  600 IF (N.LT.3) GO TO 999
      NA=N
      DO 601 I=1,N
        X(I)=DBLE(I)/DBLE(N+1)
        X(I)=X(I)*(X(I)-1.0D 0)
  601 CONTINUE
      RETURN
  610 CONTINUE
      IF (N.LT.5) GO TO 999
      IF (MOD(N-5,3).NE.0) N=N-MOD(N-5,3)
      DO 611 I=1,N
        X(I)=-1.0D 0
  611 CONTINUE
      NA=7*((N-5)/3+1)
      RETURN
  620 IF (N.LT.3) GO TO 999
      DO 621 I=1,N
        IF(MOD(I,2).EQ.1) THEN
          X(I)=-1.2D 0
        ELSE
          X(I)= 1.0D 0
        ENDIF
  621 CONTINUE
      NA=2*(N-1)
      RETURN
  630 IF (N.LT.14) GO TO 999
      N=N/2
      DO 631 I=1,N
        X(I)=5.0D 0
  631 CONTINUE
      Y(1)=SIN(1.0D 0)
      NA=13*(N-6)
      RETURN
  720 IF(N.GE.5) THEN
      N=N-MOD(N,2)
      NA=N
      DO 721 I=1,N
      IF (MOD(I,8).EQ.1) X(I)=1.0D-1
      IF (MOD(I,8).EQ.2.OR.MOD(I,8).EQ.0) X(I)=2.0D-1
      IF (MOD(I,8).EQ.3.OR.MOD(I,8).EQ.7) X(I)=3.0D-1
      IF (MOD(I,8).EQ.4.OR.MOD(I,8).EQ.6) X(I)=4.0D-1
      IF (MOD(I,8).EQ.5) X(I)=5.0D-1
  721 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  740 IF(N.GE.3) THEN
      DO 741 I=1,N
      X(I)=0.0D 0
  741 CONTINUE
      NA=N
      ELSE
      IERR=1
      ENDIF
      RETURN
  750 IF(N.GE.3) THEN
      IF (MOD(N,2).NE.1) N=N-1
      NA=N
      DO 751 I=1,N
      X(I)=1.0D 0
  751 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  760 IF(N.GE.3) THEN
      DO 761 I=1,N
      X(I)=-1.0D 0
  761 CONTINUE
      NA=N
      ELSE
      IERR=1
      ENDIF
      RETURN
  780 IF(N.GE.5) THEN
      NA=N
      DO 781 I=1,N
      X(I)=-2.0D 0
  781 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  790 IF(N.GE.7) THEN
      NA=N
      DO 791 I=1,N
      X(I)=-3.0D 0
  791 CONTINUE
      XMAX=1.0D 1
      ELSE
      IERR=1
      ENDIF
      RETURN
  810 IF(N.GE.2) THEN
      N=N-MOD(N,2)
      NA=N
      DO 811 I=1,N
      IF(MOD(I,2).EQ.1) THEN
      X(I)=9.0D 1
      ELSE
      X(I)=6.0D 1
      ENDIF
  811 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  830 IF(N.GE.4) THEN
      N=N-MOD(N,4)
      NA=N
      DO 831 I=1,N
      IF (MOD(I,4).EQ.1) THEN
      X(I)=1.0D 0
      ELSE IF(MOD(I,4).EQ.2) THEN
      X(I)=2.0D 0
      ELSEIF(MOD(I,4).EQ.3) THEN
      X(I)=2.0D 0
      ELSE
      X(I)=2.0D 0
      ENDIF
  831 CONTINUE
      XMAX=1.0D 1
      ELSE
      IERR=1
      ENDIF
      RETURN
  840 IF(N.GE.3) THEN
      DO 841 I=1,N
      X(I)=-1.0D 0
  841 CONTINUE
      NA=N
      ELSE
      IERR=1
      ENDIF
      RETURN
  860 IF(N.GE.2) THEN
      N=N-MOD(N,2)
      NA=N
      DO 861 I=1,N
      IF (MOD(I,2).EQ.1) THEN
      X(I)=0.0D 0
      ELSE
      X(I)=1.0D 0
      ENDIF
  861 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  870 IF (N.GE.4) THEN
      N=N-MOD(N,4)
      NA=N
      DO 871 I=2,N,2
      X(I-1)=-3.0D 0
      X(I)=-1.0D 0
  871 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  880 IF (N.GE.3) THEN
      DO 881 I=1,N
      X(I)=1.5D 0
  881 CONTINUE
      XMAX=1.0D 0
      NA=N
      ELSE
        IERR=1
      ENDIF
      RETURN
  900 IF(N.GE.3) THEN
      DO 901 I=1,N
      X(I)=1.0D 1
  901 CONTINUE
      NA=N
      ELSE
      IERR=1
      ENDIF
      RETURN
  910 IF (N.GE.3) THEN
      DO 911 I=1,N
      X(I)=1.0D 0
  911 CONTINUE
      PAR=1.0D 1
      NA=N
      ELSE
        IERR=1
      ENDIF
      RETURN
  920 IF (N.GE.5) THEN
      PAR=5.0D 2/DBLE(N+2)
      NA=N
      DO 921 I=1,N
      X(I)=((DBLE(I)+0.5D 0)/DBLE(N+2)-0.5D 0)**2
  921 CONTINUE
      ELSE
        IERR=1
      ENDIF
      RETURN
  930 IF (N.GE.10) THEN
      N=N-MOD(N,2)
      M=N/2
      PAR=5.0D 2
      NA=N
      DO 931 I=1,M
      X(I)=(DBLE(I)/DBLE(M+1)-0.5D 0)**2
  931 CONTINUE
      DO 932 I=M+1,N
      K=I-M
      X(I)=DBLE(K)/DBLE(M+1)-0.5D 0
  932 CONTINUE
      ELSE
        IERR=1
      ENDIF
      RETURN
  940 IF (N.GE.16) THEN
      M=INT(SQRT(DBLE(N)))
      PAR=6.8D 0/DBLE(M+1)**2
      N=M*M
      K=0
      DO 942 J=1,M
      DO 941 I=1,M
      K=K+1
      X(K)=0.0D 0
  941 CONTINUE
  942 CONTINUE
      NA=N
      ELSE
      IERR=1
      ENDIF
      RETURN
  950 IF (N.GE.16) THEN
      M=INT(SQRT(DBLE(N)))
      PAR=1.0D 0/DBLE(M+1)**2
      N=M*M
      K=0
      DO 952 J=1,M
      DO 951 I=1,M
      K=K+1
      X(K)=-1.0D 0
  951 CONTINUE
  952 CONTINUE
      NA=N
      ELSE
      IERR=1
      ENDIF
      RETURN
  960 IF (N.GE.16) THEN
      M=INT(SQRT(DBLE(N)))
      PAR=1.0D 0/DBLE(M+1)**2
      N=M*M
      K=0
      DO 962 J=1,M
      DO 961 I=1,M
      K=K+1
      X(K)=0.0D 0
  961 CONTINUE
  962 CONTINUE
      NA=N
      ELSE
      IERR=1
      ENDIF
      RETURN
  970 IF (N.GE.16) THEN
      M=INT(SQRT(DBLE(N)))
      PAR=5.0D 1/DBLE(M+1)
      N=M*M
      K=0
      DO 972 J=1,M
      DO 971 I=1,M
      K=K+1
      X(K)=1.0D 0-DBLE(I)*DBLE(J)/DBLE(M+1)**2
  971 CONTINUE
  972 CONTINUE
      NA=N
      ELSE
      IERR=1
      ENDIF
      RETURN
  980 IF (N.GE.16) THEN
      M=INT(SQRT(DBLE(N)))
      PAR=1.0D 0/DBLE(M+1)
      N=M*M
      K=0
      DO 982 J=1,M
      DO 981 I=1,M
      K=K+1
      X(K)=0.0D 0
  981 CONTINUE
  982 CONTINUE
      NA=N
      ELSE
      IERR=1
      ENDIF
      RETURN
  990 IF (N.GE.16) THEN
      M=INT(SQRT(DBLE(N)))
      N=M*M
      PAR=500.0D 0/DBLE(M+2)**4
      K=0
      DO 992 J=1,M
      DO 991 I=1,M
      K=K+1
      X(K)=0.0D 0
  991 CONTINUE
  992 CONTINUE
      NA=N
      ELSE
        IERR=1
      ENDIF
      RETURN
  800 IF (N.GE.16) THEN
      M=INT(SQRT(DBLE(N)))
      N=M*M
      PAR=500.0D 0
      K=0
      DO 802 J=1,M
      DO 801 I=1,M
      K=K+1
      X(K)=0.0D 0
  801 CONTINUE
  802 CONTINUE
      NA=N
      ELSE
        IERR=1
      ENDIF
      RETURN
  240 IF(N.GE.2) THEN
      NA=N
      DO 241 I=1,N
      IF (MOD(I,2).EQ.1) THEN
      X(I)=1.0D 0
      ELSE
      X(I)=3.0D 0
      ENDIF
  241 CONTINUE
      ELSE
      IERR=1
      ENDIF
      RETURN
  410 N1=N-1
      DO 411 I=1,N1
        X(I)=-1.2D 0
  411 CONTINUE
      X(N)=-1.0D 0
      RETURN
  420 DO 421 I=1,N
        X(I)=2.0D 0
  421 CONTINUE
      RETURN
  650 DO 651 I=1,N
        X(I)=1.5D 0
  651 CONTINUE
      RETURN
  660 DO 661 I=1,N
        X(I)=0.0D 0
  661 CONTINUE
      RETURN
  670 DO 671 I=1,N
        X(I)=-1.0D 0
  671 CONTINUE
      RETURN
  680 DO 681 I=1,N
        X(I)=-1.0D 0
  681 CONTINUE
      RETURN
  690 DO 691 I=1,N
        X(I)=0.5D 0
  691 CONTINUE
      RETURN
  340 DO 341 I=1,N
        X(I)=0.5D 0
  341 CONTINUE
      RETURN
  360 DO 361 I=1,N
        X(I)=1.0D 0
  361 CONTINUE
      RETURN
  380 DO 381 I=1,N
        X(I)=-1.0D 0
  381 CONTINUE
      RETURN
  430 ALF=5
      BET=14
      GAM=3
      DO 431 I=1,N
        X(I)=0.0D 0
  431 CONTINUE
      DO 433 I=1,N
      F=DBLE(BET*N)*X(I)+(DBLE(I)-DBLE(N)/2.0D 0)**GAM
      DO 432 J=1,N
        IF(J.NE.I) THEN
        T=SQRT(X(J)**2+DBLE(I)/DBLE(J))
        S1=LOG(T)
        F=F+T*(SIN(S1)**ALF+COS(S1)**ALF)
      ENDIF
  432 CONTINUE
        Z(I)=-F
  433 CONTINUE
      S=DBLE(BET*N)/DBLE(BET**2*N**2-(ALF+1)**2*(N-1)**2)
      DO 434 I=1,N
        X(I)=S*Z(I)
  434 CONTINUE
      RETURN
  440 DO 441 I=1,N
        X(I)=1.0D 0
  441 CONTINUE
      RETURN
  270 DO 271 I=1,N
        X(I)=1.0D 0
  271 CONTINUE
      RETURN
  280 DO 281 I=1,N
        X(I)=1.0D 0
  281 CONTINUE
      RETURN
  290 DO 291 I=1,N
        X(I)=1.0D 0
  291 CONTINUE
      RETURN
  300 T=DBLE(2)/DBLE(N+2)
      N1=N/2
      DO 301 I=1,N1
        S=DBLE(I)*T
        X(I)=S*(1.0D 0-S)
        X(N1+I)=X(I)
  301 CONTINUE
      RETURN
  710 N1=INT(SQRT(DBLE(N)))
      N=N1*N1
      DO 711 I=1,N
        X(I)=1.0D 0
  711 CONTINUE
      RETURN
  820 N1=INT(SQRT(DBLE(N)))
      N=N1*N1
      DO 821 I=1,N
        X(I)=1.0D 0
  821 CONTINUE
      RETURN
  999 IERR=1
      RETURN
      END
* SUBROUTINE TFBU28                ALL SYSTEMS                92/12/01
C PORTABILITY : ALL SYSTEMS
C 92/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  VALUES OF MODEL FUNCTIONS FOR UNCONSTRAINED MINIMIZATION.
*  UNIVERSAL VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  NA  NUMBER OF APPROXIMATED FUNCTIONS.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  F  VALUE OF THE MODEL FUNCTION.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TFBU28(N,X,F,G,NEXT)
      use matrix_routines
      INTEGER N,NA,NEXT
      DOUBLE PRECISION X(N),G(N),F
      DOUBLE PRECISION FA
      DOUBLE PRECISION A,B,C,D,E,P,Q,R,U,V,W,ALFA,SX(1000)
      DOUBLE PRECISION A1,A2,A3,A4,EX,D1S,D2S,H,PI
      DATA PI /3.14159265358979323D 0/
      DOUBLE PRECISION GA1(2),GA2(2),GA3(6),GA4(6)
      DOUBLE PRECISION AL,AL1,AL2,ALF,BE,BE1,BE2,BET,CA,CB,FF,FG,GA,
     & GAM,H2,S,S1,S2,S3,T,T1
      INTEGER I,J,K,L,M,IA,IB,IC,I1,I2,J1,J2,KA,LA,N1,ND
      DOUBLE PRECISION Y(20),PAR
      COMMON /EMPR28/ Y,PAR,NA,M
      F=0.0D 0
      CALL MXVSET(N,0.0D 0,G)
      GO TO (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,
     &  170,180,190,200,210,220,230,250,310,320,330,350,370,390,400,
     &  450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,
     &  600,610,620,630,720,740,750,760,780,790,810,830,840,860,870,
     &  880,900,910,920,930,940,950,960,970,980,990,800,240,410,420,
     &  650,660,670,680,690,340,360,380,430,440,270,280,290,300,710,
     &  820),NEXT
   10 DO 11 J=2,N
      A=X(J-1)**2-X(J)
      B=X(J-1)-1.0D 0
      F=F+1.0D 2*A**2+B**2
      G(J-1)=G(J-1)+4.0D 2*X(J-1)*A+2.0D 0*B
      G(J)=G(J)-2.0D 2*A
   11 CONTINUE
      RETURN
   20 DO 21 J=2,N-2,2
      A=X(J-1)**2-X(J)
      B=X(J-1)-1.0D 0
      C=X(J+1)**2-X(J+2)
      D=X(J+1)-1.0D 0
      U=X(J)+X(J+2)-2.0D 0
      V=X(J)-X(J+2)
      F=F+1.0D 2*A**2+B**2+9.0D 1*C**2+D**2+
     &  1.0D 1*U**2+0.1D 0*V**2
      G(J-1)=G(J-1)+4.0D 2*X(J-1)*A+2.0D 0*B
      G(J)=G(J)-2.0D 2*A+2.0D 1*U+0.2D 0*V
      G(J+1)=G(J+1)+3.6D 2*X(J+1)*C+2.0D 0*D
      G(J+2)=G(J+2)-1.8D 2*C+2.0D 1*U-0.2D 0*V
   21 CONTINUE
      RETURN
   30 DO 31 J=2,N-2,2
      A=X(J-1)+1.0D 1*X(J)
      B=X(J+1)-X(J+2)
      C=X(J)-2.0D 0*X(J+1)
      D=X(J-1)-X(J+2)
      F=F+A**2+5.0D 0*B**2+C**4+1.0D 1*D**4
      G(J-1)=G(J-1)+2.0D 0*A+4.0D 1*D**3
      G(J)=G(J)+2.0D 1*A+4.0D 0*C**3
      G(J+1)=G(J+1)-8.0D 0*C**3+1.0D 1*B
      G(J+2)=G(J+2)-4.0D 1*D**3-1.0D 1*B
   31 CONTINUE
      RETURN
   40 DO 41 J=2,N-2,2
      A=EXP(X(J-1))
      B=A-X(J)
      D=X(J)-X(J+1)
      P=X(J+1)-X(J+2)
      C=COS(P)
      Q=SIN(P)/COS(P)
      U=X(J-1)
      V=X(J+2)-1.0D 0
      F=F+B**4+1.0D 2*D**6+Q**4+U**8+V**2
      B=4.0D 0*B**3
      D=6.0D 2*D**5
      Q=4.0D 0*Q**3/C**2
      G(J-1)=G(J-1)+A*B+8.0D 0*U**7
      G(J)=G(J)+D-B
      G(J+1)=G(J+1)+Q-D
      G(J+2)=G(J+2)+2.0D 0*V-Q
   41 CONTINUE
      RETURN
   50 P=7.0D 0/3.0D 0
      DO 51 J=1,N
      A=(3.0D 0-2.0D 0*X(J))*X(J)+1.0D 0
      IF (J.GT.1) A=A-X(J-1)
      IF (J.LT.N) A=A-X(J+1)
      F=F+ABS(A)**P
      B=P*ABS(A)**(P-1.0D 0)*SIGN(1.0D 0,A)
      G(J)=G(J)+B*(3.0D 0-4.0D 0*X(J))
      IF (J.GT.1) G(J-1)=G(J-1)-B
      IF (J.LT.N) G(J+1)=G(J+1)-B
   51 CONTINUE
      RETURN
   60 P=7.0D 0/3.0D 0
      DO 63 J=1,N
      A=(2.0D 0+5.0D 0*X(J)**2)*X(J)+1.0D 0
      DO 61 I=MAX(1,J-5),MIN(N,J+1)
      A=A+X(I)*(1.0D 0+X(I))
   61 CONTINUE
      B=P*ABS(A)**(P-1.0D 0)*SIGN(1.0D 0,A)
      F=F+ABS(A)**P
      G(J)=G(J)+B*(2.0D 0+1.5D 1*X(J)**2)
      DO 62 I=MAX(1,J-5),MIN(N,J+1)
      G(I)=G(I)+B*(1.0D 0+2.0D 0*X(I))
   62 CONTINUE
   63 CONTINUE
      RETURN
   70 P=7.0D 0/3.0D 0
      K=N/2
      DO 71 J=1,N
      A=(3.0D 0-2.0D 0*X(J))*X(J)+1.0D 0
      IF (J.GT.1) A=A-X(J-1)
      IF (J.LT.N) A=A-X(J+1)
      B=P*ABS(A)**(P-1.0D 0)*SIGN(1.0D 0,A)
      F=F+ABS(A)**P
      G(J)=G(J)+B*(3.0D 0-4.0D 0*X(J))
      IF (J.GT.1) G(J-1)=G(J-1)-B
      IF (J.LT.N) G(J+1)=G(J+1)-B
      IF (J.LE.K) THEN
      F=F+ABS(X(J)+X(J+K))**P
      A=X(J)+X(J+K)
      B=P*ABS(A)**(P-1.0D 0)*SIGN(1.0D 0,A)
      G(J)=G(J)+B
      G(J+K)=G(J+K)+B
      ENDIF
   71 CONTINUE
      RETURN
   80 K=N/2
      DO 83 J=1,N
      P=0.0D 0
      DO 81 I=J-2,J+2
      IF (I.LT.1.OR.I.GT.N) GO TO 81
      A=5.0D 0*(1.0D 0+MOD(I,5)+MOD(J,5))
      B=DBLE(I+J)/1.0D 1
      P=P+A*SIN(X(I))+B*COS(X(I))
   81 CONTINUE
      IF (J.GT.K) THEN
      I=J-K
      A=5.0D 0*(1.0D 0+MOD(I,5)+MOD(J,5))
      B=DBLE(I+J)/1.0D 1
      P=P+A*SIN(X(I))+B*COS(X(I))
      ELSE
      I=J+K
      A=5.0D 0*(1.0D 0+MOD(I,5)+MOD(J,5))
      B=DBLE(I+J)/1.0D 1
      P=P+A*SIN(X(I))+B*COS(X(I))
      ENDIF
      F=F+(DBLE(N+J)-P)**2/DBLE(N)
      P=2.0D 0*(DBLE(N+J)-P)/DBLE(N)
      DO 82 I=J-2,J+2
      IF (I.LT.1.OR.I.GT.N) GO TO 82
      A=5.0D 0*(1.0D 0+MOD(I,5)+MOD(J,5))
      B=DBLE(I+J)/1.0D 1
      G(I)=G(I)-P*(A*COS(X(I))-B*SIN(X(I)))
   82 CONTINUE
      IF (J.GT.K) THEN
      I=J-K
      A=5.0D 0*(1.0D 0+MOD(I,5)+MOD(J,5))
      B=DBLE(I+J)/1.0D 1
      G(I)=G(I)-P*(A*COS(X(I))-B*SIN(X(I)))
      ELSE
      I=J+K
      A=5.0D 0*(1.0D 0+MOD(I,5)+MOD(J,5))
      B=DBLE(I+J)/1.0D 1
      G(I)=G(I)-P*(A*COS(X(I))-B*SIN(X(I)))
      ENDIF
   83 CONTINUE
      RETURN
   90 K=N/2
      Q=1.0D 0/DBLE(N)
      DO 92 J=1,N
      P=0.0D 0
      DO 91 I=J-2,J+2
      IF (I.LT.1.OR.I.GT.N) GO TO 91
      A=5.0D 0*(1.0D 0+MOD(I,5)+MOD(J,5))
      B=DBLE(I+J)/1.0D 1
      P=P+A*SIN(X(I))+B*COS(X(I))
      G(I)=G(I)+Q*(A*COS(X(I))-B*SIN(X(I)))
   91 CONTINUE
      IF (J.GT.K) THEN
      I=J-K
      A=5.0D 0*(1.0D 0+MOD(I,5)+MOD(J,5))
      B=DBLE(I+J)/1.0D 1
      P=P+A*SIN(X(I))+B*COS(X(I))
      G(I)=G(I)+Q*(A*COS(X(I))-B*SIN(X(I)))
      ELSE
      I=J+K
      A=5.0D 0*(1.0D 0+MOD(I,5)+MOD(J,5))
      B=DBLE(I+J)/1.0D 1
      P=P+A*SIN(X(I))+B*COS(X(I))
      G(I)=G(I)+Q*(A*COS(X(I))-B*SIN(X(I)))
      ENDIF
      F=F+(P+DBLE(J)*(1.0D 0-COS(X(J))))*Q
      G(J)=G(J)+Q*DBLE(J)*SIN(X(J))
   92 CONTINUE
      RETURN
  100 K=N/2
      DO 102 J=1,N
      P=0.0D 0
      Q=1.0D 0+DBLE(J)/1.0D 1
      DO 101 I=J-2,J+2
      IF (I.LT.1.OR.I.GT.N) GO TO 101
      A=5.0D 0*(1.0D 0+MOD(I,5)+MOD(J,5))
      B=1.0D 0+DBLE(I)/1.0D 1
      C=DBLE(I+J)/1.0D 1
      P=P+A*SIN(Q*X(J)+B*X(I)+C)
      R=A*COS(Q*X(J)+B*X(I)+C)/DBLE(N)
      G(J)=G(J)+R*Q
      G(I)=G(I)+R*B
  101 CONTINUE
      IF (J.GT.K) THEN
      I=J-K
      A=5.0D 0*(1.0D 0+MOD(I,5)+MOD(J,5))
      B=1.0D 0+DBLE(I)/1.0D 1
      C=DBLE(I+J)/1.0D 1
      P=P+A*SIN(Q*X(J)+B*X(I)+C)
      R=A*COS(Q*X(J)+B*X(I)+C)/DBLE(N)
      G(J)=G(J)+R*Q
      G(I)=G(I)+R*B
      ELSE
      I=J+K
      A=5.0D 0*(1.0D 0+MOD(I,5)+MOD(J,5))
      B=1.0D 0+DBLE(I)/1.0D 1
      C=DBLE(I+J)/1.0D 1
      P=P+A*SIN(Q*X(J)+B*X(I)+C)
      R=A*COS(Q*X(J)+B*X(I)+C)/DBLE(N)
      G(J)=G(J)+R*Q
      G(I)=G(I)+R*B
      ENDIF
      F=F+P
  102 CONTINUE
      F=F/DBLE(N)
      RETURN
  110 P=-0.2008D-2
      Q=-0.1900D-2
      R=-0.0261D-2
      DO 112 I=0,N-5,5
      A=1.0D 0
      B=0.0D 0
      DO 111 J=1,5
      A=A*X(I+J)
      B=B+X(I+J)**2
  111 CONTINUE
      W=EXP(A)
      A=A*W
      B=B-1.0D 1-P
      C=X(I+2)*X(I+3)-5.0D 0*X(I+4)*X(I+5)-Q
      D=X(I+1)**3+X(I+2)**3+1.0D 0-R
      F=F+W+1.0D 1*(B**2+C**2+D**2)
      G(I+1)=G(I+1)+A/X(I+1)+2.0D 1*(2.0D 0*B*X(I+1)+
     &  3.0D 0*D*X(I+1)**2)
      G(I+2)=G(I+2)+A/X(I+2)+2.0D 1*(2.0D 0*B*X(I+2)+
     &  C*X(I+3)+3.0D 0*D*X(I+2)**2)
      G(I+3)=G(I+3)+A/X(I+3)+2.0D 1*(2.0D 0*B*X(I+3)+C*X(I+2))
      G(I+4)=G(I+4)+A/X(I+4)+2.0D 1*(2.0D 0*B*X(I+4)-
     &  5.0D 0*C*X(I+5))
      G(I+5)=G(I+5)+A/X(I+5)+2.0D 1*(2.0D 0*B*X(I+5)-5.0D 0*C*X(I+4))
  112 CONTINUE
      RETURN
  120 C=0.0D 0
      DO 121 J=2,N,2
      A=X(J-1)-3.0D 0
      B=X(J-1)-X(J)
      C=C+A
      F=F+1.0D-4*A**2-B+EXP(2.0D 1*B)
      G(J-1)=G(J-1)+2.0D-4*A-1.0D 0+2.0D 1*EXP(2.0D 1*B)
      G(J)=G(J)+1.0D 0-2.0D 1*EXP(2.0D 1*B)
  121 CONTINUE
      F=F+C**2
      DO 122 J=2,N,2
      G(J-1)=G(J-1)+2.0D 0*C
  122 CONTINUE
      RETURN
  130 DO 131 J=2,N,2
      A=X(J)**2
      IF (A.EQ.0.0D 0) A=1.0D-60
      B=X(J-1)**2
      IF (B.EQ.0.0D 0) B=1.0D-60
      C=A+1.0D 0
      D=B+1.0D 0
      F=F+B**C+A**D
      P=0.0D 0
      IF (A.GT.P) P=LOG(A)
      Q=0.0D 0
      IF (B.GT.Q) Q=LOG(B)
      G(J-1)=G(J-1)+2.0D 0*X(J-1)*(C*B**A+P*A**D)
      G(J)=G(J)+2.0D 0*X(J)*(D*A**B+Q*B**C)
  131 CONTINUE
      RETURN
  140 P=1.0D 0/DBLE(N+1)
      Q=0.5D 0*P**2
      DO 141 J=1,N
      A=2.0D 0*X(J)+Q*(X(J)+DBLE(J)*P+1.0D 0)**3
      IF(J.GT.1) A=A-X(J-1)
      IF(J.LT.N) A=A-X(J+1)
      F=F+A**2
      G(J)=G(J)+A*(4.0D 0+6.0D 0*Q*(X(J)+DBLE(J)*P+1.0D 0)**2.0D 0)
      IF(J.GT.1) G(J-1)=G(J-1)-2.0D 0*A
      IF(J.LT.N) G(J+1)=G(J+1)-2.0D 0*A
  141 CONTINUE
      RETURN
  150 P=1.0D 0/DBLE(N+1)
      Q=2.0D 0/P
      R=2.0D 0*P
      DO 151 J=2,N
      A=X(J-1)-X(J)
      F=F+Q*X(J-1)*A
      G(J-1)=G(J-1)+Q*(2.0D 0*X(J-1)-X(J))
      G(J)=G(J)-Q*X(J-1)
      IF (ABS(A).LE.1.0D-6) THEN
      F=F+R*EXP(X(J))*(1.0D 0+A/2.0D 0*(1.0D 0+A/3.0D 0*(1.0D 0+
     &  A/4.0D 0)))
      G(J-1)=G(J-1)+R*EXP(X(J))*(1.0D 0/2.0D 0+A*(1.0D 0/3.0D 0+
     &  A/8.0D 0))
      G(J)=G(J)+R*EXP(X(J))*(1.0D 0/2.0D 0+A*(1.0D 0/6.0D 0+
     &  A/24.0D 0))
      ELSE
      B=EXP(X(J-1))-EXP(X(J))
      F=F+R*B/A
      G(J-1)=G(J-1)+R*(EXP(X(J-1))*A-B)/A**2
      G(J)=G(J)-R*(EXP(X(J))*A-B)/A**2
      ENDIF
  151 CONTINUE
      F=F+Q*X(N)**2+R*(EXP(X(1))-1.0D 0)/X(1)
     &             +R*(EXP(X(N))-1.0D 0)/X(N)
      G(1)=G(1)+R*(EXP(X(1))*(X(1)-1.0D 0)+1.0D 0)/X(1)**2
      G(N)=G(N)+2.0D 0*Q*X(N)
     &         +R*(EXP(X(N))*(X(N)-1.0D 0)+1.0D 0)/X(N)**2
      RETURN
  160 DO 161 J=1,N
      A=DBLE(J)*(1.0D 0-COS(X(J)))
      IF(J.GT.1) A=A+DBLE(J)*SIN(X(J-1))
      IF(J.LT.N) A=A-DBLE(J)*SIN(X(J+1))
      F=F+A
      A=DBLE(J)*SIN(X(J))
      G(J)=G(J)+A
      IF(J.GT.1) G(J-1)=G(J-1)+DBLE(J)*COS(X(J-1))
      IF(J.LT.N) G(J+1)=G(J+1)-DBLE(J)*COS(X(J+1))
  161 CONTINUE
      RETURN
  170 P=1.0D 0/DBLE(N+1)
      DO 171 J=1,N
      IF (J.EQ.1) THEN
      F=F+0.25D 0*X(J)**2/P+1.25D-1*X(J+1)**2/P+
     & P*(EXP(X(J))-1.0D 0)
      G(J)=G(J)+0.5D 0*X(J)/P+P*EXP(X(J))
      G(J+1)=G(J+1)+0.25D 0*X(J+1)/P
      ELSE IF (J.EQ.N) THEN
      F=F+0.25D 0*X(J)**2/P+1.25D-1*X(J-1)**2/P+
     & P*(EXP(X(J))-1.0D 0)
      G(J)=G(J)+0.5D 0*X(J)/P+P*EXP(X(J))
      G(J-1)=G(J-1)+0.25D 0*X(J-1)/P
      ELSE
      F=F+1.25D-1*(X(J+1)-X(J-1))**2/P+P*(EXP(X(J))-1.0D 0)
      A=0.25D 0*(X(J+1)-X(J-1))/P
      G(J)=G(J)+P*EXP(X(J))
      G(J-1)=G(J-1)-A
      G(J+1)=G(J+1)+A
      ENDIF
  171 CONTINUE
      RETURN
  180 P=1.0D 0/DBLE(N+1)
      DO 181 J=1,N
      Q=DBLE(J)*P
      IF (J.EQ.1) THEN
      F=F+0.5D 0*X(J)**2/P+0.25D 0*X(J+1)**2/P-
     & P*(X(J)**2+2.0D 0*X(J)*Q)
      G(J)=G(J)+X(J)/P-2.0D 0*P*(X(J)+Q)
      G(J+1)=G(J+1)+0.5D 0*X(J+1)/P
      ELSE IF (J.EQ.N) THEN
      F=F+0.5D 0*X(J)**2/P+0.25D 0*X(J-1)**2/P-
     & P*(X(J)**2+2.0D 0*X(J)*Q)
      G(J)=G(J)+X(J)/P-2.0D 0*P*(X(J)+Q)
      G(J-1)=G(J-1)+0.5D 0*X(J-1)/P
      ELSE
      F=F+2.5D-1*(X(J+1)-X(J-1))**2/P-P*(X(J)**2+
     & 2.0D 0*X(J)*Q)
      A=0.5D 0*(X(J+1)-X(J-1))/P
      G(J)=G(J)-2.0D 0*P*(X(J)+Q)
      G(J-1)=G(J-1)-A
      G(J+1)=G(J+1)+A
      ENDIF
  181 CONTINUE
      RETURN
  190 P=1.0D 0/DBLE(N+1)
      DO 191 J=1,N
      Q=EXP(2.0D 0*DBLE(J)*P)
      IF (J.EQ.1) THEN
      R=1.0D 0/3.0D 0
      F=F+0.5D 0*(X(J)-R)**2/P+7.0D 0*R**2+
     & 2.5D-1*(X(J+1)-R)**2/P+P*(X(J)**2+2.0D 0*X(J)*Q)
      A=0.5D 0*(X(J+1)-R)/P
      G(J)=G(J)+2.0D 0*P*(X(J)+Q)+(X(J)-R)/P
      G(J+1)=G(J+1)+A
      ELSE IF (J.EQ.N) THEN
      R=EXP(2.0D 0)/3.0D 0
      F=F+0.5D 0*(X(J)-R)**2/P+7.0D 0*R**2+
     & 2.5D-1*(X(J-1)-R)**2/P+P*(X(J)**2+2.0D 0*X(J)*Q)
      A=0.5D 0*(X(J-1)-R)/P
      G(J)=G(J)+2.0D 0*P*(X(J)+Q)+(X(J)-R)/P
      G(J-1)=G(J-1)+A
      ELSE
      F=F+2.5D-1*(X(J+1)-X(J-1))**2/P+P*(X(J)**2+
     & 2.0D 0*X(J)*Q)
      A=0.5D 0*(X(J+1)-X(J-1))/P
      G(J)=G(J)+2.0D 0*P*(X(J)+Q)
      G(J-1)=G(J-1)-A
      G(J+1)=G(J+1)+A
      ENDIF
  191 CONTINUE
      RETURN
  200 P=1.0D 0/DBLE(N+1)
      DO 201 J=1,N
      A=EXP(-2.0D 0*X(J)**2)
      IF (J.EQ.1) THEN
      F=F+(0.5D 0*X(J)**2/P-P)+
     & (2.5D-1*X(J+1)**2/P-P)*A
      B=0.5D 0*X(J+1)/P
      G(J)=G(J)+X(J)/P-4.0D 0*X(J)*A*P*(B**2-1.0D 0)
      G(J+1)=G(J+1)+A*B
      ELSE IF (J.EQ.N) THEN
      F=F+(0.5D 0*X(J)**2/P-P)*EXP(-2.0D 0)+
     & (2.5D-1*X(J-1)**2/P-P)*A
      B=0.5D 0*X(J-1)/P
      G(J)=G(J)+X(J)/P*EXP(-2.0D 0)-4.0D 0*X(J)*A*P*(B**2-1.0D 0)
      G(J-1)=G(J-1)+A*B
      ELSE
      F=F+(2.5D-1*(X(J+1)-X(J-1))**2/P-P)*A
      B=0.5D 0*(X(J+1)-X(J-1))/P
      G(J)=G(J)-4.0D 0*X(J)*A*P*(B**2-1.0D 0)
      G(J-1)=G(J-1)-A*B
      G(J+1)=G(J+1)+A*B
      ENDIF
  201 CONTINUE
      RETURN
  210 P=1.0D 0/DBLE(N+1)
      DO 211 J=1,N
      IF (J.EQ.1) THEN
      A=0.5D 0*(X(J+1)-1.0D 0)/P
      B=(X(J)-1.0D 0)/P
      U=ATAN(A)
      V=ATAN(B)
      F=F+P*(X(J)**2+A*U-LOG(SQRT(1.0D 0+A**2)))+
     &   0.5D 0*P*(1.0D 0+B*V-LOG(SQRT(1.0D 0+B**2)))
      G(J)=G(J)+2.0D 0*P*X(J)+0.5D 0*V
      G(J+1)=G(J+1)+0.5D 0*U
      ELSE IF (J.EQ.N) THEN
      A=0.5D 0*(2.0D 0-X(J-1))/P
      B=(2.0D 0-X(J))/P
      U=ATAN(A)
      V=ATAN(B)
      F=F+P*(X(J)**2+A*U-LOG(SQRT(1.0D 0+A**2)))+
     &   0.5D 0*P*(4.0D 0+B*V-LOG(SQRT(1.0D 0+B**2)))
      G(J)=G(J)+2.0D 0*P*X(J)-0.5D 0*V
      G(J-1)=G(J-1)-0.5D 0*U
      ELSE
      A=0.5D 0*(X(J+1)-X(J-1))/P
      U=ATAN(A)
      F=F+P*(X(J)**2+A*U-LOG(SQRT(1.0D 0+A**2)))
      G(J)=G(J)+2.0D 0*P*X(J)
      G(J-1)=G(J-1)-0.5D 0*U
      G(J+1)=G(J+1)+0.5D 0*U
      ENDIF
  211 CONTINUE
      RETURN
  220 P=1.0D 0/DBLE(N+1)
      DO 221 J=1,N
      IF (J.EQ.1) THEN
      A= 0.5D 0*X(J+1)/P
      B= X(J)/P
      F=F+P*(1.0D 2*(X(J)-A**2)**2+(1.0D 0-A)**2)+
     &   0.5D 0*P*(1.0D 2*B**4+(1.0D 0-B)**2)
      G(J)=G(J)+2.0D 2*P*(X(J)-A**2)+2.0D 2*B**3-(1.0D 0-B)
      G(J+1)=G(J+1)-2.0D 2*(X(J)-A**2)*A-(1.0D 0-A)
      ELSE IF (J.EQ.N) THEN
      A=-0.5D 0*X(J-1)/P
      B=-X(J)/P
      F=F+P*(1.0D 2*(X(J)-A**2)**2+(1.0D 0-A)**2)+
     &   0.5D 0*P*(1.0D 2*B**4+(1.0D 0-B)**2)
      G(J)=G(J)+2.0D 2*P*(X(J)-A**2)-2.0D 2*B**3+(1.0D 0-B)
      G(J-1)=G(J-1)+2.0D 2*(X(J)-A**2)*A+(1.0D 0-A)
      ELSE
      A=0.5D 0*(X(J+1)-X(J-1))/P
      F=F+P*(1.0D 2*(X(J)-A**2)**2+(1.0D 0-A)**2)
      G(J)=G(J)+2.0D 2*P*(X(J)-A**2)
      G(J-1)=G(J-1)+2.0D 2*(X(J)-A**2)*A+(1.0D 0-A)
      G(J+1)=G(J+1)-2.0D 2*(X(J)-A**2)*A-(1.0D 0-A)
      ENDIF
  221 CONTINUE
      RETURN
  230 A=1.0D 0
      B=1.0D-3
      C=0.0D 0
      D=0.0D 0
      DO 231 J=1,N
      C=C+(X(J)-1.0D 0)**2
      D=D+X(J)**2
  231 CONTINUE
      F=A*C+B*(D-2.5D-1)**2
      DO 232 J=1,N
      G(J)=2.0D 0*A*(X(J)-1.0D 0)+4.0D 0*B*(D-2.5D-1)*X(J)
  232 CONTINUE
      RETURN
  250 A=1.0D 0
      B=0.0D 0
      C=0.0D 0
      D=0.0D 0
      F=0.0D 0
      U=EXP(X(N))
      V=EXP(X(N-1))
      DO 251 J=1,N
      IF (J.LE.N/2) F=F+(X(J)-1.0D 0)**2
      IF (J.LE.N-2) THEN
      B=B+(X(J)+2.0D 0*X(J+1)+1.0D 1*X(J+2)-1.0D 0)**2
      C=C+(2.0D 0*X(J)+X(J+1)-3.0D 0)**2
      ENDIF
      D=D+X(J)**2-DBLE(N)
  251 CONTINUE
      F=F+A*(1.0D 0+U*B+B*C+V*C)+D**2
      DO 252 J=1,N
      IF (J.LE.N/2) G(J)=G(J)+2.0D 0*(X(J)-1.0D 0)
      IF (J.LE.N-2) THEN
      P=A*(U+C)*(X(J)+2.0D 0*X(J+1)+1.0D 1*X(J+2)-1.0D 0)
      Q=A*(V+B)*(2.0D 0*X(J)+X(J+1)-3.0D 0)
      G(J)=G(J)+2.0D 0*P+4.0D 0*Q
      G(J+1)=G(J+1)+4.0D 0*P+2.0D 0*Q
      G(J+2)=G(J+2)+2.0D 1*P
      ENDIF
      G(J)=G(J)+4.0D 0*D*X(J)
  252 CONTINUE
      G(N-1)=G(N-1)+A*V*C
      G(N)=G(N)+A*U*B
      RETURN
  310 DO 311 KA=1,NA
      IF(MOD(KA,2).EQ.1) THEN
      FA=1.0D 1*(X(KA+1)-X(KA)**2)
      G(KA)=G(KA)-2.0D 1*X(KA)*FA
      G(KA+1)=G(KA+1)+1.0D 1*FA
      ELSE
      FA=1.0D 0-X(KA-1)
      G(KA-1)=G(KA-1)-FA
      ENDIF
      F=F+FA**2
  311 CONTINUE
      F=0.5D 0*F
      RETURN
  320 DO 321 KA=1,NA
      IF(MOD(KA,4).EQ.1) THEN
      FA=X(KA)+1.0D 1*X(KA+1)
      G(KA)=G(KA)+FA
      G(KA+1)=G(KA+1)+1.0D 1*FA
      ELSEIF(MOD(KA,4).EQ.2) THEN
      FA=2.23606797749979D 0*(X(KA+1)-X(KA+2))
      G(KA+1)=G(KA+1)+2.23606797749979D 0*FA
      G(KA+2)=G(KA+2)-2.23606797749979D 0*FA
      ELSEIF(MOD(KA,4).EQ.3) THEN
      A=X(KA-1)-2.0D 0*X(KA)
      FA=A**2
      G(KA-1)=G(KA-1)+2.0D 0*A*FA
      G(KA)=G(KA)-4.0D 0*A*FA
      ELSE
      FA=3.16227766016838D 0*(X(KA-3)-X(KA))**2
      A=2.0D 0*(X(KA-3)-X(KA))
      G(KA-3)=G(KA-3)+3.16227766016838D 0*A*FA
      G(KA)=G(KA)-3.16227766016838D 0*A*FA
      ENDIF
      F=F+FA**2
  321 CONTINUE
      F=0.5D 0*F
      RETURN
  330 DO 333 KA=1,NA
      IF(KA.LE.N) THEN
      FA=(X(KA)-1.0D 0)/3.16227766016838D 0**5
      G(KA)=G(KA)+1.0D 0/3.16227766016838D 0**5*FA
      ELSE
      FA=-0.25D 0
      DO 331 J=1,N
      FA=FA+X(J)**2
  331 CONTINUE
      DO 332 J=1,N
      G(J)=G(J)+2.0D 0*X(J)*FA
  332 CONTINUE
      ENDIF
      F=F+FA**2
  333 CONTINUE
      F=0.5D 0*F
      RETURN
  350 DO 354 KA=1,NA
      IF(KA.LE.N) THEN
      FA=X(KA)-1.0D 0
      G(KA)=G(KA)+FA
      ELSE
      FA=0.0D 0
      DO 351 J=1,N
      FA=FA+DBLE(J)*(X(J)-1.0D 0)
  351 CONTINUE
      IF(KA.EQ.N+1) THEN
      DO 352 J=1,N
      G(J)=G(J)+DBLE(J)*FA
  352 CONTINUE
      ELSE IF(KA.EQ.N+2) THEN
      DO 353 J=1,N
      G(J)=G(J)+2.0D 0*DBLE(J)*FA**3
  353 CONTINUE
      FA=FA**2
      ENDIF
      ENDIF
      F=F+FA**2
  354 CONTINUE
      F=0.5D 0*F
      RETURN
  370 DO 376 KA=1,NA
      IF(KA.LT.N) THEN
      A=0.0D 0
      DO 371 J=1,N
      A=A+X(J)
  371 CONTINUE
      FA=X(KA)+A-DBLE(N+1)
      DO 372 J=1,N
      G(J)=G(J)+FA
  372 CONTINUE
      G(KA)=G(KA)+FA
      ELSE
      A=1.0D 0
      DO 373 J=1,N
      A=A*X(J)
  373 CONTINUE
      FA=A-1.0D 0
      I=0
      DO 374 J=1,N
      B=X(J)
      IF(B.EQ.0.0D 0.AND.I.EQ.0) I=J
  374 CONTINUE
      IF(I.NE.J) A=A*B
      IF(I.EQ.0) THEN
      DO 375 J=1,N
      G(J)=G(J)+A/X(J)*FA
  375 CONTINUE
      ELSE
      G(I)=G(I)+A*FA
      ENDIF
      ENDIF
      F=F+FA**2
  376 CONTINUE
      F=0.5D 0*F
      RETURN
  390 DO 393 KA=1,NA
      U=1.0D 0/DBLE(N+1)
      V=DBLE(KA)*U
      A=0.0D 0
      B=0.0D 0
      DO 391 J=1,N
      W=DBLE(J)*U
      IF(J.LE.KA) THEN
      A=A+W*(X(J)+W+1.0D 0)**3
      ELSE
      B=B+(1.0D 0-W)*(X(J)+W+1.0D 0)**3
      ENDIF
  391 CONTINUE
      FA=X(KA)+U*((1.0D 0-V)*A+V*B)/2.0D 0
      F=F+FA**2
      DO 392 J=1,N
      W=DBLE(J)*U
      A=(X(J)+W+1.0D 0)**2
      IF(J.LE.KA) THEN
      G(J)=G(J)+1.5D 0*U*(1.0D 0-V)*W*A*FA
      ELSE
      G(J)=G(J)+1.5D 0*U*(1.0D 0-W)*V*A*FA
      ENDIF
  392 CONTINUE
      G(KA)=G(KA)+FA
  393 CONTINUE
      F=0.5D 0*F
      RETURN
  400 DO 401 KA=1,NA
      FA=(3.0D 0-2.0D 0*X(KA))*X(KA)+1.0D 0
      IF(KA.GT.1) FA=FA-X(KA-1)
      IF(KA.LT.N) FA=FA-2.0D 0*X(KA+1)
      F=F+FA**2
      G(KA)=G(KA)+(3.0D 0-4.0D 0*X(KA))*FA
      IF(KA.GT.1) G(KA-1)=G(KA-1)-FA
      IF(KA.LT.N) G(KA+1)=G(KA+1)-2.0D 0*FA
  401 CONTINUE
      F=0.5D 0*F
      RETURN
  450 DO 451 KA=1,NA
      I=KA
      FA=(3.0D 0-2.0D 0*X(I))*X(I)+1.0D 0
      IF (I.GT.1) FA=FA-X(I-1)
      IF (I.LT.N) FA=FA-X(I+1)
      F=F+FA**2
      G(I)=G(I)+(3.0D 0-4.0D 0*X(I))*FA
      IF (I.GT.1) G(I-1)=G(I-1)-FA
      IF (I.LT.N) G(I+1)=G(I+1)-FA
  451 CONTINUE
      F=0.5D 0*F
      RETURN
  460 DO 463 KA=1,NA
      I=KA
      FA=(2.0D 0+5.0D 0*X(I)**2)*X(I)+1.0D 0
      DO 461 J=MAX(1,I-5),MIN(N,I+1)
      FA=FA+X(J)*(1.0D 0+X(J))
  461 CONTINUE
      F=F+FA**2
      DO 462 J=MAX(1,I-5),MIN(N,I+1)
      G(J)=G(J)+(1.0D 0+2.0D 0*X(J))*FA
  462 CONTINUE
      G(I)=G(I)+(2.0D 0+1.5D 1*X(I)**2)*FA
  463 CONTINUE
      F=0.5D 0*F
      RETURN
  470 DO 471 KA=1,NA
      I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
      FA=X(I)+X(I+1)*((5.0D 0-X(I+1))*X(I+1)-2.0D 0)-1.3D 1
      G(I)=G(I)+FA
      G(I+1)=G(I+1)+(1.0D 1*X(I+1)-3.0D 0*X(I+1)**2-2.0D 0)*FA
      ELSE
      FA=X(I)+X(I+1)*((1.0D 0+X(I+1))*X(I+1)-1.4D 1)-2.9D 1
      G(I)=G(I)+FA
      G(I+1)=G(I+1)+(2.0D 0*X(I+1)+3.0D 0*X(I+1)**2-1.4D 1)*FA
      ENDIF
      F=F+FA**2
  471 CONTINUE
      F=0.5D 0*F
      RETURN
  480 DO 481 KA=1,NA
      I=MOD(KA,N/2)+1
      J=I+N/2
      M=5*N
      IF (KA.LE.M/2) THEN
      IA=1
      ELSE
      IA=2
      ENDIF
      IB=5-KA/(M/4)
      IC=MOD(KA,5)+1
      FA=(X(I)**IA-X(J)**IB)**IC
      F=F+FA**2
      A=DBLE(IA)
      B=DBLE(IB)
      C=DBLE(IC)
      D=X(I)**IA-X(J)**IB
      IF (D.NE.0.0D 0) THEN
      E=C*D**(IC-1)
      IF (X(I).EQ.0.0D 0.AND.IA.LE.1) THEN
      ELSE
      G(I)=G(I)+E*A*X(I)**(IA-1)*FA
      ENDIF
      IF (X(J).EQ.0.0D 0.AND.IB.LE.1) THEN
      ELSE
      G(J)=G(J)-E*B*X(J)**(IB-1)*FA
      ENDIF
      ENDIF
  481 CONTINUE
      F=0.5D 0*F
      RETURN
  490 DO 491 KA=1,NA
      I=2*((KA+5)/6)-1
      IF (MOD(KA,6).EQ.1) THEN
      FA=X(I)+3.0D 0*X(I+1)*(X(I+2)-1.0D 0)+X(I+3)**2-1.0D 0
      G(I)=G(I)+FA
      G(I+1)=G(I+1)+3.0D 0*(X(I+2)-1.0D 0)*FA
      G(I+2)=G(I+2)+3.0D 0*X(I+1)*FA
      G(I+3)=G(I+3)+2.0D 0*X(I+3)*FA
      ELSEIF (MOD(KA,6).EQ.2) THEN
      FA=(X(I)+X(I+1))**2+(X(I+2)-1.0D 0)**2-X(I+3)-3.0D 0
      G(I)=G(I)+2.0D 0*(X(I)+X(I+1))*FA
      G(I+1)=G(I+1)+2.0D 0*(X(I)+X(I+1))*FA
      G(I+2)=G(I+2)+2.0D 0*(X(I+2)-1.0D 0)*FA
      G(I+3)=G(I+3)-FA
      ELSEIF (MOD(KA,6).EQ.3) THEN
      FA=X(I)*X(I+1)-X(I+2)*X(I+3)
      G(I)=G(I)+X(I+1)*FA
      G(I+1)=G(I+1)+X(I)*FA
      G(I+2)=G(I+2)-X(I+3)*FA
      G(I+3)=G(I+3)-X(I+2)*FA
      ELSEIF (MOD(KA,6).EQ.4) THEN
      FA=2.0D 0*X(I)*X(I+2)+X(I+1)*X(I+3)-3.0D 0
      G(I)=G(I)+2.0D 0*X(I+2)*FA
      G(I+1)=G(I+1)+X(I+3)*FA
      G(I+2)=G(I+2)+2.0D 0*X(I)*FA
      G(I+3)=G(I+3)+X(I+1)*FA
      ELSEIF (MOD(KA,6).EQ.5) THEN
      FA=(X(I)+X(I+1)+X(I+2)+X(I+3))**2+(X(I)-1.0D 0)**2
      G(I)=G(I)+(2.0D 0*(X(I)+X(I+1)+X(I+2)+X(I+3))+
     & 2.0D 0*(X(I)-1.0D 0))*FA
      G(I+1)=G(I+1)+2.0D 0*(X(I)+X(I+1)+X(I+2)+X(I+3))*FA
      G(I+2)=G(I+2)+2.0D 0*(X(I)+X(I+1)+X(I+2)+X(I+3))*FA
      G(I+3)=G(I+3)+2.0D 0*(X(I)+X(I+1)+X(I+2)+X(I+3))*FA
      ELSE
      FA=X(I)*X(I+1)*X(I+2)*X(I+3)+(X(I+3)-1.0D 0)**2-1.0D 0
      G(I)=G(I)+X(I+1)*X(I+2)*X(I+3)*FA
      G(I+1)=G(I+1)+X(I)*X(I+2)*X(I+3)*FA
      G(I+2)=G(I+2)+X(I)*X(I+1)*X(I+3)*FA
      G(I+3)=G(I+3)+(X(I)*X(I+1)*X(I+2)+2.0D 0*(X(I+3)-1.0D 0))*FA
      ENDIF
      F=F+FA**2
  491 CONTINUE
      F=0.5D 0*F
      RETURN
  500 DO 501 KA=1,NA
      I=(KA+1)/2
      J=MOD(KA,2)
      IF (J.EQ.0) THEN
      FA=6.0D 0-EXP(2.0D 0*X(I))-EXP(2.0D 0*X(I+1))
      G(I)=G(I)-2.0D 0*EXP(2.0D 0*X(I))*FA
      G(I+1)=G(I+1)-2.0D 0*EXP(2.0D 0*X(I+1))*FA
      ELSEIF (I.EQ.1) THEN
      FA=4.0D 0-EXP(X(I))-EXP(X(I+1))
      G(I)=G(I)-EXP(X(I))*FA
      G(I+1)=G(I+1)-EXP(X(I+1))*FA
      ELSEIF (I.EQ.N) THEN
      FA=8.0D 0-EXP(3.0D 0*X(I-1))-EXP(3.0D 0*X(I))
      G(I-1)=G(I-1)-3.0D 0*EXP(3.0D 0*X(I-1))*FA
      G(I)=G(I)-3.0D 0*EXP(3.0D 0*X(I))*FA
      ELSE
      FA=8.0D 0-EXP(3.0D 0*X(I-1))-EXP(3.0D 0*X(I))+
     & 4.0D 0-EXP(X(I))-EXP(X(I+1))
      G(I-1)=G(I-1)-3.0D 0*EXP(3.0D 0*X(I-1))*FA
      G(I)=G(I)-(3.0D 0*EXP(3.0D 0*X(I))+EXP(X(I)))*FA
      G(I+1)=G(I+1)-EXP(X(I+1))*FA
      ENDIF
      F=F+FA**2
  501 CONTINUE
      F=0.5D 0*F
      RETURN
  510 DO 511 KA=1,NA
      I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
      FA=1.0D 1*(2.0D 0*X(I)/(1.0D 0+X(I)**2)-X(I+1))
      G(I)=G(I)+2.0D 1*(1.0D 0-X(I)**2)/(1.0D 0+X(I)**2)**2*FA
      G(I+1)=G(I+1)-1.0D 1*FA
      ELSE
      FA=X(I)-1.0D 0
      G(I)=G(I)+FA
      ENDIF
      F=F+FA**2
  511 CONTINUE
      F=0.5D 0*F
      RETURN
  520 DO 521 KA=1,NA
      I=3*((KA+5)/6)-2
      IF (MOD(KA,6).EQ.1) THEN
      FA=1.0D 1*(X(I)**2-X(I+1))
      G(I)=G(I)+2.0D 1*X(I)*FA
      G(I+1)=G(I+1)-1.0D 1*FA
      ELSEIF (MOD(KA,6).EQ.2) THEN
      FA=X(I+2)-1.0D 0
      G(I+2)=G(I+2)+FA
      ELSEIF (MOD(KA,6).EQ.3) THEN
      FA=(X(I+3)-1.0D 0)**2
      G(I+3)=G(I+3)+2.0D 0*(X(I+3)-1.0D 0)*FA
      ELSEIF (MOD(KA,6).EQ.4) THEN
      FA=(X(I+4)-1.0D 0)**3
      G(I+4)=G(I+4)+3.0D 0*(X(I+4)-1.0D 0)**2*FA
      ELSEIF (MOD(KA,6).EQ.5) THEN
      FA=X(I)**2*X(I+3)+SIN(X(I+3)-X(I+4))-1.0D 1
      G(I)=G(I)+2.0D 0*X(I)*X(I+3)*FA
      G(I+3)=G(I+3)+(X(I)**2+COS(X(I+3)-X(I+4)))*FA
      G(I+4)=G(I+4)-COS(X(I+3)-X(I+4))*FA
      ELSE
      FA=X(I+1)+(X(I+2)**2*X(I+3))**2-2.0D 1
      G(I+1)=G(I+1)+FA
      G(I+2)=G(I+2)+4.0D 0*X(I+2)*(X(I+2)*X(I+3))**2*FA
      G(I+3)=G(I+3)+2.0D 0*X(I+2)**4*X(I+3)*FA
      ENDIF
      F=F+FA**2
  521 CONTINUE
      F=0.5D 0*F
      RETURN
  530 DO 531 KA=1,NA
      I=3*((KA+6)/7)-2
      IF (MOD(KA,7).EQ.1) THEN
      FA=1.0D 1*(X(I)**2-X(I+1))
      G(I)=G(I)+2.0D 1*X(I)*FA
      G(I+1)=G(I+1)-1.0D 1*FA
      ELSEIF (MOD(KA,7).EQ.2) THEN
      FA=1.0D 1*(X(I+1)**2-X(I+2))
      G(I+1)=G(I+1)+2.0D 1*X(I+1)*FA
      G(I+2)=G(I+2)-1.0D 1*FA
      ELSEIF (MOD(KA,7).EQ.3) THEN
      FA=(X(I+2)-X(I+3))**2
      G(I+2)=G(I+2)+2.0D 0*(X(I+2)-X(I+3))*FA
      G(I+3)=G(I+3)-2.0D 0*(X(I+2)-X(I+3))*FA
      ELSEIF (MOD(KA,7).EQ.4) THEN
      FA=(X(I+3)-X(I+4))**2
      G(I+3)=G(I+3)+2.0D 0*(X(I+3)-X(I+4))*FA
      G(I+4)=G(I+4)-2.0D 0*(X(I+3)-X(I+4))*FA
      ELSEIF (MOD(KA,7).EQ.5) THEN
      FA=X(I)+X(I+1)**2+X(I+2)-3.0D 1
      G(I)=G(I)+FA
      G(I+1)=G(I+1)+2.0D 0*X(I+1)*FA
      G(I+2)=G(I+2)+FA
      ELSEIF (MOD(KA,7).EQ.6) THEN
      FA=X(I+1)-X(I+2)**2+X(I+3)-1.0D 1
      G(I+1)=G(I+1)+FA
      G(I+2)=G(I+2)-2.0D 0*X(I+2)*FA
      G(I+3)=G(I+3)+FA
      ELSE
      FA=X(I)*X(I+4)-1.0D 1
      G(I)=G(I)+X(I+4)*FA
      G(I+4)=G(I+4)+X(I)*FA
      ENDIF
      F=F+FA**2
  531 CONTINUE
      F=0.5D 0*F
      RETURN
  540 DO 546 KA=1,NA
      I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      FA=-Y(L)
      DO 542 K=1,3
      A=DBLE(K*K)/DBLE(L)
      DO 541 J=1,4
      IF (X(I+J).EQ.0) X(I+J)=1.0D-16
      A=A*SIGN(1.0D 0,X(I+J))*ABS(X(I+J))**(DBLE(J)/DBLE(K*L))
  541 CONTINUE
      FA=FA+A
  542 CONTINUE
      F=F+FA**2
      DO 545 K=1,3
      A=DBLE(K*K)/DBLE(L)
      DO 543 J=1,4
      A=A*SIGN(1.0D 0,X(I+J))*ABS(X(I+J))**(DBLE(J)/DBLE(K*L))
  543 CONTINUE
      DO 544 J=1,4
      G(I+J)=G(I+J)+(DBLE(J)/DBLE(K*L))*A/X(I+J)*FA
  544 CONTINUE
  545 CONTINUE
  546 CONTINUE
      F=0.5D 0*F
      RETURN
  550 DO 556 KA=1,NA
      I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      FA=-Y(L)
      DO 552 K=1,3
      A=0.0D 0
      DO 551 J=1,4
      A=A+X(I+J)*(DBLE(J)/DBLE(K*L))
  551 CONTINUE
      FA=FA+EXP(A)*DBLE(K*K)/DBLE(L)
  552 CONTINUE
      F=F+FA**2
      DO 555 K=1,3
      A=0.0D 0
      DO 553 J=1,4
      A=A+X(I+J)*(DBLE(J)/DBLE(K*L))
  553 CONTINUE
      A=EXP(A)*DBLE(K*K)/DBLE(L)
      DO 554 J=1,4
      G(I+J)=G(I+J)+A*(DBLE(J)/DBLE(K*L))*FA
  554 CONTINUE
  555 CONTINUE
  556 CONTINUE
      F=0.5D 0*F
      RETURN
  560 DO 563 KA=1,NA
      I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      FA=-Y(L)
      DO 561 J=1,4
      FA=FA+DBLE((1-2*MOD(J,2))*L*J*J)*SIN(X(I+J))+
     & DBLE(L*L*J)*COS(X(I+J))
  561 CONTINUE
      F=F+FA**2
      DO 562 J=1,4
      G(I+J)=G(I+J)+(DBLE((1-2*MOD(J,2))*L*J*J)*COS(X(I+J))-
     & DBLE(L*L*J)*SIN(X(I+J)))*FA
  562 CONTINUE
  563 CONTINUE
      F=0.5D 0*F
      RETURN
  570 DO 571 KA=1,NA
      ALFA=0.5D 0
      IF (KA.EQ.1) THEN
      FA=ALFA-(1.0D 0-ALFA)*X(3)-X(1)*(1.0D 0+4.0D 0*X(2))
      G(1)=G(1)-(1.0D 0+4.0D 0*X(2))*FA
      G(2)=G(2)-4.0D 0*X(1)*FA
      G(3)=G(3)+(ALFA-1.0D 0)*FA
      ELSEIF(KA.EQ.2) THEN
      FA=-(2.0D 0-ALFA)*X(4)-X(2)*(1.0D 0+4.0D 0*X(1))
      G(1)=G(1)-4.0D 0*X(2)*FA
      G(2)=G(2)-(1.0D 0+4.0D 0*X(1))*FA
      G(4)=G(4)+(ALFA-2.0D 0)*FA
      ELSEIF(KA.EQ.N-1) THEN
      FA=ALFA*X(N-3)-X(N-1)*(1.0D 0+4.0D 0*X(N))
      G(N-3)=G(N-3)+ALFA*FA
      G(N-1)=G(N-1)-(1.0D 0+4.0D 0*X(N))*FA
      G(N)=G(N)-4.0D 0*X(N-1)*FA
      ELSEIF (KA.EQ.N) THEN
      FA=ALFA*X(N-2)-(2.0D 0-ALFA)-X(N)*(1.0D 0+4.0D 0*X(N-1))
      G(N-2)=G(N-2)+ALFA*FA
      G(N-1)=G(N-1)-4.0D 0*X(N)*FA
      G(N)=G(N)-(1.0D 0+4.0D 0*X(N-1))*FA
      ELSEIF (MOD(KA,2).EQ.1) THEN
      FA=ALFA*X(KA-2)-(1.0D 0-ALFA)*X(KA+2)-
     & X(KA)*(1.0D 0+4.0D 0*X(KA+1))
      G(KA-2)=G(KA-2)+ALFA*FA
      G(KA)=G(KA)-(1.0D 0+4.0D 0*X(KA+1))*FA
      G(KA+1)=G(KA+1)-4.0D 0*X(KA)*FA
      G(KA+2)=G(KA+2)+(ALFA-1.0D 0)*FA
      ELSE
      FA=ALFA*X(KA-2)-(2.0D 0-ALFA)*X(KA+2)-
     & X(KA)*(1.0D 0+4.0D 0*X(KA-1))
      G(KA-2)=G(KA-2)+ALFA*FA
      G(KA-1)=G(KA-1)-4.0D 0*X(KA)*FA
      G(KA)=G(KA)-(1.0D 0+4.0D 0*X(KA-1))*FA
      G(KA+2)=G(KA+2)+(ALFA-2.0D 0)*FA
      ENDIF
      F=F+FA**2
  571 CONTINUE
      F=0.5D 0*F
      RETURN
  580 DO 581 KA=1,NA
      IF (KA.LT.2) THEN
      FA=4.0D 0*(X(KA)-X(KA+1)**2)
      G(KA)=G(KA)+4.0D 0*FA
      G(KA+1)=G(KA+1)-8.0D 0*X(KA+1)*FA
      ELSEIF (KA.LT.N) THEN
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))+
     & 4.0D 0*(X(KA)-X(KA+1)**2)
      G(KA-1)=G(KA-1)-8.0D 0*X(KA)*FA
      G(KA)=G(KA)+(24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+6.0D 0)*FA
      G(KA+1)=G(KA+1)-8.0D 0*X(KA+1)*FA
      ELSE
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))
      G(KA-1)=G(KA-1)-8.0D 0*X(KA)*FA
      G(KA)=G(KA)+(24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+2.0D 0)*FA
      ENDIF
      F=F+FA**2
  581 CONTINUE
      F=0.5D 0*F
      RETURN
  590 DO 591 KA=1,NA
      IF (KA.EQ.1) THEN
      FA=-2.0D 0*X(KA)**2+3.0D 0*X(KA)-2.0D 0*X(KA+1)+
     & 3.0D 0*X(N-4)-X(N-3)-X(N-2)+0.5D 0*X(N-1)-X(N)+1.0D 0
      G(N-4)=G(N-4)+3.0D 0*FA
      G(N-3)=G(N-3)-FA
      G(N-2)=G(N-2)-FA
      G(N-1)=G(N-1)+0.50D 0*FA
      G(N)=G(N)-FA
      G(KA)=G(KA)-(4.0D 0*X(KA)-3.0D 0)*FA
      G(KA+1)=G(KA+1)-2.0D 0*FA
      ELSEIF (KA.LE.N-1) THEN
      FA=-2.0D 0*X(KA)**2+3.0D 0*X(KA)-X(KA-1)-2.0D 0*X(KA+1)+
     & 3.0D 0*X(N-4)-X(N-3)-X(N-2)+0.5D 0*X(N-1)-X(N)+1.0D 0
      G(N-4)=G(N-4)+3.0D 0*FA
      G(N-3)=G(N-3)-FA
      G(N-2)=G(N-2)-FA
      G(N-1)=G(N-1)+0.50D 0*FA
      G(N)=G(N)-FA
      G(KA-1)=G(KA-1)-FA
      G(KA)=G(KA)-(4.0D 0*X(KA)-3.0D 0)*FA
      G(KA+1)=G(KA+1)-2.0D 0*FA
      ELSE
      FA=-2.0D 0*X(N)**2+3.0D 0*X(N)-X(N-1)+
     & 3.0D 0*X(N-4)-X(N-3)-X(N-2)+0.5D 0*X(N-1)-X(N)+1.0D 0
      G(N-4)=G(N-4)+3.0D 0*FA
      G(N-3)=G(N-3)-FA
      G(N-2)=G(N-2)-FA
      G(N-1)=G(N-1)+0.50D 0*FA
      G(N)=G(N)-(4.0D 0*X(N)-2.0D 0)*FA
      ENDIF
      F=F+FA**2
  591 CONTINUE
      F=0.5D 0*F
      RETURN
  600 DO 601 KA=1,NA
      U=1.0D 0/DBLE(N+1)
      V=DBLE(KA)*U
      FA=2.0D 0*X(KA)+0.5D 0*U*U*(X(KA)+V+1.0D 0)**3+1.0D 0
      IF(KA.GT.1) FA=FA-X(KA-1)
      IF(KA.LT.N) FA=FA-X(KA+1)
      F=F+FA**2
      G(KA)=G(KA)+(2.0D 0+1.5D 0*U**2*(X(KA)+V+1.0D 0)**2)*FA
      IF(KA.GT.1) G(KA-1)=G(KA-1)-FA
      IF(KA.LT.N) G(KA+1)=G(KA+1)-FA
  601 CONTINUE
      F=0.5D 0*F
      RETURN
  610 DO 611 KA=1,NA
      I=3*((KA+6)/7)-2
      IF (MOD(KA,7).EQ.1) THEN
      FA=1.0D 1*(X(I)**2-X(I+1))
      G(I)=G(I)+2.0D 1*X(I)*FA
      G(I+1)=G(I+1)-1.0D 1*FA
      ELSEIF (MOD(KA,7).EQ.2) THEN
      FA=X(I+1)+X(I+2)-2.0D 0
      G(I+1)=G(I+1)+FA
      G(I+2)=G(I+2)+FA
      ELSEIF (MOD(KA,7).EQ.3) THEN
      FA=X(I+3)-1.0D 0
      G(I+3)=G(I+3)+FA
      ELSEIF (MOD(KA,7).EQ.4) THEN
      FA=X(I+4)-1.0D 0
      G(I+4)=G(I+4)+FA
      ELSEIF (MOD(KA,7).EQ.5) THEN
      FA=X(I)+3.0D 0*X(I+1)
      G(I)=G(I)+FA
      G(I+1)=G(I+1)+3.0D 0*FA
      ELSEIF (MOD(KA,7).EQ.6) THEN
      FA=X(I+2)+X(I+3)-2.0D 0*X(I+4)
      G(I+2)=G(I+2)+FA
      G(I+3)=G(I+3)+FA
      G(I+4)=G(I+4)-2.0D 0*FA
      ELSE
      FA=1.0D 1*(X(I+1)**2-X(I+4))
      G(I+1)=G(I+1)+2.0D 1*X(I+1)*FA
      G(I+4)=G(I+4)-1.0D 1*FA
      ENDIF
      F=F+FA**2
  611 CONTINUE
      F=0.5D 0*F
      RETURN
  620 DO 621 KA=1,NA
      I=KA/2
      IF (KA.EQ.1) THEN
      FA=X(KA)-1.0D 0
      G(KA)=G(KA)+FA
      ELSE IF (MOD(KA,2).EQ.0) THEN
      FA=1.0D 1*(X(I)**2-X(I+1))
      G(I)=G(I)+2.0D 1*X(I)*FA
      G(I+1)=G(I+1)-1.0D 1*FA
      ELSE
      A=2.0D 0*EXP(-(X(I)-X(I+1))**2)
      B=EXP(-2.0D 0*(X(I+1)-X(I+2))**2)
      FA=A+B
      G(I)=G(I)-2.0D 0*(X(I)-X(I+1))*A*FA
      G(I+1)=G(I+1)+(2.0D 0*(X(I)-X(I+1))*A-4.0D 0*(X(I+1)-X(I+2))*B)*FA
      G(I+2)=G(I+2)+4.0D 0*(X(I+1)-X(I+2))*B*FA
      ENDIF
      F=F+FA**2
  621 CONTINUE
      F=0.5D 0*F
      RETURN
  630 DO 633 KA=1,NA
      IA=MIN(MAX(MOD(KA,13)-2,1),7)
      IB=(KA+12)/13
      I=IA+IB-1
      IF (IA.EQ.7) THEN
      J=IB
      ELSE
      J=IA+IB
      ENDIF
      C=3.0D 0*DBLE(IA)/1.0D 1
      A=0.0D 0
      DO 631 L=0,6
      IF (IB+L.NE.I.AND.IB+L.NE.J) A=A+SIN(X(IB+L))-Y(1)
  631 CONTINUE
      FA=(1.0D 0+COS(C))*(X(J)-SIN(X(I))-1.0D 0)+
     $ 5.0D 0*(X(I)-2.0D 0)*EXP(SIN(C*X(J)))+0.5D 0*A
      F=F+FA**2
      A=COS(C)
      B=EXP(SIN(C*X(J)))
      G(I)=G(I)-(COS(X(I))*(1.0D 0+A)-5.0D 0*B)*FA
      G(J)=G(J)+((1.0D 0+A)+5.0D 0*(X(I)-2.0D 0)*C*COS(C*X(J))*B)*FA
      DO 632 L=0,6
      IF (IB+L.NE.I.AND.IB+L.NE.J)
     & G(IB+L)=G(IB+L)+0.5D 0*COS(X(IB+L))*FA
  632 CONTINUE
  633 CONTINUE
      F=0.5D 0*F
      RETURN
  720 DO 721 KA=1,NA
      A1=0.414214D 0
      IF (KA.EQ.1) THEN
      FA=X(1)-(1.0D 0-X(1))*X(3)-A1*(1.0D 0+4.0D 0*X(2))
      G(1)=G(1)+(1.0D 0+X(3))*FA
      G(2)=G(2)-4.0D 0*A1*FA
      G(3)=G(3)-(1.0D 0-X(1))*FA
      ELSEIF (KA.EQ.2) THEN
      FA=-(1.0D 0-X(1))*X(4)-A1*(1.0D 0+4.0D 0*X(2))
      G(1)=G(1)+X(4)*FA
      G(2)=G(2)-4.0D 0*A1*FA
      G(4)=G(4)-(1.0D 0-X(1))*FA
      ELSEIF (KA.EQ.3) THEN
      FA=A1*X(1)-(1.0D 0-X(1))*X(5)-X(3)*(1.0D 0+4.0D 0*X(2))
      G(1)=G(1)+(A1+X(5))*FA
      G(2)=G(2)-4.0D 0*X(3)*FA
      G(3)=G(3)-(1.0D 0+4.0D 0*X(2))*FA
      G(5)=G(5)-(1.0D 0-X(1))*FA
      ELSEIF (KA.LE.N-2) THEN
      FA=X(1)*X(KA-2)-(1.0D 0-X(1))*X(KA+2)-
     & X(KA)*(1.0D 0+4.0D 0*X(KA-1))
      G(1)=G(1)+(X(KA-2)+X(KA+2))*FA
      G(KA-2)=G(KA-2)+X(1)*FA
      G(KA-1)=G(KA-1)-4.0D 0*X(KA)*FA
      G(KA)=G(KA)-(1.0D 0+4.0D 0*X(KA-1))*FA
      G(KA+2)=G(KA+2)-(1.0D 0-X(1))*FA
      ELSEIF (KA.EQ.N-1) THEN
      FA=X(1)*X(N-3)-X(N-1)*(1.0D 0+4.0D 0*X(N-2))
      G(1)=G(1)+X(N-3)*FA
      G(N-3)=G(N-3)+X(1)*FA
      G(N-2)=G(N-2)-4.0D 0*X(N-1)*FA
      G(N-1)=G(N-1)-(1.0D 0+4.0D 0*X(N-2))*FA
      ELSE
      FA=X(1)*X(N-2)-(1.0D 0-X(1))-X(N)*(1.0D 0+4.0D 0*X(N-1))
      G(1)=G(1)+(X(N-2)+1.0D 0)*FA
      G(N-2)=G(N-2)+X(1)*FA
      G(N-1)=G(N-1)-4.0D 0*X(N)*FA
      G(N)=G(N)-(1.0D 0+4.0D 0*X(N-1))*FA
      ENDIF
      F=F+FA**2
  721 CONTINUE
      F=0.5D 0*F
      RETURN
  740 DO 741 KA=1,NA
      IF (KA.LT.2) THEN
      FA=3.0D 0*X(KA)**3+2.0D 0*X(KA+1)-5.0D 0+
     & SIN(X(KA)-X(KA+1))*SIN(X(KA)+X(KA+1))
      D1S=COS(X(KA)-X(KA+1))*SIN(X(KA)+X(KA+1))
      D2S=SIN(X(KA)-X(KA+1))*COS(X(KA)+X(KA+1))
      G(KA)=G(KA)+(9.0D 0*X(KA)**2+D1S+D2S)*FA
      G(KA+1)=G(KA+1)+(2.0D 0-D1S+D2S)*FA
      ELSEIF (KA.LT.N) THEN
      FA=3.0D 0*X(KA)**3+2.0D 0*X(KA+1)-5.0D 0+
     & SIN(X(KA)-X(KA+1))*SIN(X(KA)+X(KA+1))+4.0D 0*X(KA)-
     & X(KA-1)*EXP(X(KA-1)-X(KA))-3.0D 0
      D1S=COS(X(KA)-X(KA+1))*SIN(X(KA)+X(KA+1))
      D2S=SIN(X(KA)-X(KA+1))*COS(X(KA)+X(KA+1))
      EX=EXP(X(KA-1)-X(KA))
      G(KA-1)=G(KA-1)-(EX+X(KA-1)*EX)*FA
      G(KA)=G(KA)+(9.0D 0*X(KA)**2+D1S+D2S+4.0D 0+X(KA-1)*EX)*FA
      G(KA+1)=G(KA+1)+(2.0D 0-D1S+D2S)*FA
      ELSE
      FA=4.0D 0*X(KA)-X(KA-1)*EXP(X(KA-1)-X(KA))-3.0D 0
      EX=EXP(X(KA-1)-X(KA))
      G(KA-1)=G(KA-1)-(EX+X(KA-1)*EX)*FA
      G(KA)=G(KA)+(4.0D 0+X(KA-1)*EX)*FA
      ENDIF
      F=F+FA**2
  741 CONTINUE
      F=0.5D 0*F
      RETURN
  750 DO 751 KA=1,NA
      IF (MOD(KA,2).EQ.1) THEN
      FA=0.0D 0
      IF (KA.NE.1) FA=FA-6.0D 0*(X(KA-2)-X(KA))**3+1.0D 1-
     & 4.0D 0*X(KA-1)-2.0D 0*SIN(X(KA-2)-X(KA-1)-X(KA))*
     & SIN(X(KA-2)+X(KA-1)-X(KA))
      IF (KA.NE.N) FA=FA+3.0D 0*(X(KA)-X(KA+2))**3-5.0D 0+
     & 2.0D 0*X(KA+1)+SIN(X(KA)-X(KA+1)-X(KA+2))*
     & SIN(X(KA)+X(KA+1)-X(KA+2))
      IF (KA.NE.1) THEN
      D1S=COS(X(KA-2)-X(KA-1)-X(KA))*SIN(X(KA-2)+X(KA-1)-X(KA))
      D2S=SIN(X(KA-2)-X(KA-1)-X(KA))*COS(X(KA-2)+X(KA-1)-X(KA))
      G(KA-2)=G(KA-2)-(18.0D 0*(X(KA-2)-X(KA))**2+2.0D 0*(D1S+D2S))*FA
      G(KA-1)=G(KA-1)-(4.0D 0-2.0D 0*(D1S-D2S))*FA
      G(KA)=G(KA)+(18.0D 0*(X(KA-2)-X(KA))**2+2.0D 0*(D1S+D2S))*FA
      ENDIF
      IF (KA.NE.N) THEN
      D1S=COS(X(KA)-X(KA+1)-X(KA+2))*SIN(X(KA)+X(KA+1)-X(KA+2))
      D2S=SIN(X(KA)-X(KA+1)-X(KA+2))*COS(X(KA)+X(KA+1)-X(KA+2))
      G(KA)=G(KA)+(9.0D 0*(X(KA)-X(KA+2))**2+D1S+D2S)*FA
      G(KA+1)=G(KA+1)+(2.0D 0-D1S+D2S)*FA
      G(KA+2)=G(KA+2)-(9.0D 0*(X(KA)-X(KA+2))**2+D1S+D2S)*FA
      ENDIF
      ELSE
      EX=EXP(X(KA-1)-X(KA)-X(KA+1))
      FA=4.0D 0*X(KA)-(X(KA-1)-X(KA+1))*EX-3.0D 0
      W=X(KA-1)-X(KA+1)
      G(KA-1)=G(KA-1)-(EX+W*EX)*FA
      G(KA)=G(KA)+(4.0D 0+W*EX)*FA
      G(KA+1)=G(KA+1)+(EX+W*EX)*FA
      ENDIF
      F=F+FA**2
  751 CONTINUE
      F=0.5D 0*F
      RETURN
  760 DO 761 KA=1,NA
      H=2.0D 0
      IF (KA.EQ.1) THEN
      FA=((3.0D 0-H*X(1))*X(1)-2.0D 0*X(2)+1.0D 0)**2
      G(1)=G(1)+2.0D 0*((3.0D 0-H*X(1))*X(1)-2.0D 0*X(2)+1.0D 0)*
     & (3.0D 0-2.0D 0*H*X(1))*FA
      G(2)=G(2)-4.0D 0*((3.0D 0-H*X(1))*X(1)-2.0D 0*X(2)+1.0D 0)*FA
      ELSEIF (KA.LE.N-1) THEN
      FA=((3.0D 0-H*X(KA))*X(KA)-X(KA-1)-2.0D 0*X(KA+1)+1.0D 0)**2
      G(KA-1)=G(KA-1)-2.0D 0*
     & ((3.0D 0-H*X(KA))*X(KA)-X(KA-1)-2.0D 0*X(KA+1)+1.0D 0)*FA
      G(KA)=G(KA)+2.0D 0*
     & ((3.0D 0-H*X(KA))*X(KA)-X(KA-1)-2.0D 0*X(KA+1)+1.0D 0)*
     & (3.0D 0-2.0D 0*H*X(KA))*FA
      G(KA+1)=G(KA+1)-4.0D 0*
     & ((3.0D 0-H*X(KA))*X(KA)-X(KA-1)-2.0D 0*X(KA+1)+1.0D 0)*FA
      ELSE
      FA=((3.0D 0-H*X(N))*X(N)-X(N-1)+1.0D 0)**2
      G(N-1)=G(N-1)-2.0D 0*((3.0D 0-H*X(N))*X(N)-X(N-1)+1.0D 0)*FA
      G(N)=G(N)+2.0D 0*((3.0D 0-H*X(N))*X(N)-X(N-1)+1.0D 0)*
     & (3.0D 0-2.0D 0*H*X(N))*FA
      ENDIF
      F=F+FA**2
  761 CONTINUE
      F=0.5D 0*F
      RETURN
  780 DO 781 KA=1,NA
      IF (KA.LT.2) THEN
      FA=4.0D 0*(X(KA)-X(KA+1)**2)+X(KA+1)-X(KA+2)**2
      G(KA)=G(KA)+4.0D 0*FA
      G(KA+1)=G(KA+1)-(8.0D 0*X(KA+1)-1.0D 0)*FA
      G(KA+2)=G(KA+2)-2.0D 0*X(KA+2)*FA
      ELSEIF (KA.LT.3) THEN
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))+
     & 4.0D 0*(X(KA)-X(KA+1)**2)+X(KA+1)-X(KA+2)**2
      G(KA-1)=G(KA-1)-8.0D 0*X(KA)*FA
      G(KA)=G(KA)+(24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+6.0D 0)*FA
      G(KA+1)=G(KA+1)-(8.0D 0*X(KA+1)-1.0D 0)*FA
      G(KA+2)=G(KA+2)-2.0D 0*X(KA+2)*FA
      ELSEIF (KA.LT.N-1) THEN
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))+
     & 4.0D 0*(X(KA)-X(KA+1)**2)+X(KA-1)**2-X(KA-2)+X(KA+1)-
     & X(KA+2)**2
      G(KA-2)=G(KA-2)-FA
      G(KA-1)=G(KA-1)-(8.0D 0*X(KA)-2.0D 0*X(KA-1))*FA
      G(KA)=G(KA)+(24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+6.0D 0)*FA
      G(KA+1)=G(KA+1)-(8.0D 0*X(KA+1)-1.0D 0)*FA
      G(KA+2)=G(KA+2)-2.0D 0*X(KA+2)*FA
      ELSEIF (KA.LT.N) THEN
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))+
     & 4.0D 0*(X(KA)-X(KA+1)**2)+X(KA-1)**2-X(KA-2)
      G(KA-2)=G(KA-2)-FA
      G(KA-1)=G(KA-1)-(8.0D 0*X(KA)-2.0D 0*X(KA-1))*FA
      G(KA)=G(KA)+(24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+6.0D 0)*FA
      G(KA+1)=G(KA+1)-8.0D 0*X(KA+1)*FA
      ELSE
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))
     & +X(KA-1)**2-X(KA-2)
      G(KA-2)=G(KA-2)-FA
      G(KA-1)=G(KA-1)-(8.0D 0*X(KA)-2.0D 0*X(KA-1))*FA
      G(KA)=G(KA)+(24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+2.0D 0)*FA
      ENDIF
      F=F+FA**2
  781 CONTINUE
      F=0.5D 0*F
      RETURN
  790 DO 791 KA=1,NA
      IF (KA.LT.2) THEN
      FA=4.0D 0*(X(KA)-X(KA+1)**2)+X(KA+1)-X(KA+2)**2+
     & X(KA+2)-X(KA+3)**2
      G(KA)=G(KA)+4.0D 0*FA
      G(KA+1)=G(KA+1)-(8.0D 0*X(KA+1)-1.0D 0)*FA
      G(KA+2)=G(KA+2)-(2.0D 0*X(KA+2)-1.0D 0)*FA
      G(KA+3)=G(KA+3)-2.0D 0*X(KA+3)*FA
      ELSEIF (KA.LT.3) THEN
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))+
     & 4.0D 0*(X(KA)-X(KA+1)**2)+X(KA-1)**2+X(KA+1)-X(KA+2)**2+
     & X(KA+2)-X(KA+3)**2
      G(KA-1)=G(KA-1)-(8.0D 0*X(KA)-2.0D 0*X(KA-1))*FA
      G(KA)=G(KA)+(24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+6.0D 0)*FA
      G(KA+1)=G(KA+1)-(8.0D 0*X(KA+1)-1.0D 0)*FA
      G(KA+2)=G(KA+2)-(2.0D 0*X(KA+2)-1.0D 0)*FA
      G(KA+3)=G(KA+3)-2.0D 0*X(KA+3)*FA
      ELSEIF (KA.LT.4) THEN
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))+
     & 4.0D 0*(X(KA)-X(KA+1)**2)+X(KA-1)**2-X(KA-2)+
     & X(KA+1)-X(KA+2)**2+X(KA-2)**2+X(KA+2)-X(KA+3)**2
      G(KA-2)=G(KA-2)+(2.0D 0*X(KA-2)-1.0D 0)*FA
      G(KA-1)=G(KA-1)-(8.0D 0*X(KA)-2.0D 0*X(KA-1))*FA
      G(KA)=G(KA)+(24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+6.0D 0)*FA
      G(KA+1)=G(KA+1)-(8.0D 0*X(KA+1)-1.0D 0)*FA
      G(KA+2)=G(KA+2)-(2.0D 0*X(KA+2)-1.0D 0)*FA
      G(KA+3)=G(KA+3)-2.0D 0*X(KA+3)*FA
      ELSEIF (KA.LT.N-2) THEN
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))+
     & 4.0D 0*(X(KA)-X(KA+1)**2)+X(KA-1)**2-X(KA-2)+
     & X(KA+1)-X(KA+2)**2+X(KA-2)**2-X(KA-3)+X(KA+2)-X(KA+3)**2
      G(KA-3)=G(KA-3)-FA
      G(KA-2)=G(KA-2)+(2.0D 0*X(KA-2)-1.0D 0)*FA
      G(KA-1)=G(KA-1)-(8.0D 0*X(KA)-2.0D 0*X(KA-1))*FA
      G(KA)=G(KA)+(24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+6.0D 0)*FA
      G(KA+1)=G(KA+1)-(8.0D 0*X(KA+1)-1.0D 0)*FA
      G(KA+2)=G(KA+2)-(2.0D 0*X(KA+2)-1.0D 0)*FA
      G(KA+3)=G(KA+3)-2.0D 0*X(KA+3)*FA
      ELSEIF (KA.LT.N-1) THEN
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))+
     & 4.0D 0*(X(KA)-X(KA+1)**2)+X(KA-1)**2-X(KA-2)+
     & X(KA+1)-X(KA+2)**2+X(KA-2)**2-X(KA-3)+X(KA+2)
      G(KA-3)=G(KA-3)-FA
      G(KA-2)=G(KA-2)+(2.0D 0*X(KA-2)-1.0D 0)*FA
      G(KA-1)=G(KA-1)-(8.0D 0*X(KA)-2.0D 0*X(KA-1))*FA
      G(KA)=G(KA)+(24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+6.0D 0)*FA
      G(KA+1)=G(KA+1)-(8.0D 0*X(KA+1)-1.0D 0)*FA
      G(KA+2)=G(KA+2)-(2.0D 0*X(KA+2)-1.0D 0)*FA
      ELSEIF (KA.LT.N) THEN
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))+
     & 4.0D 0*(X(KA)-X(KA+1)**2)+X(KA-1)**2-X(KA-2)+
     & X(KA+1)+X(KA-2)**2-X(KA-3)
      G(KA-3)=G(KA-3)-FA
      G(KA-2)=G(KA-2)+(2.0D 0*X(KA-2)-1.0D 0)*FA
      G(KA-1)=G(KA-1)-(8.0D 0*X(KA)-2.0D 0*X(KA-1))*FA
      G(KA)=G(KA)+(24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+6.0D 0)*FA
      G(KA+1)=G(KA+1)-(8.0D 0*X(KA+1)-1.0D 0)*FA
      ELSE
      FA=8.0D 0*X(KA)*(X(KA)**2-X(KA-1))-2.0D 0*(1.0D 0-X(KA))
     & +X(KA-1)**2-X(KA-2)+X(KA-2)**2-X(KA-3)
      G(KA-3)=G(KA-3)-FA
      G(KA-2)=G(KA-2)+(2.0D 0*X(KA-2)-1.0D 0)*FA
      G(KA-1)=G(KA-1)-(8.0D 0*X(KA)-2.0D 0*X(KA-1))*FA
      G(KA)=G(KA)+(24.0D 0*X(KA)**2-8.0D 0*X(KA-1)+2.0D 0)*FA
      ENDIF
      F=F+FA**2
  791 CONTINUE
      F=0.5D 0*F
      RETURN
  810 DO 811 KA=1,NA
      IF(MOD(KA,2).EQ.1) THEN
      FA=X(KA)+((5.0D 0-X(KA+1))*X(KA+1)-2.0D 0)*X(KA+1)-1.3D 1
      G(KA)=G(KA)+FA
      G(KA+1)=G(KA+1)+(10.0D 0*X(KA+1)-3.0D 0*X(KA+1)**2-2.0D 0)*FA
      ELSE
      FA=X(KA-1)+((X(KA)+1.0D 0)*X(KA)-1.4D 1)*X(KA)-2.9D 1
      G(KA-1)=G(KA-1)+FA
      G(KA)=G(KA)+(3.0D 0*X(KA)**2+2.0D 0*X(KA)-1.4D 1)*FA
      ENDIF
      F=F+FA**2
  811 CONTINUE
      F=0.5D 0*F
      RETURN
  830 DO 831 KA=1,NA
      IF(MOD(KA,4).EQ.1) THEN
      A=EXP(X(KA))-X(KA+1)
      FA=A**2
      G(KA)=G(KA)+2.0D 0*A*EXP(X(KA))*FA
      G(KA+1)=G(KA+1)-2.0D 0*A*FA
      ELSEIF (MOD(KA,4).EQ.2) THEN
      FA=1.0D 1*(X(KA)-X(KA+1))**3
      A=3.0D 1*(X(KA)-X(KA+1))**2*FA
      G(KA)=G(KA)+A
      G(KA+1)=G(KA+1)-A
      ELSEIF (MOD(KA,4).EQ.3) THEN
      A=X(KA)-X(KA+1)
      FA=(SIN(A)/COS(A))**2
      B=2.0D 0*SIN(A)/(COS(A))**3*FA
      G(KA)=G(KA)+B
      G(KA+1)=G(KA+1)-B
      ELSE
      FA=X(KA)-1.0D 0
      G(KA)=G(KA)+FA
      ENDIF
      F=F+FA**2
  831 CONTINUE
      F=0.5D 0*F
      RETURN
  840 DO 841 KA=1,NA
      IF(KA.LT.2) THEN
      FA=X(KA)*(0.5D 0*X(KA)-3.0D 0)-1.0D 0+2.0D 0*X(KA+1)
      G(KA)=G(KA)+(X(KA)-3.0D 0)*FA
      G(KA+1)=G(KA+1)+2.0D 0*FA
      ELSEIF (KA.LT.N) THEN
      FA=X(KA-1)+X(KA)*(0.5D 0*X(KA)-3.0D 0)-1.0D 0+2.0D 0*X(KA+1)
      G(KA-1)=G(KA-1)+FA
      G(KA)=G(KA)+(X(KA)-3.0D 0)*FA
      G(KA+1)=G(KA+1)+2.0D 0*FA
      ELSE
      FA=X(KA-1)+X(KA)*(0.5D 0*X(KA)-3.0D 0)-1.0D 0
      G(KA-1)=G(KA-1)+FA
      G(KA)=G(KA)+(X(KA)-3.0D 0)*FA
      ENDIF
      F=F+FA**2
  841 CONTINUE
      F=0.5D 0*F
      RETURN
  860 DO 861 KA=1,NA
      IF(MOD(KA,2).EQ.1) THEN
      FA=1.0D 4*X(KA)*X(KA+1)-1.0D 0
      G(KA)=G(KA)+1.0D 4*X(KA+1)*FA
      G(KA+1)=G(KA+1)+1.0D 4*X(KA)*FA
      ELSE
      FA=EXP(-X(KA-1))+EXP(-X(KA))-1.0001D 0
      G(KA-1)=G(KA-1)-EXP(-X(KA-1))*FA
      G(KA)=G(KA)-EXP(-X(KA))*FA
      ENDIF
      F=F+FA**2
  861 CONTINUE
      F=0.5D 0*F
      RETURN
  870 DO 871 KA=1,NA
      IF(MOD(KA,4).EQ.1) THEN
      FA=-2.0D 2*X(KA)*(X(KA+1)-X(KA)**2)-(1.0D 0-X(KA))
      G(KA)=G(KA)-(2.0D 2*(X(KA+1)-3.0D 0*X(KA)**2)-1.0D 0)*FA
      G(KA+1)=G(KA+1)-2.0D 2*X(KA)*FA
      ELSEIF(MOD(KA,4).EQ.2) THEN
      FA=2.0D 2*(X(KA)-X(KA-1)**2)+2.02D 1*(X(KA)-1.0D 0)+
     &  1.98D 1*(X(KA+2)-1.0D 0)
      G(KA-1)=G(KA-1)-4.0D 2*X(KA-1)*FA
      G(KA)=G(KA)+2.202D 2*FA
      G(KA+2)=G(KA+2)+1.98D 1*FA
      ELSEIF(MOD(KA,4).EQ.3) THEN
      FA=-1.8D 2*X(KA)*(X(KA+1)-X(KA)**2)-(1.0D 0-X(KA))
      G(KA)=G(KA)-(1.8D 2*(X(KA+1)-3.0D 0*X(KA)**2)-1.0D 0)*FA
      G(KA+1)=G(KA+1)-1.8D 2*X(KA)*FA
      ELSE
      FA=1.8D 2*(X(KA)-X(KA-1)**2)+2.02D 1*(X(KA)-1.0D 0)+
     & 1.98D 1*(X(KA-2)-1.0D 0)
      G(KA-2)=G(KA-2)+1.98D 1*FA
      G(KA-1)=G(KA-1)-3.6D 2*X(KA-1)*FA
      G(KA)=G(KA)+2.002D 2*FA
      ENDIF
      F=F+FA**2
  871 CONTINUE
      F=0.5D 0*F
      RETURN
  880 DO 881 KA=1,NA
      IF (KA.LT.2) THEN
      A=EXP(COS(DBLE(KA)*(X(KA)+X(KA+1))))
      B=A*DBLE(KA)*SIN(DBLE(KA)*(X(KA)+X(KA+1)))
      FA=X(KA)-A
      G(KA+1)=G(KA+1)+B*FA
      G(KA)=G(KA)+(B+1.0D 0)*FA
      ELSEIF (KA.LT.N) THEN
      A=EXP(COS(DBLE(KA)*(X(KA-1)+X(KA)+X(KA+1))))
      B=A*SIN(DBLE(KA)*(X(KA-1)+X(KA)+X(KA+1)))*DBLE(KA)
      FA=X(KA)-A
      G(KA-1)=G(KA-1)+B*FA
      G(KA+1)=G(KA+1)+B*FA
      G(KA)=G(KA)+(B+1.0D 0)*FA
      ELSE
      A=EXP(COS(DBLE(KA)*(X(KA-1)+X(KA))))
      B=A*SIN(DBLE(KA)*(X(KA-1)+X(KA)))*DBLE(KA)
      FA=X(KA)-A
      G(KA-1)=G(KA-1)+B*FA
      G(KA)=G(KA)+(B+1.0D 0)*FA
      ENDIF
      F=F+FA**2
  881 CONTINUE
      F=0.5D 0*F
      RETURN
  900 DO 901 KA=1,NA
      IF(KA.EQ.1) THEN
      FA=3.0D 0*X(KA)*(X(KA+1)-2.0D 0*X(KA))+0.25D 0*X(KA+1)**2
      G(KA)=G(KA)+3.0D 0*(X(KA+1)-4.0D 0*X(KA))*FA
      G(KA+1)=G(KA+1)+(3.0D 0*X(KA)+0.5D 0*X(KA+1))*FA
      ELSEIF(KA.EQ.N) THEN
      FA=3.0D 0*X(KA)*(2.0D 1-2.0D 0*X(KA)+X(KA-1))+
     & 0.25D 0*(2.0D 1-X(KA-1))**2
      G(KA-1)=G(KA-1)+(3.0D 0*X(KA)-0.5D 0*(2.0D 1-X(KA-1)))*FA
      G(KA)=G(KA)+3.0D 0*(2.0D 1-4.0D 0*X(KA)+X(KA-1))*FA
      ELSE
      FA=3.0D 0*X(KA)*(X(KA+1)-2.0D 0*X(KA)+X(KA-1))+
     & 0.25D 0*(X(KA+1)-X(KA-1))**2
      G(KA-1)=G(KA-1)+(3.0D 0*X(KA)-0.5D 0*(X(KA+1)-X(KA-1)))*FA
      G(KA)=G(KA)+3.0D 0*(X(KA+1)-4.0D 0*X(KA)+X(KA-1))*FA
      G(KA+1)=G(KA+1)+(3.0D 0*X(KA)+0.5D 0*(X(KA+1)-X(KA-1)))*FA
      ENDIF
      F=F+FA**2
  901 CONTINUE
      F=0.5D 0*F
      RETURN
  910 DO 911 KA=1,NA
      H=1.0D 0/DBLE(N+1)
      IF (KA.LT.2) THEN
      FA=2.0D 0*X(KA)+PAR*H**2*SINH(PAR*X(KA))-X(KA+1)
      G(KA)=G(KA)+(2.0D 0+PAR**2*H**2*COSH(PAR*X(KA)))*FA
      G(KA+1)=G(KA+1)-FA
      ELSE IF (KA.LT.N) THEN
      FA=2.0D 0*X(KA)+PAR*H**2*SINH(PAR*X(KA))-X(KA-1)-X(KA+1)
      G(KA-1)=G(KA-1)-FA
      G(KA)=G(KA)+(2.0D 0+PAR**2*H**2*COSH(PAR*X(KA)))*FA
      G(KA+1)=G(KA+1)-FA
      ELSE
      FA=2.0D 0*X(KA)+PAR*H**2*SINH(PAR*X(KA))-X(KA-1)-1.0D 0
      G(KA)=G(KA)+(2.0D 0+PAR**2*H**2*COSH(PAR*X(KA)))*FA
      G(KA-1)=G(KA-1)-FA
      ENDIF
      F=F+FA**2
  911 CONTINUE
      F=0.5D 0*F
      RETURN
  920 DO 921 KA=1,NA
      FA=6.0D 0*X(KA)
      A1=0.0D 0
      A2=0.0D 0
      A3=0.0D 0
      IF (KA.GT.1) THEN
      FA=FA-4.0D 0*X(KA-1)
      A1=A1-X(KA-1)
      A2=A2+X(KA-1)
      A3=A3+2.0D 0*X(KA-1)
      ENDIF
      IF (KA.GT.2) THEN
      FA=FA+X(KA-2)
      A3=A3-X(KA-2)
      ENDIF
      IF (KA.LT.N-1) THEN
      FA=FA+X(KA+2)
      A3=A3+X(KA+2)
      ENDIF
      IF (KA.LT.N) THEN
      FA=FA-4.0D 0*X(KA+1)
      A1=A1+X(KA+1)
      A2=A2+X(KA+1)
      A3=A3-2.0D 0*X(KA+1)
      ENDIF
      IF (KA.GE.N-1) THEN
      FA=FA+1.0D 0
      A3=A3+1.0D 0
      ENDIF
      IF (KA.GE.N) THEN
      FA=FA-4.0D 0
      A1=A1+1.0D 0
      A2=A2+1.0D 0
      A3=A3-2.0D 0
      ENDIF
      FA=FA-0.5D 0*PAR*(A1*A2-X(KA)*A3)
      F=F+FA**2
      G(KA)=G(KA)+6.0D 0*FA
      GA1(1)=0.0D 0
      GA1(2)=0.0D 0
      GA2(1)=0.0D 0
      GA2(2)=0.0D 0
      IF (KA.GT.1) THEN
      G(KA-1)=G(KA-1)-(4.0D 0-PAR*X(KA))*FA
      GA1(1)=-1.0D 0
      GA2(1)= 1.0D 0
      ENDIF
      IF (KA.GT.2) THEN
      G(KA-2)=G(KA-2)+(1.0D 0-0.5D 0*PAR*X(KA))*FA
      ENDIF
      IF (KA.LT.N-1) THEN
      G(KA+2)=G(KA+2)+(1.0D 0+0.5D 0*PAR*X(KA))*FA
      ENDIF
      IF (KA.LT.N) THEN
      G(KA+1)=G(KA+1)-(4.0D 0+PAR*X(KA))*FA
      GA1(2)= 1.0D 0
      GA2(2)= 1.0D 0
      ENDIF
      G(KA)=G(KA)+0.5D 0*PAR*A3*FA
      IF (KA.GT.1)
     & G(KA-1)=G(KA-1)-0.5D 0*PAR*(GA1(1)*A2+A1*GA2(1))*FA
      IF (KA.LT.N)
     & G(KA+1)=G(KA+1)-0.5D 0*PAR*(GA1(2)*A2+A1*GA2(2))*FA
  921 CONTINUE
      F=0.5D 0*F
      RETURN
  930 DO 931 KA=1,NA
      H=1.0D 0/DBLE(M+1)
      IF(KA.LE.M) THEN
      J=KA+M
      FA=6.0D 0*X(KA)
      A1=0.0D 0
      A2=0.0D 0
      IF (KA.EQ.1) THEN
      A1=A1+1.0D 0
      ENDIF
      IF (KA.GT.1) THEN
      FA=FA-4.0D 0*X(KA-1)
      A1=A1-X(J-1)
      A2=A2+2.0D 0*X(KA-1)
      ENDIF
      IF (KA.GT.2) THEN
      FA=FA+X(KA-2)
      A2=A2-X(KA-2)
      ENDIF
      IF (KA.LT.M-1) THEN
      FA=FA+X(KA+2)
      A2=A2+X(KA+2)
      ENDIF
      IF (KA.LT.M) THEN
      FA=FA-4.0D 0*X(KA+1)
      A1=A1+X(J+1)
      A2=A2-2.0D 0*X(KA+1)
      ENDIF
      IF (KA.EQ.M) THEN
      A1=A1+1.0D 0
      ENDIF
      FA=FA+0.5D 0*PAR*H*(X(KA)*A2+X(J)*A1*H**2)
      ELSE
      J=KA-M
      FA=-2.0D 0*X(KA)
      A1=0.0D 0
      A2=0.0D 0
      IF (J.EQ.1) THEN
      A2=A2+1.0D 0
      ENDIF
      IF (J.GT.1) THEN
      FA=FA+X(KA-1)
      A1=A1-X(J-1)
      A2=A2-X(KA-1)
      ENDIF
      IF (J.LT.M) THEN
      FA=FA+X(KA+1)
      A1=A1+X(J+1)
      A2=A2+X(KA+1)
      ENDIF
      IF (J.EQ.M) THEN
      A2=A2+1.0D 0
      ENDIF
      FA=FA+0.5D 0*PAR*H*(X(KA)*A1+X(J)*A2)
      ENDIF
      F=F+FA**2
      IF(KA.LE.M) THEN
      G(KA)=G(KA)+6.0D 0*FA
      IF (KA.GT.1) THEN
      G(KA-1)=G(KA-1)-(4.0D 0-PAR*H*X(KA))*FA
      G(J-1)=G(J-1)-0.5D 0*PAR*H**3*X(J)*FA
      ENDIF
      IF (KA.GT.2) THEN
      G(KA-2)=G(KA-2)+(1.0D 0-0.5D 0*PAR*H*X(KA))*FA
      ENDIF
      IF (KA.LT.M-1) THEN
      G(KA+2)=G(KA+2)+(1.0D 0+0.5D 0*PAR*H*X(KA))*FA
      ENDIF
      IF (KA.LT.M) THEN
      G(KA+1)=G(KA+1)-(4.0D 0+PAR*H*X(KA))*FA
      G(J+1)=G(J+1)+0.5D 0*PAR*H**3*X(J)*FA
      ENDIF
      G(KA)=G(KA)+0.5D 0*PAR*H*A2*FA
      G(J)=G(J)+0.5D 0*PAR*H**3*A1*FA
      ELSE
      G(KA)=G(KA)-2.0D 0*FA
      IF (J.GT.1) THEN
      G(KA-1)=G(KA-1)+(1.0D 0-0.5D 0*PAR*H*X(J))*FA
      G(J-1)=G(J-1)-0.5D 0*PAR*H*X(KA)*FA
      ENDIF
      IF (J.LT.M) THEN
      G(KA+1)=G(KA+1)+(1.0D 0+0.5D 0*PAR*H*X(J))*FA
      G(J+1)=G(J+1)+0.5D 0*PAR*H*X(KA)*FA
      ENDIF
      G(KA)=G(KA)+0.5D 0*PAR*H*A1*FA
      G(J)=G(J)+0.5D 0*PAR*H*A2*FA
      ENDIF
  931 CONTINUE
      F=0.5D 0*F
      RETURN
  940 DO 941 KA=1,NA
      FA=4.0D 0*X(KA)-PAR*EXP(X(KA))
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      IF(I.GT.1) FA=FA-X(KA-1)
      IF(I.LT.M) FA=FA-X(KA+1)
      IF(J.GT.1) FA=FA-X(KA-M)
      IF(J.LT.M) FA=FA-X(KA+M)
      F=F+FA**2
      G(KA)=G(KA)+(4.0D 0-PAR*EXP(X(KA)))*FA
      IF(J.GT.1) G(KA-M)=G(KA-M)-FA
      IF(I.GT.1) G(KA-1)=G(KA-1)-FA
      IF(I.LT.M) G(KA+1)=G(KA+1)-FA
      IF(J.LT.M) G(KA+M)=G(KA+M)-FA
  941 CONTINUE
      F=0.5D 0*F
      RETURN
  950 DO 951 KA=1,NA
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      FA=4.0D 0*X(KA)+PAR*X(KA)**3/(1.0D 0+PAR*DBLE(I)**2+
     & PAR*DBLE(J)**2)
      IF(I.EQ.1) FA=FA-1.0D 0
      IF(I.GT.1) FA=FA-X(KA-1)
      IF(I.LT.M) FA=FA-X(KA+1)
      IF(I.EQ.M) FA=FA-2.0D 0+EXP(DBLE(J)/DBLE(M+1))
      IF(J.EQ.1) FA=FA-1.0D 0
      IF(J.GT.1) FA=FA-X(KA-M)
      IF(J.LT.M) FA=FA-X(KA+M)
      IF(J.EQ.M) FA=FA-2.0D 0+EXP(DBLE(I)/DBLE(M+1))
      F=F+FA**2
      G(KA)=G(KA)+(4.0D 0+3.0D 0*PAR*X(KA)**2/(1.0D 0+PAR*DBLE(I)**2+
     & PAR*DBLE(J)**2))*FA
      IF(J.GT.1) G(KA-M)=G(KA-M)-FA
      IF(I.GT.1) G(KA-1)=G(KA-1)-FA
      IF(I.LT.M) G(KA+1)=G(KA+1)-FA
      IF(J.LT.M) G(KA+M)=G(KA+M)-FA
  951 CONTINUE
      F=0.5D 0*F
      RETURN
  960 DO 961 KA=1,NA
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      A1=DBLE(I)/DBLE(M+1)
      A2=DBLE(J)/DBLE(M+1)
      FA=4.0D 0*X(KA)-PAR*SIN(2.0D 0*PI*X(KA))-
     & 1.0D 4*((A1-0.25D 0)**2+(A2-0.75D 0)**2)*PAR
      IF(I.EQ.1) FA=FA-X(KA+1)-PAR*SIN(PI*X(KA+1)*DBLE(M+1))
      IF(I.GT.1.AND.I.LT.M) FA=FA-X(KA+1)-X(KA-1)-
     & PAR*SIN(PI*(X(KA+1)-X(KA-1))*DBLE(M+1))
      IF(I.EQ.M) FA=FA-X(KA-1)+PAR*SIN(PI*X(KA-1)*DBLE(M+1))
      IF(J.EQ.1) FA=FA-X(KA+M)-PAR*SIN(PI*X(KA+M)*DBLE(M+1))
      IF(J.GT.1.AND.J.LT.M) FA=FA-X(KA+M)-X(KA-M)-
     & PAR*SIN(PI*(X(KA+M)-X(KA-M))*DBLE(M+1))
      IF(J.EQ.M) FA=FA-X(KA-M)+PAR*SIN(PI*X(KA-M)*DBLE(M+1))
      F=F+FA**2
      G(KA)=G(KA)+(4.0D 0-2.0D 0*PI*PAR*COS(2.0D 0*PI*X(KA)))*FA
      IF(I.EQ.1) G(KA+1)=G(KA+1)-
     & (1.0D 0+PI*DBLE(M+1)*PAR*COS(PI*X(KA+1)*DBLE(M+1)))*FA
      IF(I.GT.1.AND.I.LT.M) THEN
      G(KA-1)=G(KA-1)-
     & (1.0D 0-PI*DBLE(M+1)*PAR*COS(PI*(X(KA+1)-X(KA-1))*DBLE(M+1)))*FA
      G(KA+1)=G(KA+1)-
     & (1.0D 0+PI*DBLE(M+1)*PAR*COS(PI*(X(KA+1)-X(KA-1))*DBLE(M+1)))*FA
      ENDIF
      IF(I.EQ.M) G(KA-1)=G(KA-1)-
     & (1.0D 0-PI*DBLE(M+1)*PAR*COS(PI*X(KA-1)*DBLE(M+1)))*FA
      IF(J.EQ.1) G(KA+M)=G(KA+M)-
     & (1.0D 0+PI*DBLE(M+1)*PAR*COS(PI*X(KA+M)*DBLE(M+1)))*FA
      IF(J.GT.1.AND.J.LT.M) THEN
      G(KA-M)=G(KA-M)-
     & (1.0D 0-PI*DBLE(M+1)*PAR*COS(PI*(X(KA+M)-X(KA-M))*DBLE(M+1)))*FA
      G(KA+M)=G(KA+M)-
     & (1.0D 0+PI*DBLE(M+1)*PAR*COS(PI*(X(KA+M)-X(KA-M))*DBLE(M+1)))*FA
      ENDIF
      IF(J.EQ.M) G(KA-M)=G(KA-M)-
     & (1.0D 0-PI*DBLE(M+1)*PAR*COS(PI*X(KA-M)*DBLE(M+1)))*FA
  961 CONTINUE
      F=0.5D 0*F
      RETURN
  970 DO 971 KA=1,NA
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      FA=8.0D 0*X(KA)**2
      IF(I.EQ.1) FA=FA-2.0D 0*X(KA)*(X(KA+1)+1.0D 0)-
     & 0.5D 0*(X(KA+1)-1.0D 0)**2-
     & 1.5D 0*X(KA)**2*(X(KA+1)-1.0D 0)*PAR
      IF(I.GT.1.AND.I.LT.M) FA=FA-2.0D 0*X(KA)*(X(KA+1)+X(KA-1))-
     & 0.5D 0*(X(KA+1)-X(KA-1))**2-
     & 1.5D 0*X(KA)**2*(X(KA+1)-X(KA-1))*PAR
      IF(I.EQ.M) FA=FA-2.0D 0*X(KA)*X(KA-1)-
     & 0.5D 0*X(KA-1)**2+
     & 1.5D 0*X(KA)**2*X(KA-1)*PAR
      IF(J.EQ.1) FA=FA-2.0D 0*X(KA)*(X(KA+M)+1.0D 0)-
     & 0.5D 0*(X(KA+M)-1.0D 0)**2
      IF(J.GT.1.AND.J.LT.M) FA=FA-2.0D 0*X(KA)*(X(KA+M)+X(KA-M))-
     & 0.5D 0*(X(KA+M)-X(KA-M))**2
      IF(J.EQ.M) FA=FA-2.0D 0*X(KA)*X(KA-M)-
     & 0.5D 0*X(KA-M)**2
      IF (I.EQ.1.AND.J.EQ.1) FA=FA-PAR/DBLE(M+1)
      F=F+FA**2
      G(KA)=G(KA)+1.6D 1*X(KA)*FA
      IF(I.EQ.1) THEN
      G(KA)=G(KA)-(2.0D 0*(X(KA+1)+1.0D 0)+3.0D 0*X(KA)*
     & (X(KA+1)-1.0D 0)*PAR)*FA
      G(KA+1)=G(KA+1)-(2.0D 0*X(KA)+(X(KA+1)-1.0D 0)+
     & 1.5D 0*X(KA)**2*PAR)*FA
      ENDIF
      IF(I.GT.1.AND.I.LT.M) THEN
      G(KA)=G(KA)-(2.0D 0*(X(KA+1)+X(KA-1))+3.0D 0*X(KA)*
     & (X(KA+1)-X(KA-1))*PAR)*FA
      G(KA-1)=G(KA-1)-(2.0D 0*X(KA)-(X(KA+1)-X(KA-1))-
     & 1.5D 0*X(KA)**2*PAR)*FA
      G(KA+1)=G(KA+1)-(2.0D 0*X(KA)+(X(KA+1)-X(KA-1))+
     & 1.5D 0*X(KA)**2*PAR)*FA
      ENDIF
      IF(I.EQ.M) THEN
      G(KA)=G(KA)-(2.0D 0*X(KA-1)-3.0D 0*X(KA)*X(KA-1)*PAR)*FA
      G(KA-1)=G(KA-1)-(2.0D 0*X(KA)+X(KA-1)-1.5D 0*X(KA)**2*PAR)*FA
      ENDIF
      IF(J.EQ.1) THEN
      G(KA)=G(KA)-2.0D 0*(X(KA+M)+1.0D 0)*FA
      G(KA+M)=G(KA+M)-(2.0D 0*X(KA)+(X(KA+M)-1.0D 0))*FA
      ENDIF
      IF(J.GT.1.AND.J.LT.M) THEN
      G(KA)=G(KA)-2.0D 0*(X(KA+M)+X(KA-M))*FA
      G(KA-M)=G(KA-M)-(2.0D 0*X(KA)-(X(KA+M)-X(KA-M)))*FA
      G(KA+M)=G(KA+M)-(2.0D 0*X(KA)+(X(KA+M)-X(KA-M)))*FA
      ENDIF
      IF(J.EQ.M) THEN
      G(KA)=G(KA)-2.0D 0*X(KA-M)*FA
      G(KA-M)=G(KA-M)-(2.0D 0*X(KA)+X(KA-M))*FA
      ENDIF
  971 CONTINUE
      F=0.5D 0*F
      RETURN
  980 DO 981 KA=1,NA
      A3=0.0D 0
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      A1=PAR*DBLE(I)
      A2=PAR*DBLE(J)
      FA=4.0D 0*X(KA)-2.0D 3*A1*A2*(1.0D 0-A1)*(1.0D 0-A2)*PAR**2
      IF(I.GT.1) THEN
      FA=FA-X(KA-1)
      A3=A3-X(KA-1)
      ENDIF
      IF(I.LT.M) THEN
      FA=FA-X(KA+1)
      A3=A3+X(KA+1)
      ENDIF
      IF(J.GT.1) THEN
      FA=FA-X(KA-M)
      A3=A3-X(KA-M)
      ENDIF
      IF(J.LT.M) THEN
      FA=FA-X(KA+M)
      A3=A3+X(KA+M)
      ENDIF
      FA=FA+2.0D 1*PAR*A3*X(KA)
      F=F+FA**2
      G(KA)=G(KA)+4.0D 0*FA
      IF(I.GT.1) THEN
      G(KA-1)=G(KA-1)-(1.0D 0+2.0D 1*PAR*X(KA))*FA
      ENDIF
      IF(I.LT.M) THEN
      G(KA+1)=G(KA+1)-(1.0D 0-2.0D 1*PAR*X(KA))*FA
      ENDIF
      IF(J.GT.1) THEN
      G(KA-M)=G(KA-M)-(1.0D 0+2.0D 1*PAR*X(KA))*FA
      ENDIF
      IF(J.LT.M) THEN
      G(KA+M)=G(KA+M)-(1.0D 0-2.0D 1*PAR*X(KA))*FA
      ENDIF
      G(KA)=G(KA)+2.0D 1*PAR*A3*FA
  981 CONTINUE
      F=0.5D 0*F
      RETURN
  990 DO 991 KA=1,NA
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      FA=2.0D 1*X(KA)-PAR*MAX(0.0D 0,X(KA))-
     & SIGN(PAR,(DBLE(I)/DBLE(M+2)-0.5D 0))
      IF (J.GT.2) THEN
        FA=FA+X(KA-M-M)
      ENDIF
      IF (J.GT.1) THEN
        IF (I.GT.1) THEN
          FA=FA+2.0D 0*X(KA-M-1)
        ENDIF
        FA=FA-8.0D 0*X(KA-M)
        IF (I.LT.M) THEN
          FA=FA+2.0D 0*X(KA-M+1)
        ENDIF
      ENDIF
      IF (I.GT.1) THEN
        IF (I.GT.2) THEN
          FA=FA+X(KA-2)
        ENDIF
        FA=FA-8.0D 0*X(KA-1)
      ENDIF
      IF (I.LT.M) THEN
        FA=FA-8.0D 0*X(KA+1)
        IF (I.LT.M-1) THEN
          FA=FA+X(KA+2)
        ENDIF
      ENDIF
      IF (J.LT.M) THEN
        IF (I.GT.1) THEN
          FA=FA+2.0D 0*X(KA+M-1)
        ENDIF
        FA=FA-8.0D 0*X(KA+M)
        IF (I.LT.M) THEN
          FA=FA+2.0D 0*X(KA+M+1)
        ENDIF
      ENDIF
      IF (J.LT.M-1) THEN
        FA=FA+X(KA+M+M)
      ENDIF
      F=F+FA**2
      G(KA)=G(KA)+(2.0D 1-PAR)*FA
      IF (J.GT.2) THEN
        G(KA-M-M)=G(KA-M-M)+FA
      ENDIF
      IF (J.GT.1) THEN
        IF (I.GT.1) THEN
          G(KA-M-1)=G(KA-M-1)+2.0D 0*FA
        ENDIF
        G(KA-M)=G(KA-M)-8.0D 0*FA
        IF (I.LT.M) THEN
          G(KA-M+1)=G(KA-M+1)+2.0D 0*FA
        ENDIF
      ENDIF
      IF (I.GT.1) THEN
        IF (I.GT.2) THEN
          G(KA-2)=G(KA-2)+FA
        ENDIF
        G(KA-1)=G(KA-1)-8.0D 0*FA
      ENDIF
      IF (I.LT.M) THEN
        G(KA+1)=G(KA+1)-8.0D 0*FA
        IF (I.LT.M-1) THEN
          G(KA+2)=G(KA+2)+FA
        ENDIF
      ENDIF
      IF (J.LT.M) THEN
        IF (I.GT.1) THEN
          G(KA+M-1)=G(KA+M-1)+2.0D 0*FA
        ENDIF
        G(KA+M)=G(KA+M)-8.0D 0*FA
        IF (I.LT.M) THEN
          G(KA+M+1)=G(KA+M+1)+2.0D 0*FA
        ENDIF
      ENDIF
      IF (J.LT.M-1) THEN
        G(KA+M+M)=G(KA+M+M)+FA
      ENDIF
  991 CONTINUE
      F=0.5D 0*F
      RETURN
  800 DO 802 KA=1,NA
      H=0.5D 0/DBLE(M+2)
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      FA=2.0D 1*X(KA)
      A1=0.0D 0
      A2=0.0D 0
      A3=0.0D 0
      A4=0.0D 0
      IF (J.GT.2) THEN
        FA=FA+X(KA-M-M)
        A4=A4+X(KA-M-M)
      ENDIF
      IF (J.GT.1) THEN
        IF (I.GT.1) THEN
          FA=FA+2.0D 0*X(KA-M-1)
          A3=A3+X(KA-M-1)
          A4=A4+X(KA-M-1)
        ENDIF
        FA=FA-8.0D 0*X(KA-M)
        A1=A1-X(KA-M)
        A4=A4-4.0D 0*X(KA-M)
        IF (I.LT.M) THEN
          FA=FA+2.0D 0*X(KA-M+1)
          A3=A3-X(KA-M+1)
          A4=A4+X(KA-M+1)
        ENDIF
      ENDIF
      IF (I.GT.1) THEN
        IF (I.GT.2) THEN
          FA=FA+X(KA-2)
          A3=A3+X(KA-2)
        ENDIF
        FA=FA-8.0D 0*X(KA-1)
        A2=A2-X(KA-1)
        A3=A3-4.0D 0*X(KA-1)
      ENDIF
      IF (I.LT.M) THEN
        FA=FA-8.0D 0*X(KA+1)
        A2=A2+X(KA+1)
        A3=A3+4.0D 0*X(KA+1)
        IF (I.LT.M-1) THEN
          FA=FA+X(KA+2)
          A3=A3-X(KA+2)
        ENDIF
      ENDIF
      IF (J.LT.M) THEN
        IF (I.GT.1) THEN
          FA=FA+2.0D 0*X(KA+M-1)
          A3=A3+X(KA+M-1)
          A4=A4-X(KA+M-1)
        ENDIF
        FA=FA-8.0D 0*X(KA+M)
        A1=A1+X(KA+M)
        A4=A4+4.0D 0*X(KA+M)
        IF (I.LT.M) THEN
          FA=FA+2.0D 0*X(KA+M+1)
          A3=A3-X(KA+M+1)
          A4=A4-X(KA+M+1)
        ENDIF
      ENDIF
      IF (J.LT.M-1) THEN
        FA=FA+X(KA+M+M)
        A4=A4-X(KA+M+M)
      ENDIF
      IF (J.EQ.M) THEN
        IF (I.GT.1) THEN
          FA=FA-H-H
          A3=A3-H
          A4=A4+H
        ENDIF
        FA=FA+8.0D 0*H
        A1=A1-H
        A4=A4-4.0D 0*H
        IF (I.LT.M) THEN
          FA=FA-2.0D 0*H
          A3=A3+H
          A4=A4+H
        ENDIF
        FA=FA+H
        A4=A4-H
      ENDIF
      IF (J.EQ.M-1) THEN
        FA=FA-H
        A4=A4+H
      ENDIF
      FA=FA+0.25D 0*PAR*(A1*A3-A2*A4)
      F=F+FA**2
      G(KA)=G(KA)+2.0D 1*FA
      A1=0.0D 0
      A2=0.0D 0
      A3=0.0D 0
      A4=0.0D 0
      GA1(1)=0.0D 0
      GA1(2)=0.0D 0
      GA2(1)=0.0D 0
      GA2(2)=0.0D 0
      DO 801 K=1,6
      GA3(K)=0.0D 0
      GA4(K)=0.0D 0
  801 CONTINUE
      IF (J.GT.2) THEN
        G(KA-M-M)=G(KA-M-M)+FA
        GA4(1)=GA4(1)+1.0D 0
        A4=A4+X(KA-M-M)
      ENDIF
      IF (J.GT.1) THEN
        IF (I.GT.1) THEN
          G(KA-M-1)=G(KA-M-1)+2.0D 0*FA
          GA3(1)=GA3(1)+1.0D 0
          GA4(2)=GA4(2)+1.0D 0
          A3=A3+X(KA-M-1)
          A4=A4+X(KA-M-1)
        ENDIF
        G(KA-M)=G(KA-M)-8.0D 0*FA
        GA1(1)=GA1(1)-1.0D 0
        A1=A1-X(KA-M)
        IF (I.LT.M) THEN
          G(KA-M+1)=G(KA-M+1)+2.0D 0*FA
          GA3(2)=GA3(2)-1.0D 0
          GA4(3)=GA4(3)+1.0D 0
          A3=A3-X(KA-M+1)
          A4=A4+X(KA-M+1)
        ENDIF
      ENDIF
      IF (I.GT.1) THEN
        IF (I.GT.2) THEN
          G(KA-2)=G(KA-2)+FA
          GA3(3)=GA3(3)+1.0D 0
          A3=A3+X(KA-2)
        ENDIF
        G(KA-1)=G(KA-1)-8.0D 0*FA
        GA2(1)=GA2(1)-1.0D 0
        A2=A2-X(KA-1)
      ENDIF
      IF (I.LT.M) THEN
        G(KA+1)=G(KA+1)-8.0D 0*FA
        GA2(2)=GA2(2)+1.0D 0
        A2=A2+X(KA+1)
        IF (I.LT.M-1) THEN
          G(KA+2)=G(KA+2)+FA
          GA3(4)=GA3(4)-1.0D 0
          A3=A3-X(KA+2)
        ENDIF
      ENDIF
      IF (J.LT.M) THEN
        IF (I.GT.1) THEN
          G(KA+M-1)=G(KA+M-1)+2.0D 0*FA
          GA3(5)=GA3(5)+1.0D 0
          GA4(4)=GA4(4)-1.0D 0
          A3=A3+X(KA+M-1)
          A4=A4-X(KA+M-1)
        ENDIF
        G(KA+M)=G(KA+M)-8.0D 0*FA
        GA1(2)=GA1(2)+1.0D 0
        A1=A1+X(KA+M)
        IF (I.LT.M) THEN
          G(KA+M+1)=G(KA+M+1)+2.0D 0*FA
          GA3(6)=GA3(6)-1.0D 0
          GA4(5)=GA4(5)-1.0D 0
          A3=A3-X(KA+M+1)
          A4=A4-X(KA+M+1)
        ENDIF
      ENDIF
      IF (J.LT.M-1) THEN
        G(KA+M+M)=G(KA+M+M)+FA
        GA4(6)=GA4(6)-1.0D 0
        A4=A4-X(KA+M+M)
      ENDIF
      IF (J.EQ.M) THEN
        IF (I.GT.1) THEN
          A3=A3-H
          A4=A4+H
        ENDIF
        A1=A1-H
        IF (I.LT.M) THEN
          A3=A3+H
          A4=A4+H
        ENDIF
        A4=A4-H
      ENDIF
      IF (J.EQ.M-1) THEN
        A4=A4+H
      ENDIF
      IF (KA.GT.M+M)
     & G(KA-M-M)=G(KA-M-M)+0.25D 0*PAR*(-A2*GA4(1))*FA
      IF (KA.GT.M+1)
     & G(KA-M-1)=G(KA-M-1)+0.25D 0*PAR*(+A1*GA3(1)-A2*GA4(2))*FA
      IF (KA.GT.M)
     & G(KA-M)=G(KA-M)+0.25D 0*PAR*(GA1(1)*A3)*FA
      IF (KA.GT.M-1)
     & G(KA-M+1)=G(KA-M+1)+0.25D 0*PAR*(+A1*GA3(2)-A2*GA4(3))*FA
      IF (KA.GT.2)
     & G(KA-2)=G(KA-2)+0.25D 0*PAR*(+A1*GA3(3))*FA
      IF (KA.GT.1)
     & G(KA-1)=G(KA-1)+0.25D 0*PAR*(-GA2(1)*A4)*FA
      IF (KA.LE.N-1)
     & G(KA+1)=G(KA+1)+0.25D 0*PAR*(-GA2(2)*A4)*FA
      IF (KA.LE.N-2)
     & G(KA+2)=G(KA+2)+0.25D 0*PAR*(+A1*GA3(4))*FA
      IF (KA.LE.N-M+1)
     & G(KA+M-1)=G(KA+M-1)+0.25D 0*PAR*(+A1*GA3(5)-A2*GA4(4))*FA
      IF (KA.LE.N-M)
     & G(KA+M)=G(KA+M)+0.25D 0*PAR*(GA1(2)*A3)*FA
      IF (KA.LE.N-M-1)
     & G(KA+M+1)=G(KA+M+1)+0.25D 0*PAR*(+A1*GA3(6)-A2*GA4(5))*FA
      IF (KA.LE.N-M-M)
     & G(KA+M+M)=G(KA+M+M)+0.25D 0*PAR*(-A2*GA4(6))*FA
  802 CONTINUE
      F=0.5D 0*F
      RETURN
  240 DO 243 KA=1,NA
      W=0.0D 0
      DO 241 I=1,N-1
      W=W+(DBLE(KA)/DBLE(KA+I))*X(I)
  241 CONTINUE
      FA=X(KA)-(1.0D 0+(0.4D 0/DBLE(N))*X(KA)*(0.5D 0+W+
     & 0.5D 0*(DBLE(KA)/DBLE(KA+N))*X(N)))
      F=F+FA**2
      W=W+0.5D 0+0.5D 0*DBLE(KA)/DBLE(KA+N)*X(N)
      DO 242 I=1,N-1
      G(I)=G(I)-0.4D 0/DBLE(N)*X(KA)*DBLE(KA)/DBLE(KA+I)*FA
  242 CONTINUE
      G(N)=G(N)-0.2D 0/DBLE(N)*X(KA)*DBLE(KA)/DBLE(KA+N)*FA
      G(KA)=G(KA)+(1.0D 0-0.4D 0*W/DBLE(N))*FA
  243 CONTINUE
      F=0.5D 0*F
      RETURN
  410 DO 411 KA=1,NA
      IF(KA.EQ.1) THEN
        FA=1.0D 0-X(1)
      G(1)=G(1)-FA
      ELSE
        FA=10.0D 0*DBLE(KA-1)*(X(KA)-X(KA-1))**2
        G(KA)=G(KA)+20.0D 0*DBLE(KA-1)*(X(KA)-X(KA-1))*FA
        G(KA-1)=G(KA-1)-20.0D 0*DBLE(KA-1)*(X(KA)-X(KA-1))*FA
      ENDIF
      F=F+FA**2
  411 CONTINUE
      F=0.5D 0*F
      RETURN
  420 DO 421 KA=1,NA
      IF(KA.EQ.N) THEN
        FA=X(KA)-0.1D 0*X(1)**2
        G(1)=G(1)-0.20D 0*X(1)*FA
        G(N)=G(N)+FA
      ELSE
        FA=X(KA)-0.1D 0*X(KA+1)**2
        G(KA)=G(KA)+FA
        G(KA+1)=G(KA+1)-0.2D 0*X(KA+1)*FA
      ENDIF
      F=F+FA**2
  421 CONTINUE
      F=0.5D 0*F
      RETURN
  650 DO 653 KA=1,NA
      S=0.0D 0
      DO 651 J=1,N
        S=S+X(J)**3
  651 CONTINUE
      FA=X(KA)-1.0D 0/DBLE(2*N)*(S+DBLE(KA))
      F=F+FA**2
      DO 652 J=1,N
        IF(J.EQ.KA) THEN
           G(J)=G(J)+(1.0D 0-3.0D 0*X(J)**2/(2.0D 0*DBLE(N)))*FA
        ELSE
           G(J)=G(J)-3.0D 0*X(J)**2/(2.0D 0*DBLE(N))*FA
        ENDIF
  652 CONTINUE
  653 CONTINUE
      F=0.5D 0*F
      RETURN
  660 DO 661 KA=1,NA
      S=(1.0D 0/DBLE(N+1))**2*EXP(X(KA))
      IF(N.EQ.1) THEN
        FA=-2.0D 0*X(KA)-S
        G(KA)=G(KA)-(2.0D 0+S)*FA
      ELSE IF(KA.EQ.1) THEN
        FA=-2.0D 0*X(KA)+X(KA+1)-S
        G(KA)=G(KA)-(2.0D 0+S)*FA
        G(KA+1)=G(KA+1)+FA
      ELSE IF(KA.EQ.N) THEN
        FA=X(KA-1)-2.0D 0*X(KA)-S
        G(KA)=G(KA)-(2.0D 0+S)*FA
        G(KA-1)=G(KA-1)+FA
      ELSE
        FA=X(KA-1)-2.0D 0*X(KA)+X(KA+1)-S
        G(KA)=G(KA)-(2.0D 0+S)*FA
        G(KA-1)=G(KA-1)+FA
        G(KA+1)=G(KA+1)+FA
      ENDIF
      F=F+FA**2
  661 CONTINUE
      F=0.5D 0*F
      RETURN
  670 DO 671 KA=1,NA
      S=0.1D 0
      IF(N.EQ.1) THEN
        FA=(3.0D 0-S*X(KA))*X(KA)+1.0D 0
        G(KA)=G(KA)+(3.0D 0-2.0D 0*S*X(KA))*FA
      ELSE IF(KA.EQ.1) THEN
        FA=(3.0D 0-S*X(KA))*X(KA)+1.0D 0-2.0D 0*X(KA+1)
        G(KA)=G(KA)+(3.0D 0-2.0D 0*S*X(KA))*FA
        G(KA+1)=G(KA+1)-2.0D 0*FA
      ELSE IF(KA.EQ.N) THEN
        FA=(3.0D 0-S*X(KA))*X(KA)+1.0D 0-X(KA-1)
        G(KA)=G(KA)+(3.0D 0-2.0D 0*S*X(KA))*FA
        G(KA-1)=G(KA-1)-FA
      ELSE
        FA=(3.0D 0-S*X(KA))*X(KA)+1.0D 0-X(KA-1)-2.0D 0*X(KA+1)
        G(KA)=G(KA)+(3.0D 0-2.0D 0*S*X(KA))*FA
        G(KA-1)=G(KA-1)-FA
        G(KA+1)=G(KA+1)-2.0D 0*FA
      ENDIF
      F=F+FA**2
  671 CONTINUE
      F=0.5D 0*F
      RETURN
  680 DO 683 KA=1,NA
      S1=1.0D 0
      S2=1.0D 0
      S3=1.0D 0
      J1=3
      J2=3
      IF(KA-J1.GT.1) THEN
        I1=KA-J1
      ELSE
        I1=1
      ENDIF
      IF(KA+J2.LT.N) THEN
        I2=KA+J2
      ELSE
        I2=N
      ENDIF
      S=0.0D 0
      DO 681 J=I1,I2
        IF(J.NE.KA) S=S+X(J)+X(J)**2
  681 CONTINUE
      FA=(S1+S2*X(KA)**2)*X(KA)+1.D 0-S3*S
      G(KA)=G(KA)+(S1+3.0D 0*S2*X(KA)**2)*FA
      DO 682 J=I1,I2
        G(J)=G(J)-S3*(1.0D 0+2.0D 0*X(J))*FA
  682 CONTINUE
      F=F+FA**2
  683 CONTINUE
      F=0.5D 0*F
      RETURN
  690 DO 691 KA=1,NA
      IF(KA.EQ.1) THEN
        FA=X(1)**2-1.0D 0
        G(1)=G(1)+2.0D 0*X(1)*FA
      ELSE
        FA=X(KA-1)**2+LOG(X(KA))-1.0D 0
        G(KA-1)=G(KA-1)+2.0D 0*X(KA-1)*FA
        G(KA)=G(KA)+1.0D 0/X(KA)*FA
      ENDIF
      F=F+FA**2
  691 CONTINUE
      F=0.5D 0*F
      RETURN
  340 DO 341 KA=1,NA
      IF(KA.EQ.1) THEN
        FA=X(1)
        G(1)=G(1)+FA
      ELSE
        FA=COS(X(KA-1))+X(KA)-1.0D 0
        G(KA)=G(KA)+FA
        G(KA-1)=G(KA-1)-SIN(X(KA-1))*FA
      ENDIF
      F=F+FA**2
  341 CONTINUE
      F=0.5D 0*F
      RETURN
  360 DO 361 KA=1,NA
      S=(1.0D 0/DBLE(N+1))**2
      IF(N.EQ.1) THEN
        FA=2.0D 0*X(KA)-1.0D 0+S*(X(KA)+SIN(X(KA)))
        G(KA)=G(KA)+(2.0D 0+S*(1.0D 0+COS(X(KA))))*FA
      ELSE IF(KA.EQ.1) THEN
        FA=2.0D 0*X(KA)-X(KA+1)+S*(X(KA)+SIN(X(KA)))
        G(KA)=G(KA)+(2.0D 0+S*(1.0D 0+COS(X(KA))))*FA
        G(KA+1)=G(KA+1)-FA
      ELSE IF(KA.EQ.N) THEN
        FA=-X(KA-1)+2.0D 0*X(KA)-1.0D 0+S*(X(KA)+SIN(X(KA)))
        G(KA)=G(KA)+(2.0D 0+S*(1.0D 0+COS(X(KA))))*FA
        G(KA-1)=G(KA-1)-FA
      ELSE
        FA=-X(KA-1)+2.0D 0*X(KA)-X(KA+1)+S*(X(KA)+SIN(X(KA)))
        G(KA)=G(KA)+(2.0D 0+S*(1.0D 0+COS(X(KA))))*FA
        G(KA-1)=G(KA-1)-FA
        G(KA+1)=G(KA+1)-FA
      ENDIF
      F=F+FA**2
  361 CONTINUE
      F=0.5D 0*F
      RETURN
  380 DO 383 KA=1,NA
      IF(KA-5.GT.1) THEN
        I1=KA-5
      ELSE
        I1=1
      ENDIF
      IF(KA+1.LT.N) THEN
        I2=KA+1
      ELSE
        I2=N
      ENDIF
      S=0.0D 0
      DO 381 J=I1,I2
        IF(J.NE.KA) S=S+X(J)*(1.0D 0+X(J))
  381 CONTINUE
      FA=X(KA)*(2.0D 0+5.0D 0*X(KA)**2)+1.0D 0-S
      F=F+FA**2
      G(KA)=G(KA)+(2.0D 0+15.0D 0*X(KA)**2)*FA
      DO 382 J=I1,I2
      IF(J.NE.KA) G(J)=G(J)-(1.0D 0+2.0D 0*X(J))*FA
  382 CONTINUE
  383 CONTINUE
      F=0.5D 0*F
      RETURN
  430 DO 433 KA=1,NA
      ALF=5
      BET=14
      GAM=3
      FA=DBLE(BET*N)*X(KA)+(DBLE(KA)-DBLE(N)/2.0D 0)**GAM
      DO 431 J=1,N
        IF(J.NE.KA) THEN
        T=SQRT(X(J)**2+DBLE(KA)/DBLE(J))
        S1=LOG(T)
        FA=FA+T*(SIN(S1)**ALF+COS(S1)**ALF)
      ENDIF
  431 CONTINUE
      F=F+FA**2
      DO 432 J=1,N
        IF(J.NE.KA) THEN
        T=SQRT(X(J)**2+DBLE(KA)/DBLE(J))
        S1=LOG(T)
        G(J)=G(J)+(X(J)*(SIN(S1)**ALF+COS(S1)**ALF+
     +     ALF*SIN(S1)**(ALF-1)*COS(S1)-
     +     ALF*SIN(S1)*COS(S1)**(ALF-1))/T)*FA
        ELSE
          G(J)=G(J)+DBLE(BET*N)*FA
        ENDIF
  432 CONTINUE
  433 CONTINUE
      F=0.5D 0*F
      RETURN
  440 DO 445 KA=1,NA
      C=0.5D 0
      H=1.0D 0/DBLE(N)
      FA=(1.0D 0-C*H/4.0D 0)
      DO 441 J=1,N
        S=C*H*DBLE(KA)/DBLE(2*(KA+J))
        IF(J.EQ.N) S=S/2.0D 0
        FA=FA-S*X(J)
  441 CONTINUE
      FA=-1.0D 0+X(KA)*FA
      F=F+FA**2
      DO 442 J=1,N
        SX(J)=C*H*DBLE(KA)/DBLE(2*(KA+J))
  442 CONTINUE
      SX(N)=0.5D 0*SX(N)
      DO 444 J=1,N
        IF(KA.NE.J) THEN
           G(J)=G(J)-SX(J)*X(KA)*FA
        ELSE
           T=1.0D 0-C*H/4.0D 0
           DO 443 L=1,N
             IF(L.EQ.KA) THEN
               T=T-2.0D 0*SX(KA)*X(KA)
             ELSE
               T=T-SX(L)*X(L)
             ENDIF
  443      CONTINUE
           G(J)=G(J)+T*FA
        ENDIF
  444 CONTINUE
  445 CONTINUE
      F=0.5D 0*F
      RETURN
  270 DO 271 KA=1,NA
      S=0.5D 0
      H=1.0D 0/DBLE(N+1)
      T=2.0D 0*H**2
      IF(KA.EQ.1) THEN
      FA=2.0D 0*X(KA)-X(KA+1)-T*X(KA)**2+H*X(KA+1)
      G(KA)=G(KA)+2.0D 0*(1.0D 0-H**2*X(KA)/S)*FA
      G(KA+1)=G(KA+1)-(1.0D 0-H)*FA
      ELSE IF(1.LT.KA.AND.KA.LT.N) THEN
      FA=-X(KA-1)+2.0D 0*X(KA)-X(KA+1)-T*X(KA)**2+H*(X(KA+1)-X(KA-1))
      G(KA)=G(KA)+2.0D 0*(1.0D 0-T*X(KA))*FA
      G(KA-1)=G(KA-1)-(1.0D 0+H)*FA
      G(KA+1)=G(KA+1)-(1.0D 0-H)*FA
      ELSE IF(KA.EQ.N) THEN
      FA=-X(KA-1)+2.0D 0*X(KA)-0.5D 0-T*X(KA)**2+H*(0.5D 0-X(KA-1))
      G(KA)=G(KA)+2.0D 0*(1.0D 0-H**2*X(KA)/S)*FA
      G(KA-1)=G(KA-1)-(1.0D 0+H)*FA
      ENDIF
      F=F+FA**2
  271 CONTINUE
      F=0.5D 0*F
      RETURN
  280 DO 288 KA=1,NA
      S=0.5D 0
      H=1.0D 0/DBLE(N+1)
      T=H**2/S
      T1=2.0D 0*H
      AL=0.0D 0
      BE=0.5D 0
      S1=0.0D 0
      DO 281 J=1,KA
        IF(J.EQ.1) THEN
           S1=S1+DBLE(J)*(X(J)**2+(X(J+1)-AL)/T1)
        ENDIF
        IF(1.LT.J.AND.J.LT.N) THEN
           S1=S1+DBLE(J)*(X(J)**2+(X(J+1)-X(J-1))/T1)
        ENDIF
        IF(J.EQ.N) THEN
           S1=S1+DBLE(J)*(X(J)**2+(BE-X(J-1))/T1)
        ENDIF
  281 CONTINUE
      S1=(1.0D 0-DBLE(KA)*H)*S1
      IF(KA.EQ.N) GO TO 283
      S2=0.0D 0
      DO 282 J=KA+1,N
        IF(J.LT.N) THEN
           S2=S2+(1.0D 0-DBLE(J)*H)*(X(J)**2+(X(J+1)-X(J-1))/T1)
        ELSE
           S2=S2+(1.0D 0-DBLE(J)*H)*(X(J)**2+(BE-X(J-1))/T1)
        ENDIF
  282 CONTINUE
      S1=S1+DBLE(KA)*S2
  283 FA=X(KA)-0.5D 0*DBLE(KA)*H-T*S1
      F=F+FA**2
      S1=H**2/S
      S2=1.0D 0-DBLE(KA)*H
      DO 284 J=1,KA
        SX(J)=DBLE(J)*S2
  284 CONTINUE
      IF(KA.EQ.N) GO TO 286
      DO 285 J=KA+1,N
        SX(J)=DBLE(KA)*(1.0D 0-DBLE(J)*H)
  285 CONTINUE
  286 G(1)=G(1)-S1*(SX(1)*2.0D 0*X(1)-SX(2)/T1)*FA
      G(N)=G(N)-S1*(SX(N-1)/T1+SX(N)*2.0D 0*X(N))*FA
      DO 287 J=2,N-1
        G(J)=G(J)-S1*((SX(J-1)-SX(J+1))/T1+SX(J)*2.0D 0*X(J))*FA
  287 CONTINUE
      G(KA)=G(KA)+FA
      RETURN
  288 CONTINUE
      F=0.5D 0*F
      RETURN
  290 DO 291 KA=1,NA
      A=-9.0D-3
      B=1.0D-3
      AL=0.0D 0
      BE=25.0D 0
      GA=20.0D 0
      CA=0.3D 0
      CB=0.3D 0
      H=(B-A)/DBLE(N+1)
      T=A+DBLE(KA)*H
      H=H**2
      S=DBLE(KA)/DBLE(N+1)
      U=AL*(1.0D 0-S)+BE*S+X(KA)
      FF=CB*EXP(GA*(U-BE))-CA*EXP(GA*(AL-U))
      FG=FF*GA
      IF(T.LE.0) THEN
        FF=FF+CA
      ELSE
        FF=FF-CB
      ENDIF
      IF(N.EQ.1) THEN
        FA=-AL+2.0D 0*X(KA)-BE+H*FF
        G(KA)=G(KA)+(2.0D 0+H*FG)*FA
      ELSEIF(KA.EQ.1) THEN
        FA=-AL+2.0D 0*X(KA)-X(KA+1)+H*FF
        G(KA)=G(KA)+(2.0D 0+H*FG)*FA
        G(KA+1)=G(KA+1)-FA
      ELSEIF(KA.LT.N) THEN
        FA=-X(KA-1)+2.0D 0*X(KA)-X(KA+1)+H*FF
        G(KA)=G(KA)+(2.0D 0+H*FG)*FA
        G(KA-1)=G(KA-1)-FA
        G(KA+1)=G(KA+1)-FA
      ELSE
        FA=-X(KA-1)+2.0D 0*X(KA)-BE+H*FF
        G(KA)=G(KA)+(2.0D 0+H*FG)*FA
        G(KA-1)=G(KA-1)-FA
      ENDIF
      F=F+FA**2
  291 CONTINUE
      F=0.5D 0*F
      RETURN
  300 DO 301 KA=1,NA
      AL1=0.0D 0
      AL2=0.0D 0
      BE1=0.0D 0
      BE2=0.0D 0
      N1=N/2
      H=1.0D 0/DBLE(N1+1)
      T=DBLE(KA)*H
      IF(KA.EQ.1) THEN
        S1=2.0D 0*X(KA)-X(KA+1)
        B=AL1
      ELSE IF(KA.EQ.N1+1) THEN
        S1=2.0D 0*X(KA)-X(KA+1)
        B=AL2
      ELSE IF(KA.EQ.N1) THEN
        S1=-X(KA-1)+2.0D 0*X(KA)
        B=BE1
      ELSE IF(KA.EQ.N) THEN
        S1=-X(KA-1)+2.0D 0*X(KA)
        B=BE2
      ELSE
        S1=-X(KA-1)+2.0D 0*X(KA)-X(KA+1)
        B=0.0D 0
      ENDIF
      IF(KA.LE.N1) THEN
        S2=X(KA)**2+X(KA)+0.1D 0*X(N1+KA)**2-1.2D 0
      ELSE
        S2=0.2D 0*X(KA-N1)**2+X(KA)**2+2.0D 0*X(KA)-0.6D 0
      ENDIF
      FA=S1+H**2*S2-B
      F=F+FA**2
      H=1.0D 0/DBLE(N1+1)**2
      IF(KA.EQ.1) THEN
        G(KA)=G(KA)+(2.0D 0+H*(2.0D 0*X(KA)+1.0D 0))*FA
        G(KA+1)=G(KA+1)-FA
        G(N1+KA)=G(N1+KA)+H*0.2D 0*X(N1+KA)*FA
      ELSE IF(KA.EQ.N1+1) THEN
        G(1)=G(1)+H*0.4D 0*X(1)*FA
        G(KA)=G(KA)+(2.0D 0+H*(2.0D 0*X(KA)+2.0D 0))*FA
        G(KA+1)=G(KA+1)-FA
      ELSE IF(KA.EQ.N1) THEN
        G(KA-1)=G(KA-1)-FA
        G(KA)=G(KA)+(2.0D 0+H*(2.0D 0*X(KA)+1.0D 0))*FA
        G(N1+KA)=G(N1+KA)+H*0.2D 0*X(N1+KA)*FA
      ELSE IF(KA.EQ.N) THEN
        G(N1)=G(N1)+H*0.4D 0*X(N1)*FA
        G(KA-1)=G(KA-1)-FA
        G(KA)=G(KA)+(2.0D 0+H*(2.0D 0*X(KA)+2.0D 0))*FA
      ELSE IF(KA.LT.N1) THEN
        G(KA-1)=G(KA-1)-FA
        G(KA)=G(KA)+(2.0D 0+H*(2.0D 0*X(KA)+1.0D 0))*FA
        G(KA+1)=G(KA+1)-FA
        G(N1+KA)=G(N1+KA)+H*0.2D 0*X(N1+KA)*FA
      ELSE
         G(KA-N1)=G(KA-N1)+H*0.4D 0*X(KA-N1)*FA
         G(KA-1)=G(KA-1)-FA
         G(KA)=G(KA)+(2.0D 0+H*(2.0D 0*X(KA)+2.0D 0))*FA
         G(KA+1)=G(KA+1)-FA
      ENDIF
  301 CONTINUE
      F=0.5D 0*F
      RETURN
  710 DO 711 KA=1,NA
      ND=INT(SQRT(DBLE(N)))
      L=MOD(KA,ND)
      IF(L.EQ.0) THEN
         K=KA/ND
         L=ND
      ELSE
         K=INT(KA/ND)+1
      ENDIF
      LA=1.0D 0
      H=1.0D 0/DBLE(ND+1)
      H2=LA*H*H
      IF(L.EQ.1.AND.K.EQ.1) THEN
         FA=4.0D 0*X(1)-X(2)-X(ND+1)+H2*EXP(X(1))
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.1) THEN
         FA=4.0D 0*X(L)-X(L-1)-X(L+1)-X(L+ND)+H2*EXP(X(L))
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.1) THEN
         FA=4.0D 0*X(ND)-X(ND-1)-X(ND+ND)+H2*EXP(X(ND))
      ENDIF
      IF(L.EQ.1.AND.1.LT.K.AND.K.LT.ND) THEN
         FA=4.0D 0*X(KA)-X(KA+1)-X(KA-ND)-X(KA+ND)+H2*EXP(X(KA))
      ENDIF
      IF(L.EQ.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         FA=4.0D 0*X(KA)-X(KA-ND)-X(KA-1)-X(KA+ND)+H2*EXP(X(KA))
      ENDIF
      IF(L.EQ.1.AND.K.EQ.ND) THEN
         FA=4.0D 0*X(KA)-X(KA+1)-X(KA-ND)+H2*EXP(X(KA))
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.ND) THEN
         FA=4.0D 0*X(KA)-X(KA-1)-X(KA+1)-X(KA-ND)+H2*EXP(X(KA))
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.ND) THEN
        FA=4.0D 0*X(KA)-X(KA-1)-X(KA-ND)+H2*EXP(X(KA))
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.1.LT.K.AND.K.LT.ND) THEN
      FA=4.0D 0*X(KA)-X(KA-1)-X(KA+1)-X(KA-ND)-X(KA+ND)+H2*EXP(X(KA))
      ENDIF
      F=F+FA**2
      IF(L.EQ.1.AND.K.EQ.1) THEN
         G(1)=G(1)+(4.0D 0+H2*EXP(X(1)))*FA
         G(2)=G(2)-FA
         G(ND+1)=G(ND+1)-FA
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.1) THEN
         G(L)=G(L)+(4.0D 0+H2*EXP(X(L)))*FA
         G(L-1)=G(L-1)-FA
         G(L+1)=G(L+1)-FA
         G(L+ND)=G(L+ND)-FA
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.1) THEN
         G(ND)=G(ND)+(4.0D 0+H2*EXP(X(ND)))*FA
         G(ND-1)=G(ND-1)-FA
         G(ND+ND)=G(ND+ND)-FA
      ENDIF
      IF(L.EQ.1.AND.1.LT.K.AND.K.LT.ND) THEN
         G(KA)=G(KA)+(4.0D 0+H2*EXP(X(KA)))*FA
         G(KA-ND)=G(KA-ND)-FA
         G(KA+1)=G(KA+1)-FA
         G(KA+ND)=G(KA+ND)-FA
      ENDIF
      IF(L.EQ.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         G(KA)=G(KA)+(4.0D 0+H2*EXP(X(KA)))*FA
         G(KA-ND)=G(KA-ND)-FA
         G(KA-1)=G(KA-1)-FA
         G(KA+ND)=G(KA+ND)-FA
      ENDIF
      IF(L.EQ.1.AND.K.EQ.ND) THEN
         G(KA)=G(KA)+(4.0D 0+H2*EXP(X(KA)))*FA
         G(KA-ND)=G(KA-ND)-FA
         G(KA+1)=G(KA+1)-FA
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.ND) THEN
         G(KA)=G(KA)+(4.0D 0+H2*EXP(X(KA)))*FA
         G(KA-ND)=G(KA-ND)-FA
         G(KA-1)=G(KA-1)-FA
         G(KA+1)=G(KA+1)-FA
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.ND) THEN
         G(KA)=G(KA)+(4.0D 0+H2*EXP(X(KA)))*FA
         G(KA-ND)=G(KA-ND)-FA
         G(KA-1)=G(KA-1)-FA
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         G(KA)=G(KA)+(4.0D 0+H2*EXP(X(KA)))*FA
         G(KA-ND)=G(KA-ND)-FA
         G(KA-1)=G(KA-1)-FA
         G(KA+1)=G(KA+1)-FA
         G(KA+ND)=G(KA+ND)-FA
      ENDIF
  711 CONTINUE
      F=0.5D 0*F
      RETURN
  820 DO 821 KA=1,NA
      ND=INT(SQRT(DBLE(N)))
      L=MOD(KA,ND)
      IF(L.EQ.0) THEN
         K=KA/ND
         L=ND
      ELSE
         K=INT(KA/ND)+1
      ENDIF
      H=1.0D 0/DBLE(ND+1)
      H2=H*H
      IF(L.EQ.1.AND.K.EQ.1) THEN
         FA=4.0D 0*X(1)-X(2)-X(ND+1)+H2*X(1)**2-24.0D 0/(H+1.0D 0)**2
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.1) THEN
         FA=4.0D 0*X(L)-X(L-1)-X(L+1)-X(L+ND)+H2*X(L)**2
     *-12.0D 0/(DBLE(L)*H+1.0D 0)**2
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.1) THEN
         FA=4.0D 0*X(ND)-X(ND-1)-X(ND+ND)+H2*X(ND)**2
     *-12.0D 0/(DBLE(ND)*H+1.0D 0)**2-12.0D 0/(H+2.0D 0)**2
      ENDIF
      IF(L.EQ.1.AND.1.LT.K.AND.K.LT.ND) THEN
         FA=4.0D 0*X(KA)-X(KA+1)-X(KA-ND)-X(KA+ND)+H2*X(KA)**2
     *-12.0D 0/(DBLE(K)*H+1.0D 0)**2
      ENDIF
      IF(L.EQ.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         FA=4.0D 0*X(KA)-X(KA-ND)-X(KA-1)-X(KA+ND)+H2*X(KA)**2
     *-12.0D 0/(DBLE(K)*H+2.0D 0)**2
      ENDIF
      IF(L.EQ.1.AND.K.EQ.ND) THEN
         FA=4.0D 0*X(KA)-X(KA+1)-X(KA-ND)+H2*X(KA)**2
     *-12.0D 0/(DBLE(ND)*H+1.0D 0)**2-12.0D 0/(H+2.0D 0)**2
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.ND) THEN
         FA=4.0D 0*X(KA)-X(KA-1)-X(KA+1)-X(KA-ND)+H2*X(KA)**2
     *-12.0D 0/(DBLE(L)*H+2.0D 0)**2
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.ND) THEN
        FA=4.0D 0*X(KA)-X(KA-1)-X(KA-ND)+H2*X(KA)**2
     *-24.0D 0/(DBLE(ND)*H+2.0D 0)**2
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         FA=4.0D 0*X(KA)-X(KA-1)-X(KA+1)-X(KA-ND)-X(KA+ND)+H2*X(KA)**2
      ENDIF
      F=F+FA**2
      IF(L.EQ.1.AND.K.EQ.1) THEN
         G(1)=G(1)+(4.0D 0+H2*X(1)*2.0D 0)*FA
         G(2)=G(2)-FA
         G(ND+1)=G(ND+1)-FA
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.1) THEN
         G(L)=G(L)+(4.0D 0+H2*X(L)*2.0D 0)*FA
         G(L-1)=G(L-1)-FA
         G(L+1)=G(L+1)-FA
         G(L+ND)=G(L+ND)-FA
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.1) THEN
         G(ND)=G(ND)+(4.0D 0+H2*X(ND)*2.0D 0)*FA
         G(ND-1)=G(ND-1)-FA
         G(ND+ND)=G(ND+ND)-FA
      ENDIF
      IF(L.EQ.1.AND.1.LT.K.AND.K.LT.ND) THEN
         G(KA)=G(KA)+(4.0D 0+H2*X(KA)*2.0D 0)*FA
         G(KA-ND)=G(KA-ND)-FA
         G(KA+1)=G(KA+1)-FA
         G(KA+ND)=G(KA+ND)-FA
      ENDIF
      IF(L.EQ.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         G(KA)=G(KA)+(4.0D 0+H2*X(KA)*2.0D 0)*FA
         G(KA-ND)=G(KA-ND)-FA
         G(KA-1)=G(KA-1)-FA
         G(KA+ND)=G(KA+ND)-FA
      ENDIF
      IF(L.EQ.1.AND.K.EQ.ND) THEN
         G(KA)=G(KA)+(4.0D 0+H2*X(KA)*2.0D 0)*FA
         G(KA-ND)=G(KA-ND)-FA
         G(KA+1)=G(KA+1)-FA
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.K.EQ.ND) THEN
         G(KA)=G(KA)+(4.0D 0+H2*X(KA)*2.0D 0)*FA
         G(KA-ND)=G(KA-ND)-FA
         G(KA-1)=G(KA-1)-FA
         G(KA+1)=G(KA+1)-FA
      ENDIF
      IF(L.EQ.ND.AND.K.EQ.ND) THEN
         G(KA)=G(KA)+(4.0D 0+H2*X(KA)*2.0D 0)*FA
         G(KA-ND)=G(KA-ND)-FA
         G(KA-1)=G(KA-1)-FA
      ENDIF
      IF(1.LT.L.AND.L.LT.ND.AND.1.LT.K.AND.K.LT.ND) THEN
         G(KA)=G(KA)+(4.0D 0+H2*X(KA)*2.0D 0)*FA
         G(KA-ND)=G(KA-ND)-FA
         G(KA-1)=G(KA-1)-FA
         G(KA+1)=G(KA+1)-FA
         G(KA+ND)=G(KA+ND)-FA
      ENDIF
  821 CONTINUE
      F=0.5D 0*F
      RETURN
      END
* SUBROUTINE TILD22             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 94/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  INITIATION OF VARIABLES FOR NONLINEAR MINIMAX APPROXIMATION.
*  LINEARLY CONSTRAINED DENSE VERSION.
*
* PARAMETERS :
*  IO  N  NUMBER OF VARIABLES.
*  IO  NA  NUMBER OF PARTIAL FUNCTIONS.
*  IO  NB  NUMBER OF BOX CONSTRAINTS.
*  IO  NC  NUMBER OF GENERAL LINEAR CONSTRAINTS.
*  RO  X(N)  VECTOR OF VARIABLES.
*  IO  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  RO  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
*  RO  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
*  IO  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  RO  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RO  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RO  CG(NF*NC) MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  RO  FMIN  LOWER BOUND FOR VALUE OF THE OBJECTIVE FUNCTION.
*  RO  XMAX  MAXIMUM STEPSIZE.
*  IO  NEXT  NUMBER OF THE TEST PROBLEM.
*  IO  IEXT  TYPE OF OBJECTIVE FUNCTION. IEXT<0-MAXIMUM OF VALUES.
*         IEXT=0-MAXIMUM OF ABSOLUTE VALUES.
*  IO  IERR  ERROR INDICATOR.
*
      SUBROUTINE TILD22(N,NA,NB,NC,X,IX,XL,XU,IC,CL,CU,CG,FMIN,XMAX,
     +                  NEXT,IEXT,IERR)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979323846D0)
      DOUBLE PRECISION FMIN,XMAX
      INTEGER IERR,IEXT,N,NA,NB,NC,NEXT
      DOUBLE PRECISION CG(N*NC),CL(NC),CU(NC),X(N),XL(N),XU(N)
      INTEGER IC(NC),IX(N)
      DOUBLE PRECISION Y(163)
      INTEGER I,J,K,L
      COMMON /EMPR22/Y
      FMIN = -1.0D60
      XMAX = 1.0D3
      IEXT = -1
      IERR = 0
      NB = 0
      GO TO (10,20,30,40,50,100,160) NEXT
   10 IF (N.GE.2 .AND. NA.GE.3) THEN
          N = 2
          NA = 3
          NC = 1
          X(1) = 1.0D0
          X(2) = 2.0D0
          IC(1) = 1
          CL(1) = 0.5D0
          CG(1) = 1.0D0
          CG(2) = 1.0D0
      ELSE
          IERR = 1
      END IF
      RETURN
   20 IF (N.GE.2 .AND. NA.GE.3) THEN
          N = 2
          NA = 3
          NC = 1
          X(1) = -2.0D0
          X(2) = -1.0D0
          IC(1) = 2
          CU(1) = -2.5D0
          CG(1) = 3.0D0
          CG(2) = 1.0D0
      ELSE
          IERR = 1
      END IF
      RETURN
   30 IF (N.GE.2 .AND. NA.GE.3) THEN
          N = 2
          NA = 3
          NB = N
          NC = 1
          X(1) = -1.0D0
          X(2) = 1.0D-2
          IX(1) = 0
          IX(2) = 1
          XL(2) = 1.0D-2
          IC(1) = 1
          CL(1) = -5.0D-1
          CG(1) = 5.0D-2
          CG(2) = -1.0D0
      ELSE
          IERR = 1
      END IF
      RETURN
   40 IF (N.GE.2 .AND. NA.GE.3) THEN
          N = 2
          NA = 3
          NB = N
          NC = 1
          X(1) = -1.0D0
          X(2) = 3.0D0
          IX(1) = 0
          IX(2) = 1
          XL(2) = 1.0D-2
          IC(1) = 1
          CL(1) = 1.0D0
          CG(1) = -9.0D-1
          CG(2) = 1.0D0
      ELSE
          IERR = 1
      END IF
      RETURN
   50 IF (N.GE.6 .AND. NA.GE.3) THEN
          N = 6
          NA = 3
          X(1) = -1.0D0
          X(2) = 0.0D0
          X(3) = 0.0D0
          X(4) = -1.0D0
          X(5) = 1.0D0
          X(6) = 1.0D0
          NC = 5*NA
          DO 60 I = 1,NC
              CU(I) = 1.0D0
              IC(I) = 2
   60     CONTINUE
          DO 70 I = 1,N*NC
              CG(I) = 0.0D0
   70     CONTINUE
          K = 1
          DO 90 I = 1,NA
              L = 2* (I-1)
              DO 80 J = 1,5
                  CG(K+L) = SIN(2.0D0*PI*DBLE(J-1)/5.0D0)
                  CG(K+L+1) = COS(2.0D0*PI*DBLE(J-1)/5.0D0)
                  K = K + N
   80         CONTINUE
   90     CONTINUE
      ELSE
          IERR = 1
      END IF
      RETURN
  100 IF (N.GE.7 .AND. NA.GE.163) THEN
          N = 7
          NA = 163
          NB = N
          DO 110 I = 1,N
              X(I) = DBLE(I)*0.5D0
              IX(I) = 0
  110     CONTINUE
          XL(1) = 0.4D0
          IX(1) = 1
          IX(7) = 5
          DO 120 I = 1,NA
              Y(I) = 2.0D0*PI*SIN(PI* (8.5D0+DBLE(I)*0.5D0)/180.0D0)
  120     CONTINUE
          NC = 7
          DO 130 I = 1,6
              CL(I) = 0.4D0
              IC(I) = 1
  130     CONTINUE
          CL(7) = 1.0D0
          CU(7) = 1.0D0
          IC(7) = 5
          DO 140 I = 1,N*NC
              CG(I) = 0.0D0
  140     CONTINUE
          K = 0
          DO 150 I = 1,6
              CG(K+I) = -1.0D0
              CG(K+I+1) = 1.0D0
              K = K + N
  150     CONTINUE
          CG(46) = -1.0D0
          CG(48) = 1.0D0
          IEXT = 0
          FMIN = 0.0D0
      ELSE
          IERR = 1
      END IF
      RETURN
  160 IF (N.GE.8 .AND. NA.GE.8) THEN
          N = 8
          NA = 8
          NB = N
          DO 170 I = 1,N
              X(I) = 0.125D0
              XL(I) = 1.0D-8
              IX(I) = 1
  170     CONTINUE
          DO 180 I = 1,40
              Y(I) = 1.0D0
              Y(I+40) = 0.1D0
  180     CONTINUE
          Y(9) = 2.0D0
          Y(10) = 0.8D0
          Y(12) = 0.5D0
          Y(18) = 1.2D0
          Y(19) = 0.8D0
          Y(20) = 1.2D0
          Y(21) = 1.6D0
          Y(22) = 2.0D0
          Y(23) = 0.6D0
          Y(24) = 0.1D0
          Y(25) = 2.0D0
          Y(26) = 0.1D0
          Y(27) = 0.6D0
          Y(28) = 2.0D0
          Y(32) = 2.0D0
          Y(33) = 1.2D0
          Y(34) = 1.2D0
          Y(35) = 0.8D0
          Y(37) = 1.2D0
          Y(38) = 0.1D0
          Y(39) = 3.0D0
          Y(40) = 4.0D0
          Y(41) = 3.0D0
          Y(42) = 1.0D0
          Y(45) = 5.0D0
          Y(48) = 6.0D0
          Y(50) = 1.0D1
          Y(53) = 5.0D0
          Y(58) = 9.0D0
          Y(59) = 1.0D1
          Y(61) = 4.0D0
          Y(63) = 7.0D0
          Y(68) = 1.0D1
          Y(70) = 3.0D0
          Y(80) = 1.1D1
          Y(81) = 0.5D0
          Y(82) = 1.2D0
          Y(83) = 0.8D0
          Y(84) = 2.0D0
          Y(85) = 1.5D0
          NC = 1
          CL(1) = 1.0D0
          CU(1) = 1.0D0
          IC(1) = 5
          DO 190 I = 1,N
              CG(I) = 1.0D0
  190     CONTINUE
      ELSE
          IERR = 1
      END IF
      RETURN
      END
* SUBROUTINE TAFU22             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 94/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  VALUES OF PARTIAL FUNCTIONS IN THE MINIMAX CRITERION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  KA  INDEX OF THE PARTIAL FUNCTION.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  FA  VALUE OF THE PARTIAL FUNCTION AT THE
*          SELECTED POINT.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TAFU22(N,KA,X,FA,NEXT)
      DOUBLE PRECISION FA
      INTEGER KA,N,NEXT
      DOUBLE PRECISION X(N)
      DOUBLE PRECISION Y(163)
      DOUBLE PRECISION A,P
      INTEGER I,J,K
      COMMON /EMPR22/Y
      GO TO (10,10,50,50,90,130,150) NEXT
   10 GO TO (20,30,40) KA
   20 FA = X(1)**2 + X(2)**2 + X(1)*X(2) - 1.0D0
      RETURN
   30 FA = SIN(X(1))
      RETURN
   40 FA = -COS(X(2))
      RETURN
   50 GO TO (60,70,80) KA
   60 FA = -EXP(X(1)-X(2))
      RETURN
   70 FA = SINH(X(1)-1.0D0) - 1.0D0
      RETURN
   80 FA = -LOG(X(2)) - 1.0D0
      RETURN
   90 GO TO (100,110,120) KA
  100 FA = -SQRT((X(1)-X(3))**2+ (X(2)-X(4))**2)
      RETURN
  110 FA = -SQRT((X(3)-X(5))**2+ (X(4)-X(6))**2)
      RETURN
  120 FA = -SQRT((X(5)-X(1))**2+ (X(6)-X(2))**2)
      RETURN
  130 A = 0.0D0
      DO 140 I = 1,N
          A = A + COS(Y(KA)*X(I))
  140 CONTINUE
      FA = (1.0D0+2.0D0*A)/1.5D1
      RETURN
  150 FA = 0.0D0
      K = 0
      DO 170 I = 1,5
          A = 0.0D0
          P = 0.0D0
          DO 160 J = 1,N
              A = A + Y(K+J)*X(J)** (1.0D0-Y(I+80))
              P = P + Y(K+J+40)*X(J)
  160     CONTINUE
          FA = FA + Y(K+KA)*P/ (X(KA)**Y(I+80)*A) - Y(K+KA+40)
          K = K + N
  170 CONTINUE
      RETURN
      END
* SUBROUTINE TAGU22             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 94/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  GRADIENTS OF PARTIAL FUNCTIONS IN THE MINIMAX CRITERION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  KA  INDEX OF THE PARTIAL FUNCTION.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  GA(N)  GRADIENT OF THE PARTIAL FUNCTION AT THE
*          SELECTED POINT.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TAGU22(N,KA,X,GA,NEXT)
      INTEGER KA,N,NEXT
      DOUBLE PRECISION GA(N),X(N)
      DOUBLE PRECISION Y(163)
      DOUBLE PRECISION A,B,C,P
      INTEGER I,J,K
      COMMON /EMPR22/Y
      GO TO (10,10,50,50,90,130,150) NEXT
   10 GO TO (20,30,40) KA
   20 GA(1) = 2.0D0*X(1) + X(2)
      GA(2) = 2.0D0*X(2) + X(1)
      RETURN
   30 GA(1) = COS(X(1))
      GA(2) = 0.0D0
      RETURN
   40 GA(1) = 0.0D0
      GA(2) = SIN(X(2))
      RETURN
   50 GO TO (60,70,80) KA
   60 GA(1) = -EXP(X(1)-X(2))
      GA(2) = EXP(X(1)-X(2))
      RETURN
   70 GA(1) = COSH(X(1)-1.0D0)
      GA(2) = 0.0D0
      RETURN
   80 GA(1) = 0.0D0
      GA(2) = -1.0D0/X(2)
      RETURN
   90 GO TO (100,110,120) KA
  100 A = SQRT((X(1)-X(3))**2+ (X(2)-X(4))**2)
      GA(1) = - (X(1)-X(3))/A
      GA(2) = - (X(2)-X(4))/A
      GA(3) = -GA(1)
      GA(4) = -GA(2)
      GA(5) = 0.0D0
      GA(6) = 0.0D0
      RETURN
  110 A = SQRT((X(3)-X(5))**2+ (X(4)-X(6))**2)
      GA(1) = 0.0D0
      GA(2) = 0.0D0
      GA(3) = - (X(3)-X(5))/A
      GA(4) = - (X(4)-X(6))/A
      GA(5) = -GA(3)
      GA(6) = -GA(4)
      RETURN
  120 A = SQRT((X(5)-X(1))**2+ (X(6)-X(2))**2)
      GA(1) = (X(5)-X(1))/A
      GA(2) = (X(6)-X(2))/A
      GA(3) = 0.0D0
      GA(4) = 0.0D0
      GA(5) = -GA(1)
      GA(6) = -GA(2)
      RETURN
  130 DO 140 I = 1,N
          GA(I) = -2.0D0*Y(KA)*SIN(Y(KA)*X(I))/1.5D1
  140 CONTINUE
      RETURN
  150 DO 160 I = 1,N
          GA(I) = 0.0D0
  160 CONTINUE
      K = 0
      DO 190 I = 1,5
          A = 0.0D0
          P = 0.0D0
          DO 170 J = 1,N
              A = A + Y(K+J)*X(J)** (1.0D0-Y(I+80))
              P = P + Y(K+J+40)*X(J)
  170     CONTINUE
          B = Y(K+KA)/ (X(KA)**Y(I+80)*A)
          DO 180 J = 1,N
              C = Y(K+J)* (1.0D0-Y(I+80))/ (X(J)**Y(I+80)*A)
              GA(J) = GA(J) + B* (Y(K+J+40)-C*P)
  180     CONTINUE
          GA(KA) = GA(KA) - B*Y(I+80)*P/X(KA)
          K = K + N
  190 CONTINUE
      RETURN
      END
* SUBROUTINE TYTIM1                MS DOS                     91/12/01
C PORTABILITY : MS DOS / MS FORTRAN v.5.0
C 91/12/01 SI : ORIGINAL VERSION
*
* PURPOSE :
*  GET TIME IN 100TH OF SEC.
*
      SUBROUTINE TYTIM1(ITIME)
      INTEGER ITIME
      REAL*4 TIME
      CALL CPU_TIME(TIME)
      ITIME=1.0D2*TIME
      END
* SUBROUTINE TYTIM2                ALL SYSTEMS                91/12/01
C PORTABILITY : ALL SYSTEMS
C 91/12/01 SI : ORIGINAL VERSION
*
* PURPOSE :
*  PRINT TIME ELAPSED.
*
      SUBROUTINE TYTIM2(ITIME)
      INTEGER ITIME
      INTEGER IHR,IT,IMIN,ISEC
      CALL TYTIM1(IT)
      IT=IT-ITIME
      IHR=IT/(60*60*100)
      IT=IT-IHR*60*60*100
      IMIN=IT/(60*100)
      IT=IT-IMIN*60*100
      ISEC=IT/100
      IT=IT-ISEC*100
      WRITE(6,10) IHR,IMIN,ISEC,IT
   10 FORMAT(' TIME=',I2,':',I2.2,':',I2.2,'.',I2.2)
      END
