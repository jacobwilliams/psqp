
PSQP - A SEQUENTIAL QUADRATIC PROGRAMMING ALGORITHM
FOR GENERAL NONLINEAR PROGRAMMING PROBLEMS.

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![CI Status](https://github.com/jacobwilliams/psqp/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/psqp/actions)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/psqp)](https://github.com/jacobwilliams/psqp/commits/master)


This is work in progress to refactor the original [PSQP](http://www.cs.cas.cz/~luksan/subroutines.html) (by Ladislav Luksan) into Modern Fortran. The latest API documentation can be found [here](https://jacobwilliams.github.io/psqp/index.html).

# Original Readme

## Introduction

The double-precision FORTRAN 77 basic subroutine PSQP is designed
to find a close approximation to a local minimum of a nonlinear
objective function F(X) with simple bounds on variables and general
nonlinear constraints. Here X is a vector of N variables and F(X), is
a smooth function. Simple bounds are assumed in the form
```
               X(I) unbounded if  IX(I) = 0,
      XL(I) <= X(I)           if  IX(I) = 1,
               X(I) <= XU(I)  if  IX(I) = 2,
      XL(I) <= X(I) <= XU(I)  if  IX(I) = 3,
      XL(I)  = X(I)  = XU(I)  if  IX(I) = 5,
```
where 1 <= I <= N. General nonlinear constraints are assumed in the form
```
               C_I(X) unbounded if  IC(I) = 0,
      CL(I) <= C_I(X)           if  IC(I) = 1,
               C_I(X) <= CU(I)  if  IC(I) = 2,
      CL(I) <= C_I(X) <= CU(I)  if  IC(I) = 3,
      CL(I)  = C_I(X)  = CU(I)  if  IC(I) = 5,
```
where C_I(X), 1 <= I <= NC, are twice continuously differentiable
functions.

To simplify user's work, an additional easy to use subroutine
PSQPN is added. It calls the basic general subroutine PSQP. All
subroutines contain a description of formal parameters and extensive
comments. Furthermore, test program TSQPN is included, which contains
several test problems. This test programs serve as an example for
using the subroutine PSQPN, verify its correctness and demonstrate
its efficiency.

In this short guide, we describe all subroutines which can be
called from the user's program. In the description of formal parameters,
we introduce a type of the argument that specifies whether the argument
must have a value defined on entry to the subroutine (I), whether it
is a value which will be returned (O), or both (U), or whether it is
an auxiliary value (A). Note that the arguments of the type I can be
changed on output under some circumstances, especially if improper
input values were given. Besides formal parameters, we can use a
COMMON /STAT/ block containing statistical information. This block,
used in each subroutine has the following form:
```
      COMMON /STAT/ NRES,NDEC,NREM,NADD,NIT,NFV,NFG,NFH
```
The arguments have the following meaning:
```
 Argument  Type Significance
 ----------------------------------------------------------------------
  NRES      O   Positive INTEGER variable that indicates the number of
                restarts.
  NDEC      O   Positive INTEGER variable that indicates the number of
                matrix decompositions.
  NRED      O   Positive INTEGER variable that indicates the number of
                reductions.
  NREM      O   Positive INTEGER variable that indicates the number of
                constraint deletions during the QP solutions.
  NADD      O   Positive INTEGER variable that indicates the number of
                constraint additions during the QP solutions.
  NIT       O   Positive INTEGER variable that indicates the number of
                iterations.
  NFV       O   Positive INTEGER variable that indicates the number of
                function evaluations.
  NFG       O   Positive INTEGER variable that specifies the number of
                gradient evaluations.
  NFH       O   Positive INTEGER variable that specifies the number of
                Hessian evaluations.
```

## Subroutine PSQN

The calling sequence is
```
      CALL PSQPN(NF,NB,NC,X,IX,XL,XU,CF,IC,CL,CU,IPAR,RPAR,F,GMAX,
     &           CMAX,IPRNT,ITERM)
```
The arguments have the following meaning.
```
 Argument  Type Significance
 ----------------------------------------------------------------------
  NF        I   Positive INTEGER variable that specifies the number of
                variables of the objective function.
  NB        I   Nonnegative INTEGER variable that specifies whether the
                simple bounds are suppressed (NB=0) or accepted (NB>0).
  NC        I   Nonnegative INTEGER variable that specifies the number
                of general nonlinear constraints; if NC=0 the
                general nonlinear constraints are suppressed.
  X(NF)     U   On input, DOUBLE PRECISION vector with the initial
                estimate to the solution. On output, the approximation
                to the minimum.
  IX(NF)    I   On input (significant only if NB>0) INTEGER vector
                containing the simple bounds types:
                   IX(I)=0 - the variable X(I) is unbounded,
                   IX(I)=1 - the lower bound X(I) >= XL(I),
                   IX(I)=2 - the upper bound X(I) <= XU(I),
                   IX(I)=3 - the two side bound XL(I) <= X(I) <= XU(I),
                   IX(I)=5 - the variable X(I) is fixed (given by its
                             initial estimate).
  XL(NF)    I   DOUBLE PRECISION vector with lower bounds for variables
                (significant only if NB>0).
  XU(NF)    I   DOUBLE PRECISION vector with upper bounds for variables
                (significant only if NB>0).
  CF(NC)    A   DOUBLE PRECISION vector which contains values of
                constraint functions (only if NC>0).
  IC(NC)    I   On input (significant only if NC>0) INTEGER vector which
                contains constraint types:
                  IC(K)=0 - the constraint CF(K) is not used,
                  IC(K)=1 - the lower constraint CF(K) >= CL(K),
                  IC(K)=2 - the upper constraint CF(K) <= CU(K),
                  IC(K)=3 - the two side constraint
                            CL(K) <= CF(K) <= CU(K),
                  IC(K)=5 - the equality constraint CF(K) = CL(K).
  CL(NC)    I   DOUBLE PRECISION vector with lower bounds for constraint
                functions (significant only if NC>0).
  CU(NC)    I   DOUBLE PRECISION vector with upper bounds for constraint
                functions (significant only if NC>0).
  IPAR(6)   A   INTEGER parameters:
                  IPAR(1)=MIT,  IPAR(2)=MFV, IPAR(3)-NONE,
                  IPAR(4)-NONE, IPAR(5)=MET, IPAR(5)=MEC.
                Parameters MIT, MFV, MET, MEC are described in Section 3
                together with other parameters of the subroutine PSQP.
  RPAR(5)   A   DOUBLE PRECISION parameters:
                  RPAR(1)=XMAX,  RPAR(2)=TOLX,  RPAR(3)=TOLC,
                  RPAR(4)=TOLG,  RPAR(5)=RPF.
                Parameters XMAX, TOLX, TOLC, TOLG, RPF are described
                in Section 3 together with other parameters of the
                subroutine PSQP.
  F         O   DOUBLE PRECISION value of the objective function at the
                solution X.
  GMAX      O   DOUBLE PRECISION maximum absolute value of a partial
                derivative of the Lagrangian function.
  CMAX      O   maximum constraint violation.
  IPRNT     I   INTEGER variable that specifies PRINT:
                  IPRNT= 0 - print is suppressed,
                  IPRNT= 1 - basic print of final results,
                  IPRNT=-1 - extended print of final results,
                  IPRNT= 2 - basic print of intermediate and final
                             results,
                  IPRNT=-2 - extended print of intermediate and final
                             results.
  ITERM     O   INTEGER variable that indicates the cause of termination:
                  ITERM= 1 - if |X - XO| was less than or equal to TOLX
                             in MTESX subsequent iterations,
                  ITERM= 2 - if |F - FO| was less than or equal to TOLF
                             in MTESF subsequent iterations,
                  ITERM= 3 - if F is less than or equal to TOLB,
                  ITERM= 4 - if GMAX is less than or equal to TOLG,
                  ITERM=11 - if NIT exceeded MIT,
                  ITERM=12 - if NFV exceeded MFV,
                  ITERM< 0 - if the method failed. If ITERM=-6, then the
                             termination criterion has not been satisfied,
                             but the point obtained is usually acceptable.
```
The subroutine PSQPN requires the user supplied subroutines OBJ,
DOBJ that define the objective function and its gradient and CON, DCON
that define constraint functions and their gradients. These subroutines
have the form
```
      SUBROUTINE  OBJ(NF,X,F)
      SUBROUTINE DOBJ(NF,X,G)
      SUBROUTINE  CON(NF,KC,X,FC)
      SUBROUTINE DCON(NF,KC,X,GC)
```
The arguments of the user supplied subroutine have the following meaning.
```
 Argument  Type Significance
 ----------------------------------------------------------------------
  NF        I   Positive INTEGER variable that specifies the number of
                variables of the objective function.
  X(NF)     I   DOUBLE PRECISION an estimate to the solution.
  F         O   DOUBLE PRECISION value of the objective function at the
                point X.
  G(NF)     O   DOUBLE PRECISION gradient of the objective function
                at the point X.
  KC        I   INTEGER index of the partial function.
  FC        O   DOUBLE PRECISION value of the KC-th partial function at
                the point X.
  GC(NF)    O   DOUBLE PRECISION gradient of the KC-th partial function
                at the point X.
```

## Subroutine PSQP

This general subroutine is called from all the subroutines described
in Section 2. The calling sequence is
```
      CALL PSQP(NF,NB,NC,X,IX,XL,XU,CF,IC,CL,CU,CG,CFO,CFD,GC,ICA,CR,
     & CZ,CP,GF,G,H,S,XO,GO,XMAX,TOLX,TOLC,TOLG,RPF,CMAX,GMAX,F,MIT,
     & MFV,MEC,IPRNT,ITERM).
```
The arguments NF, NB, NC, X, IX, XL, XU, CF, IC, CL, CU, CMAX, GMAX,
F, IPRNT, ITERM, have the same meaning as in Section 2. Other arguments
have the following meaning:
```
 Argument  Type Significance
 ----------------------------------------------------------------------
  CG(NF*NC) A   DOUBLE PRECISION elements of the constraint Jacobian
                matrix.
  CFO(NC+1) A   DOUBLE PRECISION vector which contains old values of
                 constraint functions.
  CFD(NC)   A   DOUBLE PRECISION vector of constraint function
                increments.
  GC(NF)    A   DOUBLE PRECISION gradient of the constraint function.
  ICA(NC)   A   INTEGER vector containing indices of active constraints.
  CR(NCR)   A   DOUBLE PRECISION matrix containing triangular
                decomposition of the orthogonal projection kernel
                (NCR is equal to NF*(NF+1)/2).
  CZ(NF)    A   DOUBLE PRECISION vector of Lagrange multipliers.
                of linear manifold defined by active constraints.
  CP(NF)    A   DOUBLE PRECISION auxiliary array.
  GF(NF)    A   DOUBLE PRECISION gradient of the objective function.
  G(NF)     A   DOUBLE PRECISION gradient of the Lagrangian function.
  H(NH)     A   DOUBLE PRECISION variable metric approximation of the
                Hessian matrix of the Lagrangian function.
  S(NF)     A   DOUBLE PRECISION direction vector.
  XO(NF)    A   DOUBLE PRECISION vector which contains increments of
                variables.
  GO(NF)    A   DOUBLE PRECISION vector which contains increments of
                gradients.
  XMAX      I   DOUBLE PRECISION maximum stepsize; the choice XMAX=0
                causes that the default value 1.0D+3 will be taken.
  TOLX      I   DOUBLE PRECISION tolerance for the change of the
                coordinate vector X; the choice TOLX=0 causes that the
                default value TOLX=1.0D-16 will be taken.
  TOLC      I   DOUBLE PRECISION tolerance for the constraint violation;
                the choice TOLC=0 causes that the default
                value TOLC=1.0D-6 will be taken.
  TOLG      I   DOUBLE PRECISION tolerance for the Lagrangian function
                gradient; the choice TOLG=0 causes that the default
                value TOLG=1.0D-6 will be taken.
  RPF       I   DOUBLE PRECISION value of the penalty coefficient; the
                choice RPF=0 causes that the default value 1.0D-4 will
                be taken.
  MIT       I   INTEGER variable that specifies the maximum number of
                iterations; the choice MIT=0 causes that the default
                value 1000 will be taken.
  MFV       I   INTEGER variable that specifies the maximum number of
                function evaluations; the choice |MFV|=0 causes that
                the default value 2000 will be taken.
  MET       I   INTEGER variable that specifies the variable metric
                update.
                  MET=1 - the BFGS update is used.
                  MET=2 - the Hoshino update is used.
                The choice MET=0 causes that the default value MET=1
                will be taken.
  MEC       I   INTEGER variable that specifies a correction when the
                negative curvature is detected:
                  MEC=1 - correction is not used.
                  MEC=2 - the Powell correction is used.
                The choice MEC=0 causes that the default value MEC=2
                will be taken.
```
The default velue RPF=0.0001 is relatively small. Therefore, larger
value (RPF=1 say) can sometimes be more suitable.

The subroutine PSQP requires the user supplied subroutines OBJ
DOBJ, CON, DCON,  which are described in Section 2.

## Subroutine PLQDB1

Since the dual range space method for solving quadratic programming
subproblems arising in sequential quadratic programming algorithms
can be used separately in many applications, we describe the subroutine
PLQDB1 in more details. The calling sequence is
```
      CALL PLQDB1(NF,NC,X,IX,XL,XU,CF,CFD,IC,ICA,CL,CU,CG,CR,CZ,G,GO,
     & H,S,MFP,KBF,KBC,IDECF,ETA2,ETA9,EPS7,EPS9,UMAX,GMAX,N,ITERQ)
```
The arguments NF, NC, X, IX, XL, XU, CF, IC, CL, CU, have the same
meaning as in Section 2 (only with the difference that the
arguments X and CF are of the type (I), i.e. they  must have a value
defined on entry to ULQDF1 and they are not changed). The arguments
CFD, ICA, CG, CR, CZ have the same meaning as in Section 3 (only with
the difference that the arguments CFD, ICA, CR, CZ are of the type (O),
i.e. their values can be used subsequently). Other arguments have the
following meaning:
```
 Argument  Type Significance
 ----------------------------------------------------------------------
  G(NF)     O   DOUBLE PRECISION gradient of the Lagrangian function.
  GO(NF)    A   DOUBLE PRECISION old gradient of the Lagrangian
                function.
  H(NH)     U   DOUBLE PRECISION Choleski decomposition of the
                approximate Hessian (NH is equal to NF*(NF+1)/2).
  S(NF)     O   DOUBLE PRECISION direction vector.
  MFP       I   INTEGER variable that specifies the type of the
                computed point.
                  MFP=1 - computation is terminated whenever an
                          arbitrary feasible point is found,
                  MFP=2 - computation is terminated whenever an
                          optimum feasible point is found,
                  MFP=3 - computation starts from the previously
                          reached point and is terminated whenever
                          an optimum feasible point is found.
  KBF       I   INTEGER variable that specifies simple bounds on
                variables.
                  KBF=0 - simple bounds are suppressed,
                  KBF=1 - one sided simple bounds,
                  KBF=2 - two sided simple bounds.
  KBC       I   INTEGER variable that specifies general linear
                constraints.
                  KBC=0 - linear constraints are suppressed,
                  KBC=1 - one sided linear constraints,
                  KBC=2 - two sided linear constraints.
  IDECF     U   INTEGER variable that specifies the type of matrix
                decomposition.
                  IDECF= 0 - no decomposition,
                  IDECF= 1 - Choleski decomposition,
                  IDECF= 9 - inversion,
                  IDECF=10 - diagonal matrix.
  ETA2      I   DOUBLE PRECISION tolerance for positive definiteness
                in the Choleski decomposition.
  ETA9      I   DOUBLE PRECISION maximum floating point number.
  EPS7      I   DOUBLE PRECISION tolerance for linear independence
                of constraints (the recommended value is 1.0D-10).
  EPS9      I   DOUBLE PRECISION tolerance for the definition of active
                constraints (the recommended value is 1.0D-8).
  UMAX      O   DOUBLE PRECISION maximum absolute value of the negative
                Lagrange multiplier.
  GMAX      O   DOUBLE PRECISION infinity norm of the gradient of the
                Lagrangian function.
  N         O   INTEGER dimension of a manifold defined by active
                constraints.
  ITERQ     O   INTEGER variable that indicates the type of the
                computed feasible point.
                  ITERQ= 1 - an arbitrary feasible point was found,
                  ITERQ= 2 - the optimum feasible point was found,
                  ITERQ=-1 - an arbitrary feasible point does not
                             exist,
                  ITERQ=-2 - the optimum feasible point does not
                             exist.
```

## Verification of the subroutines

Subroutine PSQPN can be verified and tested using the program
TSQPN. This program calls the subroutines TIND07 (initiation),
TFFU07 (objective function evaluation), TFGU07 (objective gradient
evaluation), TCFU07 (constraint functions evaluation) and TFGU07
(constraint gradients evaluation) containing 34 equality constrained
test problems with at most 20 variables. The results obtained
by the program TSQPN on a PC computer with Microsoft Power Station
Fortran compiler have the following form.
```
NIT=   7  NFV=   7  NFG=   7  F=-1.41421      C=0.8E-08  G=0.4E-06  ITERM=  4
NIT=  12  NFV=  12  NFG=  12  F=-1.00000      C=0.5E-10  G=0.6E-07  ITERM=  4
NIT=  10  NFV=  11  NFG=  10  F=-30.0000      C=0.3E-10  G=0.8E-08  ITERM=  4
NIT=   6  NFV=   6  NFG=   6  F= 1.28072      C=0.0E+00  G=0.4E-06  ITERM=  4
NIT=  36  NFV=  40  NFG=  36  F=0.284597E-01  C=0.0E+00  G=0.1E-06  ITERM=  4
NIT=  11  NFV=  13  NFG=  11  F= 1.00000      C=0.0E+00  G=0.4E-06  ITERM=  4
NIT=   6  NFV=   6  NFG=   6  F=-6961.81      C=0.7E-08  G=0.1E-06  ITERM=  4
NIT=   4  NFV=   4  NFG=   4  F= 40.1987      C=0.0E+00  G=0.2E-06  ITERM=  4
NIT=   9  NFV=   9  NFG=   9  F= 13.3567      C=0.3E-06  G=0.8E-07  ITERM=  4
NIT=  16  NFV=  23  NFG=  16  F=-22.6274      C=0.9E-11  G=0.3E-06  ITERM=  4
NIT=  14  NFV=  14  NFG=  14  F= 1.00000      C=0.0E+00  G=0.8E-06  ITERM=  4
NIT=  10  NFV=  13  NFG=  10  F= 6.00000      C=0.2E-12  G=0.1E-06  ITERM=  4
NIT=  54  NFV=  71  NFG=  54  F= 6299.84      C=0.4E-15  G=0.2E-07  ITERM=  4
NIT=   8  NFV=   8  NFG=   8  F=-.834032      C=0.3E-07  G=0.3E-08  ITERM=  4
NIT=  30  NFV=  30  NFG=  30  F=-1162.12      C=0.0E+00  G=0.1E-06  ITERM=  4
NIT=  74  NFV=  76  NFG=  74  F= 4.75530      C=0.0E+00  G=0.4E-06  ITERM=  4
NIT=  22  NFV=  22  NFG=  22  F= 727.679      C=0.4E-14  G=0.7E-06  ITERM=  4
NIT=  11  NFV=  14  NFG=  11  F=-44.0000      C=0.9E-12  G=0.8E-07  ITERM=  4
NIT=  29  NFV=  35  NFG=  29  F=-210.408      C=0.2E-13  G=0.1E-06  ITERM=  4
NIT=   4  NFV=   4  NFG=   4  F=-30665.0      C=0.3E-08  G=0.1E-06  ITERM=  4
NIT=  11  NFV=  30  NFG=  11  F=-.528034E+07  C=0.0E+00  G=0.1E-07  ITERM=  4
NIT=  21  NFV=  21  NFG=  21  F=-1.90516      C=0.1E-08  G=0.2E-10  ITERM=  4
NIT=  10  NFV=  12  NFG=  10  F= 5.00000      C=0.2E-11  G=0.8E-07  ITERM=  4
NIT=  45  NFV=  50  NFG=  45  F= 135.076      C=0.3E-14  G=0.3E-07  ITERM=  4
NIT=  10  NFV=  10  NFG=  10  F= 4.07125      C=0.0E+00  G=0.7E-10  ITERM=  4
NIT=  16  NFV=  19  NFG=  16  F= 680.630      C=0.1E-11  G=0.5E-06  ITERM=  4
NIT=  26  NFV=  35  NFG=  26  F= 221.878      C=0.2E-15  G=0.3E-06  ITERM=  4
NIT=  32  NFV=  33  NFG=  32  F= 3.95116      C=0.9E-11  G=0.6E-06  ITERM=  4
NIT=  21  NFV=  21  NFG=  21  F= 7049.25      C=0.4E-07  G=0.1E-06  ITERM=  4
NIT=  16  NFV=  16  NFG=  16  F=-.866025      C=0.2E-08  G=0.2E-09  ITERM=  4
NIT=  22  NFV=  24  NFG=  22  F= 24.3062      C=0.5E-12  G=0.3E-06  ITERM=  4
NIT=  22  NFV=  22  NFG=  22  F= 32.3487      C=0.8E-12  G=0.1E-06  ITERM=  4
NIT=  77  NFV=  77  NFG=  77  F= 174.787      C=0.6E-11  G=0.4E-10  ITERM=  4
NIT=  32  NFV=  35  NFG=  32  F= 133.728      C=0.9E-13  G=0.8E-06  ITERM=  4
NITER =  734    NFVAL =  823    NSUCC =   34
TIME= 0:00:00.03
```

The rows corresponding to individual test problems contain the number of
iterations NIT, the number of function evaluations NFV, the number of
gradient evaluations NFG, the final value of the objective function F,
the constraint violation C, the norm of the Lagrangian function gradient
G and the cause of termination ITERM.
