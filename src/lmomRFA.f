C
C  Fortran code for the R package "lmomRFA"
C  Based on the LMOMENTS Fortran package, version 3.04.
C
C  The following routines are called from R functions:
C
C    QKAP
C    REGTST
C
C  The following routines are called from other Fortran routines.
C  Those marked (=) are the same as in the lmom package.
C
C    DERF    (=)
C    DIGAMD  (=)
C    DLGAMA  (=)
C    LMRGEV  (=)
C    LMRGLO  (=)
C    LMRGNO  (=)
C    LMRGPA  (=)
C    LMRPE3  (=)
C    PELGEV  (=)
C    PELGLO  (=)
C    PELGNO  (=)
C    PELGPA  (=)
C    PELKAP  (=)
C    PELPE3  (=)
C    PELWAK
C    QKAP
C    SAMLMU  (=)
C    SORT
C
C  Note that routine PELWAK in this package is not the same as in the
C  lmom R package.  The version in this package retains the behaviour
C  of PELWAK in the LMOMENTS Fortran package: if a fit using all five
C  parameters fails, an attempt is made to fit the distribution with
C  four parameters and lower bound zero.
C
C===================================================== REGTST.FOR
      SUBROUTINE REGTST(NSITES,LEN,XMOM,NSIM,RMOM,D,VOBS,VBAR,VSD,
     *  H,Z,PARA,RPARA,T4FIT,WORK,X,MAXREC)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmomRFA"                       *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Based on the routine of the same name in the LMOMENTS Fortran      *
C*  package, version 3.04.  See IBM Research Report RC20525, 'Fortran  *
C*  routines for use with the method of L-moments, version 3', and its *
C*  updates.                                                           *
C*                                                                     *
C*  Version 2.6    January 2014                                        *
C*                                                                     *
C***********************************************************************
C
C  Summary of differences from the Fortran version
C  * Omit arguments NAMES,NPROB,PROB,KPRINT,KOUT,A,B,SEED
C  * Add arguments RPARA,T4FIT,WORK,X,MAXREC
C  * No printing
C  * Return (as *output* arguments) all numbers that used to be printed
C  * No use of PARAMETERs
C  * Calls to PELxxx routines are to the versions adapted for R,
C    so gain an IFAIL parameter
C  * Compute sample L-moments via SAMLMU, not SAMLMR
C  * Don't stop when unable to invert SS matrix
C  * Fit the candidate (and Wakeby) distributions even when NSIM.EQ.0
C  * Compute T4FIT more accurately
C  * Don't compute quantiles of candidate distributions
C  * Use R's random-number generator, via C routines RANGET,CURAND,RANPUT
C
C  CALCULATES THREE STATISTICS USEFUL IN REGIONAL FREQUENCY ANALYSIS
C
C  DISCORDANCY MEASURE, D(I), FOR INDIVIDUAL SITES IN A REGION.
C      LARGE VALUES MIGHT BE USED AS A FLAG TO INDICATE POTENTIAL ERRORS
C      IN THE DATA AT THE SITE.  "LARGE" MIGHT BE 3 FOR REGIONS WITH 15
C      OR MORE SITES, BUT LESS (EXACT VALUES IN ARRAY DC1) FOR SMALLER
C      REGIONS.
C
C  HETEROGENEITY MEASURES, H(J), FOR A REGION BASED UPON EITHER:-
C      J=1: THE WEIGHTED S.D. OF THE L-CVS OR
C      J=2: THE AVERAGE DISTANCE FROM THE SITE TO THE REGIONAL AVERAGE
C           ON A GRAPH OF L-CV VS. L-SKEWNESS
C      J=3: THE AVERAGE DISTANCE FROM THE SITE TO THE REGIONAL AVERAGE
C           ON A GRAPH OF L-SKEWNESS VS. L-KURTOSIS
C
C      IN PRACTICE H(1) IS PROBABLY SUFFICIENT.  A VALUE GREATER THAN
C      (SAY) 1.0 SUGGESTS THAT FURTHER SUBDIVISION OF THE REGION SHOULD
C      BE CONSIDERED AS IT MIGHT IMPROVE QUANTILE ESTIMATES.
C
C  GOODNESS-OF-FIT MEASURES, Z(K), FOR 5 CANDIDATE DISTRIBUTIONS:
C      K=1: GENERALIZED LOGISTIC
C      K=2: GENERALIZED EXTREME VALUE
C      K=3: GENERALIZED NORMAL (LOGNORMAL)
C      K=4: PEARSON TYPE III (3-PARAMETER GAMMA)
C      K=5: GENERALIZED PARETO
C
C      PROVIDED THAT THE REGION IS ACCEPTABLY CLOSE TO HOMOGENEOUS,
C      THE FIT MAY BE JUDGED ACCEPTABLE AT 10% SIGNIFICANCE LEVEL
C      IF Z(K) IS LESS THAN 1.645 IN ABSOLUTE VALUE.
C
C  FOR FURTHER DETAILS SEE J.R.M. HOSKING AND J.R. WALLIS (1997),
C  "REGIONAL FREQUENCY ANALYSIS: AN APPROACH BASED ON L-MOMENTS",
C  CAMBRIDGE UNIVERSITY PRESS, CHAPTERS 3-5.
C
C  PARAMETERS OF ROUTINE:
C  NSITES * INPUT* NUMBER OF SITES IN REGION
C  LEN    * INPUT* ARRAY OF LENGTH NSITES. RECORD LENGTHS AT EACH SITE.
C  XMOM   * INPUT* ARRAY OF DIMENSION (5,NSITES). ARRAY CONTAINING
C                  THE FIRST 5 SAMPLE L-MOMENTS FOR EACH SITE, IN THE
C                  ORDER MEAN, L-CV, L-SKEWNESS, L-KURTOSIS, T-5, I.E
C                  XMOM(I,J) CONTAINS THE I'TH L-MOMENT FOR SITE J.
C                    N.B. XMOM(2,.) CONTAINS L-CV, NOT THE USUAL L-2!
C  NSIM   * INPUT* NUMBER OF SIMULATED WORLDS FOR HETEROGENEITY AND
C                  GOODNESS-OF-FIT TESTS.
C                    NOTE: NSIM=0 WILL FORCE RETURN AT COMPLETION OF
C                  OUTLIER TEST.  NSIM=1 WILL SUPPRESS CALCULATION OF
C                  H AND Z STATISTICS, BUT PARAMETER AND QUANTILE
C                  ESTIMATES WILL BE FOUND.
C  RMOM   *OUTPUT* ARRAY OF LENGTH 5. ON EXIT, CONTAINS THE REGIONAL
C                  WEIGHTED AVERAGE L-MOMENT RATIOS.
C  D      *OUTPUT* ARRAY OF LENGTH NSITES. ON EXIT, CONTAINS THE
C                  DISCORDANCY MEASURE (D STATISTIC) FOR EACH SITE.
C  VOBS   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE REGIONAL
C                  OBSERVED VALUES OF 3 HETEROGENEITY STATISTICS:
C                  (1) WEIGHTED S.D. OF L-CVS;
C                  (2) AVERAGE OF L-CV/L-SKEW DISTANCES;
C                  (3) AVERAGE OF L-SKEW/L-KURTOSIS DISTANCES.
C  VBAR   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE MEAN OF THE
C                  SIMULATED VALUES OF THE 3 HETEROGENEITY STATISTICS.
C  VSD    *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE S.D. OF THE
C                  SIMULATED VALUES OF THE 3 HETEROGENEITY STATISTICS.
C  H      *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS HETEROGENEITY
C                  MEASURES (H STATISTICS), I.E. H=(VOBS-VBAR)/VSD.
C  Z      *OUTPUT* ARRAY OF LENGTH 5. ON EXIT, CONTAINS GOODNESS-OF-FIT
C                  MEASURES (Z STATISTICS) FOR 5 DISTRIBUTIONS:
C                  (1) GEN. LOGISTIC, (2) GEN. EXTREME VALUE,
C                  (3) GEN. NORMAL, (4) PEARSON TYPE III,
C                  (5) GEN. PARETO.
C  PARA   *OUTPUT* ARRAY OF DIMENSION (5,6). ON EXIT, IF NSIM.GE.1,
C                  CONTAINS PARAMETERS OF GROWTH CURVES FITTED BY THE
C                  ABOVE 5 DISTRIBUTIONS, PLUS WAKEBY.
C  RPARA  *output* Array of length 4.  On exit, contains parameters of
C                  kappa distribution fitted to regional L-moments.
C  T4FIT  *output* Array of length 5.  On exit, contains the L-kurtosis
C                  values of the  5 fitted candidate distributions.
C  WORK   * local* Array of length 3*NSITES.  Used as workspace.
C  X      * local* Array of length MAXREC.  Used as workspace.
C  MAXREC * input* Must be at least as large as the largest of the
C                  record lengths.
C
C  OTHER FORTRAN ROUTINES USED: DERF,DIGAMD,DLGAMA,LMRGEV,LMRGLO,LMRGNO,
C    LMRGPA,LMRPE3,PELGEV,PELGLO,PELGNO,PELGPA,PELKAP,PELPE3,PELWAK,
C    QKAP,SAMLMU,SORT
C
C  C ROUTINES USED: CURAND,RANGET,RANPUT
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION D(NSITES),H(3),PARA(5,6),
     *  RMOM(5),RPARA(4),SMAT(3,3),TMOM(4),T4FIT(5),
     *  VBAR(3),VOBS(3),VSD(3),WORK(NSITES,3),X(MAXREC),XMOM(5,NSITES),
     *  Z(5)
      INTEGER LEN(NSITES)
      DATA ZERO/0D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/
      EXTERNAL DERF
C
C         INITIALIZE ARRAYS
C
      NMAX=0
      SUMLEN=0
      DO 10 I=1,NSITES
      NREC=LEN(I)
      IF(NREC.GT.NMAX)NMAX=NREC
      SUMLEN=SUMLEN+NREC
   10 D(I)=ZERO
      DO 20 K=1,3
      VOBS(K)=ZERO
      VBAR(K)=ZERO
      VSD(K)=ZERO
      H(K)=ZERO
   20 CONTINUE
      DO 30 IDIST=1,5
   30 Z(IDIST)=ZERO
      DO 40 IPARA=1,5
      DO 40 IDIST=1,6
   40 PARA(IPARA,IDIST)=ZERO
C
C         CALCULATE THE WEIGHTED MEAN OF L-CV, L-SKEW, L-KURTOSIS
C
      DO 60 K=2,5
      RMOM(K)=ZERO
      DO 50 I=1,NSITES
   50 RMOM(K)=RMOM(K)+LEN(I)*XMOM(K,I)
   60 RMOM(K)=RMOM(K)/SUMLEN
      RMOM(1)=ONE
C
C         CALCULATE SUM OF SQUARES MATRIX
C
      IF(NSITES.LE.3)GOTO 135
      SUM2=ZERO
      SUM3=ZERO
      SUM4=ZERO
      DO 70 I=1,NSITES
      SUM2=SUM2+XMOM(2,I)
      SUM3=SUM3+XMOM(3,I)
      SUM4=SUM4+XMOM(4,I)
   70 CONTINUE
      SUM2=SUM2/NSITES
      SUM3=SUM3/NSITES
      SUM4=SUM4/NSITES
      DO 80 I=1,NSITES
      WORK(I,1)=XMOM(2,I)-SUM2
      WORK(I,2)=XMOM(3,I)-SUM3
      WORK(I,3)=XMOM(4,I)-SUM4
   80 CONTINUE
      DO 100 J=1,3
      DO 100 K=J,3
      SMAT(J,K)=ZERO
      DO 90 I=1,NSITES
   90 SMAT(J,K)=SMAT(J,K)+WORK(I,J)*WORK(I,K)
  100 CONTINUE
C
C         INVERT SUM OF SQUARES MATRIX
C
      DO 110 K=1,3
      IF(SMAT(1,1).LE.ZERO)GOTO 140
      TEMP0=ONE/SMAT(1,1)
      TEMP1=-SMAT(1,2)*TEMP0
      TEMP2=-SMAT(1,3)*TEMP0
      IF(K.GT.2)TEMP1=-TEMP1
      IF(K.GT.1)TEMP2=-TEMP2
      SMAT(1,1)=SMAT(2,2)+TEMP1*SMAT(1,2)
      SMAT(1,2)=SMAT(2,3)+TEMP1*SMAT(1,3)
      SMAT(2,2)=SMAT(3,3)+TEMP2*SMAT(1,3)
      SMAT(1,3)=TEMP1
      SMAT(2,3)=TEMP2
      SMAT(3,3)=TEMP0
  110 CONTINUE
      SMAT(2,1)=SMAT(1,2)
      SMAT(3,1)=SMAT(1,3)
      SMAT(3,2)=SMAT(2,3)
C
C         CALCULATE DISCORDANCY MEASURES (D STATISTICS)
C
      FACTOR=NSITES/THREE
      DO 130 I=1,NSITES
      DO 120 J=1,3
      DO 120 K=1,3
  120 D(I)=D(I)+WORK(I,J)*WORK(I,K)*SMAT(J,K)
      D(I)=D(I)*FACTOR
      WORK(I,1)=D(I)
  130 CONTINUE
      CALL SORT(WORK(1,1),NSITES)
      GOTO 140
  135 DO 138 I=1,NSITES
  138 D(I)=ONE
  140 CONTINUE
C
C         FIT DISTRIBUTIONS
C
      CALL PELGLO(RMOM,PARA(1,1),IFAIL)
      CALL PELGEV(RMOM,PARA(1,2),IFAIL)
      CALL PELGNO(RMOM,PARA(1,3),IFAIL)
      CALL PELPE3(RMOM,PARA(1,4),IFAIL)
      CALL PELGPA(RMOM,PARA(1,5),IFAIL)
      CALL PELWAK(RMOM,PARA(1,6),IFAIL)
C
      IF(NSIM.LE.1)RETURN
C
C         FIT KAPPA DISTRIBUTION TO REGIONAL L-MOMENTS
C
      CALL PELKAP(RMOM,RPARA,IFAIL)
      IF(IFAIL.EQ.0)GOTO 180
      CALL PELGLO(RMOM,RPARA,IFAIL)
      RPARA(4)=-ONE
  180 CONTINUE
C
C         INITIALIZE R's RANDOM NUMBER GENERATOR
C
      CALL RANGET()
C
C         START THE NSIM REPETITIONS
C
      T4BAR=ZERO
      T4SD=ZERO
      DO 220 ISIM=1,NSIM
      SUM2=ZERO
      SUM3=ZERO
      SUM4=ZERO
C
C         Check for user interrupt from R
C
      CALL RCHKUSR()
C
C         START OF LOOP OVER SITES
C
      DO 200 I=1,NSITES
      NREC=LEN(I)
C
C         GET VECTOR OF UNIFORM RANDOM NUMBERS
C
      CALL CURAND(NREC,X)
C
C         TRANSFORM FROM UNIFORM TO KAPPA
C
      CALL QKAP(X,NREC,RPARA)
C
C         FIND L-MOMENTS OF SIMULATED DATA
C
      CALL SORT(X,NREC)
      CALL SAMLMU(X,NREC,TMOM,4,IFAIL)
      CV=TMOM(2)/TMOM(1)
      WORK(I,1)=CV
      WORK(I,2)=TMOM(3)
      WORK(I,3)=TMOM(4)
      SUM2=SUM2+NREC*CV
      SUM3=SUM3+NREC*TMOM(3)
      SUM4=SUM4+NREC*TMOM(4)
C
C         END OF LOOP OVER SITES
C
  200 CONTINUE
C
      SUM2=SUM2/SUMLEN
      SUM3=SUM3/SUMLEN
      SUM4=SUM4/SUMLEN
      T4BAR=T4BAR+SUM4
      T4SD=T4SD+SUM4**2
C
C         CALCULATE HETEROGENEITY V-STATISTICS FOR SIMULATED DATA
C
      IF(NSITES.EQ.1)GOTO 215
      V1=ZERO
      V2=ZERO
      V3=ZERO
      DO 210 I=1,NSITES
      NREC=LEN(I)
      TEMP2=(WORK(I,1)-SUM2)**2
      TEMP3=(WORK(I,2)-SUM3)**2
      TEMP4=(WORK(I,3)-SUM4)**2
      V1=V1+NREC*TEMP2
      V2=V2+NREC*DSQRT(TEMP2+TEMP3)
      V3=V3+NREC*DSQRT(TEMP3+TEMP4)
  210 CONTINUE
      V1=DSQRT(V1/SUMLEN)
      V2=V2/SUMLEN
      V3=V3/SUMLEN
      VBAR(1)=VBAR(1)+V1
      VBAR(2)=VBAR(2)+V2
      VBAR(3)=VBAR(3)+V3
      VSD(1)=VSD(1)+V1**2
      VSD(2)=VSD(2)+V2**2
      VSD(3)=VSD(3)+V3**2
  215 CONTINUE
C
C         END OF THE NSIM REPETITIONS
C
  220 CONTINUE
C
C         FINALIZE R'S RANDOM NUMBER GENERATOR (WON'T NEED IT AGAIN)
C
      CALL RANPUT()
C
C         CALCULATE HETEROGENEITY V-STATISTICS FOR OBSERVED DATA
C
      IF(NSITES.EQ.1)GOTO 235
      V1=ZERO
      V2=ZERO
      V3=ZERO
      DO 225 I=1,NSITES
      NREC=LEN(I)
      TEMP2=(XMOM(2,I)-RMOM(2))**2
      TEMP3=(XMOM(3,I)-RMOM(3))**2
      TEMP4=(XMOM(4,I)-RMOM(4))**2
      V1=V1+NREC*TEMP2
      V2=V2+NREC*DSQRT(TEMP2+TEMP3)
      V3=V3+NREC*DSQRT(TEMP3+TEMP4)
  225 CONTINUE
      VOBS(1)=DSQRT(V1/SUMLEN)
      VOBS(2)=V2/SUMLEN
      VOBS(3)=V3/SUMLEN
C
C         CALCULATE HETEROGENEITY MEASURES (H STATISTICS)
C
      DO 230 J=1,3
      VBAR(J)=VBAR(J)/NSIM
      VSD(J)=DSQRT((VSD(J)-NSIM*VBAR(J)**2)/(NSIM-ONE))
      H(J)=(VOBS(J)-VBAR(J))/VSD(J)
  230 CONTINUE
  235 CONTINUE
C
C         FIND TAU-4 VALUES OF EACH CANDIDATE DISTRIBUTION
C
      CALL LMRGLO(PARA(1,1),TMOM,4,IFAIL)
      T4FIT(1)=TMOM(4)
      CALL LMRGEV(PARA(1,2),TMOM,4,IFAIL)
      T4FIT(2)=TMOM(4)
      CALL LMRGNO(PARA(1,3),TMOM,4,IFAIL)
      T4FIT(3)=TMOM(4)
      CALL LMRPE3(PARA(1,4),TMOM,4,IFAIL)
      T4FIT(4)=TMOM(4)
      CALL LMRGPA(PARA(1,5),TMOM,4,IFAIL)
      T4FIT(5)=TMOM(4)
C
C         CALCULATE GOODNESS-OF-FIT MEASURES (Z STATISTICS)
C
      T4BAR=T4BAR/NSIM
      T4SD=DSQRT((T4SD-NSIM*T4BAR**2)/(NSIM-ONE))
      DO 240 IDIST=1,5
      Z(IDIST)=(T4FIT(IDIST)+T4BAR-TWO*RMOM(4))/T4SD
  240 CONTINUE
C
      RETURN
C
      END
C===================================================== LMRGEV.FOR
      SUBROUTINE LMRGEV(PARA,XMOM,NMOM,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    July 2008                                           *
C*                                                                     *
C***********************************************************************
C
C  L-MOMENT RATIOS FOR THE GENERALIZED EXTREME-VALUE DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  PARA   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE PARAMETERS OF THE
C                  DISTRIBUTION, IN THE ORDER XI, ALPHA, K (LOCATION,
C                  SCALE, SHAPE).
C  XMOM   *OUTPUT* ARRAY OF LENGTH NMOM. ON EXIT, CONTAINS THE L-MOMENTS
C                  LAMBDA-1, LAMBDA-2, TAU-3, TAU-4, ... .
C  NMOM   * INPUT* NUMBER OF L-MOMENTS TO BE FOUND. AT MOST 20.
C
C  OTHER ROUTINES USED: DLGAMA
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(3),XMOM(NMOM),ZMOM(20)
      DATA ZERO/0D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/,FOUR/4D0/,SIX/6D0/
C
C         ARRAY ZMOM CONTAINS THE L-MOMENT RATIOS OF THE STANDARD
C         GUMBEL DISTRIBUTION (XI=0, ALPHA=1).
C         ZMOM(1) IS EULER'S CONSTANT, ZMOM(2) IS LOG(2).
C
      DATA ZMOM/
     *  0.57721 56649 01532 861D 0,  0.69314 71805 59945 309D 0,
     *  0.16992 50014 42312 363D 0,  0.15037 49927 88438 185D 0,
     *  0.55868 35005 77583 138D-1,  0.58110 02399 99710 876D-1,
     *  0.27624 25842 97309 125D-1,  0.30556 37665 79053 126D-1,
     *  0.16465 02822 58328 802D-1,  0.18784 66242 98170 912D-1,
     *  0.10932 82150 63027 148D-1,  0.12697 31266 76329 530D-1,
     *  0.77898 28180 57231 804D-2,  0.91483 61796 21999 726D-2,
     *  0.58333 23893 28363 588D-2,  0.69010 42875 90348 154D-2,
     *  0.45326 79701 80679 549D-2,  0.53891 68113 26595 459D-2,
     *  0.36240 77677 72368 790D-2,  0.43238 76086 05538 096D-2/
C
C         SMALL IS USED TO TEST WHETHER K IS EFFECTIVELY ZERO
C
      DATA SMALL/1D-6/
C
      IFAIL=0
      U=PARA(1)
      A=PARA(2)
      G=PARA(3)
      IF(A.LE.ZERO.OR.G.LE.-ONE)GOTO 1000
      IF(NMOM.GT.20)GOTO 1010
C
C         TEST FOR K=0
C
      IF(DABS(G).GT.SMALL)GOTO 20
      XMOM(1)=U+A*ZMOM(1)
      IF(NMOM.EQ.1)RETURN
      XMOM(2)=A*ZMOM(2)
      IF(NMOM.EQ.2)RETURN
      DO 10 I=3,NMOM
   10 XMOM(I)=ZMOM(I)
      RETURN
   20 CONTINUE
C
C         FIRST 2 MOMENTS
C
      GAM=DEXP(DLGAMA(ONE+G))
      XMOM(1)=U+A*(ONE-GAM)/G
      IF(NMOM.EQ.1)RETURN
      XX2=ONE-TWO**(-G)
      XMOM(2)=A*XX2*GAM/G
      IF(NMOM.EQ.2)RETURN
C
C         HIGHER MOMENTS
C
      Z0=ONE
      DO 50 J=3,NMOM
      DJ=J
      BETA=(ONE-DJ**(-G))/XX2
      Z0=Z0*(FOUR*DJ-SIX)/DJ
      Z=Z0*THREE*(DJ-ONE)/(DJ+ONE)
      SUM=Z0*BETA-Z
      IF(J.EQ.3)GOTO 40
      DO 30 I=2,J-2
      DI=I
      Z=Z*(DI+DI+ONE)*(DJ-DI)/((DI+DI-ONE)*(DJ+DI))
      SUM=SUM-Z*XMOM(I+1)
   30 CONTINUE
   40 XMOM(J)=SUM
   50 CONTINUE
      RETURN
C
 1000 IFAIL=7000
      RETURN
 1010 IFAIL=7010
      RETURN
C
C7000 FORMAT(' *** ERROR *** ROUTINE LMRGEV : PARAMETERS INVALID')
C7010 FORMAT(' *** ERROR *** ROUTINE LMRGEV : PARAMETER NMOM TOO LARGE')
      END
C===================================================== LMRGLO.FOR
      SUBROUTINE LMRGLO(PARA,XMOM,NMOM,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    July 2008                                           *
C*                                                                     *
C***********************************************************************
C
C  L-MOMENT RATIOS FOR THE GENERALIZED LOGISTIC DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  PARA   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE PARAMETERS OF THE
C                  DISTRIBUTION, IN THE ORDER XI, ALPHA, K (LOCATION,
C                  SCALE, SHAPE).
C  XMOM   *OUTPUT* ARRAY OF LENGTH NMOM. ON EXIT, CONTAINS THE L-MOMENTS
C                  LAMBDA-1, LAMBDA-2, TAU-3, TAU-4, ... .
C  NMOM   * INPUT* NUMBER OF L-MOMENTS TO BE FOUND. AT MOST 20.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(3),XMOM(NMOM),Z(10,20)
      DATA ZERO/0D0/,ONE/1D0/
      DATA PI/3.141592653589793238D0/
C
C         SMALL IS USED TO DECIDE WHETHER TO APPROXIMATE THE FIRST 2
C         L-MOMENTS BY A POWER-SERIES EXPANSION WHEN G IS NEAR ZERO.
C         C1,C2 ARE COEFFICIENTS OF THIS POWER-SERIES EXPANSION.
C         C1 IS PI**2/6, C2 IS 7*PI**4/360.
C
      DATA SMALL/1D-4/
      DATA C1,C2/
     *  0.16449 34066 84822 644D 1,  0.18940 65658 99449 184D 1/
C
C         Z-ARRAY CONTAINS COEFFICIENTS OF THE REPRESENTATIONS OF
C         L-MOMENT RATIOS AS POLYNOMIALS IN THE SHAPE PARAMETER K
C
      DATA Z(1,3)/1D0/
      DATA (Z(I, 4),I=1, 2)/
     *  0.16666 66666 66666 667D 0,  0.83333 33333 33333 333D 0/
      DATA (Z(I, 5),I=1, 2)/
     *  0.41666 66666 66666 667D 0,  0.58333 33333 33333 333D 0/
      DATA (Z(I, 6),I=1, 3)/
     *  0.66666 66666 66666 667D-1,  0.58333 33333 33333 333D 0,
     *  0.35000 00000 00000 000D 0/
      DATA (Z(I, 7),I=1, 3)/
     *  0.23333 33333 33333 333D 0,  0.58333 33333 33333 333D 0,
     *  0.18333 33333 33333 333D 0/
      DATA (Z(I, 8),I=1, 4)/
     *  0.35714 28571 42857 143D-1,  0.42083 33333 33333 333D 0,
     *  0.45833 33333 33333 333D 0,  0.85119 04761 90476 190D-1/
      DATA (Z(I, 9),I=1, 4)/
     *  0.15099 20634 92063 492D 0,  0.51562 50000 00000 000D 0,
     *  0.29791 66666 66666 667D 0,  0.35466 26984 12698 413D-1/
      DATA (Z(I,10),I=1, 5)/
     *  0.22222 22222 22222 222D-1,  0.31889 32980 59964 727D 0,
     *  0.47997 68518 51851 852D 0,  0.16550 92592 59259 259D 0,
     *  0.13398 36860 67019 400D-1/
      DATA (Z(I,11),I=1, 5)/
     *  0.10650 79365 07936 508D 0,  0.44766 31393 29805 996D 0,
     *  0.36081 01851 85185 185D 0,  0.80390 21164 02116 402D-1,
     *  0.46285 27336 86067 019D-2/
      DATA (Z(I,12),I=1, 6)/
     *  0.15151 51515 15151 515D-1,  0.25131 61375 66137 566D 0,
     *  0.46969 52160 49382 716D 0,  0.22765 04629 62962 963D 0,
     *  0.34713 95502 64550 265D-1,  0.14727 13243 54657 688D-2/
      DATA (Z(I,13),I=1, 6)/
     *  0.79569 50456 95045 695D-1,  0.38976 59465 02057 613D 0,
     *  0.39291 73096 70781 893D 0,  0.12381 31062 61022 928D 0,
     *  0.13499 87139 91769 547D-1,  0.43426 15974 56041 900D-3/
      DATA (Z(I,14),I=1, 7)/
     *  0.10989 01098 90109 890D-1,  0.20413 29966 32996 633D 0,
     *  0.44773 66255 14403 292D 0,  0.27305 34428 27748 383D 0,
     *  0.59191 74382 71604 938D-1,  0.47768 77572 01646 091D-2,
     *  0.11930 26366 63747 775D-3/
      DATA (Z(I,15),I=1, 7)/
     *  0.61934 52050 59490 774D-1,  0.34203 17593 92870 504D 0,
     *  0.40701 37051 73427 396D 0,  0.16218 91928 06752 331D 0,
     *  0.25249 21002 35155 791D-1,  0.15509 34276 62872 107D-2,
     *  0.30677 82085 63922 850D-4/
      DATA (Z(I,16),I=1, 8)/
     *  0.83333 33333 33333 333D-2,  0.16976 83649 02293 474D 0,
     *  0.42219 12828 68366 202D 0,  0.30542 71728 94620 811D 0,
     *  0.84082 79399 72285 210D-1,  0.97243 57914 46208 113D-2,
     *  0.46528 02829 88616 322D-3,  0.74138 06706 96146 887D-5/
      DATA (Z(I,17),I=1, 8)/
     *  0.49716 60284 16028 416D-1,  0.30276 58385 89871 328D 0,
     *  0.41047 33000 89185 506D 0,  0.19483 90265 03251 764D 0,
     *  0.38659 80637 04648 526D-1,  0.34139 94076 42897 226D-2,
     *  0.12974 16173 71825 705D-3,  0.16899 11822 91033 482D-5/
      DATA (Z(I,18),I=1, 9)/
     *  0.65359 47712 41830 065D-2,  0.14387 48475 95085 690D 0,
     *  0.39643 28537 10259 464D 0,  0.32808 41807 20899 471D 0,
     *  0.10797 13931 65194 318D 0,  0.15965 33699 32077 769D-1,
     *  0.11012 77375 69143 819D-2,  0.33798 23645 82066 963D-4,
     *  0.36449 07853 33601 627D-6/
      DATA (Z(I,19),I=1, 9)/
     *  0.40878 45705 49276 431D-1,  0.27024 42907 25441 519D 0,
     *  0.40759 95245 14551 521D 0,  0.22211 14264 89320 008D 0,
     *  0.52846 38846 29533 398D-1,  0.59829 82392 72872 761D-2,
     *  0.32859 39655 65898 436D-3,  0.82617 91134 22830 354D-5,
     *  0.74603 37711 50646 605D-7/
      DATA (Z(I,20),I=1,10)/
     *  0.52631 57894 73684 211D-2,  0.12381 76557 53054 913D 0,
     *  0.37185 92914 44794 917D 0,  0.34356 87476 70189 607D 0,
     *  0.13019 86628 12524 058D 0,  0.23147 43648 99477 023D-1,
     *  0.20519 25194 79869 981D-2,  0.91205 82581 07571 930D-4,
     *  0.19023 86116 43414 884D-5,  0.14528 02606 97757 497D-7/
C
      IFAIL=0
      U=PARA(1)
      A=PARA(2)
      G=PARA(3)
      IF(A.LE.ZERO.OR.DABS(G).GE.ONE)GOTO 1000
      IF(NMOM.GT.20)GOTO 1010
C
C         FIRST 2 MOMENTS
C
      GG=G*G
      ALAM1=-G*(C1+GG*C2)
      ALAM2=ONE+GG*(C1+GG*C2)
      IF(DABS(G).GT.SMALL)ALAM2=G*PI/DSIN(G*PI)
      IF(DABS(G).GT.SMALL)ALAM1=(ONE-ALAM2)/G
      XMOM(1)=U+A*ALAM1
      IF(NMOM.EQ.1)RETURN
      XMOM(2)=A*ALAM2
      IF(NMOM.EQ.2)RETURN
C
C         HIGHER MOMENTS
C
      DO 20 M=3,NMOM
      KMAX=M/2
      SUM=Z(KMAX,M)
      DO 10 K=KMAX-1,1,-1
   10 SUM=SUM*GG+Z(K,M)
      IF(M.NE.M/2*2)SUM=-G*SUM
      XMOM(M)=SUM
   20 CONTINUE
      RETURN
C
 1000 IFAIL=7000
      RETURN
 1010 IFAIL=7010
      RETURN
C
C7000 FORMAT(' *** ERROR *** ROUTINE LMRGLO : PARAMETERS INVALID')
C7010 FORMAT(' *** ERROR *** ROUTINE LMRGLO : PARAMETER NMOM TOO LARGE')
      END
C===================================================== LMRGNO.FOR
      SUBROUTINE LMRGNO(PARA,XMOM,NMOM,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    July 2008                                           *
C*                                                                     *
C***********************************************************************
C
C  L-MOMENT RATIOS FOR THE GENERALIZED NORMAL DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  PARA   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE PARAMETERS OF THE
C                  DISTRIBUTION, IN THE ORDER XI, ALPHA, K (LOCATION,
C                  SCALE, SHAPE).
C  XMOM   *OUTPUT* ARRAY OF LENGTH NMOM. ON EXIT, CONTAINS THE L-MOMENTS
C                  LAMBDA-1, LAMBDA-2, TAU-3, TAU-4, ... .
C  NMOM   * INPUT* NUMBER OF L-MOMENTS TO BE FOUND. AT MOST 20.
C
C  OTHER ROUTINES USED: DERF
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(3),XMOM(NMOM),EST(20),ESTX(20),SUM(20),
     *  ZMOM(20)
      EXTERNAL DERF
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/
C
C         ARRAY ZMOM CONTAINS L-MOMENTS OF THE STANDARD NORMAL DIST.
C
      DATA ZMOM/
     *  0D0,   0.56418 95835 47756 287D 0,
     *  0D0,   0.12260 17195 40890 947D 0,
     *  0D0,   0.43661 15389 50024 944D-1,
     *  0D0,   0.21843 13603 32508 776D-1,
     *  0D0,   0.12963 50158 01507 746D-1,
     *  0D0,   0.85296 21241 91705 402D-2,
     *  0D0,   0.60138 90151 79323 333D-2,
     *  0D0,   0.44555 82586 47650 150D-2,
     *  0D0,   0.34264 32435 78076 985D-2,
     *  0D0,   0.27126 79630 48139 365D-2/
C
C         RRT2 IS 1/SQRT(2), RRTPI IS 1/SQRT(PI)
C
      DATA RRT2 /0.70710 67811 86547 524D0/
      DATA RRTPI/0.56418 95835 47756 287D0/
C
C         RANGE,EPS,MAXIT CONTROL THE ITERATIVE PROCEDURE FOR NUMERICAL
C         INTEGRATION
C
      DATA RANGE/5D0/,EPS/1D-8/,MAXIT/10/
C
C----
C Next lines to assuage compilers that warn that variables
C 'might be used uninitialized'
      NOTCGD=0
C----
      IFAIL=0
      U=PARA(1)
      A=PARA(2)
      G=PARA(3)
      IF(A.LE.ZERO)GOTO 1000
      IF(NMOM.GT.20)GOTO 1010
C
C         TEST FOR K=0
C
      IF(DABS(G).GT.EPS)GOTO 5
      XMOM(1)=U
      IF(NMOM.EQ.1)RETURN
      XMOM(2)=A*ZMOM(2)
      IF(NMOM.EQ.2)RETURN
      DO 2 I=3,NMOM
    2 XMOM(I)=ZMOM(I)
      RETURN
    5 CONTINUE
C
C         LAMBDA-1
C
      EGG=DEXP(HALF*G*G)
      ALAM1=(ONE-EGG)/G
      XMOM(1)=U+A*ALAM1
      IF(NMOM.EQ.1)RETURN
C
C         LAMBDA-2
C
      ALAM2=EGG*DERF(HALF*G)/G
      XMOM(2)=A*ALAM2
      IF(NMOM.EQ.2)RETURN
C
C         HIGHER MOMENTS. THE INTEGRAL DEFINING LAMBDA-R IS EVALUATED
C         BY ITERATIVE APPLICATION OF THE TRAPEZIUM RULE.
C
C         - INITIAL ESTIMATE, USING 16 ORDINATES  (THE 'DO 20' LOOP
C           CALCULATES LEGENDRE POLYNOMIALS RECURSIVELY)
C
      CC=-G*RRT2
      XMIN=CC-RANGE
      XMAX=CC+RANGE
      DO 10 M=3,NMOM
   10 SUM(M)=ZERO
      N=16
      XINC=(XMAX-XMIN)/N
      DO 30 I=1,N-1
      X=XMIN+I*XINC
      E=DEXP(-((X-CC)**2))
      D=DERF(X)
      P1=ONE
      P=D
      DO 20 M=3,NMOM
      C1=M+M-3
      C2=M-2
      C3=M-1
      P2=P1
      P1=P
      P=(C1*D*P1-C2*P2)/C3
   20 SUM(M)=SUM(M)+E*P
   30 CONTINUE
      DO 40 M=3,NMOM
   40 EST(M)=SUM(M)*XINC
C
C         - DOUBLE THE NUMBER OF ORDINATES UNTIL CONVERGED
C
      DO 90 IT=1,MAXIT
      DO 50 M=3,NMOM
   50 ESTX(M)=EST(M)
      N=N*2
      XINC=(XMAX-XMIN)/N
      DO 70 I=1,N-1,2
      X=XMIN+I*XINC
      E=DEXP(-((X-CC)**2))
      D=DERF(X)
      P1=ONE
      P=D
      DO 60 M=3,NMOM
      C1=M+M-3
      C2=M-2
      C3=M-1
      P2=P1
      P1=P
      P=(C1*D*P1-C2*P2)/C3
   60 SUM(M)=SUM(M)+E*P
   70 CONTINUE
C
C         --- TEST FOR CONVERGENCE
C
      NOTCGD=0
      DO 80 M=NMOM,3,-1
      EST(M)=SUM(M)*XINC
      IF(DABS(EST(M)-ESTX(M)).GT.EPS*DABS(EST(M)))NOTCGD=M
   80 CONTINUE
      IF(NOTCGD.EQ.0)GOTO 100
   90 CONTINUE
C
      IFAIL=7100+(NOTCGD-1)
  100 CONTINUE
      CONST=-DEXP(CC*CC)*RRTPI/(ALAM2*G)
      DO 110 M=3,NMOM
  110 XMOM(M)=CONST*EST(M)
      RETURN
C
 1000 IFAIL=7000
      RETURN
 1010 IFAIL=7010
      RETURN
C
C7000 FORMAT(' *** ERROR *** ROUTINE LMRGNO : PARAMETERS INVALID')
C7010 FORMAT(' *** ERROR *** ROUTINE LMRGNO : PARAMETER NMOM TOO LARGE')
C7100 FORMAT(' ** WARNING ** ROUTINE LMRGNO :',
C    *  ' ITERATION HAS NOT CONVERGED. ONLY THE FIRST',I3,
C    *  ' L-MOMENTS ARE RELIABLE.')
      END
C===================================================== LMRGPA.FOR
      SUBROUTINE LMRGPA(PARA,XMOM,NMOM,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    July 2008                                           *
C*                                                                     *
C***********************************************************************
C
C  L-MOMENT RATIOS FOR THE GENERALIZED PARETO DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  PARA   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE PARAMETERS OF THE
C                  DISTRIBUTION, IN THE ORDER XI, ALPHA, K (LOCATION,
C                  SCALE, SHAPE).
C  XMOM   *OUTPUT* ARRAY OF LENGTH NMOM. ON EXIT, CONTAINS THE L-MOMENTS
C                  LAMBDA-1, LAMBDA-2, TAU-3, TAU-4, ... .
C  NMOM   * INPUT* NUMBER OF L-MOMENTS TO BE FOUND. AT MOST 20.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(5),XMOM(NMOM)
      DATA ZERO/0D0/,ONE/1D0/,TWO/2D0/
C
      IFAIL=0
      U=PARA(1)
      A=PARA(2)
      G=PARA(3)
      IF(A.LE.ZERO.OR.G.LT.-ONE)GOTO 1000
      IF(NMOM.GT.20)GOTO 1010
C
C         LAMBDA-1
C
      Y=ONE/(ONE+G)
      XMOM(1)=U+A*Y
      IF(NMOM.EQ.1)RETURN
C
C         LAMBDA-2
C
      Y=Y/(TWO+G)
      XMOM(2)=A*Y
      IF(NMOM.EQ.2)RETURN
C
C         HIGHER MOMENTS
C
      Y=ONE
      DO 10 M=3,NMOM
      AM=M-TWO
      Y=Y*(AM-G)/(M+G)
      XMOM(M)=Y
   10 CONTINUE
      RETURN
C
 1000 IFAIL=7000
      RETURN
 1010 IFAIL=7010
      RETURN
C
C7000 FORMAT(' *** ERROR *** ROUTINE LMRGPA : PARAMETERS INVALID')
C7010 FORMAT(' *** ERROR *** ROUTINE LMRGPA : PARAMETER NMOM TOO LARGE')
      END
C===================================================== LMRPE3.FOR
      SUBROUTINE LMRPE3(PARA,XMOM,NMOM,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    July 2008                                           *
C*                                                                     *
C***********************************************************************
C
C  L-MOMENT RATIOS FOR THE PEARSON TYPE 3 DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  PARA   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE PARAMETERS OF THE
C                  DISTRIBUTION, IN THE ORDER MU, SIGMA, GAMMA (MEAN,
C                  S.D., SKEWNESS).
C  XMOM   *OUTPUT* ARRAY OF LENGTH NMOM. ON EXIT, CONTAINS UP TO 4 OF
C                  THE L-MOMENTS LAMBDA-1, LAMBDA-2, TAU-3, TAU-4.
C  NMOM   * INPUT* NUMBER OF L-MOMENTS TO BE FOUND. AT MOST 4.
C
C  OTHER ROUTINES USED: DLGAMA
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(3),XMOM(NMOM)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,FOUR/4D0/
C
C         SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
C
      DATA SMALL/1D-6/
C
C         CONST IS 1/SQRT(PI)
C
      DATA CONST/0.56418 95835 47756 287D0/
C
C         COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATIONS
C         A0 IS 1/SQRT(3*PI)
C         C0 IS TAU-4 FOR THE NORMAL DISTRIBUTION
C
      DATA A0      / 0.32573501D+00/
      DATA A1,A2,A3/ 0.16869150D+00, 0.78327243D-01,-0.29120539D-02/
      DATA B1,B2   / 0.46697102D+00, 0.24255406D+00/
      DATA C0      / 0.12260172D+00/
      DATA C1,C2,C3/ 0.53730130D-01, 0.43384378D-01, 0.11101277D-01/
      DATA D1,D2   / 0.18324466D+00, 0.20166036D+00/
      DATA E1,E2,E3/ 0.23807576D+01, 0.15931792D+01, 0.11618371D+00/
      DATA F1,F2,F3/ 0.51533299D+01, 0.71425260D+01, 0.19745056D+01/
      DATA G1,G2,G3/ 0.21235833D+01, 0.41670213D+01, 0.31925299D+01/
      DATA H1,H2,H3/ 0.90551443D+01, 0.26649995D+02, 0.26193668D+02/
C
      IFAIL=0
      SD=PARA(2)
      IF(SD.LE.ZERO)GOTO 1000
      IF(NMOM.GT.4)GOTO 1010
C
C         LAMBDA-1
C
      XMOM(1)=PARA(1)
      IF(NMOM.EQ.1)RETURN
C
C         LAMBDA-2
C
      GAMMA=PARA(3)
      IF(DABS(GAMMA).LT.SMALL)GOTO 20
      ALPHA=FOUR/(GAMMA*GAMMA)
      BETA=DABS(HALF*SD*GAMMA)
      ALAM2=CONST*DEXP(DLGAMA(ALPHA+HALF)-DLGAMA(ALPHA))
      XMOM(2)=ALAM2*BETA
      IF(NMOM.EQ.2)RETURN
C
C         HIGHER MOMENTS
C
      IF(ALPHA.LT.ONE)GOTO 10
      Z=ONE/ALPHA
      XMOM(3)=DSQRT(Z)*(((A3*Z+A2)*Z+A1)*Z+A0)/((B2*Z+B1)*Z+ONE)
      IF(GAMMA.LT.ZERO)XMOM(3)=-XMOM(3)
      IF(NMOM.EQ.3)RETURN
      XMOM(4)=(((C3*Z+C2)*Z+C1)*Z+C0)/((D2*Z+D1)*Z+ONE)
      RETURN
C
   10 Z=ALPHA
      XMOM(3)=(((E3*Z+E2)*Z+E1)*Z+ONE)/(((F3*Z+F2)*Z+F1)*Z+ONE)
      IF(GAMMA.LT.ZERO)XMOM(3)=-XMOM(3)
      IF(NMOM.EQ.3)RETURN
      XMOM(4)=(((G3*Z+G2)*Z+G1)*Z+ONE)/(((H3*Z+H2)*Z+H1)*Z+ONE)
      IF(NMOM.GT.4)IFAIL=7010
      RETURN
C
C         CASE OF ZERO SKEWNESS
C
   20 XMOM(1)=PARA(1)
      IF(NMOM.EQ.1)RETURN
      XMOM(2)=CONST*PARA(2)
      IF(NMOM.EQ.2)RETURN
      XMOM(3)=0
      IF(NMOM.EQ.3)RETURN
      XMOM(4)=C0
      IF(NMOM.GT.4)IFAIL=7010
      RETURN
C
 1000 IFAIL=7000
      RETURN
 1010 IFAIL=7010
      RETURN
C
C7000 FORMAT(' *** ERROR *** ROUTINE LMRPE3 : PARAMETERS INVALID')
C7010 FORMAT(' *** ERROR *** ROUTINE LMRPE3 : PARAMETER NMOM TOO LARGE')
      END
C===================================================== PELGEV.FOR
      SUBROUTINE PELGEV(XMOM,PARA,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    July 2008                                           *
C*                                                                     *
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE GENERALIZED EXTREME-VALUE
C  DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2, TAU-3.
C  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA, K (LOCATION, SCALE, SHAPE).
C
C  OTHER ROUTINES USED: DLGAMA
C
C  METHOD: FOR  -0.8 LE TAU3 LT 1,  K IS APPROXIMATED BY RATIONAL
C  FUNCTIONS AS IN DONALDSON (1996, COMMUN. STATIST. SIMUL. COMPUT.).
C  IF TAU3 IS OUTSIDE THIS RANGE, NEWTON-RAPHSON ITERATION IS USED.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(3),PARA(3)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/
      DATA P8/0.8D0/,P97/0.97D0/
C
C         SMALL IS USED TO TEST WHETHER K IS EFFECTIVELY ZERO
C         EPS,MAXIT CONTROL THE TEST FOR CONVERGENCE OF N-R ITERATION
C
      DATA SMALL/1D-5/,EPS/1D-6/,MAXIT/20/
C
C         EU IS EULER'S CONSTANT
C         DL2 IS LOG(2), DL3 IS LOG(3)
C
      DATA EU/0.57721566D0/,DL2/0.69314718D0/,DL3/1.0986123D0/
C
C         COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATIONS FOR K
C
      DATA A0,A1,A2/ 0.28377530D0,-1.21096399D0,-2.50728214D0/
      DATA A3,A4   /-1.13455566D0,-0.07138022D0/
      DATA B1,B2,B3/ 2.06189696D0, 1.31912239D0, 0.25077104D0/
      DATA C1,C2,C3/ 1.59921491D0,-0.48832213D0, 0.01573152D0/
      DATA D1,D2   /-0.64363929D0, 0.08985247D0/
C
      IFAIL=0
      T3=XMOM(3)
      IF(XMOM(2).LE.ZERO)GOTO 1000
      IF(DABS(T3).GE.ONE)GOTO 1000
      IF(T3.LE.ZERO)GOTO 10
C
C         RATIONAL-FUNCTION APPROXIMATION FOR TAU3 BETWEEN 0 AND 1
C
      Z=ONE-T3
      G=(-ONE+Z*(C1+Z*(C2+Z*C3)))/(ONE+Z*(D1+Z*D2))
      IF(DABS(G).LT.SMALL)GOTO 50
      GOTO 40
C
C         RATIONAL-FUNCTION APPROXIMATION FOR TAU3 BETWEEN -0.8 AND 0
C
   10 G=(A0+T3*(A1+T3*(A2+T3*(A3+T3*A4))))/(ONE+T3*(B1+T3*(B2+T3*B3)))
      IF(T3.GE.-P8)GOTO 40
C
C         NEWTON-RAPHSON ITERATION FOR TAU3 LESS THAN -0.8
C
      IF(T3.LE.-P97)G=ONE-DLOG(ONE+T3)/DL2
      T0=(T3+THREE)*HALF
      DO 20 IT=1,MAXIT
      X2=TWO**(-G)
      X3=THREE**(-G)
      XX2=ONE-X2
      XX3=ONE-X3
      T=XX3/XX2
      DERIV=(XX2*X3*DL3-XX3*X2*DL2)/(XX2*XX2)
      GOLD=G
      G=G-(T-T0)/DERIV
      IF(DABS(G-GOLD).LE.EPS*G)GOTO 30
   20 CONTINUE
      IFAIL=7020
   30 CONTINUE
C
C         ESTIMATE ALPHA,XI
C
   40 PARA(3)=G
      GAM=DEXP(DLGAMA(ONE+G))
      PARA(2)=XMOM(2)*G/(GAM*(ONE-TWO**(-G)))
      PARA(1)=XMOM(1)-PARA(2)*(ONE-GAM)/G
      RETURN
C
C         ESTIMATED K EFFECTIVELY ZERO
C
   50 PARA(3)=ZERO
      PARA(2)=XMOM(2)/DL2
      PARA(1)=XMOM(1)-EU*PARA(2)
      RETURN
C
 1000 IFAIL=7000
      RETURN
C
C 7000 FORMAT(' *** ERROR *** ROUTINE PELGEV : L-MOMENTS INVALID')
C 7020 FORMAT(' ** WARNING ** ROUTINE PELGEV :',
C    *  ' ITERATION HAS NOT CONVERGED. RESULTS MAY BE UNRELIABLE.')
      END
C===================================================== PELGLO.FOR
      SUBROUTINE PELGLO(XMOM,PARA,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    July 2008                                           *
C*                                                                     *
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE GENERALIZED LOGISTIC
C  DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2, TAU-3.
C  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA, K (LOCATION, SCALE, SHAPE).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(3),PARA(3)
      DATA ZERO/0D0/,ONE/1D0/
      DATA PI/3.141592653589793238D0/
C
C         SMALL IS USED TO TEST WHETHER K IS EFFECTIVELY ZERO
C
      DATA SMALL/1D-6/
C
C         ESTIMATE K
C
      IFAIL=0
      G=-XMOM(3)
      IF(XMOM(2).LE.ZERO.OR.DABS(G).GE.ONE)GOTO 1000
      IF(DABS(G).LE.SMALL)GOTO 10
C
C         ESTIMATE ALPHA, XI
C
      GG=G*PI/DSIN(G*PI)
      A=XMOM(2)/GG
      PARA(1)=XMOM(1)-A*(ONE-GG)/G
      PARA(2)=A
      PARA(3)=G
      RETURN
C
C         ESTIMATED K EFFECTIVELY ZERO
C
   10 PARA(3)=ZERO
      PARA(2)=XMOM(2)
      PARA(1)=XMOM(1)
      RETURN
C
 1000 IFAIL=7000
      RETURN
C
C 7000 FORMAT(' *** ERROR *** ROUTINE PELGLO : L-MOMENTS INVALID')
      END
C===================================================== PELGNO.FOR
      SUBROUTINE PELGNO(XMOM,PARA,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    July 2008                                           *
C*                                                                     *
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE GENERALIZED NORMAL
C  DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2, TAU-3. ABS(TAU3) MAY NOT EXCEED 0.95.
C  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA, K (LOCATION, SCALE, SHAPE).
C
C  OTHER ROUTINES USED: DERF
C
C  METHOD: RATIONAL-FUNCTION APPROXIMATION OF K IN TERMS OF TAU-3
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(3),PARA(3)
      EXTERNAL DERF
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/
      DATA P95/0.95D0/
      DATA ROOTPI/1.772453850905516027D0/
C
C         COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATION
C         A0 IS 0.5*SQRT(3/PI)
C
      DATA A0,A1,A2,A3/
     *  0.20466534D+01,-0.36544371D+01,0.18396733D+01,-0.20360244D+00/
      DATA B1,B2,B3/-0.20182173D+01,0.12420401D+01,-0.21741801D+00/
C
C         SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
C
      DATA SMALL/1D-8/
C
      IFAIL=0
      T3=XMOM(3)
      IF(XMOM(2).LE.ZERO.OR.DABS(T3).GE.ONE)GOTO 1000
      IF(DABS(T3).GE.P95)GOTO 1010
      IF(DABS(T3).LE.SMALL)GOTO 30
C
      TT=T3*T3
      G=-T3*(A0+TT*(A1+TT*(A2+TT*A3)))/(ONE+TT*(B1+TT*(B2+TT*B3)))
      E=DEXP(HALF*G*G)
      A=XMOM(2)*G/(E*DERF(HALF*G))
      U=XMOM(1)+A*(E-ONE)/G
      PARA(1)=U
      PARA(2)=A
      PARA(3)=G
      RETURN
C
   30 PARA(1)=XMOM(1)
      PARA(2)=XMOM(2)*ROOTPI
      PARA(3)=ZERO
      RETURN
C
 1000 IFAIL=7000
      RETURN
 1010 IFAIL=7010
      PARA(1)=ZERO
      PARA(2)=-ONE
      PARA(3)=ZERO
      RETURN
C
C 7000 FORMAT(' *** ERROR *** ROUTINE PELGNO : L-MOMENTS INVALID')
C 7010 FORMAT(' *** ERROR *** ROUTINE PELGNO :',
C    *  ' TAU-3 TOO LARGE FOR ROUTINE')
      END
C===================================================== PELGPA.FOR
      SUBROUTINE PELGPA(XMOM,PARA,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    July 2008                                           *
C*                                                                     *
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR  THE GENERALIZED PARETO
C  DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2, TAU-3.
C  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA, K (LOCATION, SCALE, SHAPE).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(3),PARA(3)
      DATA ZERO/0D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/
C
      IFAIL=0
      T3=XMOM(3)
      IF(XMOM(2).LE.ZERO)GOTO 1000
      IF(DABS(T3).GE.ONE)GOTO 1000
      G=(ONE-THREE*T3)/(ONE+T3)
      PARA(3)=G
      PARA(2)=(ONE+G)*(TWO+G)*XMOM(2)
      PARA(1)=XMOM(1)-PARA(2)/(ONE+G)
      RETURN
C
 1000 IFAIL=7000
      RETURN
C
C 7000 FORMAT(' *** ERROR *** ROUTINE PELGPA : L-MOMENTS INVALID')
      END
C===================================================== PELKAP.FOR
      SUBROUTINE PELKAP(XMOM,PARA,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    July 2008                                           *
C*                                                                     *
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE KAPPA DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 4. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2, TAU-3, TAU-4.
C  PARA   *OUTPUT* ARRAY OF LENGTH 4. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA, K, H.
C  IFAIL  *OUTPUT* FAIL FLAG. ON EXIT, IT IS SET AS FOLLOWS.
C                  0  SUCCESSFUL EXIT
C                  7000  L-MOMENTS INVALID
C                  2  (TAU-3, TAU-4) LIES ABOVE THE GENERALIZED-LOGISTIC
C                     LINE (SUGGESTS THAT L-MOMENTS ARE NOT CONSISTENT
C                     WITH ANY KAPPA DISTRIBUTION WITH H.GT.-1)
C                  3  ITERATION FAILED TO CONVERGE
C                  4  UNABLE TO MAKE PROGRESS FROM CURRENT POINT IN
C                     ITERATION
C                  5  ITERATION ENCOUNTERED NUMERICAL DIFFICULTIES -
C                     OVERFLOW WOULD HAVE BEEN LIKELY TO OCCUR
C                  6  ITERATION FOR H AND K CONVERGED, BUT OVERFLOW
C                     WOULD HAVE OCCURRED WHEN CALCULATING XI AND ALPHA
C
C  N.B.  PARAMETERS ARE SOMETIMES NOT UNIQUELY DEFINED BY THE FIRST 4
C  L-MOMENTS. IN SUCH CASES THE ROUTINE RETURNS THE SOLUTION FOR WHICH
C  THE H PARAMETER IS LARGEST.
C
C  OTHER ROUTINES USED: DLGAMA,DIGAMD
C
C  THE SHAPE PARAMETERS K AND H ARE ESTIMATED USING NEWTON-RAPHSON
C  ITERATION ON THE RELATIONSHIP BETWEEN (TAU-3,TAU-4) AND (K,H).
C  THE CONVERGENCE CRITERION IS THAT TAU-3 AND TAU-4 CALCULATED FROM
C  THE ESTIMATED VALUES OF K AND H SHOULD DIFFER BY LESS THAN 'EPS'
C  FROM THE VALUES SUPPLIED IN ARRAY XMOM.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(4),PARA(4)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/,FOUR/4D0/
      DATA FIVE/5D0/,SIX/6D0/,TWELVE/12D0/,TWENTY/20D0/,THIRTY/30D0/
      DATA P725/0.725D0/,P8/0.8D0/
C
C         EPS,MAXIT CONTROL THE TEST FOR CONVERGENCE OF N-R ITERATION
C         MAXSR IS THE MAX. NO. OF STEPLENGTH REDUCTIONS PER ITERATION
C         HSTART IS THE STARTING VALUE FOR H
C         BIG IS USED TO INITIALIZE THE CRITERION FUNCTION
C         OFLEXP IS SUCH THAT DEXP(OFLEXP) JUST DOES NOT CAUSE OVERFLOW
C         OFLGAM IS SUCH THAT DEXP(DLGAMA(OFLGAM)) JUST DOES NOT CAUSE
C           OVERFLOW
C
      DATA EPS/1D-6/,MAXIT/20/,MAXSR/10/,HSTART/1.001D0/,BIG/10D0/
      DATA OFLEXP/170D0/,OFLGAM/53D0/
C
C----
C Next lines to assuage compilers that warn that variables
C 'might be used uninitialized'
      DEL1=ZERO
      DEL2=ZERO
      XG=ZERO
      XH=ZERO
C----
      T3=XMOM(3)
      T4=XMOM(4)
      DO 10 I=1,4
   10 PARA(I)=ZERO
C
C         TEST FOR FEASIBILITY
C
      IF(XMOM(2).LE.ZERO)GOTO 1000
      IF(DABS(T3).GE.ONE.OR.DABS(T4).GE.ONE)GOTO 1000
      IF(T4.LE.(FIVE*T3*T3-ONE)/FOUR)GOTO 1000
      IF(T4.GE.(FIVE*T3*T3+ONE)/SIX )GOTO 1010
C
C         SET STARTING VALUES FOR N-R ITERATION:
C         G IS CHOSEN TO GIVE THE CORRECT VALUE OF TAU-3 ON THE
C         ASSUMPTION THAT H=1 (I.E. A GENERALIZED PARETO FIT) -
C         BUT H IS ACTUALLY SET TO 1.001 TO AVOID NUMERICAL
C         DIFFICULTIES WHICH CAN SOMETIMES ARISE WHEN H=1 EXACTLY
C
      G=(ONE-THREE*T3)/(ONE+T3)
      H=HSTART
      Z=G+H*P725
      XDIST=BIG
C
C         START OF NEWTON-RAPHSON ITERATION
C
      DO 100 IT=1,MAXIT
C
C         REDUCE STEPLENGTH UNTIL WE ARE NEARER TO THE REQUIRED
C         VALUES OF TAU-3 AND TAU-4 THAN WE WERE AT THE PREVIOUS STEP
C
      DO 40 I=1,MAXSR
C
C         - CALCULATE CURRENT TAU-3 AND TAU-4
C
C           NOTATION:
C           U.    - RATIOS OF GAMMA FUNCTIONS WHICH OCCUR IN THE PWM'S
C                   BETA-SUB-R
C           ALAM. - L-MOMENTS (APART FROM A LOCATION AND SCALE SHIFT)
C           TAU.  - L-MOMENT RATIOS
C
      IF(G.GT.OFLGAM)GOTO 1020
      IF(H.GT.ZERO)GOTO 20
      U1=DEXP(DLGAMA(  -ONE/H-G)-DLGAMA(  -ONE/H+ONE))
      U2=DEXP(DLGAMA(  -TWO/H-G)-DLGAMA(  -TWO/H+ONE))
      U3=DEXP(DLGAMA(-THREE/H-G)-DLGAMA(-THREE/H+ONE))
      U4=DEXP(DLGAMA( -FOUR/H-G)-DLGAMA( -FOUR/H+ONE))
      GOTO 30
   20 U1=DEXP(DLGAMA(  ONE/H)-DLGAMA(  ONE/H+ONE+G))
      U2=DEXP(DLGAMA(  TWO/H)-DLGAMA(  TWO/H+ONE+G))
      U3=DEXP(DLGAMA(THREE/H)-DLGAMA(THREE/H+ONE+G))
      U4=DEXP(DLGAMA( FOUR/H)-DLGAMA( FOUR/H+ONE+G))
   30 CONTINUE
      ALAM2=U1-TWO*U2
      ALAM3=-U1+SIX*U2-SIX*U3
      ALAM4=U1-TWELVE*U2+THIRTY*U3-TWENTY*U4
      IF(ALAM2.EQ.ZERO)GOTO 1020
      TAU3=ALAM3/ALAM2
      TAU4=ALAM4/ALAM2
      E1=TAU3-T3
      E2=TAU4-T4
C
C         - IF NEARER THAN BEFORE, EXIT THIS LOOP
C
      DIST=DMAX1(DABS(E1),DABS(E2))
      IF(DIST.LT.XDIST)GOTO 50
C
C         - OTHERWISE, HALVE THE STEPLENGTH AND TRY AGAIN
C
      DEL1=HALF*DEL1
      DEL2=HALF*DEL2
      G=XG-DEL1
      H=XH-DEL2
   40 CONTINUE
C
C         TOO MANY STEPLENGTH REDUCTIONS
C
      IFAIL=4
      RETURN
C
C         TEST FOR CONVERGENCE
C
   50 CONTINUE
      IF(DIST.LT.EPS)GOTO 110
C
C         NOT CONVERGED: CALCULATE NEXT STEP
C
C         NOTATION:
C         U1G  - DERIVATIVE OF U1 W.R.T. G
C         DL2G - DERIVATIVE OF ALAM2 W.R.T. G
C         D..  - MATRIX OF DERIVATIVES OF TAU-3 AND TAU-4 W.R.T. G AND H
C         H..  - INVERSE OF DERIVATIVE MATRIX
C         DEL. - STEPLENGTH
C
      XG=G
      XH=H
      XZ=Z
      XDIST=DIST
      RHH=ONE/(H*H)
      IF(H.GT.ZERO)GOTO 60
      U1G=-U1*DIGAMD(  -ONE/H-G)
      U2G=-U2*DIGAMD(  -TWO/H-G)
      U3G=-U3*DIGAMD(-THREE/H-G)
      U4G=-U4*DIGAMD( -FOUR/H-G)
      U1H=      RHH*(-U1G-U1*DIGAMD(  -ONE/H+ONE))
      U2H=  TWO*RHH*(-U2G-U2*DIGAMD(  -TWO/H+ONE))
      U3H=THREE*RHH*(-U3G-U3*DIGAMD(-THREE/H+ONE))
      U4H= FOUR*RHH*(-U4G-U4*DIGAMD( -FOUR/H+ONE))
      GOTO 70
   60 U1G=-U1*DIGAMD(  ONE/H+ONE+G)
      U2G=-U2*DIGAMD(  TWO/H+ONE+G)
      U3G=-U3*DIGAMD(THREE/H+ONE+G)
      U4G=-U4*DIGAMD( FOUR/H+ONE+G)
      U1H=      RHH*(-U1G-U1*DIGAMD(  ONE/H))
      U2H=  TWO*RHH*(-U2G-U2*DIGAMD(  TWO/H))
      U3H=THREE*RHH*(-U3G-U3*DIGAMD(THREE/H))
      U4H= FOUR*RHH*(-U4G-U4*DIGAMD( FOUR/H))
   70 CONTINUE
      DL2G=U1G-TWO*U2G
      DL2H=U1H-TWO*U2H
      DL3G=-U1G+SIX*U2G-SIX*U3G
      DL3H=-U1H+SIX*U2H-SIX*U3H
      DL4G=U1G-TWELVE*U2G+THIRTY*U3G-TWENTY*U4G
      DL4H=U1H-TWELVE*U2H+THIRTY*U3H-TWENTY*U4H
      D11=(DL3G-TAU3*DL2G)/ALAM2
      D12=(DL3H-TAU3*DL2H)/ALAM2
      D21=(DL4G-TAU4*DL2G)/ALAM2
      D22=(DL4H-TAU4*DL2H)/ALAM2
      DET=D11*D22-D12*D21
      H11= D22/DET
      H12=-D12/DET
      H21=-D21/DET
      H22= D11/DET
      DEL1=E1*H11+E2*H12
      DEL2=E1*H21+E2*H22
C
C         TAKE NEXT N-R STEP
C
      G=XG-DEL1
      H=XH-DEL2
      Z=G+H*P725
C
C         REDUCE STEP IF G AND H ARE OUTSIDE THE PARAMETER SPACE
C
      FACTOR=ONE
      IF(G.LE.-ONE)FACTOR=P8*(XG+ONE)/DEL1
      IF(H.LE.-ONE)FACTOR=DMIN1(FACTOR,P8*(XH+ONE)/DEL2)
      IF(Z.LE.-ONE)FACTOR=DMIN1(FACTOR,P8*(XZ+ONE)/(XZ-Z))
      IF(H.LE.ZERO.AND.G*H.LE.-ONE)
     *  FACTOR=DMIN1(FACTOR,P8*(XG*XH+ONE)/(XG*XH-G*H))
      IF(FACTOR.EQ.ONE)GOTO 80
      DEL1=DEL1*FACTOR
      DEL2=DEL2*FACTOR
      G=XG-DEL1
      H=XH-DEL2
      Z=G+H*P725
   80 CONTINUE
C
C         END OF NEWTON-RAPHSON ITERATION
C
  100 CONTINUE
C
C         NOT CONVERGED
C
      IFAIL=3
      RETURN
C
C         CONVERGED
C
  110 IFAIL=0
      PARA(4)=H
      PARA(3)=G
      TEMP=DLGAMA(ONE+G)
      IF(TEMP.GT.OFLEXP)GOTO 1030
      GAM=DEXP(TEMP)
      TEMP=(ONE+G)*DLOG(DABS(H))
      IF(TEMP.GT.OFLEXP)GOTO 1030
      HH=DEXP(TEMP)
      PARA(2)=XMOM(2)*G*HH/(ALAM2*GAM)
      PARA(1)=XMOM(1)-PARA(2)/G*(ONE-GAM*U1/HH)
      RETURN
C
 1000 IFAIL=7000
      RETURN
 1010 IFAIL=2
      RETURN
 1020 IFAIL=5
      RETURN
 1030 IFAIL=6
      RETURN
C
      END
C===================================================== PELPE3.FOR
      SUBROUTINE PELPE3(XMOM,PARA,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    July 2008                                           *
C*                                                                     *
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE PEARSON TYPE 3 DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2 AND TAU-3.
C  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER MU, SIGMA, GAMMA (MEAN, S.D., SKEWNESS).
C
C  OTHER ROUTINES USED: DLGAMA
C
C  METHOD: RATIONAL APPROXIMATION IS USED TO EXPRESS ALPHA, THE SHAPE
C  PARAMETER OF THE GAMMA DISTRIBUTION, AS A FUNCTION OF TAU-3.
C  RELATIVE ACCURACY OF THE APPROXIMATION IS BETTER THAN 3E-5.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(3),PARA(3)
      DATA ZERO/0D0/,THIRD/0.33333333D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/
C
C         SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
C
      DATA SMALL/1D-6/
C
C         CONSTANTS USED IN MINIMAX APPROXIMATIONS
C
      DATA C1,C2,C3/ 0.2906D0,  0.1882D0,  0.0442D0/
      DATA D1,D2,D3/ 0.36067D0,-0.59567D0, 0.25361D0/
      DATA D4,D5,D6/-2.78861D0, 2.56096D0,-0.77045D0/
      DATA PI3,ROOTPI/9.4247780D0,1.7724539D0/
C
      IFAIL=0
      T3=DABS(XMOM(3))
      IF(XMOM(2).LE.ZERO.OR.T3.GE.ONE)GOTO 1000
      IF(T3.LE.SMALL)GOTO 100
      IF(T3.GE.THIRD)GOTO 10
      T=PI3*T3*T3
      ALPHA=(ONE+C1*T)/(T*(ONE+T*(C2+T*C3)))
      GOTO 20
   10 CONTINUE
      T=ONE-T3
      ALPHA=T*(D1+T*(D2+T*D3))/(ONE+T*(D4+T*(D5+T*D6)))
   20 CONTINUE
      RTALPH=DSQRT(ALPHA)
      BETA=ROOTPI*XMOM(2)*DEXP(DLGAMA(ALPHA)-DLGAMA(ALPHA+HALF))
      PARA(1)=XMOM(1)
      PARA(2)=BETA*RTALPH
      PARA(3)=TWO/RTALPH
      IF(XMOM(3).LT.ZERO)PARA(3)=-PARA(3)
      RETURN
C
C         ZERO SKEWNESS
C
  100 CONTINUE
      PARA(1)=XMOM(1)
      PARA(2)=XMOM(2)*ROOTPI
      PARA(3)=ZERO
      RETURN
C
 1000 IFAIL=7000
      DO 1010 I=1,3
 1010 PARA(I)=ZERO
      RETURN
C
C 7000 FORMAT(' *** ERROR *** ROUTINE PELPE3 : L-MOMENTS INVALID')
      END
C===================================================== PELWAK.FOR
      SUBROUTINE PELWAK(XMOM,PARA,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmomRFA"                       *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 0.5    June 2007                                           *
C*                                                                     *
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE WAKEBY DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 5. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2, TAU-3, TAU-4, TAU-5.
C  PARA   *OUTPUT* ARRAY OF LENGTH 5. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA, BETA, GAMMA, DELTA.
C  IFAIL  *OUTPUT* FAIL FLAG. ON EXIT, IT IS SET AS FOLLOWS.
C                  0 SUCCESSFUL EXIT
C                  1 ESTIMATES COULD ONLY BE OBTAINED BY SETTING XI=0
C                  2 ESTIMATES COULD ONLY BE OBTAINED BY FITTING A
C                    GENERALIZED PARETO DISTRIBUTION
C                  7000 L-MOMENTS INVALID
C
C  PROCEDURE:
C  1. LOOK FOR A SOLUTION WITH XI UNCONSTRAINED;
C  2. IF NONE FOUND, LOOK FOR A SOLUTION WITH XI=0;
C  3. IF NONE FOUND, FIT A GENERALIZED PARETO DISTRIBUTION TO THE
C     FIRST 3 L-MOMENTS.
C  ESTIMATES ARE CALCULATED USING THE FORMULAS GIVEN BY GREENWOOD ET AL.
C  (1979, WATER RESOUR. RES., TABLE 5), BUT EXPRESSED IN TERMS OF
C  L-MOMENTS RATHER THAN PROBABILITY WEIGHTED MOMENTS.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(5),PARA(5)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/,FOUR/4D0/
      DATA X2/2D0/,X3/3D0/,X4/4D0/,X5/5D0/,X7/7D0/,X8/8D0/,X9/9D0/,
     *  X10/10D0/,X11/11D0/,X16/16D0/,X25/25D0/,X29/29D0/,X32/32D0/,
     *  X35/35D0/,X85/85D0/,X125/125D0/,X203/203D0/
C
      IF(XMOM(2).LE.ZERO)GOTO 1000
      IF(DABS(XMOM(3)).GE.ONE)GOTO 1000
      IF(DABS(XMOM(4)).GE.ONE)GOTO 1000
      IF(DABS(XMOM(5)).GE.ONE)GOTO 1000
      IFAIL=0
C
C         CALCULATE THE L-MOMENTS (LAMBDA'S)
C
      ALAM1=XMOM(1)
      ALAM2=XMOM(2)
      ALAM3=XMOM(3)*ALAM2
      ALAM4=XMOM(4)*ALAM2
      ALAM5=XMOM(5)*ALAM2
C
C         ESTIMATE N1,N2,N3,C1,C2,C3 WHEN XI.NE.0
C
      XN1= X3*ALAM2-X25*ALAM3 +X32*ALAM4
      XN2=-X3*ALAM2 +X5*ALAM3  +X8*ALAM4
      XN3= X3*ALAM2 +X5*ALAM3  +X2*ALAM4
      XC1= X7*ALAM2-X85*ALAM3+X203*ALAM4-X125*ALAM5
      XC2=-X7*ALAM2+X25*ALAM3  +X7*ALAM4 -X25*ALAM5
      XC3= X7*ALAM2 +X5*ALAM3  -X7*ALAM4  -X5*ALAM5
C
C         ESTIMATE B AND D
C
      XA=XN2*XC3-XC2*XN3
      XB=XN1*XC3-XC1*XN3
      XC=XN1*XC2-XC1*XN2
      DISC=XB*XB-FOUR*XA*XC
      IF(DISC.LT.ZERO)GOTO 10
      DISC=DSQRT(DISC)
      ROOT1=HALF*(-XB+DISC)/XA
      ROOT2=HALF*(-XB-DISC)/XA
      B= DMAX1(ROOT1,ROOT2)
      D=-DMIN1(ROOT1,ROOT2)
      IF(D.GE.ONE)GOTO 10
C
C         ESTIMATE A, C AND XI
C
      A=(ONE+B)*(TWO+B)*(THREE+B)/
     *  (FOUR*(B+D))*((ONE+D)*ALAM2-(THREE-D)*ALAM3)
      C=-(ONE-D)*(TWO-D)*(THREE-D)/
     *  (FOUR*(B+D))*((ONE-B)*ALAM2-(THREE+B)*ALAM3)
      XI=ALAM1-A/(ONE+B)-C/(ONE-D)
C
C         CHECK FOR VALID PARAMETERS
C
      IF(C.GE.ZERO.AND.A+C.GE.ZERO)GOTO 30
C
C         CAN'T FIND VALID ESTIMATES FOR XI UNRESTRICTED, SO TRY XI=0
C
C         ESTIMATE B AND D FOR XI=0
C
   10 IFAIL=1
      XI=ZERO
      ZN1=X4*ALAM1-X11*ALAM2+X9*ALAM3
      ZN2=-ALAM2+X3*ALAM3
      ZN3=ALAM2+ALAM3
      ZC1=X10*ALAM1-X29*ALAM2+X35*ALAM3-X16*ALAM4
      ZC2=-ALAM2+X5*ALAM3-X4*ALAM4
      ZC3=ALAM2-ALAM4
      ZA=ZN2*ZC3-ZC2*ZN3
      ZB=ZN1*ZC3-ZC1*ZN3
      ZC=ZN1*ZC2-ZC1*ZN2
      DISC=ZB*ZB-FOUR*ZA*ZC
      IF(DISC.LT.ZERO)GOTO 20
      DISC=DSQRT(DISC)
      ROOT1=HALF*(-ZB+DISC)/ZA
      ROOT2=HALF*(-ZB-DISC)/ZA
      B= DMAX1(ROOT1,ROOT2)
      D=-DMIN1(ROOT1,ROOT2)
      IF(D.GE.ONE)GOTO 20
C
C         ESTIMATE A AND C
C
      A= (ONE+B)*(TWO+B)/(B+D)*(ALAM1-(TWO-D)*ALAM2)
      C=-(ONE-D)*(TWO-D)/(B+D)*(ALAM1-(TWO+B)*ALAM2)
      IF(C.GE.ZERO.AND.A+C.GE.ZERO)GOTO 30
C
C         CAN'T FIND VALID ESTIMATES EVEN WITH XI=0 -
C         FIT GENERALIZED PARETO DISTRIBUTION INSTEAD
C
   20 IFAIL=2
      D=-(ONE-THREE*XMOM(3))/(ONE+XMOM(3))
      C=(ONE-D)*(TWO-D)*XMOM(2)
      B=ZERO
      A=ZERO
      XI=XMOM(1)-C/(ONE-D)
      IF(D.GT.ZERO)GOTO 30
      A=C
      B=-D
      C=ZERO
      D=ZERO
C
C         COPY RESULTS INTO ARRAY PARA
C
   30 PARA(1)=XI
      PARA(2)=A
      PARA(3)=B
      PARA(4)=C
      PARA(5)=D
      RETURN
C
 1000 IFAIL=7000
      DO 1010 I=1,5
 1010 PARA(I)=ZERO
      END
C===================================================== QKAP.FOR
      SUBROUTINE QKAP(X,N,PARA)
C
C  QUANTILE FUNCTION OF THE KAPPA DISTRIBUTION
C  Implemented as a subroutine for a vector of inputs
C  Does not test validity of parameters
C
C  Arguments
C  X      *in/out* Array of length N. On input, contains the
C                  probabilities (arguments of the quantile function).
C                  On exit, contains the quantiles.
C  N      * input* Size of array X.
C  PARA   * input* Array of length 4. Contains the parameters of the
C                  distribution.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N),PARA(4)
      DATA ZERO/0D0/,ONE/1D0/
      IF (PARA(4).EQ.ZERO) THEN
        DO 10 I=1,N
        X(I)=PARA(1)+PARA(2)/PARA(3)*(ONE-(-LOG(X(I)))**PARA(3))
   10   CONTINUE
      ELSE IF (PARA(4).EQ.-ONE) THEN
        DO 20 I=1,N
        X(I)=PARA(1)+PARA(2)/PARA(3)*(ONE-(ONE/X(I)-ONE)**PARA(3))
   20   CONTINUE
      ELSE
        DO 30 I=1,N
        X(I)=PARA(1)+PARA(2)/PARA(3)*
     *    (ONE-((ONE-X(I)**PARA(4))/PARA(4))**PARA(3))
   30   CONTINUE
      END IF
      RETURN
      END
C===================================================== SAMLMU.FOR
      SUBROUTINE SAMLMU(X,N,XMOM,NMOM,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmom"                          *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    July 2008                                           *
C*                                                                     *
C***********************************************************************
C
C  SAMPLE L-MOMENTS OF A DATA ARRAY
C
C  PARAMETERS OF ROUTINE:
C  X      * INPUT* ARRAY OF LENGTH N. CONTAINS THE DATA, IN ASCENDING
C                  ORDER.
C  N      * INPUT* NUMBER OF DATA VALUES
C  XMOM   *OUTPUT* ARRAY OF LENGTH NMOM. CONTAINS THE SAMPLE L-MOMENTS,
C                  STORED AS DESCRIBED BELOW.
C  NMOM   * INPUT* NUMBER OF L-MOMENTS TO BE FOUND.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        Note: Fortran 77 would require dimension of COEF to be known
C        at compile time.
      DOUBLE PRECISION X(N),XMOM(NMOM),COEF(2,NMOM)
      DATA ZERO/0D0/,ONE/1D0/,TWO/2D0/
C
      DN=N
      DO 10 J=1,NMOM
   10 XMOM(J)=ZERO
      IF(NMOM.LE.2)GOTO 100
C
C         UNBIASED ESTIMATES OF L-MOMENTS -- THE 'DO 30' LOOP
C         RECURSIVELY CALCULATES DISCRETE LEGENDRE POLYNOMIALS, VIA
C         EQ.(9) OF NEUMAN AND SCHONBACH (1974, INT.J.NUM.METH.ENG.)
C
      DO 20 J=3,NMOM
      TEMP=ONE/DFLOAT((J-1)*(N-J+1))
      COEF(1,J)=DFLOAT(J+J-3)*TEMP
      COEF(2,J)=DFLOAT((J-2)*(N+J-2))*TEMP
   20 CONTINUE
      TEMP=-DN-ONE
      CONST=ONE/(DN-ONE)
      NHALF=N/2
      DO 40 I=1,NHALF
      TEMP=TEMP+TWO
      XI=X(I)
      XII=X(N+1-I)
      TERMP=XI+XII
      TERMN=XI-XII
      XMOM(1)=XMOM(1)+TERMP
      S1=ONE
      S=TEMP*CONST
      XMOM(2)=XMOM(2)+S*TERMN
      DO 30 J=3,NMOM,2
      S2=S1
      S1=S
      S=COEF(1,J)*TEMP*S1-COEF(2,J)*S2
      XMOM(J)=XMOM(J)+S*TERMP
      IF(J.EQ.NMOM)GOTO 30
      JJ=J+1
      S2=S1
      S1=S
      S=COEF(1,JJ)*TEMP*S1-COEF(2,JJ)*S2
      XMOM(JJ)=XMOM(JJ)+S*TERMN
   30 CONTINUE
   40 CONTINUE
      IF(N.EQ.NHALF+NHALF)GOTO 60
      TERM=X(NHALF+1)
      S=ONE
      XMOM(1)=XMOM(1)+TERM
      DO 50 J=3,NMOM,2
      S=-COEF(2,J)*S
      XMOM(J)=XMOM(J)+S*TERM
   50 CONTINUE
C
C         L-MOMENT RATIOS
C
   60 CONTINUE
      XMOM(1)=XMOM(1)/DN
      IF(XMOM(2).EQ.ZERO)GOTO 1010
      DO 70 J=3,NMOM
   70 XMOM(J)=XMOM(J)/XMOM(2)
      XMOM(2)=XMOM(2)/DN
      RETURN
C
C         AT MOST TWO L-MOMENTS
C
  100 CONTINUE
      SUM1=ZERO
      SUM2=ZERO
      TEMP=-DN+ONE
      DO 110 I=1,N
      SUM1=SUM1+X(I)
      SUM2=SUM2+X(I)*TEMP
      TEMP=TEMP+TWO
  110 CONTINUE
      XMOM(1)=SUM1/DN
      IF(NMOM.EQ.1)RETURN
      XMOM(2)=SUM2/(DN*(DN-ONE))
      RETURN
C
 1010 IFAIL=7010
      DO 1020 J=2,NMOM
 1020 XMOM(J)=ZERO
      RETURN
C
C7010 FORMAT(' *** ERROR *** ROUTINE SAMLMU :',
C    *  ' ALL DATA VALUES EQUAL')
      END
C===================================================== DERF.FOR
      DOUBLE PRECISION FUNCTION DERF(X)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  ERROR FUNCTION
C
C  BASED ON ALGORITHM 5666, J.F.HART ET AL. (1968) 'COMPUTER
C  APPROXIMATIONS'
C
C  ACCURATE TO 15 DECIMAL PLACES
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA ZERO/0D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/,FOUR/4D0/,P65/0.65D0/
C
C         COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATION
C
      DATA P0,P1,P2,P3,P4,P5,P6/
     *  0.22020 68679 12376 1D3,    0.22121 35961 69931 1D3,
     *  0.11207 92914 97870 9D3,    0.33912 86607 83830 0D2,
     *  0.63739 62203 53165 0D1,    0.70038 30644 43688 1D0,
     *  0.35262 49659 98910 9D-1/
      DATA Q0,Q1,Q2,Q3,Q4,Q5,Q6,Q7/
     *  0.44041 37358 24752 2D3,   0.79382 65125 19948 4D3,
     *  0.63733 36333 78831 1D3,   0.29656 42487 79673 7D3,
     *  0.86780 73220 29460 8D2,   0.16064 17757 92069 5D2,
     *  0.17556 67163 18264 2D1,   0.88388 34764 83184 4D-1/
C
C         C1 IS SQRT(2), C2 IS SQRT(2/PI)
C         BIG IS THE POINT AT WHICH DERF=1 TO MACHINE PRECISION
C
      DATA C1/1.4142 13562 37309 5D0/
      DATA C2/7.9788 45608 02865 4D-1/
      DATA BIG/6.25D0/,CRIT/5D0/
C
      DERF=ZERO
      IF(X.EQ.ZERO)RETURN
      XX=DABS(X)
      IF(XX.GT.BIG)GOTO 20
      EXPNTL=DEXP(-X*X)
      ZZ=DABS(X*C1)
      IF(XX.GT.CRIT)GOTO 10
      DERF=EXPNTL*((((((P6*ZZ+P5)*ZZ+P4)*ZZ+P3)*ZZ+P2)*ZZ+P1)*ZZ+P0)/
     *  (((((((Q7*ZZ+Q6)*ZZ+Q5)*ZZ+Q4)*ZZ+Q3)*ZZ+Q2)*ZZ+Q1)*ZZ+Q0)
      IF(X.GT.ZERO)DERF=ONE-TWO*DERF
      IF(X.LT.ZERO)DERF=TWO*DERF-ONE
      RETURN
C
   10 DERF=EXPNTL*C2/(ZZ+ONE/(ZZ+TWO/(ZZ+THREE/(ZZ+FOUR/(ZZ+P65)))))
      IF(X.GT.ZERO)DERF=ONE-DERF
      IF(X.LT.ZERO)DERF=DERF-ONE
      RETURN
C
   20 DERF=ONE
      IF(X.LT.ZERO)DERF=-ONE
      RETURN
      END
C===================================================== DIGAMD.FOR
      DOUBLE PRECISION FUNCTION DIGAMD(X)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  DIGAMMA FUNCTION (EULER'S PSI FUNCTION) - THE FIRST DERIVATIVE OF
C  LOG(GAMMA(X))
C
C  BASED ON ALGORITHM AS103, APPL. STATIST. (1976) VOL.25 NO.3
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/
      DATA SMALL/1D-9/,CRIT/13D0/
C
C         C1...C7 ARE THE COEFFTS OF THE ASYMPTOTIC EXPANSION OF DIGAMD
C         D1 IS  -(EULER'S CONSTANT)
C
      DATA C1,C2,C3,C4,C5,C6,C7,D1/
     *  0.83333 33333 33333 333D-1,  -0.83333 33333 33333 333D-2,
     *  0.39682 53968 25396 825D-2,  -0.41666 66666 66666 666D-2,
     *  0.75757 57575 75757 575D-2,  -0.21092 79609 27960 928D-1,
     *  0.83333 33333 33333 333D-1,  -0.57721 56649 01532 861D 0/
      DIGAMD=ZERO
      IF(X.LE.ZERO)GOTO 1000
C
C         USE SMALL-X APPROXIMATION IF X.LE.SMALL
C
      IF(X.GT.SMALL)GOTO 10
      DIGAMD=D1-ONE/X
      RETURN
C
C         REDUCE TO DIGAMD(X+N) WHERE X+N.GE.CRIT
C
   10 Y=X
   20 IF(Y.GE.CRIT)GOTO 30
      DIGAMD=DIGAMD-ONE/Y
      Y=Y+ONE
      GOTO 20
C
C         USE ASYMPTOTIC EXPANSION IF Y.GE.CRIT
C
   30 DIGAMD=DIGAMD+DLOG(Y)-HALF/Y
      Y=ONE/(Y*Y)
      SUM=((((((C7*Y+C6)*Y+C5)*Y+C4)*Y+C3)*Y+C2)*Y+C1)*Y
      DIGAMD=DIGAMD-SUM
      RETURN
C
 1000 CONTINUE
C     WRITE(6,7000)X
      RETURN
C
C7000 FORMAT(' *** ERROR *** ROUTINE DIGAMD :',
C    *  ' ARGUMENT OUT OF RANGE :',D24.16)
      END
C===================================================== DLGAMA.FOR
      DOUBLE PRECISION FUNCTION DLGAMA(X)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  LOGARITHM OF GAMMA FUNCTION
C
C  BASED ON ALGORITHM ACM291, COMMUN. ASSOC. COMPUT. MACH. (1966)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA SMALL,CRIT,BIG,TOOBIG/1D-7,13D0,1D9,2D36/
C
C         C0 IS 0.5*LOG(2*PI)
C         C1...C7 ARE THE COEFFTS OF THE ASYMPTOTIC EXPANSION OF DLGAMA
C
      DATA C0,C1,C2,C3,C4,C5,C6,C7/
     *   0.91893 85332 04672 742D 0,  0.83333 33333 33333 333D-1,
     *  -0.27777 77777 77777 778D-2,  0.79365 07936 50793 651D-3,
     *  -0.59523 80952 38095 238D-3,  0.84175 08417 50841 751D-3,
     *  -0.19175 26917 52691 753D-2,  0.64102 56410 25641 026D-2/
C
C         S1 IS -(EULER'S CONSTANT), S2 IS PI**2/12
C
      DATA S1/-0.57721 56649 01532 861D 0/
      DATA S2/ 0.82246 70334 24113 218D 0/
C
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/
      DLGAMA=ZERO
      IF(X.LE.ZERO)GOTO 1000
      IF(X.GT.TOOBIG)GOTO 1000
C
C         USE SMALL-X APPROXIMATION IF X IS NEAR 0, 1 OR 2
C
      IF(DABS(X-TWO).GT.SMALL)GOTO 10
      DLGAMA=DLOG(X-ONE)
      XX=X-TWO
      GOTO 20
   10 IF(DABS(X-ONE).GT.SMALL)GOTO 30
      XX=X-ONE
   20 DLGAMA=DLGAMA+XX*(S1+XX*S2)
      RETURN
   30 IF(X.GT.SMALL)GOTO 40
      DLGAMA=-DLOG(X)+S1*X
      RETURN
C
C         REDUCE TO DLGAMA(X+N) WHERE X+N.GE.CRIT
C
   40 SUM1=ZERO
      Y=X
      IF(Y.GE.CRIT)GOTO 60
      Z=ONE
   50 Z=Z*Y
      Y=Y+ONE
      IF(Y.LT.CRIT)GOTO 50
      SUM1=SUM1-DLOG(Z)
C
C         USE ASYMPTOTIC EXPANSION IF Y.GE.CRIT
C
   60 SUM1=SUM1+(Y-HALF)*DLOG(Y)-Y+C0
      SUM2=ZERO
      IF(Y.GE.BIG)GOTO 70
      Z=ONE/(Y*Y)
      SUM2=((((((C7*Z+C6)*Z+C5)*Z+C4)*Z+C3)*Z+C2)*Z+C1)/Y
   70 DLGAMA=SUM1+SUM2
      RETURN
C
 1000 CONTINUE
C     WRITE(6,7000)X
      RETURN
C
C7000 FORMAT(' *** ERROR *** ROUTINE DLGAMA :',
C    *  ' ARGUMENT OUT OF RANGE :',D24.16)
      END
C===================================================== SORT.FOR
      SUBROUTINE SORT(X,N)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  SORTS THE ARRAY X INTO ASCENDING ORDER
C
C  PARAMETERS OF ROUTINE:
C  X      *IN/OUT* ARRAY OF LENGTH N. CONTAINS THE NUMBERS TO BE SORTED.
C                  ON EXIT, CONTAINS THE SORTED NUMBERS.
C  N      * INPUT* NUMBER OF ELEMENTS TO BE SORTED
C
C  METHOD USED IS SHELL SORT WITH SEQUENCE OF INCREMENTS AS IN
C  D.F.KNUTH (1969) 'THE ART OF COMPUTER PROGRAMMING', VOL.3, P.95
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N)
      IF(N.LE.1)RETURN
      J=4
      DO 10 I=1,100
      J=3*J+1
      IF(J.GE.N)GOTO 20
   10 CONTINUE
   20 CONTINUE
      M=(J/3)
      DO 60 MM=1,100
      M=M/3
      IF(M.EQ.0)RETURN
      DO 50 I=M+1,N
      TEST=X(I)
      J=I
      DO 30 JJ=1,100
      J=J-M
      IF(J.LE.0)GOTO 40
      IF(TEST.GE.X(J))GOTO 40
   30 X(J+M)=X(J)
   40 CONTINUE
   50 X(J+M)=TEST
   60 CONTINUE
      END
