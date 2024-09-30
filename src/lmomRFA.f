C
C  Fortran code for the R package "lmomRFA".
C  Based on the LMOMENTS Fortran package, version 3.04.
C
C  The following routines are called from R functions:
C
C    QKAP
C    REGTST
C
C  The following routines are called from other Fortran routines.
C
C    PELWAK
C    QKAP
C
C
C  File lmomrfa-lmom.f contains the following routines, copied from
C  package "lmom". All are called from other Fortran routines.
C
C    LMRGEV
C    LMRGLO
C    LMRGNO
C    LMRGPA
C    LMRPE3
C    PELGEV
C    PELGLO
C    PELGNO
C    PELGPA
C    PELKAP
C    PELPE3
C    DIGAMD
C    XERF
C    XLGAMA
C    SAMLM
C
C  Note that routine PELWAK in this package is not the same as in R
C  package "lmom".  The version in this package retains the behaviour
C  of PELWAK in the LMOMENTS Fortran package: if a fit using all five
C  parameters fails, an attempt is made to fit the distribution with
C  four parameters and lower bound zero.
C
C===================================================== pelwak.f
      SUBROUTINE PELWAK(XMOM,PARA,IFAIL)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmomRFA"                       *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 1.0    March 2009                                          *
C*                                                                     *
C*  Version 3.0    August 2023                                         *
C*  * Code cleanup:                                                    *
C*    - Specific names of intrinsic functions changed to generic.      *
C*    - All DO loops now end with CONTINUE.                            *
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
      IF(ABS(XMOM(3)).GE.ONE)GOTO 1000
      IF(ABS(XMOM(4)).GE.ONE)GOTO 1000
      IF(ABS(XMOM(5)).GE.ONE)GOTO 1000
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
      DISC=SQRT(DISC)
      ROOT1=HALF*(-XB+DISC)/XA
      ROOT2=HALF*(-XB-DISC)/XA
      B= MAX(ROOT1,ROOT2)
      D=-MIN(ROOT1,ROOT2)
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
      DISC=SQRT(DISC)
      ROOT1=HALF*(-ZB+DISC)/ZA
      ROOT2=HALF*(-ZB-DISC)/ZA
      B= MAX(ROOT1,ROOT2)
      D=-MIN(ROOT1,ROOT2)
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
      PARA(I)=ZERO
 1010 CONTINUE
      END
C===================================================== qkap.f
      SUBROUTINE QKAP(X,N,PARA)
C***********************************************************************
C*                                                                     *
C*  Fortran code written for R package "lmomRFA"                       *
C*                                                                     *
C*  J. R. M. Hosking <jrmhosking@gmail.com>                            *
C*                                                                     *
C*  Version 2.5    June 2013                                           *
C*                                                                     *
C***********************************************************************
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
C===================================================== regtst.f
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
C*  * Summary of differences from the Fortran version:                 *
C*    - Omit arguments NAMES,NPROB,PROB,KPRINT,KOUT,A,B,SEED.          *
C*    - Add arguments RPARA,T4FIT,WORK,X,MAXREC.                       *
C*    - No printing.                                                   *
C*    - Return (as *output* arguments) all numbers that used to be     *
C*      printed.                                                       *
C*    - No use of PARAMETERs.                                          *
C*    - Calls to PELxxx routines are to the versions adapted for R,    *
C*      so gain an IFAIL parameter.                                    *
C*    - Compute sample L-moments via SAMLMU, not SAMLMR.               *
C*    - Don't stop when unable to invert SS matrix.                    *
C*    - Fit the candidate (and Wakeby) distributions even when         *
C*      NSIM.EQ.0.                                                     *
C*    - Compute T4FIT more accurately.                                 *
C*    - Don't compute quantiles of candidate distributions.            *
C*    - Use R's random-number generator, via C routines RANGET,        *
C*      CURAND, RANPUT.                                                *
C*                                                                     *
C*  Version 3.5    February 2023                                       *
C*  * Use routine SAMLM from current package "lmom" rather than        *
C*    routine SAMLMU from an earlier version.                          *
C*  * No need to call routine SORT.                                    *
C*                                                                     *
C*  Version 3.6    August 2023                                         *
C*  * Code cleanup:                                                    *
C*    - Specific names of intrinsic functions changed to generic.      *
C*    - All DO loops now end with CONTINUE.                            *
C*    - No shared termination statements for DO loops.                 *
C*    - Some IF constructs changed to IF-THEN(-ELSE)-ENDIF, sometimes  *
C*      with rearrangement of code.                                    *
C*                                                                     *
C***********************************************************************
C
C  Calculates three statistics useful in regional frequency analysis
C
C  Discordancy measure, D(I), for individual sites in a region.
C      Large values might be used as a flag to indicate potential errors
C      in the data at the site. "Large" might be 3 for regions with 15
C      or more sites, but less (exact values in array DC1) for smaller
C      regions.
C
C  Heterogeneity measures, H(J), for a region based upon either:-
C      J=1: the weighted s.d. of the L-CVs or
C      J=2: the average distance from the site to the regional average
C           on a graph of L-CV vs. L-skewness
C      J=3: the average distance from the site to the regional average
C           on a graph of L-skewness vs. L-kurtosis
C
C      In practice H(1) is probably sufficient. A value greater than
C      (say) 1.0 suggests that further subdivision of the region should
C      be considered as it might improve quantile estimates.
C
C  Goodness-of-fit measures, Z(K), for 5 candidate distributions:
C      K=1: generalized logistic
C      K=2: generalized extreme value
C      K=3: generalized normal (lognormal)
C      K=4: Pearson type III (3-parameter gamma)
C      K=5: generalized Pareto
C
C      Provided that the region is acceptably close to homogeneous,
C      the fit may be judged acceptable at 10% significance level
C      if Z(K) is less than 1.645 in absolute value.
C
C  For further details see J.R.M. Hosking and J.R. Wallis (1997),
C  "Regional Frequency Analysis: an approach based on L-moments",
C  Cambridge University Press, chapters 3-5.
C
C  Parameters of routine:
C  NSITES * input* Number of sites in region
C  LEN    * input* Array of length NSITES. Record lengths at each site.
C  XMOM   * input* Array of dimension (5,NSITES). Array containing
C                  the first 5 sample L-moments for each site, in the
C                  order mean, L-CV, L-skewness, L-kurtosis, t-5, i.e
C                  XMOM(I,J) contains the I'th L-moment for site J.
C                    n.b. xmom(2,.) contains l-cv, not the usual l-2!
C  NSIM   * input* Number of simulated worlds for heterogeneity and
C                  goodness-of-fit tests.
C                    Note: NSIM=0 will force return at completion of
C                  outlier test. NSIM=1 will suppress calculation of
C                  H and Z statistics, but parameter and quantile
C                  estimates will be found.
C  RMOM   *output* Array of length 5. On exit, contains the regional
C                  weighted average L-moment ratios.
C  D      *output* Array of length nsites. on exit, contains the
C                  discordancy measure (D statistic) for each site.
C  VOBS   *output* Array of length 3. On exit, contains the regional
C                  observed values of three heterogeneity statistics:
C                  (1) weighted s.d. of L-CVs;
C                  (2) average of L-CV/L-skew distances;
C                  (3) average of L-skew/L-kurtosis distances.
C  VBAR   *output* Array of length 3. On exit, contains the mean of the
C                  simulated values of the 3 heterogeneity statistics.
C  VSD    *output* Array of length 3. On exit, contains the s.d. of the
C                  simulated values of the 3 heterogeneity statistics.
C  H      *output* Array of length 3. On exit, contains heterogeneity
C                  measures (H statistics), i.e. H=(VOBS-VBAR)/VSD.
C  Z      *output* Array of length 5. On exit, contains goodness-of-fit
C                  measures (Z statistics) for 5 distributions:
C                  (1) gen. logistic, (2) gen. extreme value,
C                  (3) gen. normal, (4) Pearson type III,
C                  (5) gen. Pareto.
C  PARA   *output* Array of dimension (5,6). On exit, if NSIM.ge.1,
C                  contains parameters of growth curves fitted by the
C                  above 5 distributions, plus Wakeby.
C  RPARA  *output* Array of length 4. On exit, contains parameters of
C                  kappa distribution fitted to regional L-moments.
C  T4FIT  *output* Array of length 5. On exit, contains the L-kurtosis
C                  values of the  5 fitted candidate distributions.
C  WORK   * local* Array of length 3*NSITES. Used as workspace.
C  X      * local* Array of length MAXREC. Used as workspace.
C  MAXREC * input* Must be at least as large as the largest of the
C                  record lengths.
C
C  Other Fortran routines used: DIGAMD,LMRGEV,LMRGLO,LMRGNO,LMRGPA,
C    LMRPE3,PELGEV,PELGLO,PELGNO,PELGPA,PELKAP,PELPE3,PELWAK,QKAP,
C    QSORT3(from R's API),SAMLM,XERF,XLGAMA
C
C  C routines used: CURAND,RANGET,RANPUT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION D(NSITES),H(3),PARA(5,6),
     *  RMOM(5),RPARA(4),SMAT(3,3),TMOM(4),T4FIT(5),
     *  VBAR(3),VOBS(3),VSD(3),WORK(NSITES,3),X(MAXREC),XMOM(5,NSITES),
     *  Z(5)
      INTEGER LEN(NSITES)
      DATA ZERO/0D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/
C
C         INITIALIZE ARRAYS
C
      NMAX=0
      SUMLEN=0
      DO 10 I=1,NSITES
      NREC=LEN(I)
      IF(NREC.GT.NMAX)NMAX=NREC
      SUMLEN=SUMLEN+NREC
      D(I)=ZERO
   10 CONTINUE
      DO 20 K=1,3
      VOBS(K)=ZERO
      VBAR(K)=ZERO
      VSD(K)=ZERO
      H(K)=ZERO
   20 CONTINUE
      DO 30 IDIST=1,5
      Z(IDIST)=ZERO
   30 CONTINUE
      DO 41 IPARA=1,5
      DO 40 IDIST=1,6
      PARA(IPARA,IDIST)=ZERO
   40 CONTINUE
   41 CONTINUE
C
C         CALCULATE THE WEIGHTED MEAN OF L-CV, L-SKEW, L-KURTOSIS
C
      DO 60 K=2,5
      RMOM(K)=ZERO
      DO 50 I=1,NSITES
      RMOM(K)=RMOM(K)+LEN(I)*XMOM(K,I)
   50 CONTINUE
      RMOM(K)=RMOM(K)/SUMLEN
   60 CONTINUE
      RMOM(1)=ONE
C
C         CALCULATE SUM OF SQUARES MATRIX
C
      IF (NSITES.LE.3) THEN
        DO 65 I=1,NSITES
        D(I)=ONE
   65   CONTINUE
      ELSE
        SUM2=ZERO
        SUM3=ZERO
        SUM4=ZERO
        DO 70 I=1,NSITES
        SUM2=SUM2+XMOM(2,I)
        SUM3=SUM3+XMOM(3,I)
        SUM4=SUM4+XMOM(4,I)
   70   CONTINUE
        SUM2=SUM2/NSITES
        SUM3=SUM3/NSITES
        SUM4=SUM4/NSITES
        DO 80 I=1,NSITES
        WORK(I,1)=XMOM(2,I)-SUM2
        WORK(I,2)=XMOM(3,I)-SUM3
        WORK(I,3)=XMOM(4,I)-SUM4
   80   CONTINUE
        DO 101 J=1,3
        DO 100 K=J,3
        SMAT(J,K)=ZERO
        DO 90 I=1,NSITES
        SMAT(J,K)=SMAT(J,K)+WORK(I,J)*WORK(I,K)
   90   CONTINUE
  100   CONTINUE
  101   CONTINUE
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
  110   CONTINUE
        SMAT(2,1)=SMAT(1,2)
        SMAT(3,1)=SMAT(1,3)
        SMAT(3,2)=SMAT(2,3)
C
C         CALCULATE DISCORDANCY MEASURES (D STATISTICS)
C
        FACTOR=NSITES/THREE
        DO 130 I=1,NSITES
        DO 121 J=1,3
        DO 120 K=1,3
        D(I)=D(I)+WORK(I,J)*WORK(I,K)*SMAT(J,K)
  120   CONTINUE
  121   CONTINUE
        D(I)=D(I)*FACTOR
        WORK(I,1)=D(I)
  130   CONTINUE
      END IF
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
      IF (IFAIL.NE.0) THEN
        CALL PELGLO(RMOM,RPARA,IFAIL)
        RPARA(4)=-ONE
      END IF
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
      CALL SAMLM(X,NREC,TMOM,4,1,1)
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
      V2=V2+NREC*SQRT(TEMP2+TEMP3)
      V3=V3+NREC*SQRT(TEMP3+TEMP4)
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
      V2=V2+NREC*SQRT(TEMP2+TEMP3)
      V3=V3+NREC*SQRT(TEMP3+TEMP4)
  225 CONTINUE
      VOBS(1)=SQRT(V1/SUMLEN)
      VOBS(2)=V2/SUMLEN
      VOBS(3)=V3/SUMLEN
C
C         CALCULATE HETEROGENEITY MEASURES (H STATISTICS)
C
      DO 230 J=1,3
      VBAR(J)=VBAR(J)/NSIM
      VSD(J)=SQRT((VSD(J)-NSIM*VBAR(J)**2)/(NSIM-ONE))
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
      T4SD=SQRT((T4SD-NSIM*T4BAR**2)/(NSIM-ONE))
      DO 240 IDIST=1,5
      Z(IDIST)=(T4FIT(IDIST)+T4BAR-TWO*RMOM(4))/T4SD
  240 CONTINUE
C
      RETURN
C
      END
