CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 22 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE RADIAL INTEGRALS INVOLVED IN THE TWO-   * C
C * CENTER ONE-ELECTRON INTEGRALS OTHER THAN OVERLAP INTEGRALS.         * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LMDL12 : IS THE HIGHEST ORDER OF THE MODIFIED BESSEL FUNCTIONS    * C
C *                                                                     * C
C *   N1,L1, N2,L2, ZETA1, ZETA2 : QUANTUM NUMBERS AND SLATER EXPONENTS * C
C *    OF THE ORBITALS BEING OVERLAPED                                  * C
C *                                                                     * C
C *   AB : IS THE DISTANCE SEPARATING THE TWO-CENTERS                   * C
C *                                                                     * C
C *   RLEG01(1 .. L1_3), RLEG02(1 .. L4_6) : CONTAINS THE ROOTS OF THE  * C
C *    GAUSS-LEGENDRE QUADRATURE $z_i$ USED TO COMPUTE THE INTEGRAND    * C
C *                                                                     * C
C *   BESF01 : IS AN ARRAY CONTAINING THE FIRST VALUES THAT SHOULD BE   * C
C *    USED IN THE ASCENDING OR DESCENDING RECURRENCE RELATIONS OF THE  * C
C *    MODIFIED BESSL FUNCTIONS $I_{\l+\hl} (z)$                        * C
C *                                                                     * C
C *   X4NL : IS AN ARRAY CONTAINING THE POWERS OF THE GAUSS-LEGENDRE    * C
C *    ROOTS AFTER A SUITABLE COORDINATE TRANSFORMATION                 * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   RADAR : IS THE ARRAY CONTAINING THE VALUES OF THE RADIAL INTEGRAL * C
C *    WHICH DEPENDS ON THE QUANTUM NUMBERS                             * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RTCOEI(LMDL12, N1,L1, N2,L2, ZETA1,ZETA2, AB,
     $                              RLEG01, RLEG02, BESF01, X4NL, RADAR)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION RLEG01(*), RLEG02(*), BESF01(*), X4NL(*), RADAR(*)
      COMMON/STOP08/ZETA1P, N1P, LP2, LP12, ISTR
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      COMMON/BESHER/HEREX(LTOT6*LMHERM), BESS01(LTOT6*(LDEV+1))
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))
      EXTERNAL GLROOT, RAD10


C.... PASSAGE TO A COMMON

      N1P    = N1
      ZETA1P = ZETA1

C.... COLLECTING THE ROOTS OF THE GAUSS-LAGUERRE QUADRATURE USED WITHIN THE
C.... RANGE [AB, +\INFINITY)

      ISP  = LEG15
c      XXXX = DGLGQ(GLROOT, AB, ZETA1+ZETA2, LAG10)
      XXXX = DGLGQ(GLROOT, AB, ZETA1+ZETA2, 8)

C.... COMPUTATION OF THE COEFFICIENTS INVOLVED IN THE GENERALIZED GEGENBAUER
C.... ADDITION THEOREM FOR THE ROOTS OF THE GAUSS-LAGUERRE QUADRATURE

      IN = LEA10-LEG15+1

      CALL HERMIT(ZETA2, AB, N2-L2, RLEG01,L1_3, RLEG02,L4_6,
     $                     RLEAG0,IN,LEA10, BESF01, X4NL, HEREX, BESS01)

      DO 5 LMD=0, LMDL12-1
       IND1 = 1
      DO 5 I=IN, LEA10
       R    = RLEAG0(I)
       IND2 = LMD+1
       CALL BAD01(IND1, IND2, ZETA2, AB, R, N2-L2, LMD, XYINT)
       INDEX         = LMD*LEA10 + I
       GEGIN0(INDEX) = XYINT
 5    CONTINUE

C.... COMPUTATION OF THE RADIAL INTEGRALS

      INDEX = 1

      DO 10 LP2=0, L2
      DO 10 LP12=IABS(L1-LP2), L1+LP2, 2
       ISTR         = 0
       XINT1        = DGLNQ(RAD10, 0.D0, AB, LEG15)
       ISTR         = ISTR + LEG15
c       XINT2        = DGLGQ(RAD10, AB, ZETA1+ZETA2, LAG10)
       XINT2        = DGLGQ(RAD10, AB, ZETA1+ZETA2, 8)
       RADAR(INDEX) = XINT1 + XINT2
       INDEX        = INDEX + 1
 10   CONTINUE

      RETURN
      END




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C * THIS SUBROUTINE RETURNS THE VALUE OF THE RADIAL INTEGRAL WITHIN THE * C
C * RANGE [0, AB]                                                       * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RAD10(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/STOP08/ZETA1, N1, LP2, LP12, ISTR
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))

      INDO = LP12*LEA10 + ISTR

      DO 5 I=1, N
       R = T(I)
       Y(I) = 1.D0/DSQRT(R) * R**(N1+LP2+1) *
     $                                   DEXP(-ZETA1*R) * GEGIN0(INDO+I)
 5    CONTINUE

      RETURN
      END
