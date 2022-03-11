CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 20 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE COMPUTES THE RADIAL INTEGRALS INVOLVED IN THE THREE-   * C
C * CENTER NUCLEAR ATTRACTION INTEGRALS. THESE VALUES ARE STORED IN     * C
C * THE ONE-DIMENSIONAL ARRAY RADAR.                                    * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LMDL12 : IS THE HIGHEST ORDER OF THE TERMS OCCURING IN THE SERIES * C
C *    EXPANSION                                                        * C
C *                                                                     * C
C *   N1,L1, N2,L2, ZETA1,ZETA2 : CONSTANTS CHRACTERIZING THE ORBITALS  * C
C *    INVOLVED IN THE INTEGRAL                                         * C
C *                                                                     * C
C *   AB, AC : ARE THE DISTANCES SUCH THAT D(C1, X2) AND D(C1, C3)      * C
C *                                                                     * C
C *   RMIN, RMAX : MIN(AB, AC) AND MAX(AB, AC) RESPECTIVELY             * C
C *                                                                     * C
C *   RLEG01, RLEG02, RLEAG0 : ARE ONE-DIMENSIONAL ARRAYS CONTAINING THE* C
C *    ROOTS OF THE GAUSS-LEGENDRE AND GAUSS-LAGUERRE QUADRATURES WITHIN* C
C *    THE RANGES [0, MIN(AB,AC)], [MIN(AB,AC), MAX(AB,AC)] AND         * C
C *    [MAX(AB,AC), +\INFTY[                                            * C
C *                                                                     * C
C *   BESF01 : IS A ONE-DIMENSIONAL ARRAY CONTAINING THE CALUES OF THE  * C
C *    BESSEL FUNCTIONS FOR EACH ROOT COLLECTED IN THE ABOVE ARRAYS     * C
C *                                                                     * C
C *   X4NL : IS AN AUXILIARY ONE-DIMENSIONAL ARRAY CONTAINING THE FOURTH* C
C *    POWERS OF A CERTAIN VARIABLE                                     * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   RADAR : IS A ONE-DIMENSIONAL ARRAY CONTAINING THE RADIAL INTEGRALS* C
C *    INVOLVED IN THE THREE-CENTER NUCLEAR ATTRACTION INTEGRAL         * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RTCNAI(LMDL12, N1,L1, N2,L2, ZETA1,ZETA2, AB,AC,
     $                                                 RMIN,RMAX, RADAR)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION RADAR(*)
      DIMENSION BESI(0:100), BESK(0:100)
      COMMON/STOP01/ZETA1P, ACP, N1P, L, LP2, LMD1, ISTR
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      COMMON/BESHER/HEREX(LTOT6*LMHERM), BESS01(LTOT6*(LDEV+1))
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))
c      COMMON/BESS0/BESF01(2*LTOT6)

      EXTERNAL GLROOT, RAD01, RAD02


C.... PASSAGE TO A COMMON

      N1P    = N1
      ZETA1P = ZETA1
      ACP    = AC

C.... COLLECTING THE ROOTS OF THE GAUSS-LAGUERRE QUADRATURE USED WITHIN
C.... THE RANGE [RMAX, +\INFTY)

      ISP  = LEG01+LEG02
      XXXX = DGLGQ(GLROOT, RMAX, ZETA1+ZETA2, LAG01)

C.... COMPUTATION OF THE COEFFICIENTS INVOLVED IN THE GENERALIZED GEGENBAUER
C.... ADDITION THEOREM FOR THE ROOTS OF THE GAUSS-LAGUERRE QUADRATURE

      IN  = LEA01-LAG01+1
      NL2 = N2-L2

      DO 5 I=IN, LEA01
       R1 = RLEAG0(I)

       RM = ZETA2*DMIN1(R1, AB)
       RP = ZETA2*DMAX1(R1, AB)

       FRAC1 = RM/RP
       FRAC2 = FRAC1*FRAC1

       CALL BESSI(-NL2, LMDL12+NL2, RM, BESI)
       CALL BESSK(LMDL12+NL2, RP, BESK)
c       CALL BESSK(RP, LMDL12+NL2, I-1, BESF01, BESK)

       CSTE = (0.5D0*RP/ZETA2)**NL2

      DO 5 LMD=0, LMDL12-1
       POWR = 1.D0
       SUMP = BESPRO(RM, NL2, LMD, BESI, BESK, 0, POWR, FRAC2)

       POWR = FRAC1
       SUMN = BESPRO(RM, NL2, LMD, BESI, BESK, 1, POWR, FRAC2)

       INDEX         = LMD*LEA01 + I
       GEGIN0(INDEX) = CSTE*(SUMP - SUMN)

 5    CONTINUE

C.... COMPUTATION OF THE RADIAL INTEGRALS

      INDEX = 1

      DO 10 L=0, LMDMAX
      DO 10 LP2=0, L2
      DO 10 LP12=IABS(L1-LP2), L1+LP2, 2
      DO 10 LMD1=IABS(L-LP12), L+LP12, 2
       ISTR  = 0
       XINT1 = DGLNQ(RAD01, 0.D0, RMIN, LEG01)
       ISTR  = ISTR + LEG01
       IF(AC .GT. RMIN)THEN
	XINT2 = DGLNQ(RAD01, RMIN, RMAX, LEG02)
       ELSE
	XINT2 = DGLNQ(RAD02, RMIN, RMAX, LEG02)
       ENDIF
       ISTR  = ISTR + LEG02
       XINT3 = DGLGQ(RAD02, RMAX, ZETA1+ZETA2, LAG01)
       RADAR(INDEX) = XINT1+XINT2+XINT3

       INDEX = INDEX + 1
 10   CONTINUE

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 20 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE VALUE OF THE RADIAL INTEGRAL WITHIN THE * C
C * RANGE [0, MIN(AB, AC)] AND [MIN(AB, AC), MAX(AB, AC)] IN THE  CASE  * C
C * MIN(AC, AB) = AC.                                                   * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RAD01(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/STOP01/ZETA1, AC, N1, L, LP2, LMD1, ISTR
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))

      INDO = LMD1*LEA01 + ISTR

      DO 5 I=1, N
       R = T(I)
       Y(I) = 1.D0/DSQRT(R) * R**(N1+LP2) * (R/AC)**(L+1) *
     $                                   DEXP(-ZETA1*R) * GEGIN0(INDO+I)
 5    CONTINUE

      RETURN
      END




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 20 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE VALUE OF THE RADIAL INTEGRAL WITHIN     * C
C * THE RANGE [MAX(AB, AC), +\INFTY) AND [MIN(AB, AC), MAX(AB, AC] IF   * C
C * AC = MAX(AC, AB).                                                   * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RAD02(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/STOP01/ZETA1, AC, N1, L, LP2, LMD1, ISTR
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))

      INDO = LMD1*LEA01 + ISTR

      DO 5 I=1, N
       R = T(I)
       Y(I) = 1.D0/DSQRT(R) * R**(N1+LP2) * (AC/R)**L *
     $                                   DEXP(-ZETA1*R) * GEGIN0(INDO+I)
 5    CONTINUE

      RETURN
      END

