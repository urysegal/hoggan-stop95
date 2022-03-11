CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 20 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE COEFFICIENTS OF THE GENERALIZED         * C
C * GEGENBAUER ADDITION THEOREM INVOLVED IN THE TRANSLATION OF A        * C
C * SLATER ORBITAL.                                                     * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LMDL12 : IS THE HIGHEST ORDER OF THE TERMS OCCURING IN THE SERIES * C
C *    EXPANSION                                                        * C
C *                                                                     * C
C *   NL2, ZETA2 : N2 - L2 AND THE CORRESPONDING SLATER EXPONENT        * C
C *                                                                     * C
C *   AB : DISTANCE SEPARATING THE TWO ORBITALS INVOLVED IN THE INTEGRAL* C
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
C *   GEGEN0 : IS A ONE-DIMENSIONAL ARRAY CONTAINING THE VALUES OF THE  * C
C *    GENERALIZED GEGENBAUER ADDITION THEOREM CORRESPONDING TO THE     * C
C *    AFORECOLLECTED ROOTS                                             * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GEGEN1(LMDL12, NL2, ZETA2, AB, RLEAG0, GEGIN0)

      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION  RLEAG0(*), GEGIN0(*)
      DIMENSION BESI(0:100), BESK(0:100)

      COMMON/BESHER/HEREX(LTOT6*LMHERM), BESS01(LTOT6*(LDEV+1))
c      COMMON/BESS0/BESF01(2*LTOT6)

C.... COMPUTATION OF THE MODIFIED BESSEL FUNCTION I AND THE HERMITE
C.... POLYNOMIALS FOR EVERY ROOT OF THE GAUSS-LEGENDRE QUADRATURE.

      DO 10 I=1, LEG01+LEG02
       R1 = RLEAG0(I)

       RM = ZETA2*DMIN1(R1, AB)
       RP = ZETA2*DMAX1(R1, AB)

       FRAC1 = RM/RP
       FRAC2 = FRAC1*FRAC1

       CALL BESSI(-NL2, LMDL12+NL2, RM, BESI)
       CALL BESSK(LMDL12+NL2, RP, BESK)
c       CALL BESSK(RP, LMDL12+NL2, I-1, BESF01, BESK)

       CSTE = (0.5D0*RP/ZETA2)**NL2

      DO 10 LMD=0, LMDL12-1
       POWR = 1.D0
       SUMP = BESPRO(RM, NL2, LMD, BESI, BESK, 0, POWR, FRAC2)

       POWR = FRAC1
       SUMN = BESPRO(RM, NL2, LMD, BESI, BESK, 1, POWR, FRAC2)

       INDEX         = LMD*LEA01 + I
       GEGIN0(INDEX) = CSTE*(SUMP - SUMN)

 10   CONTINUE


      RETURN
      END
