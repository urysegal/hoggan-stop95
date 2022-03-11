CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 22 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE COEFFICIENTS OF THE GENERALIZED         * C
C * GEGENBAUER ADDITION THEOREM INVOLVED IN THE TRANSLATION OF A SLATER * C
C * ORBITAL.                                                            * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LMDL12 : HIGHEST ORDER OF THE MODIFIED BESSEL FUNCTIONS           * C
C *                                                                     * C
C *   NL2 : IS THE ORDER OF THE DERIVATIVE WITH RESPECT TO ZETA2        * C
C *                                                                     * C
C *   ZETA2 : IS THE SLATER EXPONENT INVOLVED IN THE SLATER ORBITAL     * C
C *    BEING TRANSLATED                                                 * C
C *                                                                     * C
C *   AB : IS THE DISTANCE SEPARATING THE CENTER OF THE STO AND THE     * C
C *    ORINGIN                                                          * C
C *                                                                     * C
C *   RLEG01(1 .. L1_3), RLEG02(1 .. L4_6) : CONTAINS THE ROOTS OF THE  * C
C *    GAUSS-LEGENDRE QUADRATURE $z_i$ USED TO COMPUTE THE INTEGRAND    * C
C *                                                                     * C
C *   RLEAG0(1 .. LMHERM) : CONTAINS THE ROOTS OF THE VROOTS OF THE     * C
C *    QUADRATURES USED TO COMPUTE THE OUTER RADIAL INTEGRALS $r_I$     * C
C *                                                                     * C
C *   BESF01 : IS AN ARRAY CONTAINING THE FIRST VALUES THAT SHOULD BE   * C
C *    USED IN THE ASCENDING OR DESCENDING RECURRENCE RELATIONS OF THE  * C
C *    MODIFIED BESSL FUNCTIONS $I_{\l+\hl} (z)$                        * C
C *                                                                     * C
C *   X4NL : IS AN ARRAY CONTAINING THE POWERS OF THE GAUSS-LEGENDRE    * C
C *    ROOTS AFTER A SUITABLE COORDINATE TRANSFORMATION                 * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   GEGIN0 : IS AN ARRAY CONTAINING THE NUMERICAL VALUES OF THE       * C
C *    INTEGRAL TRANSFORM FOR EACH VALUE OF $r_i$                       * C
C *                                                                     * C
C * INTERMEDIARY INPUT/OUTPUT : HEREX AND BESS00 ARRAYS (SEE ABIN03.f)  * C
C *                                                                     * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GEGEN8(LMDL12, NL2, ZETA2, AB, RLEG01, RLEG02, RLEAG0,
     $                                             BESF01, X4NL, GEGIN0)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION RLEG01(*), RLEG02(*), RLEAG0(*),
     $                                     BESF01(*), X4NL(*), GEGIN0(*)
      COMMON/BESHER/HEREX(LTOT6*LMHERM), BESS01(LTOT6*(LDEV+1))

C.... COMPUTATION OF THE MODIFIED BESSEL FUNCTION I AND THE HERMITE
C.... POLYNOMIALS FOR EACH ROOT OF THE GAUSS-LEGENDRE QUADRATURE

      IN = 1

      CALL ABIN02(ZETA2, AB, NL2, RLEG01,L1_3, RLEG02,L4_6,
     $                     RLEAG0,IN,LEG15, BESF01, X4NL, HEREX, BESS01)

C.... COMPUTATION OF THE GEGENBAUER'S COEFFICIENTS FOR THE ROOTS
C.... CORRESPONDING TO THE RANGE [0, AB]

      DO 5 LMD=0, LMDL12-1
       IND1 = 1
      DO 5 I=1, LEG15
       R    = RLEAG0(I)
       IND2 = LMD+1
       CALL BAD01(IND1, IND2, ZETA2, AB, R, NL2, LMD, XYINT)
       INDEX         = LMD*LEA10 + I
       GEGIN0(INDEX) = XYINT
 5    CONTINUE


      RETURN
      END
