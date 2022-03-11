CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 19 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS THE DERIVATIVES OF THE PRODUCT OF THE MODIFIED * C
C * BESSEL FUNCTIONS INVOLVED IN THE GENERALIZED GEGENBAUER ADDITION    * C
C * THEOREM.                                                            * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LMDL12 :IS THE HIGHEST ORDER OF THE BESSEL FUNCTIONS BEING DERIVED* C
C *   NL2 : (N2-L2), IS THE ORDER OF THE DERIVATIVE BEING COMPUTED      * C
C *   ZETA2 : SLATER EXPONENT INVOLVED IN THE STO BEING TRANSLATED      * C
C *   AB : DISTANCE SEPARATING THE CENTERS OF THE T.C.C.D.P.            * C
C *   RLEAG0 : ONE-DIMENSIONAL ARRAY CONTAINING THE ROOTS OF THE OUTER  * C
C *    RADIAL INTEGRAL, I.E. OVER R2                                    * C
C *   ROOTS : ONE-DIMENSIONAL ARRAY CONTAINING THE ROOTS OF THE INNER   * C
C *    RADIAL INTEGRAL, I.E. OVER R1                                    * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   GEGIN0 : ONE-DIMENSIONAL ARRAY WHERE THE GEGENBAUERS' COEFFICIENTS* C
C *    ARE STORED                                                       * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GEGEN9(LMDL12, NL2, ZETA2, AB, RLEAG0, ROOTS, GEGIN0)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION RLEAG0(*), ROOTS(*), GEGIN0(*)
      DIMENSION BESI(0:100), BESK(0:100)

      DO 5 I=1, LEA11
       ISIT = (I-1)*LEA12

      DO 5 J=1, LEA12-LAG12
       R1 = ROOTS(ISIT+J)

       RM = ZETA2 * DMIN1(R1, AB)
       RP = ZETA2 * DMAX1(R1, AB)

       FRAC1 = RM/RP
       FRAC2 = FRAC1*FRAC1

       CALL BESSI(-NL2, LMDL12+NL2, RM, BESI)
       CALL BESSK(LMDL12+NL2, RP, BESK)

       CSTE = (0.5D0*RP/ZETA2)**NL2

       INDEX           = (I-1)*LEA12*LMDL12 + LMD*LEA12 + J
       GEGIN0(INDEX)   = CSTE*(SUMP - SUMN)

      DO 5 LMD=0, LMDL12-1

       POWR = 1.D0
       SUMP = BESPRO(RM, NL2, LMD, BESI, BESK, 0, POWR, FRAC2)

       POWR = FRAC1
       SUMN = BESPRO(RM, NL2, LMD, BESI, BESK, 1, POWR, FRAC2)

       INDEX           = (I-1)*LEA12*LMDL12 + LMD*LEA12 + J
       GEGIN0(INDEX)   = CSTE*(SUMP - SUMN)

 5    CONTINUE

      RETURN
      END
