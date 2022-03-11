CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 19 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE RADIAL PART OF THE T.C.C.D.P. FOR A     * C
C * GIVEN LMD, L, N AND ZETA.                                           * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   INDEX : RUNNING INDEX USED TO STORE THE NUMERICAL VALUES OF THE   * C
C *    RADIAL INTEGRALS INVOLVED IN THE T.C.C.D.P.                      * C
C *   ZETA1,ZETA2 : SLATER EXPONENTS USED FOR THE G-LAGUERRE INTEGRATION* C
C *    RULE                                                             * C
C *   AB : DISTANCE SEPARTING THE CENTERS OF THE T.C.C.D.P.             * C
C *   RMIN, RMAX : MIN(AC, AD) AND MAX(AC, AD)                          * C
C *   RLEAG0 : ONE-DIMENSIONAL ARRAY CONTAINING THE ROOTS OF THE OUTER  * C
C *    RADIAL INTEGRAL                                                  * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   VL : ONE-DIMENSIONAL ARRAY CONTAINING THE NUMERICAL VALUES OF THE * C
C *    RADIAL INTEGRALS INVOLVED IN THE T.C.C.D.P.                      * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE OUADHIAS(INDEX, ZETA1, ZETA2, AB, RLEAG0, VL)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION RLEAG0(*), VL(*)
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))
      COMMON/EXCH04/ZETA1P, R2, N1P, LP2, L, LMD2, LMDL12, IO, ISTR

      EXTERNAL TCCDP6, TCCDP7


      DO 5 IO=1, LEA11
       R2        = RLEAG0(IO)
       IF(R2 .LT. AB)THEN
	ISTR      = 0
	XINT1     = DGLNQ(TCCDP6, 0.D0, R2, LEG18)
	ISTR      = ISTR + LEG18
	XINT2     = DGLNQ(TCCDP7, R2, AB, LEG19)
	ISTR      = ISTR + LEG19
	XINT3     = DGLGQ(TCCDP7, AB, ZETA1+ZETA2, LAG12)

       ELSE
	ISTR      = 0
	XINT1     = DGLNQ(TCCDP6, 0.D0, AB, LEG18)
	ISTR      = ISTR + LEG18
	XINT2     = DGLNQ(TCCDP6, AB, R2, LEG19)
	ISTR      = ISTR + LEG19
	XINT3     = DGLGQ(TCCDP7, R2, ZETA1+ZETA2, LAG12)
       ENDIF

       VL(INDEX) = XINT1 + XINT2 + XINT3

       INDEX     = INDEX + 1

 5    CONTINUE

      RETURN
      END
