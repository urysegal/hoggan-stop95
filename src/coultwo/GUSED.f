CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE COEFFICIENTS NEEDED IN THE SOLID        * C
C * SPHERICAL HARMONICS ADDITION AND MULTIPLICATION THEOREMS.           * C
C *       EXCLUSIVELY FOR THE TWO-CENTER COULOMB INTEGRALS.             * C
C *             CONJ{Y(L3,M3) [R-AB]} Y(L4,M4) [R-AB]                   * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   DFAC : DOUBLE FACTORIALS, (2L + 1)!!                              * C
C *                                                                     * C
C *   L3,M3, L4,M4 : PARAMETERS OF THE SOLID SPHERICAL HARMONICS TO BE  * C
C *    EXPANDED AND MULTIPLIED                                          * C
C *                                                                     * C
C *   AB : DISTANCE SEPARATING THE CENTERS                              * C
C *                                                                     * C
C *   YLMAB : ONE-DIMENSIONAL ARRAY CONTAINING THE REAL PART OF THE     * C
C *    SPHERICAL HARMONICS FOR (L=0..20, M=0..20)                       * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   ARGNT1, ARGNT2 : ONE-DIMENSIONAL ARRAYS WHERE THE COEFFICIENTS    * C
C *    OF THE EXPANSION AND THE MULTIPLICATION ARE STORED               * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GUSED(DFAC, L3,M3, L4,M4, AB,YLMAB, ARGNT1, ARGNT2)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION DFAC(0:*), YLMAB(0:*),
     $                       ARGNT1(*),ARGNT2(*),ARGNTX(100),ARGNTY(100)


      CALL GAUNT(L3,M3, L4,M4, LMIN,LMAX,M34, NGAUNT, ARGNTX)
      KG1 = 1
      K1  = 1
      K2  = 1

      DO 5 L34=IABS(L3-L4), L3+L4, 2
       IF(L34 .GE. LMIN)THEN
	ARGNT1(K1) = DFAC(L34) * ARGNTX(KG1) * (-AB)**L34
	KG1        = KG1 + 1
       ELSE
	ARGNT1(K1) = 0.D0
       ENDIF
       K1 = K1 + 1

      ABLP = -AB
      DO 5 LP34=0, L34
       ABLP = -ABLP/AB
      DO 5 MP34=MAX(-LP34, -L34+LP34+M34), MIN(LP34, L34-LP34+M34)
       CALL GAUNT(LP34,MP34, L34,M34, LMIN1,LMAX1,M34P, NGAUNT, ARGNTY)
       NYL        = ((L34-LP34)*(L34-LP34+1))/2 + IABS(M34-MP34)
       MUPW       = (M34-MP34 - IABS(M34-MP34))/2
       ARGNT2(K2) = (-1)**MUPW *
     $             ABLP*YLMAB(NYL)*ARGNTY(1)/(DFAC(LP34)*DFAC(L34-LP34))
       K2        = K2 + 1
 5    CONTINUE

      RETURN
      END
