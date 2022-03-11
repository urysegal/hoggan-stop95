CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS A TWO-CENTER COULOMB INTEGRALS OVER A SLATER   * C
C * A SLATER ORBITAL BASIS SET USING A TWO-RANGE ON-CENTER EXPANSION    * C
C * METHOD.                                                             * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   N1,L1,M1, N2,L2,M2, N3,L3,M3, N4,L4,M4 : QUANTUM NUMBERS          * C
C *                                                                     * C
C *   ZETA1, ZETA2, ZETA3, ZETA4 : SLATER EXPONENETS                    * C
C *                                                                     * C
C *   AB, YLMAB, PHIAB : DESCRIPTION OF THE VECTOR AB                   * C
C *                                                                     * C
C *   ARGNT1, ARGNT2 : COEFFICIENTS OF THE MULTIPLICATION AND THE       * C
C *    ADDITION THEOREM OF THE SOLID SPHERICAL HARMONICS                * C
C *                                                                     * C
C *   RADAR : VALUES OF THE RADIAL INTEGRAL                             * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   TCCV, TCCV1 : TWO-CENTER COULOMB INTEGRALS : <A A' | B B'> AND    * C
C *    <A' A | B B'>                                                    * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TCCINS(N1,L1,M1, N2,L2,M2, N3,L3,M3, N4,L4,M4, ZETA1,
     $    ZETA2,ZETA3,ZETA4, AB,YLMAB,PHIAB, ARGNT1,ARGNT2, RADAR,
     $                                                      TCCV, TCCV1)

      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      INCLUDE "../lib95/CONST.INCL"
      DIMENSION YLMAB(0:*), ARGNT1(*), ARGNT2(*), RADAR(*),
     $       ARGNTA(L1MAX+1), ARGNTB(L1MAX+L2MAX+L3MAX+L4MAX+1), INX(7),
     $       ARGNTC(L1MAX+L2MAX+L3MAX+L4MAX+1)

      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:100)


C.... THE NORMALIZATION CONSTANTS

       A1 = (2.D0*ZETA1)**N1 * DSQRT((2.D0*ZETA1)/FACT(2*N1))
       A2 = (2.D0*ZETA2)**N2 * DSQRT((2.D0*ZETA2)/FACT(2*N2))
       A3 = (2.D0*ZETA3)**N3 * DSQRT((2.D0*ZETA3)/FACT(2*N3))
       A4 = (2.D0*ZETA4)**N4 * DSQRT((2.D0*ZETA4)/FACT(2*N4))

       CSTE = 64.D0 * PI3 * A1 * A2 * A3 * A4

C.... INITIALIZATION OF THE T.C.C.I.

       TCCV  = 0.D0
       TCCV1 = 0.d0

C.... INITIALIZATION OF THE INDEXING ARRAY

       CALL RAZ1(INX, 1, 4)

       CALL GAUNT(L2,M2, L1,M1, LMIN1,LMAX1,M12, NGAUNT,ARGNTA)

       KL = 1

       DO 10 L=IABS(L1-L2), L1+L2, 2
	INX(1) = INX(1) + INX(2) + INX(3) + INX(4)
	CALL RAZ1(INX, 2, 4)
	IF(L .GE. LMIN1)THEN
	 AUX1 = ARGNTA(KL) / DBLE(L+L+1)
	 KL   = KL + 1
	ELSE
	 AUX1 = 0.D0
	ENDIF

       KL34  = 1
       KLP34 = 1

       DO 10 L34=IABS(L3-L4), L3+L4, 2
	INX(2) = INX(2) + INX(3) + INX(4)
	CALL RAZ1(INX, 3, 4)
	AUX2   = ARGNT1(KL34)
	KL34   = KL34 + 1

       DO 10 LP34=0, L34
	INX(3) = INX(3) + INX(4)
	MAXMP  = MAX(-LP34, -L34+LP34+M4-M3)
	MINMP  = MIN( LP34,  L34-LP34+M4-M3)
	IF(MINMP .LT. MAXMP)THEN
	 INX(4) = MIN(L, LP34) + 1
	ENDIF
       DO 10 MP34=MAXMP, MINMP
	AUX3  = ARGNT2(KLP34)
	KLP34 = KLP34 + 1
	CALL RAZ1(INX, 4, 4)
	IF(AUX1 .NE. 0.D0)THEN
	 CALL GAUNT(LP34,MP34, L,M1-M2, LMIN2,LMAX2,MU, NGAUNT, ARGNTB)
	 CALL GAUNT(LP34,MP34, L,M2-M1, LMIN3,LMAX3,MU, NGAUNT, ARGNTC)
	ENDIF

	KO = 1
	KX = 1

       DO 10 LMD=IABS(L-LP34), L+LP34, 2
	INX(4) = INX(4) + 1
	INDIC  = INX(1) + INX(2) + INX(3) + INX(4)
	IF(LMD .GE. LMIN2)THEN
	 NYL  = (LMD*(LMD+1))/2 + IABS((M1-M2)-MP34)
	 MUPW = ((M1-M2)-MP34 - IABS((M1-M2)-MP34))/2
	 AUX4 = (-1)**MUPW * ARGNTB(KO) * YLMAB(NYL)
	 KO   = KO + 1
	ELSE
	 AUX4 = 0.D0
	ENDIF

	IF(LMD .GE. LMIN3)THEN
	 NYL  = (LMD*(LMD+1))/2 + IABS((M2-M1)-MP34)
	 MUPW = ((M2-M1)-MP34 - IABS((M2-M1)-MP34))/2
	 AUX5 = (-1)**MUPW * ARGNTC(KX) * YLMAB(NYL)
	 KX   = KX + 1
	ELSE
	 AUX5 = 0.D0
	ENDIF

	 TCCV1 =  TCCV1 + AUX1*AUX2*AUX3*AUX5*RADAR(INDIC)
	 TCCV  =  TCCV  + AUX1*AUX2*AUX3*AUX4*RADAR(INDIC)

 10    CONTINUE

	TCCV  =                 CSTE/DSQRT(AB) * TCCV
	TCCV1 = (-1)**(M1+M2) * CSTE/DSQRT(AB) * TCCV1

C.... INCREMENTING THE NUMBER OF TOTAL INTEGRALS

       NINT = NINT + 1

      RETURN
      END
