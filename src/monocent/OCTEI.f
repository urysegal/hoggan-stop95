CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            NOV 21 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS A SET OF ONE-CENTER TWO-ELECTRON INTEGRALS     * C
C * OVER A SLATER ORBITAL BASIS SET, USING A TWO-RANGE ONE-CENTER       * C
C * EXPANSION METHOD.                                                   * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LM123 : IS THE HIGHEST VALUE OF THE QUANTUM NUMBER L              * C
C *                                                                     * C
C *   NORB1 : NUMBER OF ORBITALS DESCRIBING THE CENTER                  * C
C *                                                                     * C
C *   NLMA, ZETAA : ARE TWO ARRAYS CONTAINING THE DESCRIPTION OF THE    * C
C *    ORBITALS, i.e. QUANTUM NUMBERS AND SLATER EXPONENTS              * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   FATIMA : IS ONE-DIMENSIONAL ARRAY CONTAINING THE TWO-ELECTRONS    * C
C *    INTEGRALS <I J | K L>. IN ORDER TO LOWER THE USED STORAGE WE HAVE* C
C *    SET I >= J AND K >= L, SO THE INDEX IS DETERMINED BY :           * C
C *                                                                     * C
C *     INDEX = NSUM*(NSUM - (N-L+1)*(N-L+2)/2) + (K-L)*NSUM +          * C
C *             NSYM - (N-J+1)*(N-J+2)/2 + (I-J+1)                      * C
C *                                                                     * C
C *    WHERE NSUM = N*(N+1)/2 WITH N : NUMBER OR ORBITALS               * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE OCTEI(NTORB,NSUM, LM123, NORB1, NLMA, ZETAA, FATIMA)
      IMPLICIT REAL*8 (A-H, O-Z)
      COMPLEX*16 OUIZA
      INCLUDE "../lib95/SIZE.INCL"
      INCLUDE "../lib95/CONST.INCL"

      DIMENSION NLMA(*), ZETAA(*), FATIMA(*)

      DIMENSION ARGNT1(L1MAX+L2MAX+1), ARGNT2(L3MAX+L4MAX+1)

      DIMENSION OUIZA(20,20,20,20)

      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:100)
      COMMON/GAUSOR/NORDER(0:99)
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      COMMON/STOP07/ZETA34, N34, L12, L, NSTRT
      COMMON/COMV00/VL(NORBM*NORBM*(L1MAX+L2MAX+1)*LEA09)
      COMMON/INTEGN/NINT
      EXTERNAL GLROOT, RAD09


C.... COMBINING THE GIVEN ATOMS

      DO 5 NAT4=1, NORB1
       LST   = 5*(NAT4-1)
       NO4   = NLMA(LST+2)
       N4    = NLMA(LST+3)
       L4    = NLMA(LST+4)
       M4    = NLMA(LST+5)
       ZETA4 = ZETAA(NAT4)

      DO 5 NAT3=NAT4, NORB1
       KST   = 5*(NAT3-1)
       NO3   = NLMA(KST+2)
       N3    = NLMA(KST+3)
       L3    = NLMA(KST+4)
       M3    = NLMA(KST+5)
       ZETA3 = ZETAA(NAT3)

       ZETA34 = ZETA3+ZETA4
       N34    = N3 + N4
       INDEX  = 0

      DO 5 NAT2=1, NORB1
       JST   = 5*(NAT2-1)
       NO2   = NLMA(JST+2)
       N2    = NLMA(JST+3)
       L2    = NLMA(JST+4)
       M2    = NLMA(JST+5)
       ZETA2 = ZETAA(NAT2)

      DO 5 NAT1=NAT2, NORB1
       IST   = 5*(NAT1-1)
       NO1   = NLMA(IST+2)
       N1    = NLMA(IST+3)
       L1    = NLMA(IST+4)
       M1    = NLMA(IST+5)
       ZETA1 = ZETAA(NAT1)

       N12    = N1+N2
       L12    = L1+L2
       ZETA12 = ZETA1+ZETA2

C.... COLLECTING THE ROOTS OF THE GAUSS-LAGUERRE QUADARTURE

       ISP = 0
       LAG09X = NORDER((N1+N2+N3+N4-1)/2 + 1)
       LEA09X = LAG09X
       XXXX = DGLGQ(GLROOT, 0.D0, ZETA12+ZETA34, LAG09X)

C.... COMPUTATION OF THE ONE-CENTER CHARGE DISTRIBUTION POTENTIAL

       INC = 0
      CALL OCCDP2(L1MAX+L2MAX,INDEX,N12,L12,ZETA12,RLEAG0,INC,LEA09X,VL)
       NSTRT = INDEX*LEA09X*(L1MAX+L2MAX+1)
       INDEX = INDEX + 1

C.... NORMALIZATION CONSTANTS

       A1 = (2.d0*ZETA1)**N1 * DSQRT((2.D0*ZETA1)/FACT(2*N1))
       A2 = (2.d0*ZETA2)**N2 * DSQRT((2.D0*ZETA2)/FACT(2*N2))
       A3 = (2.d0*ZETA3)**N3 * DSQRT((2.D0*ZETA3)/FACT(2*N3))
       A4 = (2.d0*ZETA4)**N4 * DSQRT((2.D0*ZETA4)/FACT(2*N4))

       CSTE = 4.D0 * PI * A1 * A2 * A3 * A4

       LSOM = L1+L2+L3+L4
       M13  = M1+M3
       M24  = M2+M4

       M14  = M1+M4
       M23  = M2+M3

       IF(2*(LSOM/2) .EQ. LSOM)THEN
	IF(M13 .EQ. M24)THEN
	 CALL GAUNT(L2,M2, L1,M1, LMIN1,LMAX1,M12, NGAUNT, ARGNT1)
	 CALL GAUNT(L3,M3, L4,M4, LMIN2,LMAX2,M34, NGAUNT, ARGNT2)

	 IF(LMIN1 .GE. LMIN2)THEN
	  K1 = (LMIN1-LMIN2)/2 + 1
	  K2 = 1
	 ELSE
	  K1 = 1
	  K2 = (LMIN2-LMIN1)/2 + 1
	 ENDIF

	 OCTES1 = 0.D0

	 DO 10 L=MAX(LMIN1, LMIN2), MIN(LMAX1, LMAX2), 2
	  AUX1 = ARGNT1(K1) * ARGNT2(K2) / DBLE(L+L+1)
	  K1   = K1 + 1
	  K2   = K2 + 1
	  XINT = DGLGQ(RAD09, 0.D0, ZETA12+ZETA34, LAG09X)
	  YINT = FTERM(N1+N2+L, ZETA12, N3+N4-L-1, ZETA34)
	  OCTES1 = OCTES1 + AUX1 * (XINT + YINT)
 10      CONTINUE

	 OCTES1 = CSTE * OCTES1

	 NINT = NINT + 1

	ELSE
	 OCTES1 = 0.D0
	ENDIF

	IF(M14 .EQ. M23)THEN
	 CALL GAUNT(L1,M1, L2,M2, LMIN1,LMAX1,M12, NGAUNT, ARGNT1)
	 CALL GAUNT(L3,M3, L4,M4, LMIN2,LMAX2,M34, NGAUNT, ARGNT2)

	 IF(LMIN1 .GE. LMIN2)THEN
	  K1 = (LMIN1-LMIN2)/2 + 1
	  K2 = 1
	 ELSE
	  K1 = 1
	  K2 = (LMIN2-LMIN1)/2 + 1
	 ENDIF

	 OCTES2 = 0.D0

	 DO 15 L=MAX(LMIN1, LMIN2), MIN(LMAX1, LMAX2), 2
	  AUX1 = ARGNT1(K1) * ARGNT2(K2) / DBLE(L+L+1)
	  K1   = K1 + 1
	  K2   = K2 + 1
	  YINT = FTERM(N1+N2+L, ZETA12, N3+N4-L-1, ZETA34)
	  XINT = DGLGQ(RAD09, 0.D0, ZETA12+ZETA34, LAG09X)
	  OCTES2 = OCTES2 + AUX1 * (XINT + YINT)
 15      CONTINUE

	 OCTES2 = CSTE * OCTES2

	 NINT = NINT + 1

	ELSE
	 OCTES2 = 0.D0
	ENDIF
       ELSE
	OCTES1 = 0.D0
	OCTES2 = 0.D0
       ENDIF

cccc   IF(OCTES1 .NE. 0.D0)THEN
       WRITE(*, 1)NO1,NO2,NO3,NO4, OCTES1, OCTES2
cccc   ENDIF

       OUIZA(NAT1,NAT2, NAT3,NAT4) = DCMPLX(OCTES1, 0.D0)
       OUIZA(NAT2,NAT1, NAT4,NAT3) = DCMPLX(OCTES1, 0.D0)

       OUIZA(NAT2,NAT1, NAT3,NAT4) = DCMPLX(OCTES2, 0.D0)
       OUIZA(NAT1,NAT2, NAT4,NAT3) = DCMPLX(OCTES2, 0.D0)

 5    CONTINUE

 1    FORMAT(4(I2, 1X), 2(D21.14,1X))

      CALL REAL01(NORB1, NLMA,
     $                                       NTORB, NSUM, OUIZA, FATIMA)

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C * THIS SUBROUTINE COMPUTES THE RADIAL INTEGRAL INVOLVED IN THE ONE-   * C
C * CENTER TWO-ELECTRON INTEGRALS.                                      * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RAD09(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/STOP07/ZETA34, N34, L12, L, NSTRT
      COMMON/COMV00/VL(NORBM*NORBM*(L1MAX+L2MAX+1)*LEA09)


      DO 5 I=1, N
       INDX = NSTRT + (I-1)*(L12+1) + (L+1)
       R2 = T(I)
       Y(I) = R2**N34 * DEXP(-ZETA34*R2) * VL(INDX)
 5    CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C * THIS FUNCTION COMPUTES THE FIRST TERM OF THE RADIAL INTEGRAL, WHICH * C
C * IS A LAPLACE TRANSFORM OF R^N EXP(-Z R).                            * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION FTERM(NL12, ZETA12, NL34, ZETA34)
      IMPLICIT REAL*8 (A-H, O-Z)

      U0 = 1.D0/ZETA12

      DO 5 I=0, NL12-1
       U0 = DBLE(I+1)/ZETA12 * U0
 5    CONTINUE

      V0 = 1.D0/ZETA34

      DO 10 I=0, NL34-1
       V0 = DBLE(I+1)/ZETA34 * V0
 10   CONTINUE

      FTERM = U0 * V0

      RETURN
      END
