CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 20 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS A SET OF ONE-CENTER OVERLAP INTEGRALS OVER A   * C
C * SLATER ORBITAL BASIS SET.                                           * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   NORB1 : NUMBER OF SLATER ORBITALS ATTACHED TO THE ATOM THE        * C
C *    ORBITALS OF WHICH ARE OVERLAPING                                 * C
C *                                                                     * C
C *   NLMA, ZETAA : ARE TWO ARRAYS CONTAINING THE PARAMETERS OF THE     * C
C *    ORBITALS, i.e. QUANTUM NUMBERS ARE SLATER EXPONENTS              * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   SOVL : IS A TWO-DIMENSIONAL ARRAY CONTAINING THE REAL OVERLAPS    * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE OCOVL(NORB1, NLMA, ZETAA, XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      INCLUDE "../lib95/CONST.INCL"

      DIMENSION NLMA(*), ZETAA(*), XCORE(N_ORB, N_ORB)

      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:100)


C.... COMBINING THE GIVEN SET OF ATOMIC ORBITALS

      DO 5 NAT1=1, NORB1
       IST   = 5*(NAT1-1)
       LINE  = NLMA(IST+2)
       N1    = NLMA(IST+3)
       L1    = NLMA(IST+4)
       M1    = NLMA(IST+5)
       ZETA1 = ZETAA(NAT1)

      DO 5 NAT2=NAT1, NORB1
       JST   = 5*(NAT2-1)
       ICOL  = NLMA(JST+2)
       N2    = NLMA(JST+3)
       L2    = NLMA(JST+4)
       M2    = NLMA(JST+5)
       ZETA2 = ZETAA(NAT2)

       IF(L1.EQ.L2 .AND. M1.EQ.M2)THEN
	A1 = (2.d0*ZETA1)**N1 * DSQRT((2.D0*ZETA1)/FACT(2*N1))
	A2 = (2.d0*ZETA2)**N2 * DSQRT((2.D0*ZETA2)/FACT(2*N2))
	OCOEV = A1*A2 * FACT(N1+N2)/(ZETA1+ZETA2)**(N1+N2+1)

ccc     WRITE(*, 1)LINE, ICOL, N1,L1,M1, ZETA1, N2,L2,M2, ZETA2
ccc     WRITE(*, *)OCOEV

       ELSE
	OCOEV = 0.D0
       ENDIF

C.... SINCE THE INTEGRALS <Xa | Xa> ARE ALL REAL, THEREFORE FOR OUT-OFF
C.... DIAGONAL ELEMENTS THE IMAGINARY PART IS NULL

       XCORE(LINE, ICOL) = OCOEV

       IF(LINE .NE. ICOL)THEN
	XCORE(ICOL, LINE) = 0.D0
       ENDIF

 5    CONTINUE

c1    FORMAT(2Z2, 2X, 3(I3, 2X), D12.6, 3(I3, 2X), D12.6)

      RETURN
      END
