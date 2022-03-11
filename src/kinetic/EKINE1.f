CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 23 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS A SET OF ONE-CENTER KINETIC ENERGY INTEGRALS   * C
C * OVER A SLATER ORBITAL BASIS SET.                                    * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   NORB1 : TOTAL NUMBER OF STOs ATTACHED TO THE CENTERS (1)          * C
C *                                                                     * C
C *   NLMA, ZETAA : ARRAYS CONTAINING THE DESCRIPTION OF THE ORBITALS,  * C
C *    i.e. QUANTUM NUMBERS AND SLATER EXPONENTS                        * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   XCORE : IS A TWO-DIMENSIONAL ARRAY CONTAINING THE CORE INTEGRALS  * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE EKINE1(NORB1, NLMA, ZETAA, XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"

      DIMENSION XCORE(N_ORB, N_ORB)
      DIMENSION NLMA(*), ZETAA(*)

      COMMON/FACT0/FACT(0:30)
      COMMON/INTEGN/NINT


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
	DN2 = DBLE(N2)
	DL2 = DBLE(L2)

	A1 = (2.D0*ZETA1)**N1 * DSQRT((2.D0*ZETA1)/FACT(2*N1))
	A2 = (2.D0*ZETA2)**N2 * DSQRT((2.D0*ZETA2)/FACT(2*N2))

	ZETAS  = ZETA1+ZETA2
	ZETAS1 = ZETAS**(N1+N2-1)
	ZETAS2 = ZETAS * ZETAS1
	ZETAS3 = ZETAS * ZETAS2

	CSTE = -0.5D0 * A1 * A2

	ZNRJ1 = CSTE *
     $          ((DN2+DL2)*(DN2-DL2-1.D0)*FACT(N1+N2-2)/ZETAS1 -
     $                     2.D0*DN2*ZETA2*FACT(N1+N2-1)/ZETAS2 +
     $                                   ZETA2*ZETA2*FACT(N1+N2)/ZETAS3)

	NINT       = NINT + 1
       ELSE
	ZNRJ1 = 0.D0
       ENDIF

ccc    WRITE(*, 1)LINE, ICOL, N1, L1, M1, ZETA1, N2, L2, M2, ZETA2
ccc    WRITE(*, 2)ZNRJ1

C.... THESE INTEGRALS ARE ULTIMATELY REAL, THEREFORE, XCORE(J, I) = 0

       XCORE(LINE, ICOL) = XCORE(LINE, ICOL) + ZNRJ1

 5    CONTINUE

c1    FORMAT(2Z2, 2X, 3(I3, 2X), D12.6, 3(I3, 2X), D12.6)
c2    FORMAT(D16.10)

      RETURN
      END
