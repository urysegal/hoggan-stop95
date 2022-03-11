CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS THE OVERLAP INTEGRALS GIVEN A SET OF  SLATER   * C
C * ATOMIC ORBITALS.                                                    * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LM123 : IS THE HIGHEST VALUE OF THE QUANTUM NUMBER L              * C
C *                                                                     * C
C *   NORB1 AND NORB2 : TOTAL NUMBER OF STOs ATTACHED TO THE CENTERS (1)* C
C *    AND (2)                                                          * C
C *                                                                     * C
C *   NLMA, NLMB, ZETAA AND ZETAB : ARRAYS CONTAINING THE DESCRIPTION   * C
C *    OF THE ORBITALS, i.e. QUANTUM NUMBERS AND SLATER EXPONENTS       * C
C *                                                                     * C
C *   AB, THETAB AND PHIAB : SPHERICAL COORDINATES OF THE VECTOR        * C
C *    SEPARATING THE CENTERS                                           * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   SOVL : IS A TWO-DIMENSIONAL ARRAY CONTAINING THE VALUES OF THE    * C
C *    OVERLAP INTEGRALS. THESE INTEGRALS HAVE BEEN TRANSFORMED IN      * C
C *    SUCH A WAY AS TO YIELDS ONLY REAL VALUES                         * C
C *                                                                     * C
C * DESCRIPTION OF THE COMMONS :                                        * C
C *   GLE01, GLE02 : CONTAINS THE ARRAYS WHICH CONTAIN THE ROOTS OF THE * C
C *    G-LEGENDRE QUADRATURES USED TO COMPUTE THE INTEGRAND OF OUR      * C
C *    INTEGRAL REPRESENTATION.                                         * C
C *                                                                     * C
C *   BESS0 : CONTAINS THE ARRAY BESF01 WHICH CONTAINS THE FIRST VALUES * C
C *    INVOLVED IN THE ASCENDING OR DESCENDING RECURRENCE RELATION OF   * C
C *    THE FUNCTION $\exp(-z) I_{\l + \hl} (z)$                         * C
C *                                                                     * C
C *   X4POW : CONTAINS THE ARRAY X4NL WHICH CONTAINS THE POWERS OF THE  * C
C *    G-LEGENDRE ROOTS $x^0, x^1, .. x^10$                             * C
C *                                                                     * C
C *   GLAE0 : CONTAINS THE ARRAY RLEAG0 WHICH CONTAINS THE ROOTS OF THE * C
C *    GAUSS QUADRATURES USED TO COMPTED THE WHOLE RADIAL INTEGRAL AND  * C
C *    THE INDEX ISP WHICH IS THE RANK FROM WHICH THE ROOTS ARE STORED  * C
C *                                                                     * C
C *   COMG00 : CONTAINS THE ARRAY GEGIN0 WHICH CONTAINS THE DERIVATIVES * C
C *    OF THE PRODUCT OF THE MODIFIED BESSEL FUNCTIONS COMPUTED via THE * C
C *    INTEGRAL REPRESENTATION DEFINED IN M.O.S. (P : 98)               * C
C *                                                                     * C
C *   INTEGN : CONTAINS THE INTEGER VARIABLE NINT WHICH IS USED TO COUNT* C
C *    THE NUMBER OF OVERLAP INTEGRALS COMPUTED IN THE PROCESS          * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE OVELP1(LM123, NORB1,NORB2, NLMA,NLMB, ZETAA,ZETAB,
     $                                           AB,THETAB,PHIAB, XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      INCLUDE "../lib95/CONST.INCL"

      DIMENSION XCORE(N_ORB, N_ORB)
      DIMENSION NLMA(*), NLMB(*), ZETAA(*), ZETAB(*),
     $       YLMAB(0:((LDEV+1)*(LDEV+2))/2), RADAR((L1MAX+1)*(L2MAX+1)),
     $       ARGNT1((L2MAX+1)*(L2MAX+2)), ARGNT2((L2MAX+1)**2), INX(7)

      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:100)
      COMMON/GLE01/RLEG01(L001+L002+L003)
      COMMON/GLE02/RLEG02(L004+L005+L006)
      COMMON/BESS0/BESF01(2*LTOT6)
      COMMON/X4POW/X4NL(10*LTOT6)
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))
      COMMON/INTEGN/NINT

C.... COMPUTATION OF THE REQUIRED SURFACE SPHERICAL HARMONICS

      CALL YLMV0(THETAB, YLMAB)

C.... COLLECTING THE ROOTS OF THE GAUSS-LEGENDRE QUADRATURE FOR COMPUTING
C.... THE RADIAL INTEGRAL OVER [0, AB]

      CALL DUMMY0(AB, LEG00)

C.... COMBINING THE GIVEN SET OF SLATER ORBITALS. WE FIRST BEGIN BY THE
C.... ORBITALS THAT SHOULD BE TRANSLATED IN ORDER TO LOWER THE NUMBER OF
C.... CALLS OF THE ROUTINE GEGEN0

      DO 5 NAT2=1, NORB2
       JST   = 5*(NAT2-1)
       ICOL  = NLMB(JST+2)
       N2    = NLMB(JST+3)
       L2    = NLMB(JST+4)
       M2    = NLMB(JST+5)
       ZETA2 = ZETAB(NAT2)

C.... COMPUTATION OF SOME COEFFICIENTS INVOLVED IN THE GENERALIZED
C.... GEGENBAUER ADDITION THEOREM. THESE FIRST COEFFICIENTS DEPEND ONLY
C.... THE ROOTS OF THE GAUSS-LEGENDRE QUADRATURE COLLECTED IN THE DUMMY0
C.... ROUTINE, THE QUANTUM NUMBERS N2, L2, ZETA2 AND OF COURSE AB.

       LMDL12 = LM123 + L2 + 1
       CALL GEGEN0(LMDL12, N2-L2, ZETA2, AB, RLEG01, RLEG02, RLEAG0,
     $                                             BESF01, X4NL, GEGIN0)

      DO 5 NAT1=1, NORB1
       IST   = 5*(NAT1-1)
       LINE  = NLMA(IST+2)
       N1    = NLMA(IST+3)
       L1    = NLMA(IST+4)
       M1    = NLMA(IST+5)
       ZETA1 = ZETAA(NAT1)

C.... COMPUTATION OF THE COEFFICIENTS INVOLVED IN THE SOLID SPHERICAL
C.... HARMONICS ADDITION AN MULTIPLICATION THEOREM

c       if(m1.ne.m2)then
c       if(iabs(m1-m2).eq.1)goto 100
c       endif
       CALL FUSED(DFAC, L1,M1, L2,M2, AB,YLMAB, ARGNT1, ARGNT2)

C.... COMPUTATION OF THE RADIAL INTEGRAL INVOLVED IN THE OVERLAP INTEGRALS

       CALL ROVELP(LMDL12, N1,L1, N2,L2, ZETA1,ZETA2, AB,
     $                              RLEG01, RLEG02, BESF01, X4NL, RADAR)

C.... THE NORMALIZATION CONSTANTS

       A1 = (2.D0*ZETA1)**N1 * DSQRT((2.D0*ZETA1)/FACT(2*N1))
       A2 = (2.D0*ZETA2)**N2 * DSQRT((2.D0*ZETA2)/FACT(2*N2))

       CSTE = 16.D0 * PI2 * A1 * A2 * DFAC(L2) * (-AB)**L2

C.... INITIALIZATION OF THE INDEXING ARRAY

       CALL RAZ1(INX, 1, 2)

C.... INITIALIZATION OF THE OVERLAP INTEGRALS
c  100 continue

       OVLP = 0.D0

       K1 = 0
       K2 = 0

       DO 10 LP2=0, L2
	INX(1) = INX(1) + INX(2)
       DO 10 MP2=MAX(-LP2, -L2+LP2+M2), MIN(LP2, L2-LP2+M2)
	K1   = K1 + 1
	AUX1 = ARGNT1(K1)
	CALL RAZ1(INX, 2, 2)
       DO 10 LP12=IABS(L1-LP2), L1+LP2, 2
	K2     = K2+1
	INX(2) = INX(2) + 1
	INDIC  = INX(1) + INX(2)
	NYL    = (LP12 * (LP12 + 1))/2 + IABS(MP2-M1)
	MUPW   = (MP2 - M1 - IABS(MP2 - M1))/2
	AUX2   = (-1)**MUPW * ARGNT2(K2) * YLMAB(NYL)
	OVLP   = OVLP + AUX1*AUX2*RADAR(INDIC)
 10    CONTINUE

C.... FOR NON-LINEAR MOLECULES, THE SPHERICAL HARMONICS ARE COMBINED IN
C.... ORDER TO AVOID IMAGINARY VALUES

       ANGLE  = (M2 - M1)*PHIAB
       COSINE = DCOS(ANGLE)
       SINE   = DSIN(ANGLE)
c
       IF(DABS(COSINE) .LT. 1.D-10)THEN
	COSINE = 0.D0
       ENDIF
       IF(DABS(SINE) .LT. 1.D-10)THEN
	SINE = 0.D0
       ENDIF
c
c       cste = cste/(4.0d0)
c       if(m1.ne.0.or.m2.ne.0)cste = cste/2.0d0
c       cste = cste/1.41d0
       OVLPR = COSINE * CSTE*OVLP/DSQRT(AB)
       OVLPI = SINE   * CSTE*OVLP/DSQRT(AB)
c       if(ovlpr.lt.0.0d0)ovlpr= - ovlpr

C.... WHEN THE OVERLAP INTEGRAL IS IMAGINARY, THE IMAGINARY PART CHANGES ITS
C.... SIGN SINCE SIN((M2-M1)PHIAB) = - SIN((M1-M2)PHIAB)

       XCORE(LINE, ICOL) = OVLPR
       XCORE(ICOL, LINE) = OVLPI

c      WRITE(*, 1)LINE, ICOL, N1, L1, M1, ZETA1, N2, L2, M2, ZETA2
c      WRITE(*, 2)OVLPR, OVLPI, AB

       NINT = NINT + 1

 5    CONTINUE

c 1    FORMAT(2Z2, 2X, 3(I3, 2X), D12.6, 3(I3, 2X), D12.6)
c 2    FORMAT(2(D16.10, 2X))
      RETURN
      END
