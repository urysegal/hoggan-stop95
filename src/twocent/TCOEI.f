CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE COMPUTES A STREAM OF OVERLAP INTEGRALS, GIVEN A SET OF * C
C * SLATER ORBITALS                                                     * C
C *                        < XA | 1/RA | XB >                           * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LM123 : IS THE HIGHEST VALUE OF THE QUANTUM NUMBER L              * C
C *                                                                     * C
C *   NUMAT : IS AN ARRAY CONTAINING THE ATOMIC NUMBERS OF THE ATOMS    * C
C *    INVOLVED IN THE MOLECULE                                         * C
C *                                                                     * C
C *   NORB1 AND NORB2 : TOTAL NUMBER OF STOs ATTACHED TO THE CENTERS (1)* C
C *    AND (2)                                                          * C
C *                                                                     * C
C *   NLMA, NLMB, ZETAA AND ZETAB : ARRAYS CONTAINING THE DESCRIPTION   * C
C *    OF THE ORBITALS, i.e. QUANTUM NUMBERS AND SLATER EXPONENTS       * C
C *                                                                     * C
C *   AB, THETAB, PHIAB : SPHERICAL COORDINATES OF THE VECTOR SEPARATING* C
C *    THE CENTERS                                                      * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   XCORE : IS A TWO-DIMENSIONAL ARRAY CONTAINING THE IMAGINARY CORE  * C
C *    INTEGRALS                                                        * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TCOEI(LM123, NUMAT,NORB1,NORB2,NLMA,NLMB,ZETAA,ZETAB,
     $                                           AB,THETAB,PHIAB, XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      INCLUDE "../lib95/CONST.INCL"

      DIMENSION XCORE(N_ORB,N_ORB)
      DIMENSION NLMA(*), NLMB(*), ZETAA(*), ZETAB(*),
     $       YLMAB(0:((LDEV+1)*(LDEV+2))/2), RADAR((L1MAX+1)*(L2MAX+1)),
     $       ARGNT1((L2MAX+1)*(L2MAX+1)),  ARGNT2((L2MAX+1)**3), INX(7)

      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:100)
      COMMON/GLE01/RLEG01(L001+L002+L003)
      COMMON/GLE02/RLEG02(L004+L005+L006)
      COMMON/BESS0/BESF01(2*LTOT6)
      COMMON/X4POW/X4NL(10*LTOT6)
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))
      COMMON/INTEGN/NINT


C.... COMPUTATION OF THE REQUIRED SPHERICAL HARMONICS

      CALL YLMV0(THETAB, YLMAB)

C.... COLLECTING THE ROOTS OF THE GAUSS QUADRATURES FOR COMPUTING THE RADIAL INTEGRAL

      CALL DUMMY0(AB, LEG15)

C.... COMBINING THE GIVEN SET OF SLATER ORBITALS

      DO 5 NAT2=1, NORB2
       JST   = 5*(NAT2-1)
       ICOL  = NLMB(JST+2)
       N2    = NLMB(JST+3)
       L2    = NLMB(JST+4)
       M2    = NLMB(JST+5)
       ZETA2 = ZETAB(NAT2)

C.... COMPUTATION OF SOME COEFFICIENTS INVOLVED IN THE GENERALIZED
C.... GEGENBAUER ADDITION THEOREM. THESE FIRST COEFFICIENTS DEPEND ONLY ON
C.... THE ROOTS OF THE GAUSS-LEGENDRE QUADRATURES COLLECTED IN THE DUMMY0
C.... ROUTINE + THE QUANTUM NUMBERS N2, L2 + ZETA2 AND OF COURSE AB.

       LMDL12 = LM123+L2+1
       CALL GEGEN8(LMDL12, N2-L2, ZETA2, AB, RLEG01, RLEG02, RLEAG0,
     $                                             BESF01, X4NL, GEGIN0)

      DO 5 NAT1=1, NORB1
       IST   = 5*(NAT1-1)
       LINE  = NLMA(IST+2)
       N1    = NLMA(IST+3)
       L1    = NLMA(IST+4)
       M1    = NLMA(IST+5)
       ZETA1 = ZETAA(NAT1)

C.... COMPUTATION OF THE COEFFICIENTS INVOLVED IN TH SOLID SPHERICAL
C.... HARMONICS ADDITION AND MULTIPLICATION THEOREMS

       CALL FUSED(DFAC, L1,M1, L2,M2, AB,YLMAB, ARGNT1, ARGNT2)

C.... COMPUTATION OF THE RADIAL INTEGRALS INVOLVED IN THE OVERLAP INTEGRALS

       CALL RTCOEI(LMDL12, N1-1,L1, N2,L2, ZETA1,ZETA2, AB,
     $                              RLEG01, RLEG02, BESF01, X4NL, RADAR)

C.... THE NORMALIZATION CONSTANTS

       A1 = (2.D0*ZETA1)**N1 * DSQRT((2.D0*ZETA1)/FACT(2*N1))
       A2 = (2.D0*ZETA2)**N2 * DSQRT((2.D0*ZETA2)/FACT(2*N2))

       CSTE = 16.D0 * PI2 * A1 * A2 * DFAC(L2) * (-AB)**L2

C.... INITIALIZATION OF THE INDEXING ARRAY

      CALL RAZ1(INX, 1, 2)

C.... INITIALIZATION OF THE OVERLAP INTEGRALS

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
       K2     = K2 + 1
       INX(2) = INX(2) + 1
       INDIC  = INX(1) + INX(2)
       NYL    = (LP12*(LP12+1))/2 + IABS(MP2-M1)
       MUPW   = (MP2-M1 - IABS(MP2-M1))/2
       AUX2   = (-1)**MUPW * ARGNT2(K2) * YLMAB(NYL)

       OVLP   = OVLP + AUX1*AUX2*RADAR(INDIC)

 10   CONTINUE

C.... FOR NON-LINEAR MOLECULES, THE SPHERICAL HARMONICS MUST BE COMBINED IN
C.... SUCH A WAY AS TO AVOID IMAGINARY VALUES

       ANGLE  = (M2 - M1)*PHIAB
       COSINE = DCOS(ANGLE)
       SINE   = DSIN(ANGLE)

       OVLPR  = COSINE * CSTE*OVLP/DSQRT(AB)
       OVLPI  = SINE   * CSTE*OVLP/DSQRT(AB)

C.... FOR LINE > ICOL THE TERM <ICOL | LINE> IS DEDUCED FROM <LINE | ICOL>
C.... BY CONJUGATING THE IMAGINARY PART -(-NUMAT*OVLPI)

       IF(LINE .LT. ICOL)THEN
	XCORE(LINE, ICOL) = XCORE(LINE, ICOL) - NUMAT*OVLPR
	XCORE(ICOL, LINE) = XCORE(ICOL, LINE) - NUMAT*OVLPI
       ELSE
	XCORE(ICOL, LINE) = XCORE(ICOL, LINE) - NUMAT*OVLPR
	XCORE(LINE, ICOL) = XCORE(LINE, ICOL) - (-NUMAT*OVLPI)
       ENDIF

ccc    WRITE(*, 1)LINE,ICOL, NUMAT*OVLPR, NUMAT*OVLPI

       NINT = NINT + 1

 5    CONTINUE

c1    FORMAT(2Z2, 1X, 2(D16.10, 2X))

      RETURN
      END
