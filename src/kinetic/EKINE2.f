CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS A SET OF TWO-CENTER KINETIC ENERGY INTEGRALS* C
C * OVER A SLATER ORBITAL BASIS SET, USING THE TWO-RANGE ONE-CENTER     * C
C * EXPANSION METHOD.                                                   * C
C *                                                                     * C
C * INPUT :                                                             * C
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
C *   XCORE : IS A TWO-DIMENSIONAL ARRAY CONTAINING THE CORE INTEGRALS  * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE EKINE2(NORB1, NORB2, NLMA,NLMB, ZETAA,ZETAB,
     $                                         AB,THETAB,PHIAB, XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"

      DIMENSION XCORE(N_ORB,N_ORB)
      DIMENSION NLMA(*), NLMB(*), ZETAA(*),ZETAB(*),
     $                                    YLMAB(0:((LDEV+1)*(LDEV+2))/2)
      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:100)
      COMMON/INTEGN/NINT


C.... COMPUTATION OF THE REQUIRED SPHERICAL HARMONICS

      CALL YLMV0(THETAB, YLMAB)

C.... COMBINATION OF THE ORBITALS

      DO 5 NAT1=1, NORB1
       IST   = 5*(NAT1-1)
       LINE  = NLMA(IST+2)
       N1    = NLMA(IST+3)
       L1    = NLMA(IST+4)
       M1    = NLMA(IST+5)
       ZETA1 = ZETAA(NAT1)

      DO 5 NAT2=1, NORB2
       JST   = 5*(NAT2-1)
       ICOL  = NLMB(JST+2)
       N2    = NLMB(JST+3)
       L2    = NLMB(JST+4)
       M2    = NLMB(JST+5)
       ZETA2 = ZETAB(NAT2)

       CALL OVLBID(AB,YLMAB, N2  ,L2,M2,ZETA2, N1,L1,M1,ZETA1, OVL1)
       CALL OVLBID(AB,YLMAB, N2-1,L2,M2,ZETA2, N1,L1,M1,ZETA1, OVL2)
       CALL OVLBID(AB,YLMAB, N2-2,L2,M2,ZETA2, N1,L1,M1,ZETA1, OVL3)

       A1 = (2.D0*ZETA1)**N1 * DSQRT((2.D0*ZETA1)/FACT(2*N1))
       A2 = (2.D0*ZETA2)**N2 * DSQRT((2.D0*ZETA2)/FACT(2*N2))

       CSTE = -0.5D0 * A1 * A2 * DFAC(L1) * (-AB)**L1

C.... GLOBAL AZIMUTHAL ANGLE

       ANGLE  = (M2-M1)*PHIAB
       COSINE = DCOS(ANGLE)
       SINE   = DSIN(ANGLE)

       ZNRJ2C = COSINE * CSTE *
     $  (ZETA2*ZETA2*OVL1 - 2.D0*N2*ZETA2*OVL2 + (N2+L2)*(N2-1-L2)*OVL3)

       ZNRJ2S = SINE * CSTE *
     $  (ZETA2*ZETA2*OVL1 - 2.D0*N2*ZETA2*OVL2 + (N2+L2)*(N2-1-L2)*OVL3)

cccc   WRITE(*, 1)LINE,ICOL, N1, L1, M1, ZETA1, N2, L2, M2, ZETA2
cccc   WRITE(*, 2)ZNRJ2C, ZNRJ2S

       XCORE(LINE, ICOL) = XCORE(LINE, ICOL) + ZNRJ2C
       XCORE(ICOL, LINE) = XCORE(ICOL, LINE) + ZNRJ2S

       NINT = NINT + 1

 5    CONTINUE

c1    FORMAT(2Z2, 2X, 3(I3, 2X), D12.6, 3(I3, 2X), D12.6)
c2    FORMAT(2(D16.10, 2X))

      RETURN
      END
