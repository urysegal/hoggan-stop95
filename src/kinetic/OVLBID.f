CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 23 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS AN UNORMALIZED OVERLAP INTEGRALS FOR A      * C
C * SET OF QUANTUM NUMBERS AND OF COURSE THE GEOMETRY.                  * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   AB, YLMAB : IS THE DISTANCE D(C1, C2) AND THE SPHERICAL HARMONICS * C
C *    DEPENDING ON THE ANGULAR PARAMETERS OF THE VECTOR AB             * C
C *                                                                     * C
C *   N1,L1,M1, N2,L2,M2, ZETA1,ZETA2 : CONSTANTS CHARACTERIZING THE    * C
C *    ORBITALS                                                         * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   OVL : IS THE UNNORMALIZED OVERLAP INTEGRAL VALUE                  * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE OVLBID(AB,YLMAB, N1,L1,M1, ZETA1, N2,L2,M2, ZETA2, OVL)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      INCLUDE "../lib95/CONST.INCL"

      DIMENSION YLMAB(0:*), ARGNT(100)

      COMMON/STPBID/ZETA1P, N1P, LP2, LP12, ISTR
      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:100)
      COMMON/GLE01/RLEG01(L001+L002+L003)
      COMMON/GLE02/RLEG02(L004+L005+L006)
      COMMON/BESS0/BESF01(2*LTOT6)
      COMMON/X4POW/X4NL(10*LTOT6)
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      COMMON/BESHER/HEREX(LTOT6*LMHERM), BESS01(LTOT6*(LDEV+1))
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))

      EXTERNAL RADBID

C.... PASSAGE TO COMMON

      ZETA1P = ZETA1
      N1P    = N1

      CALL DUMBID(AB, ZETA1+ZETA2)

C.... COMPUTATION OF THE MODIFIED BESSEL FUNCTIONS WITH THE HERMITE
C.... POLYNOMIALS

      IN = 1

      CALL ABIN02(ZETA2, AB, N2-L2, RLEG01,L1_3, RLEG02,L4_6,
     $                     RLEAG0,IN,LEA00, BESF01, X4NL, HEREX, BESS01)

C.... COMPUTATION OF THE GEGENBAUER'S COEFFICIENTS

      INDEX = 1

      DO 5 LMD=0, L1+L2
       IND1 = 1
      DO 5 I=1, LEA00
       R    = RLEAG0(I)
       IND2 = LMD+1
       CALL BAD01(IND1, IND2, ZETA2, AB, R, N2-L2, LMD, XYINT)
       GEGIN0(INDEX) = XYINT
       INDEX         = INDEX + 1
 5    CONTINUE

C***********************************************************************

      OVL  = 0.D0
      ABLP = -AB
      DO 10 LP2=0, L2
       ABLP = -ABLP/AB
      DO 10 MP2=MAX(-LP2, -L2+LP2+M2), MIN(LP2, L2-LP2+M2)
       CALL GAUNT(LP2,MP2, L2,M2, LMIN,LMAX,MP22, NGAUNT,ARGNT)
       NYL  = ((L2-LP2)*(L2-LP2+1))/2 + IABS(M2-MP2)
       MUPW = (M2-MP2 - IABS(M2-MP2))/2
       AUX1 = (-1)**MUPW *
     $        ABLP*YLMAB(NYL)*ARGNT(1)/(DFAC(LP2)*DFAC(L2-LP2))
       CALL GAUNT(L1,M1, LP2,MP2, LMIN1,LMAX1,M12, NGAUNT, ARGNT)
       KGAUNT = 0

      DO 10 LP12=LMIN1, LMAX1, 2
       KGAUNT = KGAUNT + 1
       NYL    = (LP12*(LP12+1))/2 + IABS(MP2-M1)
       MUPW   = (MP2-M1-IABS(MP2-M1))/2
       AUX2   = (-1)**MUPW * ARGNT(KGAUNT) * YLMAB(NYL)

       ISTR   = 0
       XINT1  = DGLNQ(RADBID, 0.D0, AB, LEG00)
       ISTR   = ISTR + LEG00
       XINT2  = DGLGQ(RADBID, AB, ZETA1+ZETA2, LAG00)

       OVL    = OVL + AUX1*AUX2*(XINT1+XINT2)

 10   CONTINUE

      OVL = 16.D0*PI2 * OVL / DSQRT(AB)

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 23 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE RADIAL INTEGRALS INVOLVED IN THE ABOVE  * C
C * OVERLAP INTEGRALS.                                                  * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RADBID(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/STPBID/ZETA1, N1, LP2, LP12, ISTR
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))

      INDO = LP12*LEA00 + ISTR

      DO 5 I=1, N
       R    = T(I)
       Y(I) = 1.D0/DSQRT(R) * R**(N1+LP2+1) *
     $                                   DEXP(-ZETA1*R) * GEGIN0(INDO+I)
 5    CONTINUE

      RETURN
      END
