CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE RADIAL INTEGRAL OCCURING IN THE TWO-    * C
C * CENTER COULOMB INTEGRALS. THESE VALUES ARE STORED IN THE ONE-       * C
C * DIMENSIONAL ARRAY 'RADAR'.                                          * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   INDIC : AN INDEX WHICH IS USED TO DETERMINE THE STARTING STORING  * C
C *    INDEX                                                            * C
C *                                                                     * C
C *   L1,L2,L3,L4 : QUANTUM NUMBERS OF THE SLATER ORBITALS              * C
C *                                                                     * C
C *   AB : DISTANCE SEPARATING THE CENTERS (1) AND (2)                  * C
C *                                                                     * C
C *   ZETA1, ZETA2 : SLATER EXPONENETS                                  * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   RADAR : ONE-DIMENSIONAL ARRAY CONTAINING THE VALUES OF THE RADIAL * C
C *    INTEGRAL INVOLVED IN THE TWO-CENTER COULOMB INTEGRALS            * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RTCCI(INDIC, L1,L2,L3,L4, AB, ZETA3,ZETA4, RADAR)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION RADAR(*)
      COMMON/STOP03/L12, L, L34, LP34, LMD, ISTR, NSTRT, MSTRT
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))
      COMMON/COMV00/VL((NORBM*(NORBM+1))/2*(L1MAX+L2MAX+1)*LEA03)
      EXTERNAL RAD04

C.... PASSAGE TO COMMON

      L12   = L1+L2
      NSTRT = INDIC*LEA03*(L1MAX+L2MAX+1)

      INDEX = 1
      LMIN  = IABS(L3-L4)
      LMAX  = L3 + L4

      DO 5 L=IABS(L1-L2), L1+L2, 2
      DO 5 L34=LMIN, LMAX, 2
       MSTRT = (L34-LMIN)/2*(L1MAX+L2MAX+L3MAX+L4MAX+1)*LEA03
      DO 5 LP34=0, L34
      DO 5 LMD=IABS(L-LP34), L+LP34, 2
       ISTR         = 0
       XINT1        = DGLNQ(RAD04, 0.D0, AB, LEG04)
c       XINT1        = DGLNQ(RAD04, 0.D0, AB, 16)
       ISTR         = ISTR + LEG04
       XINT2        = DGLGQ(RAD04, AB, ZETA3+ZETA4, LAG03)
c       XINT2        = DGLGQ(RAD04, AB, ZETA3+ZETA4, 64)
       RADAR(INDEX) = XINT1 + XINT2
       INDEX        = INDEX + 1
 5    CONTINUE

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C * THIS SUBROUTINE IS CALLED FROM THE DGLNQ AND DGLGQ. IT RETURNS THE  * C
C * VALUE OF THE RADIAL INVOLVED IN THE TWO-CENTER COULOMB INTEGRALS    * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RAD04(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/STOP03/L12, L, L34, LP34, LMD, ISTR, NSTRT, MSTRT
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))
      COMMON/COMV00/VL((NORBM*(NORBM+1))/2*(L1MAX+L2MAX+1)*LEA03)


      INDO = MSTRT + LMD*LEA03 + ISTR

      DO 5 I=1, N
       INDX = NSTRT + (ISTR+I-1)*(L12+1) + L+1
       R2   = T(I)
       Y(I) = 1.D0/DSQRT(R2) * R2**(LP34+2) * VL(INDX) * GEGIN0(INDO+I)
 5    CONTINUE

      RETURN
      END
