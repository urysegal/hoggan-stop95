CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 23 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS THE ROOTS OF THE GAUSS-LEGENDRE AND THOSE OF   * C
C * THE GAUSS-LAGUERRE QUADRATURES FOR THE COMPUTATION OF THE TWO-      * C
C * CENTER KINETIC ENERGY INTEGRALS.                                    * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   AB : IS THE DISTANCE SEPARATING THE CENTERS C1 AND C2             * C
C *                                                                     * C
C *   ZETAS : IS THE SUM ZETA1 + ZETA2                                  * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   THE OUTPUT OF THIS ROUTINE IS RETURNED VIA THE GLAE0 COMMON       * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DUMBID(AB, ZETAS)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      EXTERNAL FGLE00

C.... COMPUTATION AND STORAGE OF THE ROOTS

      ISP  = 0
      XXXX = DGLNQ(FGLE00, 0.D0, AB, LEG00)
      ISP  = ISP + LEG00
      XXXX = DGLGQ(FGLE00, AB, ZETAS, LAG00)


      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 23 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS THE ROOTS OF THE GAUSS QUADRATURES INVOLVED    * C
C * IN THE OVLBID PACKAGE.                                              * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE FGLE00(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/GLAE0/RLEAG0(LMHERM), ISP

      DO 5 I=1, N
       RLEAG0(ISP+I) = T(I)
       Y(I)          = 0.D0
 5    CONTINUE

      RETURN
      END



