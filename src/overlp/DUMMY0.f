CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS THE ROOTS OF THE GAUSS-LEGENDRE  QUADRATURE    * C
C * USED FOR THE OVERLAP INTEGRALS.                                     * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   AB : THE UPPER LIMIT OF THE RANGE WITHIN WHICH THE ROOTS ARE      * C
C *    COMPUTED, i.e. [0, AB]                                           * C
C *                                                                     * C
C *   LEGP : IS THE ORDER OF THE OF THE G-LEGENDRE QUADRATURE           * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   RLEAG0 : IS AN ARRAY CONTAINING THE ROOTS                         * C
C *                                                                     * C
C *   ISP : IS THE RANK FROM WHICH THE ROOTS ARE STORED                 * C
C *                                                                     * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DUMMY0(AB, LEGP)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      EXTERNAL GLROOT

C.... COMPUTATION AND STORAGE OF THE ROOTS

      ISP  = 0
      XXXX = DGLNQ(GLROOT, 0.D0, AB, LEGP)

      RETURN
      END
