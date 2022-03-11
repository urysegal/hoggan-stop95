CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            JAN 06 1995               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE WRITES THE INTEGRALS IN A FILE NAMED "INT.XXX"         * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   NTORB : TOTAL NUMBER OF ATOMIC ORBITALS USED TO DESCRIBE THE      * C
C *    MOLECULAR SYSTEM                                                 * C
C *                                                                     * C
C *   NSUM : (NTOB * (NTORB + 1))/2                                     * C
C *                                                                     * C
C *   FATIMA : IS A ONE-DIMENSIONAL ARRAY CONTAINING THE REAL TWO-      * C
C *    ELECTRON INTEGRALS                                               * C
C *                                                                     * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE STPOUT(NTORB, NSUM, FATIMA)
      IMPLICIT REAL*8 (A-H, O-Z)

      DIMENSION FATIMA(*)

      END
