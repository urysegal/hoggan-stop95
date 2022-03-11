C.... A INSERER DANS LA BOUCLE D'APRES, DANS TETRAC.f
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 19 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS THE GAUSS-LAGUERRE ROOTS FO THE INNERMOST      * C
C * INTEGRALS.                                                          * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DUMMY5(ZETA12, AB, RLEAG0)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION RLEAG0(*)
      COMMON/XROOTS/ROOTS(LEA11*LEA12), ISQ
      EXTERNAL FGLA10


      DO 5 I=1, LEA11
       R2   = RLEAG0(I)
       ISQ  = I*LEA12 - LAG12
       XXXX = DGLGQ(FGLA10, DMAX1(R2, AB), ZETA12, LAG12)
 5    CONTINUE

      RETURN
      END
