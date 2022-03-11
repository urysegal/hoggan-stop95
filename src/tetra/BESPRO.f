CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 19 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS FUNCTION RETURNS THE PRODUCT OF THE MODFIED BESSEL FUNCTIONS   * C
C * I AND K OCCURING IN THE EXPRESSION OF THE DERIVATIVES OF THE SUCH   * C
C * A PRODUCT.                                                          * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION BESPRO(RMIN, NLX, LMD, BESI, BESK, NINF, POWR, FRAC2)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER P, Q
      DIMENSION BESI(0:*), BESK(0:*)
      COMMON/CNP0/CNP(0:30, 0:30)

      IF(RMIN .NE. 0.D0)THEN

       BESPRO  = 0.D0
       POWFR   = POWR

       DO 5 M=NINF, NLX, 2

	CSTM  = CNP(NLX, M) * POWFR
	POWFR = POWFR * FRAC2

	SUMI = 0.D0
	DO 10 P=0, M
	 SUMI = SUMI + CNP(M, P) * BESI(NLX+LMD+M-2*P)
 10     CONTINUE

	SUMK = 0.D0
	DO 15 Q=0, NLX-M
	 INDEX = IABS(2*(LMD+NLX-M-2*Q) + 1)/2
	 SUMK  = SUMK + CNP(NLX-M, Q) * BESK(INDEX)
 15    CONTINUE

       BESPRO = BESPRO + CSTM * SUMI * SUMK

 5     CONTINUE
      ELSE
       BESPRO = 0.D0
      ENDIF


      RETURN
      END




