CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 20 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE TESTS WHETHER OR NOT A GIVEN SERIES SHOULD BE          * C
C * ACCELERATED                                                         * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   NEPSIL : IS THE ORDER OF THE EPSILON ALGORITHM                    * C
C *                                                                     * C
C *   ARVAL : IS AN ARRAY CONTAINING THE TERMS OF THE SERIES            * C
C *                                                                     * C
C *   MAXLMD : IS THE HIGHEST ORDER OF THE TERM INVOLVED IN THE SERIES  * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   RELL : IS THE RELATIVE ACCURACY                                   * C
C *                                                                     * C
C *   EPS2 : IS THE SUM OF THE MAXLMD TERMS                             * C
C *                                                                     * C
C *   EPSIL : IS THE RESULTING VALUE AFTER THE ACCLERATION PROCESS      * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ACCEL(NEPSIL, ARVAL, LIMIT, CONST, VALUE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION ARVAL(0:*), EPS((LDEV+1)*(LDEV+1))

      DATA CUTOFF/1.D-12/

C.... INITIALIZATION OF THE EPSILON-ALGORITHM WORKING MATRIX (PARTIAL SUMS)
C.... WHEN ALL THE TERMS OF THE SERIES ARE NUL, THERE IS NOTHING TO
C.... ACCELERATE

      SUM = 0.0d0
      VALUE = 0.0d0
      SUM    = ARVAL(0)
      EPS(1) = CONST * SUM
      NEPOX  = 1

      DO 15 I=1, LIMIT
       VAROX = DABS(CONST * ARVAL(I))

       IF(VAROX .GE. CUTOFF)THEN
	SUM          = SUM + ARVAL(I)
	EPS(NEPOX+1) = CONST * SUM
	NEPOX        = NEPOX + 1
       ENDIF

 15   CONTINUE

C.... IF THE REAL OR IMAGINARY PART OF THE INTEGRAL IS LOWER THAN A GIVEN
C.... VALUE IT IS SET TO ZERO AND NO ACCELARATION IS DONE

      IF(NEPOX .GE. 3)THEN
       EPS1 = EPS(NEPOX-1)
       EPS2 = EPS(NEPOX)
       RELL = DABS(EPS2-EPS1)/DABS(EPS2)

       IF(RELL .LT. 1.D-7)THEN
	VALUE = EPS2
       ELSE
	VALUE = EPS94(EPS, NEPOX, MIN(NEPSIL, NEPOX-1))
       ENDIF

      ELSE
       VALUE = EPS(NEPOX)

      ENDIF

      RETURN
      END
