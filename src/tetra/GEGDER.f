CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 19 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS THE DERIVATIVES OF THE PRODUCT OF THE MODIFIED * C
C * FUNCTIONS INVOLVED IN THE GENERALIZED GEGENBAUER ADDITION THEOREM.  * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   MAXLMD : IS HIGHEST BESSEL FUNCTIONS BEING DERIVED                * C
C *   NORBX  : TOTAL NUMBER OF ATOMIC ORBITALS ATTACHED TO THE Xth ATOM * C
C *   DX     : DISTANCE OF THE Xth ATOM WITH REPECT TO THE ROGIN        * C
C *   NLMX, ZETA : ONE DIMENSIONAL ARRAYS CONTAINING THE PARAMETERS OF  * C
C *    THE ATOMIC ORBITALS, I.E. N, L, M AND ZETA                       * C
C *   ROOTS  : ONE-DIMENSIONAL ARRAY CONTAINING THE ROOTS FOR WHICH THE * C
C *    DERIVATIVES ARE COMPUTED                                         * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   GEGIN0 : ONE DIMENSIONAL ARRAY WHERE THE DERIVATIVES ARE STORED   * C
C *                        GEGIN0(NATX, LMD, R2)                        * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GEGDER(MAXLMD, NORBX,NLMX,ZETA,DX, LEGX,ROOTS, GEGIN0)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER XST, TRANS
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION NLMX(*), ZETA(*), ROOTS(*), GEGIN0(*)


      TRANS = LEGX * (MAXLMD+1)

      DO 5 NATX=1, NORBX
       XST   = 5*(NATX-1)
       NOX   = NLMX(XST+2)
       NX    = NLMX(XST+3)
       LX    = NLMX(XST+4)
       MX    = NLMX(XST+5)
       ZETAX = ZETA(NATX)

       NSTR = (NATX-1)*TRANS
       CALL DERBES(DX, ZETAX, NX-LX, MAXLMD, LEGX, ROOTS, NSTR, GEGIN0)

 5    CONTINUE

      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 19 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE DERIVATIVES OF THE PRODUCT OF THE MODI- * C
C * FIED BESSEL FUNCTIONS OF THE SECOND KIND WITH HALF ODD INTEGER      * C
C * ORDERS. A SIMPLE ALGORITHM BASED ON THE WELL KNOWN RECURRENCE       * C
C * RELATIONS IS USED.                                                  * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   DX : IS THE DISTANCE OF THE Xth ATOM WITH RESPECT TO THE ORIGIN   * C
C *   ZETAX : IS THE ZETA SLATER EXPONENT FOR THE Xth ORBITAL           * C
C *   NLX : IS THE ORDER OF THE DERIVATIVE, I.E. NX-LX                  * C
C *   MAXLMD : IS THE HIGHEST BESSEL FUNCTIONS BEING DERIVED            * C
C *   LEGX : IS THE TOTAL NUMBER OF THE ROOTS (R) WHERE THE DERIVATIVES * C
C *    ARE COMPUTED                                                     * C
C *   ROOTS : A ONE-DIMENSIONAL ARRAY CONTAING THE ROOTS                * C
C *   INDEX : IS A RUNNING INDEX USED TO STORE THE NUMERICAL RESULT IN A* C
C *    ONE-DIMENSIONAL ARRAY                                            * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   GEGIN0 : A ONE-DIMENSIONAL ARRAY WHERE THE DERIVATIVES ARE STORED * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DERBES(DX,ZETAX,NLX, MAXLMD, LEGX, ROOTS, NSTR,GEGIN0)
      INCLUDE "../lib95/SIZE.INCL"
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION ROOTS(*), GEGIN0(*), BESI(0:100), BESK(0:100)

      DO 5 I=1, LEGX
       R = ROOTS(I)

       RM = ZETAX * DMIN1(R, DX)
       RP = ZETAX * DMAX1(R, DX)

       FRAC1 = RM/RP
       FRAC2 = FRAC1 * FRAC1

       CALL BESSI(-NLX, MAXLMD+NLX, RM, BESI)
       CALL BESSK(MAXLMD+NLX, RP, BESK)

       CSTE = (0.5D0*RP/ZETAX)**NLX

       DO 10 LMD=0, MAXLMD
	POWR = 1.D0
	SUMP = BESPRO(RM, NLX, LMD, BESI, BESK, 0, POWR, FRAC2)

	POWR = FRAC1
	SUMN = BESPRO(RM, NLX, LMD, BESI, BESK, 1, POWR, FRAC2)

	INDEX = NSTR + LMD*LEGX + I
	GEGIN0(INDEX) = CSTE * (SUMP - SUMN)
 10    CONTINUE
 5    CONTINUE

      RETURN
      END

