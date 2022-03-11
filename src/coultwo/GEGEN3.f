CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS THE WHOLE COEFFICIENTS OF THE GENERALIZED      * C
C * GEGENBAUER ADDITION THEOREM CORRESPONDING TO THE PRODUCT OF TWO     * C
C * SLATER ORBITALS CENTERED ON THE SAME POINT.                         * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LM123 : IS THE HIGHEST VALUE OF THE QUANTUM NUMBER L              * C
C *                                                                     * C
C *   N3,L3,ZETA3, N4,L4,ZETA4 : ARE THE CONSTANTS CHARACTERIZING THE   * C
C *    SLATER ORBITALS                                                  * C
C *                                                                     * C
C *   AB : IS THE DISTANCE SEPARATING THE CENTERS (1) AND (2)           * C
C *                                                                     * C
C *   RLEG01, RLEG02 : ARE ONE-DIMENSIONAL ARRAYS CONTAINING THE        * C
C *    ROOTS OF THE GAUSS-LEGENDRE QUADRATURES                          * C
C *                                                                     * C
C *   BESF01 : ONE-DIMENSIONAL ARRAY WHERE THE VALUES OF THE BESSEL     * C
C *    FUNCTION NUMERICAL VALUES ARE STORED                             * C
C *                                                                     * C
C *   X4NL : ONE-DIMENSIONAL ARRAY CONTAINING THE FOURTH POWERS OF      * C
C *    A CERTAIN VARIABLE                                               * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   GEGEN0 : ONE-DIMENSIONAL ARRAY WHERE THE VALUES OF THE GENERALIZED* C
C *    GEGENBAUER ADDITION THEOREM ARE STORED                           * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GEGEN3(LM123, N3,L3,ZETA3, N4,L4,ZETA4, AB, RLEG01,
     $                                     RLEG02, BESF01, X4NL, GEGIN0)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION  RLEG01(*), RLEG02(*), BESF01(*), X4NL(*), GEGIN0(*)
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      COMMON/BESHER/HEREX(LTOT6*LMHERM), BESS01(LTOT6*(LDEV+1))
      EXTERNAL GLROOT


      N34    = N3+N4-1
      ZETA34 = ZETA3+ZETA4

C.... COMPUTATION OF THE GAUSS-LAGUERRE ROOTS

      ISP  = LEA03-LAG03
      XXXX = DGLGQ(GLROOT, AB, ZETA34, LAG03)

C.... APPLYING THE MULTIPLICATION THEOREM OF THE SPHERICAL HARMONICS

      LMIN = IABS(L3-L4)
      LMAX = L3 + L4

      DO 5 L34=LMIN, LMAX, 2

       IN = 1
       CALL HERMIT(ZETA34, AB, N34-L34, RLEG01,L1_3, RLEG02,L4_6,
     $                     RLEAG0,IN,LEA03, BESF01, X4NL, HEREX, BESS01)

C.... NSTRT IS THE STARTING INDEX

       NSTRT = (L34-LMIN)/2*(L1MAX+L2MAX+L3MAX+L4MAX+1)*LEA03
       INDEX = 1

      DO 5 LMD=0, LM123+LM123+L3+L4
       IND1 = 1
      DO 5 I=1, LEA03
       R2   = RLEAG0(I)
       IND2 = LMD+1
       CALL BAD01(IND1, IND2, ZETA34, AB, R2, N34-L34, LMD, XYINT)
       GEGIN0(NSTRT+INDEX) = XYINT
       INDEX               = INDEX + 1
 5    CONTINUE

      RETURN
      END
