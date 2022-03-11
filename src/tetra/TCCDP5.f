CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 19 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE COMPUTES THE RADIAL PART OF THE TWO-CENTER DISTRIBU * C
C * TION POTENTIAL.                                                     * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LM123 : MAXIMUM OF ALL L'S                                        * C
C *   N1,L1, ZETA1, N2,L2, ZETA2 : CHARACTERESTICS OF THE TWO STO'S FOR * C
C *    WHICH THE T.C.C.D.P. IS COMPUTED                                 * C
C *   AB : DISTANCE BETWEEN THE CENTERS DEFINING THE T.C.C.D.P.         * C
C *   RMIN, RMAX : MIN(AC, AD) AND MAX(AC, AD)                          * C
C *   RLEAG0 : ONE-DIMENSIONAL ARRAY CONTAINING THE ROOTS OF THE OUTER  * C
C *    RADIAL INTEGRAL                                                  * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   VL : ONE-DIMENSIONAL ARRAY CONTAINING THE RADIAL PART OF THE      * C
C *    T.C.C.D.P.                                                       * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TCCDP5(KINDEX, LM123, N1,L1, N2,L2, ZETA1,ZETA2,
     $                                         AB,RMIN,RMAX, RLEAG0, VL)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION VL(*), RLEAG0(*), BESI(0:100), BESK(0:100)
      COMMON/XROOTS/ROOTS(LEA11*LEA12), ISP
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))
      COMMON/EXCH04/ZETA1P, R2, N1P, LP2, L, LMD2, LMDL12, IO, ISTR
      EXTERNAL FGLA10

C.... PASSAGE TO COMMON

      ZETA1P = ZETA1
      N1P    = N1
      LMDL12 = LMDMAX+LM123+L2+1
      NL2    = N2-L2

C.... COLLECTING THE ROOTS OF GAUSS-LAGUERRE QUADRATURE FOR INTEGRATING OVER R1

      DO 5 I=1, LEA11
       R2   = RLEAG0(I)
       ISP  = I*LEA12 - LAG12
       XXXX = DGLGQ(FGLA10, DMAX1(R2, AB), ZETA1+ZETA2, LAG12)
 5    CONTINUE

C.... COMPUTATION OF THE GEGENBAUER COEFFICIENTS FOR THE ABOVE ROOTS

      IN = LEA12-LAG12+1

      DO 10 I=1, LEA11
       ISIT = (I-1)*LEA12

      DO 10 J=IN, LEA12
       R1 = ROOTS(ISIT+J)

       RM = ZETA2*DMIN1(R1, AB)
       RP = ZETA2*DMAX1(R1, AB)

       FRAC1 = RM/RP
       FRAC2 = FRAC1*FRAC1

       CALL BESSI(-NL2, LMDL12+NL2, RM, BESI)
       CALL BESSK(LMDL12+NL2, RP, BESK)

       CSTE = (0.5D0*RP/ZETA2)**NL2

      DO 10 LMD=0, LMDL12-1

       POWR = 1.D0
       SUMP = BESPRO(RM, NL2, LMD, BESI, BESK, 0, POWR, FRAC2)

       POWR = FRAC1
       SUMN = BESPRO(RM, NL2, LMD, BESI, BESK, 1, POWR, FRAC2)

       INDEX = (I-1)*LEA12*LMDL12 + LMD*LEA12 + J
       GEGIN0(INDEX) = CSTE * (SUMP - SUMN)

 10   CONTINUE


C.... COMPUTATION OF THE RADIAL PART OF THE T.C.C.D.P

      DO 15 L=0, LMDMAX
      DO 15 LP2=0, L2
      DO 15 LP12=IABS(L1-LP2), L1+LP2, 2
      DO 15 LMD2=IABS(L-LP12), L+LP12, 2

       CALL OUADHIAS(KINDEX, ZETA1, ZETA2, AB, RLEAG0, VL)

 15   CONTINUE

      RETURN
      END
