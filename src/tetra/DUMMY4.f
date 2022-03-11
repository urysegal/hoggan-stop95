CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 19 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE ROOTS OF THE GAUSS QUADRATURES USED TO  * C
C * COMPUTE THE DOUBLE INTEGRAL OVER R1 AND R2.                         * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   ZETAMY : SCALING FACTOR FOR THE GAUSS-LAGUERRE QUADRATURE         * C
C *                                                                     * C
C *   RMIN, RMED, RMAX : RADIUS OF THE CONCENTRIC SPHERES               * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   THE ROOTS FOR R2 ARE RETURNED IN THE RLEAG0 ARRAY, WHILE THOSE OF * C
C *    R1 ARE RETURNED IN THR ROOTS ARRAY                               * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DUMMY4(ZETAMY, RMIN, RMAX, AB)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      COMMON/XROOTS/ROOTS(LEA11*LEA12), ISQ
      EXTERNAL FGLE09, FGLA09, FGLE10


C.... COMPUTATION AND STORAGE OF THE ROOTS OF THE OUTER INTEGRAL (OVER R2)

      ISP  = 0
      XXXX = DGLNQ(FGLE09, 0.D0, RMIN, LEG16)
      ISP  = ISP + LEG16
      XXXX = DGLNQ(FGLE09, RMIN, RMAX, LEG17)
      ISP  = ISP + LEG17
      XXXX = DGLGQ(FGLE09, RMAX, ZETAMY, LAG11)


C.... COMPUTATION AND STORAGE OF THE ROOTS OF THE INNERMOST INTEGRAL (R1)

      DO 5 I=1, LEA11
       R2 = RLEAG0(I)

       ISQ  = (I-1)*LEA12
       XXXX = DGLNQ(FGLE10, 0.D0, DMIN1(R2, AB), LEG18)

       ISQ  = ISQ+LEG18
       XXXX = DGLNQ(FGLE10, DMIN1(R2, AB), DMAX1(R2, AB), LEG19)

 5    CONTINUE

      RETURN
      END



