CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE RETURNS THE ROOTS OF THE GAUSS QUADRATURES USED TO   *C
C* COMPUTES THE DOUBLE INTEGRAL OVER R1 AND R2.                      *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DUMMY3(ZETAMY, RMIN, RMAX)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      COMMON/XROOTS/ROOTS(LEA07*LEA08), ISQ
      EXTERNAL FGLE07, FGLA07, FGLE08, FGLA08


C.... COMPUTATION AND STORAGE OF THE ROOTS OF THE OUTER INTEGRAL (OVER R2)

      ISP  = 0
      XXXX = DGLNQ(FGLE07, 0.D0, RMIN, LEG10)
      ISP  = ISP + LEG10
      XXXX = DGLNQ(FGLE07, RMIN, RMAX, LEG11)
      ISP  = ISP + LEG11
      XXXX = DGLGQ(FGLA07, RMAX, ZETAMY, LAG07)

C.... COMPUTATION AND STORAGE OF THE ROOTS OF THE INNERMOST INTEGRAL (R1)

      DO 5 I=1, LEG10
       R2   = RLEAG0(I)
       ISQ  = (I-1)*LEA08
       XXXX = DGLNQ(FGLE08, 0.D0, R2, LEG12)
       ISQ  = ISQ + LEG12
       XXXX = DGLNQ(FGLE08, R2, RMIN, LEG13)
       ISQ  = ISQ + LEG13
       XXXX = DGLNQ(FGLE08, RMIN, RMAX, LEG14)
 5    CONTINUE

      DO 10 I=1, LEG11
       R2   = RLEAG0(LEG10+I)
       ISQ  = (LEG10+I-1)*LEA08
       XXXX = DGLNQ(FGLE08, 0.D0, RMIN, LEG12)
       ISQ  = ISQ + LEG12
       XXXX = DGLNQ(FGLE08, RMIN, R2, LEG13)
       ISQ  = ISQ + LEG13
       XXXX = DGLNQ(FGLE08, R2, RMAX, LEG14)
 10   CONTINUE

      DO 15 I=1, LAG07
       R2   = RLEAG0(LEG10+LEG11+I)
       ISQ  = (LEG10+LEG11+I-1)*LEA08
       XXXX = DGLNQ(FGLE08, 0.D0, RMIN, LEG12)
       ISQ  = ISQ + LEG12
       XXXX = DGLNQ(FGLE08, RMIN, RMAX, LEG13)
       ISQ  = ISQ + LEG13
       XXXX = DGLNQ(FGLE08, RMAX, R2, LEG14)
 15   CONTINUE

      RETURN
      END
