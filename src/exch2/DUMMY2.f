CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE RETURNS THE ROOTS OF THE GAUSS-LEGENDRE QUADRATURES  *C
C* INVOLVED IN THE COMPUTATION OF THE T.C.C.D.P AND THE TWO-CENTER   *C
C* EXCHANGE INTEGRALS.                                               *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DUMMY2(ZETAM, AB)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      COMMON/XROOTS/ROOTS(LEA05*LEA06), ISQ
      COMMON/R1R2P/R12POW(LEA05*LEA06*(LMDMAX+1)), R2
      COMMON/R1POW/R1PUIS((L1MAX+L2MAX+2)*LEA05*LEA06)
      COMMON/R2POW/R2PUIS((L3MAX+L4MAX+3)*LEA05)
      COMMON/R1SQR/R1SQRT(LEA05*LEA06)
      COMMON/R2SQR/R2SQRT(LEA05)
      EXTERNAL FGLE05, FGLA05, FGLE06

C.... COMPUTATION AND STORAGE OF THE ROOTS OF THE G-LEGENDRE AND G-LAGUERRE OF THE T.C.EX.I

      ISP  = 0
      XXXX = DGLNQ(FGLE05, 0.D0, AB, LEG07)
      ISP  = LEG07
      XXXX = DGLGQ(FGLA05, AB, ZETAM, LAG05)

C.... COMPUTATION AND STORAGE OF THE ROOTS OF THE G-LEGENDRE NEEDED IN THE T.C.C.D.P

      DO 5 I=1, LEA05
       R2   = RLEAG0(I)
       ISQ  = (I-1) * LEA06
       XXXX = DGLNQ(FGLE06, 0.D0, DMIN1(AB, R2), LEG08)
       ISQ  = ISQ + LEG08
       XXXX = DGLNQ(FGLE06, DMIN1(AB, R2), DMAX1(AB, R2), LEG09)
 5    CONTINUE

      DO 10 LMD=1, LMDMAX
       INOW = (LMD-1)*LEA05*LEA06
       INEX = INOW + LEA05*LEA06
      DO 10 J=1, LEA05
      DO 10 K=1, LEG08+LEG09
       I = (J-1)*LEA06 + K
       R12POW(INEX+I) = R12POW(I) * R12POW(INOW+I)
 10   CONTINUE

      DO 15 LMD=1, L1MAX+L2MAX+1
       INOW = (LMD-1)*LEA05*LEA06
       INEX = INOW + LEA05*LEA06
      DO 15 J=1, LEA05
      DO 15 K=1, LEG08+LEG09
       I = (J-1)*LEA06 + K
       R1PUIS(INEX+I) = R1PUIS(I) * R1PUIS(INOW+I)
 15   CONTINUE

      DO 20 J=1, LEA05
      DO 20 K=1, LEG08+LEG09
       I = (J-1)*LEA06 + K
       R1SQRT(I) = 1.D0/DSQRT(ROOTS(I))
 20   CONTINUE

C.... POWERS OF THE ROOTS INVOLVED IN THE OUTER INTEGRAL (I.E. R2)

      DO 25 LMD=1, L3MAX+L2MAX+2
       INOW = (LMD-1)*LEA05
       INEX = INOW+LEA05
      DO 25 I=1, LEA05
       R2PUIS(INEX+I) = R2PUIS(I) * R2PUIS(INOW+I)
 25   CONTINUE


      DO 30 I=1, LEA05
       R2SQRT(I) = 1.D0/DSQRT(RLEAG0(I))
 30   CONTINUE


      RETURN
      END


      SUBROUTINE EXPOS(ZETA3, RLEAG0, LEA05, SEXP1)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION RLEAG0(*), SEXP1(*)


      DO 5 I=1, LEA05
       SEXP1(I) = DEXP(-ZETA3*RLEAG0(I))
 5    CONTINUE


      RETURN
      END
