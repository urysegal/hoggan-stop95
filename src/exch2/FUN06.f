CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE RETURNS THE ROOTS OF THE G-LEGENDRE QUADRATURE OVER    *C
C* THE RANGE [0, R2], INVOLVED IN THE T.C.C.D.P.                       *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE FGLE06(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/XROOTS/ROOTS(LEA05*LEA06), ISP
      COMMON/R1R2P/R12POW(LEA05*LEA06*(LMDMAX+1)), R2
      COMMON/R1POW/R1PUIS((L1MAX+L2MAX+2)*LEA05*LEA06)

      DO 5 I=1, N
       X             = T(I)
       ROOTS(ISP+I)  = X
       R12POW(ISP+I) = X/R2
       R1PUIS(ISP+I) = X
       Y(I)          = 0.D0
 5    CONTINUE

      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE RETURNS THE ROOTS OF THE GAUSS-LEGENDRE QUADRATURE   *C
C* WITHIN THE RANGE [R2, +INFINITY[.                                 *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE FGLA06(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/XROOTS/ROOTS(LEA05*LEA06), ISP
      COMMON/R1R2P/R12POW(LEA05*LEA06*(LMDMAX+1)), R2
      COMMON/R1POW/R1PUIS((L1MAX+L2MAX+2)*LEA05*LEA06)

      IST = LEA05*LEA06

      DO 5 I=1, N
       X                 = T(I)
       ROOTS(ISP+I)      = X
       R12POW(ISP+I)     = 1.D0
       R12POW(IST+ISP+I) = R2/X
       R1PUIS(ISP+I)     = X
       Y(I)         = 0.D0
 5    CONTINUE

      RETURN
      END


