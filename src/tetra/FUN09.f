CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE RETURNS THE ROOTS OF THE G-LEGENDRE QUADRATURE OVER    *C
C* THE RANGE [0, AB] FOR THE FOUR-CENTER INTEGRALS.                    *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE FGLE09(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/GLAE0/RLEAG0(LMHERM), ISP

      DO 5 I=1, N
       RLEAG0(I+ISP) = T(I)
       Y(I)          = 0.D0
 5    CONTINUE

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE RETURNS THE ROOTS OF THE G-LAGUERRE QUADRATURE OVER    *C
C* THE RANGE [AC, +\INFTY) FOR THE FOUR-CENTER INTEGRALS.        THEY  *C
C* ARE STORED AT THE END OF THE ABOVE ONES, THANKS TO THE INDEX ISP    *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      SUBROUTINE FGLA09(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/GLAE0/RLEAG0(LMHERM), ISP

      DO 5 I=1, N
       RLEAG0(I+ISP) = T(I)
       Y(I)          = 0.D0
 5    CONTINUE

      RETURN
      END





