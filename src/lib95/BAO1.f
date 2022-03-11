      SUBROUTINE BA01(Z, INDEX, ISTRT, BESF01, BESS00)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "SIZE.INCL"
      INCLUDE "CONST.INCL"
      DIMENSION BESS00(*), BESF01(*)


      ZINV = 1.D0 / Z

      IF(Z .LE. 10.D0)THEN

       IND01           = INDEX + LDEV

       B20             = BESF01(ISTRT+1)
       B19             = BESF01(ISTRT+2)
       BESS00(IND01)   = B20
       BESS00(IND01-1) = B19

       DO 5 I=1, LDEV-1
	II  = LDEV - I
	B18 = DBLE(2*II + 1) * ZINV * B19 + B20
	B20 = B19
	B19 = B18
	BESS00(IND01-I-1) = B18
 5     CONTINUE
      ELSE

       B0              = BESF01(ISTRT+1)
       BESS00(INDEX)   = B0
       B1              = BESF01(ISTRT+2)
       BESS00(INDEX+1) = B1
       DO 10 I=1, LDEV-1
	B2 = B0 - DBLE(2*I+1) * ZINV * B1
	B0 = B1
	B1 = B2
	BESS00(INDEX+I+1) = B2
 10    CONTINUE
      ENDIF

      INDEX = INDEX + LDEV + 1

      RETURN
      END



      FUNCTION BA02(N, ZN, Z)
      IMPLICIT REAL*8 (A-H, O-Z)

      MMAX = IDINT(Z) + 5

      Z2 = 0.25D0 * Z*Z
      SUM = ZN

      DO 5 M=0, MMAX
       ZN = Z2/(DBLE(M+1)*(N+M+1.5D0)) * ZN
       SUM = SUM + ZN
 5    CONTINUE

      BA02 = SUM

      RETURN
      END




      SUBROUTINE BA03(N, Z, PGAM19, PGAM20)
      IMPLICIT REAL*8 (A-H, O-Z)
      DATA PI, SPI/3.14159265358979324D0, 1.77245385090551603D0/

      PGAM20 = DEXP(-Z) * DSQRT(2.D0*Z)/SPI

      DO 5 I=0, N-1
       PGAM19 = PGAM20
       PGAM20 = (0.5D0*Z) / (DBLE(I) + 1.5D0) * PGAM20
 5    CONTINUE

      RETURN
      END
