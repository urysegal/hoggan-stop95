
      SUBROUTINE BESSI(MINN, MAXN, Z, T)
      IMPLICIT REAL*8 (A-H, O-Z)
      DATA SPI, SQR2/1.77245385090551603D0, 1.41421356237309505D0/

      DIMENSION T(0:*)

      DZ1 = DEXP(-Z)

      IF(Z .GE. 30.D0)THEN
       SQRZ = DSQRT(Z)
       DZ2  = DZ1 * DZ1
       T(0) = (1.D0 - DZ2)/(SPI * SQRZ * SQR2)
       T(1) = ((1.D0 - 1.D0/Z) + DZ2 * (1.D0 + 1.D0/Z))/(SQR2*SQRZ*SPI)

       DO 25 I=1, MAXN-1
	T(I+1) = T(I-1) - DBLE(2*I+1)/Z * T(I)
 25    CONTINUE

      ELSE
       INX      = IABS(MINN) + MAXN

       IF(Z .NE. 0.D0)THEN
	T(INX)   = DZ1 * BISER(MAXN, Z)
	T(INX-1) = DZ1 * BISER(MAXN-1, Z)

	DO 5 I=MAXN-1, MINN+1, -1
	 IN      = IABS(MINN) + I
	 T(IN-1) = T(IN+1) + DBLE(2*I+1)/Z*T(IN)
 5      CONTINUE
       ELSE
	DO 10 I=MINN, MAXN
	 T(INX + I) = 0.D0
 10     CONTINUE
       ENDIF
      ENDIF

      RETURN
      END




      FUNCTION BISER(N, Z)
      IMPLICIT REAL*8 (A-H, O-Z)
      LOGICAL*1 BOOL
      DATA SPI/1.77245385090551603D0/

      VL = DSQRT(2.D0*Z)/SPI

      DO 5 LMD=0, N-1
       VL = (0.5D0*Z)/(LMD+1.5D0) * VL
 5    CONTINUE

      UM  = VL
      SUM = UM

      DO  10 M=0, 40
       SOX = SUM
       UM  = (0.25D0*Z*Z)/((M+1)*(N+M+1.5D0)) * UM
       SUM = SUM + UM
       REL = DABS(SUM - SOX)/DABS(SUM)
       IF(REL .LE. 1.D-10)THEN
	GOTO 15
       ENDIF
 10   CONTINUE

 15   BISER = SUM

      RETURN
      END
