      SUBROUTINE YGAUNT(GONE, YIWAN, LMAX2)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER YIWAN
      REAL*4 GONE
      INCLUDE "./lib95/SIZE.INCL"
      DIMENSION GONE((LMDMAX+1)**2 * (LMDMAX+6+1)**2, LMDMAX+3+1)
      DIMENSION ARGNT(100), YIWAN(0:*)

      LMAX2 = (LMDMAX+6+1)*(LMDMAX+6+1)

      DO 5 L1=0, LMDMAX
       I1 = L1 * L1 * LMAX2 + 1

       I2 = -LMAX2
      DO 5 M1=-L1, L1
       I2  = I2 + LMAX2

      DO 5 L2=0, LMDMAX+6
       I3 = L2

       I4 = 0
      DO 5 M2=-L2, L2
       I4 = I4 + 1

       IT = I1 + I2 + I3 + I4

       CALL GAUNT(L2,M2, L1,M1, L3M,L3X,M3, NGAUNT, ARGNT)

       KG = 1
       KX = 1

      DO 5 L3=IABS(L1-L2), L1+L2, 2
       IF(L3 .GE. L3M)THEN
	GONE(IT, KX) = SNGL(ARGNT(KG))
c        write(*,*)"ygaunt",argnt(kg),gone(it,kx)
	KG           = KG + 1
       ELSE
	GONE(IT, KX) = 0.E0
       ENDIF

       KX = KX + 1

 5    CONTINUE

      CPOUT = MCLOCK()

ccc      write(*, *)1.d-2 * (cpout-cpin)

      YIWAN(0) = 1

      DO 10 I=1, 99
       YIWAN(I)  = - YIWAN(I-1)
 10   CONTINUE

      RETURN
      END
