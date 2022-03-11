      FUNCTION EPS94(EPS, NMAX, NINF)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION EPS(*)


      IF(NMAX .GE. 3)THEN

       DO 5 I=1, NMAX-2
	EPS(NMAX+I) = (EPS(I)*EPS(I+2) - EPS(I+1)*EPS(I+1))/
     $                               (EPS(I+2) - 2.D0*EPS(I+1) + EPS(I))
 5     CONTINUE
       I5 = NMAX + I - 1
      ENDIF

      DO 10 K=2, MIN((NMAX/2), (NINF/2))
       K1 = (K-2) * (NMAX-(K-2)+1)
       K2 = (K-1) * (NMAX-(K-1)+1)
       K3 =    K  * (NMAX-K+1)

      DO 15 M=1, NMAX-2*K
       I1 = K1 + (M+2)
       I2 = K2 + M
       I3 = K2 + (M+1)
       I4 = K2 + (M+2)
       I5 = K3 + M

       V1 = EPS(I1) - EPS(I3)
       V2 = EPS(I4) - EPS(I3)
       V3 = EPS(I2) - EPS(I3)

       DNUM = V1 * (V2 + V3) - V2 * V3
       DNOM = V1 * V2 * V3


       IF(DNOM .EQ. 0.D0)THEN
	EPS(I5) = EPS(I3)
       ELSE
	EPS(I5) = DNOM/DNUM + EPS(I3)
       ENDIF

 15   CONTINUE

 10   CONTINUE

      EPS94 = EPS(I5)

      RETURN
      END
