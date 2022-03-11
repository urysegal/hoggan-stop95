      FUNCTION HN(N, Z)
      IMPLICIT REAL*8 (A-H, O-Z)

      IF(N .EQ. 0)THEN
       HN = 1.D0
      ELSE
       IF(N .EQ. 1)THEN
	HN = 2.D0 * Z
       ELSE
	IF(N .EQ. 2)THEN
	 Z2 = Z*Z
	 HN = 4.D0*Z2 - 2.D0
	ELSE
	 IF(N .EQ. 3)THEN
	  Z2 = Z*Z
	  HN = (8.D0*Z2 - 12.D0)*Z
	 ELSE
	  IF(N .EQ. 4)THEN
	   Z2 = Z*Z
	   HN = (16.D0*Z2 - 48.D0)*Z2 + 12.D0
	  ELSE
	   IF(N .EQ. 5)THEN
	    Z2 = Z*Z
	    HN = ((32.D0*Z2-160.D0)*Z2+120.D0)*Z
	   ELSE
	    Z2 = Z*Z
	    H4 = (16.D0*Z2 - 48.D0)*Z2 + 12.D0
	    H5 = ((32.D0*Z2-160.D0)*Z2+120.D0)*Z
	    DO 5 I=5, N-1
	     H6 = 2.D0 * Z * H5 - 2.D0 * DBLE(I) * H4
	     H4 = H5
	     H5 = H6
 5          CONTINUE
	    HN = H6
	   ENDIF
	  ENDIF
	 ENDIF
	ENDIF
       ENDIF
      ENDIF
      RETURN
      END
