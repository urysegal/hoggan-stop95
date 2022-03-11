      SUBROUTINE CARTCO(XYZO, XYZE, XYZV, VMOD, THETV, PHIV)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION XYZO(3), XYZE(3), XYZV(3)
      DATA PIM /1.57079632679489662D0/


      VX   = XYZE(1) - XYZO(1)
      VY   = XYZE(2) - XYZO(2)
      VZ   = XYZE(3) - XYZO(3)

      VMOD = DSQRT(VX*VX + VY*VY + VZ*VZ)

      XYZV(1) = VX
      XYZV(2) = VY
      XYZV(3) = VZ

      THETV = DACOS(VZ/VMOD)

      IF(VX .EQ. 0.D0)THEN
       IF(VY .EQ. 0.D0)THEN
	PHIV = 0.D0
       ELSE
	PHIV = DSIGN(1.D0, VY) * PIM
       ENDIF
      ELSE
       SGNX = DSIGN(1.D0, VX)
       PHIV = DABS(SGNX - 1.D0) * PIM + DATAN(VY/VX)
      ENDIF

      RETURN
      END
