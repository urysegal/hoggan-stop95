      SUBROUTINE CARTES(RTPA, RTPB, AB, THETAB, PHIAB)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION RTPA(*), RTPB(*)
      INCLUDE "CONST.INCL"
      DATA PIM /1.57079632679489662D0/

      RA = RTPA(1)
      RB = RTPB(1)

      THETA = RADEG * RTPA(2)
      THETB = RADEG * RTPB(2)

      PHIA  = RADEG * RTPA(3)
      PHIB  = RADEG * RTPB(3)

      XA = RA * DSIN(THETA) * DCOS(PHIA)
      YA = RA * DSIN(THETA) * DSIN(PHIA)
      ZA = RA * DCOS(THETA)

      XB = RB * DSIN(THETB) * DCOS(PHIB)
      YB = RB * DSIN(THETB) * DSIN(PHIB)
      ZB = RB * DCOS(THETB)

      XV = XB - XA
      YV = YB - YA
      ZV = ZB - ZA

      AB = DSQRT(XV*XV + YV*YV + ZV*ZV)

      THETAB = DACOS(ZV/AB)

      IF(XV .EQ. 0.D0)THEN
       IF(YV .EQ. 0.D0)THEN
	PHIAB = 0.D0
       ELSE
	PHIAB = DSIGN(1.D0, YV) * PIM
       ENDIF
      ELSE
       SGNX  = DSIGN(1.D0, XV)
       PHIAB = DABS(SGNX - 1.D0) * PIM + DATAN(YV/XV)
      ENDIF

      RETURN
      END

