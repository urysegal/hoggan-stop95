      SUBROUTINE ISBOOL(M, IS, BOOL0, BOOL1)
      IMPLICIT REAL*8 (A-H, O-Z)
      LOGICAL*1 BOOL0, BOOL1

      IF(M.NE.0)THEN
       IS    = ISIGN(1, M)
       BOOL1 = M.GT.0
       BOOL0 = .FALSE.
      ELSE
       IS    = 0
       BOOL1 = .FALSE.
       BOOL0 = .TRUE.
      ENDIF

      RETURN
      END
