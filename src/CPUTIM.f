      FUNCTION CPUTIM()
      IMPLICIT REAL*8 (A-H, O-Z)

      CPUTIM = 1.D-2 * MCLOCK()

      RETURN
      END
