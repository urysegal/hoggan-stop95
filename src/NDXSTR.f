      SUBROUTINE NDXSTR(NORBAT, IC, NLMSTR, ZETSTR)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION NORBAT(*)
      INTEGER ZETSTR

      NORB = 0

      DO 5 I=1, IC-1
       NORB = NORB + NORBAT(I)
 5    CONTINUE

      NLMSTR = 5*NORB + 1
      ZETSTR = NORB + 1

      RETURN
      END

