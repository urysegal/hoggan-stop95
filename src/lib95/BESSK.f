
C VERSION ASCENDANTE
C Dand cette version, j'ai simplifie l'ancien code pour n'utiliser que
C la recurrence ascendante.
C Peut etre y faudra t il ajouter un check pour underflow, mais
C sur beaucoup de machine, cela se passe bien (mis a zero)
C
       SUBROUTINE BESSK(INDEX, Z, BESS00)
       IMPLICIT REAL*8 (A-H, O-Z)
       INCLUDE "SIZE.INCL"
       INCLUDE "CONST.INCL"
       DIMENSION BESS00(0:*)

       SQ_ROOT_PI = 0.77245385D+00

       BESS00(0) = SQ_ROOT_PI/DSQRT(2.D0 * Z) * DEXP(-Z)
       BESS00(1) = BESS00(0) * (1.D0 + 1.D0/Z)

       DO 5 I=1, INDEX-1
        BESS00(I+1) = BESS00(I-1) + DBLE(2*I + 1)/Z * BESS00(I)
  5    CONTINUE

       RETURN
       END

