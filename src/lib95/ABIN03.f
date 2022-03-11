      SUBROUTINE ABIN02(ZETA2, A, NL, RLEG01,L1_3, RLEG02,L4_6,
     $                     RLEAG0,LEAS,LEAE, BESF01, X4NL, HEREX,BESS00)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "CONST.INCL"
      DIMENSION RLEG01(*),RLEG02(*),RLEAG0(*),HEREX(*),BESS00(*),
     $          BESF01(*), X4NL(*)


      INDEX = 1

      DO 20 I=1, L1_3
       X  = RLEG01(I)
       X2 = X * X
       X4 = X2*X2
       CALL BA01(1.D0/(X4*X4), INDEX, 2*(I-1), BESF01, BESS00)
 20   CONTINUE

      DO 25 J=1, L4_6
       X  = RLEG02(J)
       X2 = X * X
       X4 = X2*X2
       CALL BA01(X4*X4, INDEX, 2*(L1_3+J-1), BESF01, BESS00)
 25   CONTINUE

      IF(NL .LT. 0)THEN
       RETURN
      ENDIF


      ENTRY HERMIT(ZETA2, A, NL, RLEG01,L1_3, RLEG02,L4_6,
     $                     RLEAG0,LEAS,LEAE, BESF01, X4NL, HEREX,BESS00)


      NDEX1 = (NL-1)*(L1_3+L4_6)
      NDEX2 = NDEX1 + L1_3
      INDEX = 1

      SQRTA = DSQRT(A) / SQRT2

      DO 5 I=LEAS, LEAE
       R         = RLEAG0(I)
       AR        = SQRTA * DSQRT(R)
       Z1        = ZETA2 * AR
       Z2        = 0.5D0*(A-R)/AR
       ARNL      = AR**NL

       DO 10 J=1, L1_3
	X            = RLEG01(J)
	X2           = X * X
	X4           = X2*X2
	ARG1         = Z1*X4
	ARG2         = Z2/X4
	HEREX(INDEX) = ARNL * HN(NL, ARG1) * DEXP(-ARG1*ARG1-ARG2*ARG2)
     $                      * X4NL(NDEX1+J)
	INDEX = INDEX + 1
 10    CONTINUE

       DO 15 K=1, L4_6
	X            = RLEG02(K)
	X2           = X * X
	X4           = X2*X2
	ARG1         = Z1/X4
	ARG2         = Z2*X4
	HEREX(INDEX) = ARNL *
     $        (HN(NL, ARG1) * DEXP(-ARG1*ARG1-ARG2*ARG2))/X4NL(NDEX2+K)
	INDEX = INDEX + 1
 15    CONTINUE
 5    CONTINUE

      RETURN
      END


