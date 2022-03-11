CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE RETURNS THE COEFFICIENTS INVOLVED IN THE GEGENBAUER  *C
C* ADDITION THEOREM FOR THE HYBRID INTEGRALS.                        *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GEGEN2(LMDL12, NL4, ZETA4, AB, RLEG01, RLEG02, RLEAG0,
     $                                             BESF01, X4NL, GEGIN0)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION RLEG01(*), RLEG02(*), RLEAG0(*),
     $                                     BESF01(*), X4NL(*), GEGIN0(*)
      COMMON/BESHER/HEREX(LTOT6*LMHERM), BESS01(LTOT6*(LDEV+1))


C.... COMPUTATION OF THE MODIFIED BESSEL FUNCTION I AND THE HERMITE
C.... POLYNOMIALS FOR EACH ROOT OF THE GAUSS-LEGENDRE QUADRATURE

      IN = 1

      CALL ABIN02(ZETA4, AB, NL4, RLEG01,L1_3, RLEG02,L4_6,
     $                     RLEAG0,IN,LEG03, BESF01, X4NL, HEREX, BESS01)

C.... COMPUTATION OF THE GEGENBAUER'S COEFFICIENTS FOR THE ROOTS
C.... CORRESPONDING TO THE RANGE [0, AB]

      DO 5 LMD=0, LMDL12-1
       IND1 = 1
      DO 5 I=IN, LEG03
       R2   = RLEAG0(I)
       IND2 = LMD+1
       CALL BAD01(IND1, IND2, ZETA4, AB, R2, NL4, LMD, XYINT)
       INDEX         = LMD*LEA02 + I
       GEGIN0(INDEX) = XYINT
 5    CONTINUE

      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE RETURNS THE COEFFICIENTS INVOLVED IN THE GEGENBAUER  *C
C* ADDITION THEOREM FOR THE HYBRID INTEGRALS.                        *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GEGOX (LMDL12, NL4, ZETA3,ZETA4, AB, RLEG01, RLEG02,
     $                                             BESF01, X4NL, GEGIN0)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION RLEG01(*), RLEG02(*), BESF01(*), X4NL(*), GEGIN0(*)

      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      COMMON/BESHER/HEREX(LTOT6*LMHERM), BESS01(LTOT6*(LDEV+1))
      EXTERNAL GLROOT


C.... COLLECTING THE ROOTS OF THE GAUSS-LAGUERRE QUADRATURE

       ISP  = LEA02-LAG02
       XXXX = DGLGQ(GLROOT, AB, ZETA3+ZETA4, LAG02)

C.... COMPUTATION OF THE MODIFIED BESSEL FUNCTION I AND THE HERMITE
C.... POLYNOMIALS FOR EACH ROOT OF THE GAUSS-LEGENDRE QUADRATURE

      IN = LEA02-LAG02+1

      CALL HERMIT(ZETA4, AB, NL4, RLEG01,L1_3, RLEG02,L4_6,
     $                     RLEAG0,IN,LEA02, BESF01, X4NL, HEREX, BESS01)

C.... COMPUTATION OF THE GEGENBAUER'S COEFFICIENTS FOR THE ROOTS
C.... CORRESPONDING TO THE RANGE [0, AB]

      DO 5 LMD=0, LMDL12-1
       IND1 = 1
      DO 5 I=IN, LEA02
       R2   = RLEAG0(I)
       IND2 = LMD+1
       CALL BAD01(IND1, IND2, ZETA4, AB, R2, NL4, LMD, XYINT)
       INDEX         = LMD*LEA02 + I
       GEGIN0(INDEX) = XYINT
 5    CONTINUE

      RETURN
      END
