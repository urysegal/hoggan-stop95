CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE RETURNS THE COEFFICIENTS OF THE GEGENBAUER ADDITION  *C
C* THEOREM INVOLVED IN THE THE TRANSLATION OF THE SLATER ORBITAL THE *C
C* LABEL OF WHICH IS (2).                                            *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GEGEN6(LMDL12, NL2, ZETA2, AB, RLEG01, RLEG02, RLEAG0,
     $                                      BESF01, ROOTS, X4NL, GEGIN0)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION RLEG01(*), RLEG02(*), RLEAG0(*), ROOTS(*),
     $                                     BESF01(*), X4NL(*), GEGIN0(*)
      COMMON/BESHER/HEREX(LTOT6*LMHERM), BESS01(LTOT6*(LDEV+1))

C.... COMPUTATION OF THE MODIFIED BESSEL FUNCTION I, ONCE AND FOR ALL

      NONE = -1
      INP  = 1

      CALL ABIN02(ZETA2, AB, NONE, RLEG01,L1_3, RLEG02,L4_6,
     $                    RLEAG0,INP,LEA07, BESF01, X4NL, HEREX, BESS01)

C.... COMPUTATION OF THE HERMIT POLYNOMIALS INVOLVED IN THE GEGENBAUER A. THEO.
C.... THE FIRST VALUES DEPENDS ONLY ON THE 2ND ORBITAL PARAMETERS AND OF
C.... COURSE ON THE ROOTS OF THE GAUSS-LEGENDRE QUADRATURES (LEG10+LEG11 ROOTS)

      IN = 1

      DO 5 I=1, LEA07
       ISIT = (I-1)*LEA08
       CALL HERMIT(ZETA2, AB, NL2, RLEG01,L1_3, RLEG02,L4_6,
     $        ROOTS(ISIT+1),IN,LEA08-LAG08, BESF01, X4NL, HEREX, BESS01)
      DO 5 LMD=0, LMDL12-1
       IND1 = 1
      DO 5 J=1, LEA08-LAG08
       R1   = ROOTS(ISIT+J)
       IND2 = LMD+1
       CALL BAD01(IND1, IND2, ZETA2, AB, R1, NL2, LMD, XYINT)
       INDEX         = (I-1)*LMDL12*LEA08 + LMD*LEA08 + J
       GEGIN0(INDEX) = XYINT
 5    CONTINUE


      RETURN
      END