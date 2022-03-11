CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE COMPUTES THE OUTER RADIAL INTEGRAL (OVER R2) NEEDED  *C
C* IN THE THREE-CENTER EXCHANGE INTEGRALS.                           *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RADP1(L1,L2, N3,L3, N4,L4, ZETA3,ZETA4, ZETAMY,
     $         AC,RMIN,RMAX, RLEG01,RLEG02, BESF01, X4NL, RLEAG0, RADAR)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION RLEG01(*), RLEG02(*), BESF01(*), X4NL(*), RLEAG0(*),
     $          RADAR(*)
      COMMON/BESHER/HEREX(LTOT6*LMHERM), BESS01(LTOT6*(LDEV+1))
      COMMON/COMV00/VL(LEA07*(LMDMAX+1)*(L1MAX+1)**3)
      COMMON/COMG01/GEGIN01(LMHERM*(LDEV+1))
      COMMON/EXCH031/ZETA3P, N3P, LP4, INDIC, LMD2, NSTR
      EXTERNAL RAD08


C.... PASSAGE TO COMMON

      ZETA3P = ZETA3
      N3P    = N3

C.... COMPUTATION OF THE GEGENBAUER COEFFICIENTS

      IN = 1
      CALL ABIN02(ZETA4, AC, N4-L4, RLEG01,L1_3, RLEG02,L4_6,
     $                     RLEAG0,IN,LEA07, BESF01, X4NL, HEREX, BESS01)

      INDEX = 1

      DO 5 LMD=0, LMDMAX+L3+L4
       IND1 = 1
      DO 5 I=1, LEA07
       R2   = RLEAG0(I)
       IND2 = LMD+1
       CALL BAD01(IND1, IND2, ZETA4, AC, R2, N4-L4, LMD, XYINT)
       GEGIN01(INDEX) = XYINT
       INDEX = INDEX + 1
 5    CONTINUE

C.... COMPUTATION OF THE RADIAL PART OF THE THREE-CENTER EXCHANGE INTEGRALS

      INDEX = 1
      INDIC = 0

      DO 10 L=0, LMDMAX
      DO 10 LP2=0, L2
      DO 10 LP12=IABS(L1-LP2), L1+LP2, 2
      DO 10 LMD1=IABS(L-LP12), L+LP12, 2
       INDIC = INDIC + 1
      DO 10 LP4=0, L4
      DO 10 LP34=IABS(L3-LP4), L3+LP4, 2
      DO 10 LMD2=IABS(L-LP34), L+LP34, 2
       NSTR  = 0
       XINT1 = DGLNQ(RAD08, 0.D0, RMIN, LEG10)
       NSTR  = NSTR + LEG10
       XINT2 = DGLNQ(RAD08, RMIN, RMAX, LEG11)
       NSTR  = NSTR + LEG11
       XINT3 = DGLGQ(RAD08, RMAX, ZETAMY, LAG07)
       RADAR(INDEX) = XINT1 + XINT2 + XINT3

       INDEX = INDEX + 1
 10   CONTINUE


      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE RETURNS THE RADIAL PART OF THE THREE-CENTER EXCHANGE *C
C* INTEGRALS OVER THE RANGE [0, +\INFTY]                             *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RAD08(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/EXCH031/ZETA3, N3, LP4, INDIC, LMD2, NSTR
      COMMON/COMV00/VL(LEA07*(LMDMAX+1)*(L1MAX+1)**3)
      COMMON/COMG01/GEGIN01(LMHERM*(LDEV+1))

      INDO = (INDIC-1)*LEA07 + NSTR
      INDP = LMD2*LEA07 + NSTR

      DO 5 I=1, N
       R2 = T(I)
       Y(I) = 1.D0/DSQRT(R2) * R2**(N3+LP4+1) * DEXP(-ZETA3*R2) *
     $                                    VL(INDO+I) * GEGIN01(INDP+I)
 5    CONTINUE

      RETURN
      END
