CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE RETURNS THE RADIAL PART INVOLVED IN THE TWO-CENTER   *C
C* CHARGE DISTRIBUTION POTENTIAL.                                    *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TCCDP2(LM123, N1,L1, N2,L2, ZETA1,ZETA2, AB,RMIN,RMAX,
     $                      RLEG01, RLEG02, BESF01, X4NL, RLEAG0, VL)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"

      DIMENSION RLEG01(*), RLEG02(*), BESF01(*), X4NL(*), RLEAG0(*),
     $          VL(*)
      COMMON/BESHER/HEREX(LTOT6*LMHERM), BESS01(LTOT6*(LDEV+1))
      COMMON/XROOTS/ROOTS(LEA07*LEA08), ISP
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))
      COMMON/EXCH03/ZETA1P, R2, N1P, LP2, L, LMD1, LMDL12, IO, ISTR
      EXTERNAL FGLE08, FGLA08

C.... PASSAGE TO COMMON

      ZETA1P = ZETA1
      N1P    = N1
      LMDL12 = LMDMAX+LM123+L2+1

C.... COLLECTING THE ROOTS OF THE GAUSS LAGUERRE QUADRATURE FOR THE T.C.C.D.P

      DO 5 I=1, LEA07
       R2   = RLEAG0(I)
       ISP  = (I-1)*LEA08 + LEA08 - LAG08
       XXXX = DGLGQ(FGLA08, DMAX1(RMAX, R2), ZETA1+ZETA2, LAG08)
 5    CONTINUE


C.... COMPUTATION OF THE HERMIT POLYNOMIALS INVOLVED IN THE GEGENBAUER A. THEO.
C.... THE FOLLOWING VALUES DEPEND ON THE PARAMETERS OF THE TWO STO's (1) & (2)
C.... AND OF COURSE ON THE ROOTS OF THE GAUSS-LAGUERRE QUADRATURE (LAG08 ROOTS)

      NONE = - 1
      INP  = 1

      CALL ABIN02(ZETA2, AB, NONE, RLEG01,L1_3, RLEG02,L4_6,
     $                    RLEAG0,INP,LEA07, BESF01, X4NL, HEREX, BESS01)

      IN = LEA08-LAG08+1

      DO 10 I=1, LEA07
       ISIT = (I-1)*LEA08
       CALL HERMIT(ZETA2, AB, N2-L2, RLEG01,L1_3, RLEG02,L4_6,
     $              ROOTS(ISIT+1),IN,LEA08, BESF01, X4NL, HEREX, BESS01)
      DO 10 LMD=0, LMDL12-1
       IND1 = 1
      DO 10 J=IN, LEA08
       R1   = ROOTS(ISIT+J)
       IND2 = LMD+1
       CALL BAD01(IND1, IND2, ZETA2, AB, R1, N2-L2, LMD, XYINT)
       INDEX = (I-1)*LMDL12*LEA08 + LMD*LEA08 + J
       GEGIN0(INDEX) = XYINT
 10   CONTINUE

C.... COMPUTATION OF THE RADIAL PART OF THE T.C.C.D.P

      INDEX = 1

      DO 15 L=0, LMDMAX
      DO 15 LP2=0, L2
      DO 15 LP12=IABS(L1-LP2), L1+LP2, 2
      DO 15 LMD1=IABS(L-LP12), L+LP12, 2
       CALL TIZIRI(INDEX, ZETA1,ZETA2, AB,RMIN,RMAX, RLEAG0, VL)
 15   CONTINUE


      RETURN
      END




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS SUBROUTINE COMPUTES THE RADIAL INTEGRAL OCCURING IN THE          *C
C* T.C.C.D.P OVER THE RANGE [0, R2].                                     *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TCCDP3(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/EXCH03/ZETA1, R2, N1, LP2, L, LMD1, LMDL12, IO, ISTR
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))
      COMMON/LEGER/NNN

      INDO = (IO-1)*LMDL12*LEA08 + LMD1*LEA08 + ISTR

      DO 5 I=1, N
       R1 = T(I)
       Y(I) = 1.D0/DSQRT(R1) * R1**(N1+LP2) * (R1/R2)**(L+1) *
     $                                  DEXP(-ZETA1*R1) * GEGIN0(INDO+I)
 5    CONTINUE

      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS SUBROUTINE COMPUTES THE RADIAL INTEGRAL OCCURING IN THE          *C
C* T.C.C.D.P OVER THE RANGE [R2, +\INFTY).                               *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TCCDP4(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/EXCH03/ZETA1, R2, N1, LP2, L, LMD1, LMDL12, IO, ISTR
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))


      INDO = (IO-1)*LMDL12*LEA08 + LMD1*LEA08 + ISTR

      DO 5 I=1, N
       R1   = T(I)
       Y(I) = 1.D0/DSQRT(R1) * R1**(N1+LP2) * (R2/R1)**L *
     $                                  DEXP(-ZETA1*R1) * GEGIN0(INDO+I)
 5    CONTINUE

      RETURN
      END

