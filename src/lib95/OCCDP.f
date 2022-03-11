      SUBROUTINE OCCDP(L12MAX, NORB1, NLMA, ZETAA, RLEAG0,INC,LEA00, VL)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION NLMA(*), ZETAA(*), RLEAG0(*), VL(*)
      COMMON/TCCI0/ZETA12, NL12
      COMMON/GAUSOR/NORDER(0:99)
      EXTERNAL POT01

      INDIC = 0

      DO 5 NAT1=1, NORB1
       IST   = 5*(NAT1-1)
       N1    = NLMA(IST+3)
       L1    = NLMA(IST+4)
       M1    = NLMA(IST+5)
       ZETA1 = ZETAA(NAT1)

      DO 5 NAT2=NAT1, NORB1
       JST   = 5*(NAT2-1)
       N2    = NLMA(JST+3)
       L2    = NLMA(JST+4)
       M2    = NLMA(JST+5)
       ZETA2 = ZETAA(NAT2)

       N12    = N1 + N2
       L12    = L1 + L2
       ZETA12 = ZETA1+ZETA2

       NSTRT = INDIC * LEA00 * (L12MAX + 1)
       INDIC = INDIC + 1

      DO 5 I=1, INC
       R2 = RLEAG0(I)

       ZETASR = ZETA12*R2
       EXPR   = DEXP(-ZETASR)
       R2N    = R2**N12
       R2L    = R2**L12

       NL12   = N12+L12
       EXPR1  = R2N * EXPR
       VOX10 = 1.D0/(R2*R2L) * DGLNQ(POT01, 0.D0, R2, 20)

       NL12  = N12-L12-1
       EXPR2 = R2N/ZETASR * EXPR
       VOX20 = R2L * DGLGQ(POT01, R2, ZETA12, NORDER(NL12/2+1))

       INDEX = NSTRT + (I-1)*(L12+1) + L12

       VL(INDEX+1) = VOX10 + VOX20

      DO 5 L=0, L12-1
       LP    = L12-1-L
       NL12  = N12 + LP + 1
       VOX11 = 1.D0/DBLE(NL12) * (ZETASR * VOX10 + EXPR1)

       NL12  = N12 - LP - 1
       VOX21 = DBLE(NL12)/ZETASR * VOX20 + EXPR2

       VL(INDEX-L) = VOX11 + VOX21

       VOX10 = VOX11
       VOX20 = VOX21

 5    CONTINUE

      RETURN
      END



      SUBROUTINE POT01(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION T(*), Y(*)
      COMMON/TCCI0/ZETA12, NL12

      DO 5 I=1, N
       R    = T(I)
       Y(I) = R**NL12 * DEXP(-ZETA12*R)
 5    CONTINUE

      RETURN
      END
