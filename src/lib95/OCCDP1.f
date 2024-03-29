      SUBROUTINE OCCDP1(L12MAX,INDIC,N12,L12,ZETAS,RLEAG0,INC,LEA00, VL)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION RLEAG0(*), VL(*)
      COMMON/TCCI0/ZETA12, NL12
      COMMON/GAUSOR/NORDER(0:99)
      EXTERNAL POT01

      ZETA12 = ZETAS

      NSTRT = INDIC * LEA00 * (L12MAX + 1)

      DO 5 I=INC+1, LEA00
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


