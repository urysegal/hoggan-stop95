CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS SUBROUTINE RETURNS THE REPULSION NUCLEUS ENERGY. NATOM IS    *C
C* THE TOTAL NUMBER OF ATOMS WHITH A CHARGE GIVEN BY NUMAT(I).       *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION ENUC(NATOM, NUMAT, XYZAT)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER AT1, AT2

      DIMENSION NUMAT(*), XYZAT(*), XYZV(3)

      ENUC = 0.D0

      WRITE(*, 1)

 1    FORMAT(27HAT1  AT2    DISTANCE (u.a) ,
     $ 18H  REPULSION ENERGY, /, 45(1H-), /)

      DO 5 AT1=1, NATOM-1
      DO 5 AT2=AT1+1, NATOM
       NAT1 = 3*(AT1-1) + 1
       NAT2 = 3*(AT2-1) + 1

       CALL CARTCO(XYZAT(NAT1), XYZAT(NAT2), XYZV, RVMOD, THETRV, PHIRV)

       ENUC = ENUC + NUMAT(AT1)*NUMAT(AT2) / RVMOD
       WRITE(*, 2)AT1, AT2, RVMOD, NUMAT(AT1)*NUMAT(AT2) / RVMOD
 2     format(I3, 2X, I3, 2X, 2(D16.10, 2X))
 5    CONTINUE

      RETURN
      END
