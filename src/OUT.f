      SUBROUTINE OUT(MOLEC, NTORB, NLMAT, ZETAAT, NATOM, XYZAT, NUMAT)
      IMPLICIT REAL*8 (A-H, O-Z)
      CHARACTER*2 ATOMS, XYZC
      CHARACTER*60 MOLEC

      DIMENSION NLMAT(*), ZETAAT(*), XYZAT(*), NUMAT(*)
      DIMENSION ATOMS(36), XYZC(3)

      DATA XYZC /'X ', 'Y ', 'Z '/
      DATA ATOMS/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     $                     'Na','Mg','Al','Si',' P',' S','Cl','Ar',
     $                     ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co',
     $                     'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr'/


C.... THE ADVERTS

c     WRITE(*, 1)
c1    FORMAT(80(1HC), /, 1HC, 78X, 1HC, /, 1HC, 16X,
c    $ 46HANOTHER FANTASTIC RUN FOR STOP  (RELEASE - 95), 16X, 1HC, /,
c    $ 1HC, 78X, 1HC, /, 1HC, 38X, 2HBY, 38X, 1HC, /, 1HC, 78X, 1HC, /,
c    $ 1HC, 26X, 26H** AHMED ** BOUFERGUENE **, 26X, 1HC, /,
c    $ 1HC, 78X, 1HC, /, 1HC, 9X,
c    $ 60HLABORATOIRE DE CATALYSE ET SPECTROCHIMIE, UNIVERSITE DE CAEN,
c    $ 9X, 1HC, /, 1HC, 78X, 1HC, /, 80(1HC))

C.... THE MOLECULAR SYSTEM

      WRITE(*, '(A60)')MOLEC

C.... THE ATOMIC ORBITALS USED

      WRITE(*, 2)
 2    FORMAT(4X, 'N', 2X, 'L', 2X, 'M', 6X, 'ZETA', /, 25(1H-))

      DO 5 I=1, NTORB
       NDEX = 5*(I-1)
       WRITE(*, 3)NLMAT(NDEX+2),NLMAT(NDEX+3),
     $                          NLMAT(NDEX+4),NLMAT(NDEX+5), ZETAAT(I)
 5    CONTINUE

 3    FORMAT(Z2, 1X, 3(I2, 1X), 1X, D12.6)

      WRITE(*, 4)
 4    FORMAT(25(1H-), //, 25H COORDINATES AND CHARGES , /, 29(1H-))

C.... OUTPUTING THE SYMBOLES OF THE ATOMS INVOLVED IN THE SYSTEM

      WRITE(*, 6)(ATOMS(NUMAT(I)), I=1, NATOM)
 6    FORMAT(5X, 12(A2, 7X))

C.... OUTPUTING THE CHARGES OF THE ATOMS

      WRITE(*, 7)(NUMAT(I), I=1, NATOM)
 7    FORMAT(5X, 12(I2, 7X))

C.... OUTPUTING THE CARTESIAN COORDINATES OF THE ATOMS (CENTERS)

      DO 10 I=1, 3
       WRITE(*, 8)XYZC(I), (XYZAT(3*(J-1)+I), J=1, NATOM)
 10   CONTINUE

 8    FORMAT(A2, 12(F8.4, 1X))

      WRITE(*, 9)
 9    FORMAT(29(1H-))
      RETURN
      END
