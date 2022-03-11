C.... ZEROTH CASE : ALL ORBITALS ARE REAL

      SUBROUTINE CMINT0(M1,M2,M3,M4, IY,JY,KY,LY, IP,JP,KP,LP,
     $                                       NTORB, NSUM, OUIZA, FATIMA)
      IMPLICIT REAL*8 (A-H, O-Z)
      COMPLEX*16 OUIZA

      DIMENSION OUIZA(20,20,20,20), FATIMA(*)


      ZRE = DREAL(OUIZA(IY,JY,KY,LY))

      JRUN1  = ISIND(IP, JP, KP, LP, NTORB, NSUM)
      KRUN1  = ISIND(KP, LP, IP, JP, NTORB, NSUM)

      FATIMA(JRUN1) = ZRE
      FATIMA(KRUN1) = ZRE

      RETURN
      END