CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS SUBROUTINE RETURNS THE TWO-ELECTRON ENERGY, WHERE NATO IS    *C
C* THE TOTAL NUMBER OF ATOMIC ORBITALS                               *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TEENER(NATO, FATIMA, PMN, E1)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER SIG
      INCLUDE "../lib95/SIZE.INCL"
c
      DIMENSION FATIMA(*), PMN(N_ORB,N_ORB)
c
      NSUM = (NATO*(NATO+1))/2
c
      E1 = 0.0D0
c
      DO 5 MU=1, NATO
      DO 5 NU=1, NATO
       IP = MAX(MU, NU)
       JP = MIN(MU, NU)
      DO 5 LMD=1, NATO
      DO 5 SIG=1, NATO
       KP = MAX(LMD, SIG)
       LP = MIN(LMD, SIG)
c
       IQ = MAX(MU, SIG)
       JQ = MIN(MU, SIG)
c
       KQ = MAX(LMD, NU)
       LQ = MIN(LMD, NU)
c
       JRUN1 = ISIND(IP,JP,KP,LP, NATO, NSUM)
       JRUN2 = ISIND(IQ,JQ,KQ,LQ, NATO, NSUM)
c
       E1 = E1 + PMN(MU,NU) * PMN(LMD,SIG) *
     $           (FATIMA(JRUN1) - 0.5D0 * FATIMA(JRUN2))
 5    CONTINUE
c
      E1 = 0.5D0 * E1
c
      RETURN
      END
