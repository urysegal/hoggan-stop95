C.... FIRST  CASE : ONLY ONE ORBITAL IS IMAGINARY

      SUBROUTINE CMINT1(M1,M2,M3,M4, BOOL1,BOOL2,BOOL3,BOOL4,
     $             IY,JY,KY,LY, IX,JX,KX,LX,
     $             IP,IQ,IR,IS, JP,JQ,JR,JS, KP,KQ,KR,KS, LP,LQ,LR,LS,
     $                                     NTORB, NSUM, OUIZA, FATIMA)
      IMPLICIT REAL*8 (A-H, O-Z)
      COMPLEX*16 OUIZA
      COMPLEX*16 Z1, Z2, ZRE, ZIM
      LOGICAL*1 BOOL1, BOOL2, BOOL3, BOOL4

      DIMENSION OUIZA(20,20,20,20), FATIMA(*)
      DIMENSION MC(0:3)
      DATA SQRT2/0.707106781186547524D0/


      IS1 = IX-IY
      IS3 = KX-KY

      MC(0) =  1
      MC(1) = -1
      MC(2) =  1
      MC(3) = -1

      Z1 = (-1)**(M1+M2+M3+M4) * OUIZA(IY, JY, KY, LY)
      Z2 =                       OUIZA(IX, JX, KX, LX)

      ZRE =                       SQRT2 * (Z1 + Z2)
      ZIM = MC(IS1) * MC(2+IS3) * SQRT2 * (Z1 - Z2)

      JRUN1 = ISIND(IP, JP, KP, LP, NTORB, NSUM)
      JRUN2 = ISIND(IS, JS, KS, LS, NTORB, NSUM)

      KRUN1 = ISIND(KP, LP, IP, JP, NTORB, NSUM)
      KRUN2 = ISIND(KS, LS, IS, JS, NTORB, NSUM)

      FATIMA(JRUN1) = DREAL(ZRE)
      FATIMA(JRUN2) = DIMAG(ZIM)

      FATIMA(KRUN1) = FATIMA(JRUN1)
      FATIMA(KRUN2) = FATIMA(JRUN2)

      RETURN
      END
