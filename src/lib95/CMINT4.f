C.... FOURTH CASE : ALL THE ORBITALS ARE IMAGINARY

      SUBROUTINE CMINT4(M1,M2,M3,M4, IY,JY,KY,LY, IX,JX,KX,LX,
     $             IP,IQ,IR,IS, JP,JQ,JR,JS, KP,KQ,KR,KS, LP,LQ,LR,LS,
     $                                       NTORB, NSUM, OUIZA, FATIMA)
      IMPLICIT REAL*8 (A-H, O-Z)
      COMPLEX*16 OUIZA, ZINT, ZRES

      DIMENSION OUIZA(20,20,20,20), ZINT(16), MC(16), ZRES(16)
      DIMENSION FATIMA(*)
      COMMON/SIGNES/NPSIG(16,16)


      ZINT(1 ) = (-1)**(M1+M2+M3+M4) * OUIZA(IY, JY, KY, LY)
      ZINT(2 ) = (-1)**(M1+M2+M3   ) * OUIZA(IY, JY, KY, LX)
      ZINT(3 ) = (-1)**(M1+M2   +M4) * OUIZA(IY, JY, KX, LY)
      ZINT(4 ) = (-1)**(M1+M2      ) * OUIZA(IY, JY, KX, LX)

      ZINT(5 ) = (-1)**(M1   +M3+M4) * OUIZA(IY, JX, KY, LY)
      ZINT(6 ) = (-1)**(M1   +M3   ) * OUIZA(IY, JX, KY, LX)
      ZINT(7 ) = (-1)**(M1      +M4) * OUIZA(IY, JX, KX, LY)
      ZINT(8 ) = (-1)**(M1         ) * OUIZA(IY, JX, KX, LX)

      ZINT(9 ) = (-1)**(   M2+M3+M4) * OUIZA(IX, JY, KY, LY)
      ZINT(10) = (-1)**(   M2+M3   ) * OUIZA(IX, JY, KY, LX)
      ZINT(11) = (-1)**(   M2   +M4) * OUIZA(IX, JY, KX, LY)
      ZINT(12) = (-1)**(   M2      ) * OUIZA(IX, JY, KX, LX)

      ZINT(13) = (-1)**(      M3+M4) * OUIZA(IX, JX, KY, LY)
      ZINT(14) = (-1)**(      M3   ) * OUIZA(IX, JX, KY, LX)
      ZINT(15) = (-1)**(         M4) * OUIZA(IX, JX, KX, LY)
      ZINT(16) =                       OUIZA(IX, JX, KX, LX)

      MC(1 ) =  1
      MC(2 ) =  1
      MC(3 ) = -1
      MC(4 ) =  1

      MC(5 ) =  1
      MC(6 ) = -1
      MC(7 ) =  1
      MC(8 ) =  1

      MC(9 ) = -1
      MC(10) =  1
      MC(11) = -1
      MC(12) = -1

      MC(13) =  1
      MC(14) =  1
      MC(15) = -1
      MC(16) =  1

      DO 5 I=1, 16
       ZRES(I) = (0.D0, 0.D0)
       DO 10 J=1, 16
	ZRES(I) = ZRES(I) + NPSIG(I, J) * ZINT(J)
 10    CONTINUE
       ZRES(I) = MC(I) * 0.25D0 * ZRES(I)
 5    CONTINUE

C.... COMPUTATION OF THE INDICES FOR THE STORAGE

      JRUN1  = ISIND(IP, JP, KP, LP, NTORB, NSUM)
      JRUN2  = ISIND(IP, JP, KQ, LQ, NTORB, NSUM)
      JRUN3  = ISIND(IP, JP, KR, LR, NTORB, NSUM)
      JRUN4  = ISIND(IP, JP, KS, LS, NTORB, NSUM)

      JRUN5  = ISIND(IQ, JQ, KP, LP, NTORB, NSUM)
      JRUN6  = ISIND(IQ, JQ, KQ, LQ, NTORB, NSUM)
      JRUN7  = ISIND(IQ, JQ, KR, LR, NTORB, NSUM)
      JRUN8  = ISIND(IQ, JQ, KS, LS, NTORB, NSUM)

      JRUN9  = ISIND(IR, JR, KP, LP, NTORB, NSUM)
      JRUN10 = ISIND(IR, JR, KQ, LQ, NTORB, NSUM)
      JRUN11 = ISIND(IR, JR, KR, LR, NTORB, NSUM)
      JRUN12 = ISIND(IR, JR, KS, LS, NTORB, NSUM)

      JRUN13 = ISIND(IS, JS, KP, LP, NTORB, NSUM)
      JRUN14 = ISIND(IS, JS, KQ, LQ, NTORB, NSUM)
      JRUN15 = ISIND(IS, JS, KR, LR, NTORB, NSUM)
      JRUN16 = ISIND(IS, JS, KS, LS, NTORB, NSUM)

      FATIMA(JRUN1 ) = DREAL(ZRES(1 ))
      FATIMA(JRUN2 ) = DIMAG(ZRES(2 ))
      FATIMA(JRUN3 ) = DIMAG(ZRES(3 ))
      FATIMA(JRUN4 ) = DREAL(ZRES(4 ))

      FATIMA(JRUN5 ) = DIMAG(ZRES(5 ))
      FATIMA(JRUN6 ) = DREAL(ZRES(6 ))
      FATIMA(JRUN7 ) = DREAL(ZRES(7 ))
      FATIMA(JRUN8 ) = DIMAG(ZRES(8 ))

      FATIMA(JRUN9 ) = DIMAG(ZRES(9 ))
      FATIMA(JRUN10) = DREAL(ZRES(10))
      FATIMA(JRUN11) = DREAL(ZRES(11))
      FATIMA(JRUN12) = DIMAG(ZRES(12))

      FATIMA(JRUN13) = DREAL(ZRES(13))
      FATIMA(JRUN14) = DIMAG(ZRES(14))
      FATIMA(JRUN15) = DIMAG(ZRES(15))
      FATIMA(JRUN16) = DREAL(ZRES(16))

      KRUN1  = ISIND(KP, LP, IP, JP, NTORB, NSUM)
      KRUN2  = ISIND(KP, LP, IQ, JQ, NTORB, NSUM)
      KRUN3  = ISIND(KP, LP, IR, JR, NTORB, NSUM)
      KRUN4  = ISIND(KP, LP, IS, JS, NTORB, NSUM)

      KRUN5  = ISIND(KQ, LQ, IP, JP, NTORB, NSUM)
      KRUN6  = ISIND(KQ, LQ, IQ, JQ, NTORB, NSUM)
      KRUN7  = ISIND(KQ, LQ, IR, JR, NTORB, NSUM)
      KRUN8  = ISIND(KQ, LQ, IS, JS, NTORB, NSUM)

      KRUN9  = ISIND(KR, LR, IP, JP, NTORB, NSUM)
      KRUN10 = ISIND(KR, LR, IQ, JQ, NTORB, NSUM)
      KRUN11 = ISIND(KR, LR, IR, JR, NTORB, NSUM)
      KRUN12 = ISIND(KR, LR, IS, JS, NTORB, NSUM)

      KRUN13 = ISIND(KS, LS, IP, JP, NTORB, NSUM)
      KRUN14 = ISIND(KS, LS, IQ, JQ, NTORB, NSUM)
      KRUN15 = ISIND(KS, LS, IR, JR, NTORB, NSUM)
      KRUN16 = ISIND(KS, LS, IS, JS, NTORB, NSUM)

      FATIMA(KRUN1 ) = FATIMA(JRUN1 )
      FATIMA(KRUN2 ) = FATIMA(JRUN5 )
      FATIMA(KRUN3 ) = FATIMA(JRUN9 )
      FATIMA(KRUN4 ) = FATIMA(JRUN13)

      FATIMA(KRUN5 ) = FATIMA(JRUN2 )
      FATIMA(KRUN6 ) = FATIMA(JRUN6 )
      FATIMA(KRUN7 ) = FATIMA(JRUN10)
      FATIMA(KRUN8 ) = FATIMA(JRUN14)

      FATIMA(KRUN9 ) = FATIMA(JRUN3 )
      FATIMA(KRUN10) = FATIMA(JRUN7 )
      FATIMA(KRUN11) = FATIMA(JRUN11)
      FATIMA(KRUN12) = FATIMA(JRUN15)

      FATIMA(KRUN13) = FATIMA(JRUN4 )
      FATIMA(KRUN14) = FATIMA(JRUN8 )
      FATIMA(KRUN15) = FATIMA(JRUN12)
      FATIMA(KRUN16) = FATIMA(JRUN16)


      RETURN
      END


      FUNCTION ISIND(I,J,K,L, NTORB, NSUM)
      IMPLICIT REAL*8 (A-H, O-Z)


      LIN = (L-1) * NTORB * NSUM
      KIN = (K-1) * NSUM
      JIN = NSUM - ((NTORB-J+1)*(NTORB-J+2))/2
      IIN = I - J + 1

      ISIND = LIN + KIN + JIN + IIN

      RETURN
      END
