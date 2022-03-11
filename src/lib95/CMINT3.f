C.... THIRD  CASE : ONLY ONE ORBITAL IS REAL

      SUBROUTINE CMINT3(M1,M2,M3,M4, BOOL1,BOOL2,BOOL3,BOOL4,
     $             IY,JY,KY,LY, IX,JX,KX,LX,
     $             IP,IQ,IR,IS, JP,JQ,JR,JS, KP,KQ,KR,KS, LP,LQ,LR,LS,
     $                                       NTORB, NSUM, OUIZA, FATIMA)
      IMPLICIT REAL*8 (A-H, O-Z)
      COMPLEX*16 OUIZA
      COMPLEX*16 Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8
      COMPLEX*16 ZRE1,ZIM1, ZRE2,ZIM2, ZRE3,ZIM3, ZRE4,ZIM4
      LOGICAL*1 BOOL1, BOOL2, BOOL3, BOOL4

      DIMENSION OUIZA(20,20,20,20), FATIMA(*)
      DATA SQRT2/0.707106781186547524D0/



      Z1 = (-1)**(M1+M2+M3+M4) * OUIZA(IY, JY, KY, LY)
      Z8 =                       OUIZA(IX, JX, KX, LX)

      IF(BOOL1)THEN
       Z2 = (-1)**(M1+M2   +M4) * OUIZA(IY, JY, KX, LY)
       Z3 = (-1)**(M1   +M3+M4) * OUIZA(IY, JX, KY, LY)
       Z4 = (-1)**(M1      +M4) * OUIZA(IY, JX, KX, LY)

       Z5 = (-1)**(   M2+M3+M4) * OUIZA(IX, JY, KY, LY)
       Z6 = (-1)**(   M2   +M4) * OUIZA(IX, JY, KX, LY)
       Z7 = (-1)**(      M3+M4) * OUIZA(IX, JX, KY, LY)

       CALL COEFFS(MC2,MC3,MC4,MC5,MC6,MC7,MC8, -1,1,1,-1,-1,1,-1)
       CALL GNDAFF(IP,JP,KP,LP, JRUN1, IP,JP,KR,LR, JRUN2,
     $             IQ,JQ,KP,LP, JRUN3, IQ,JQ,KR,LR, JRUN4,
     $             IR,JR,KP,LP, JRUN5, IR,JR,KR,LR, JRUN6,
     $             IS,JS,KP,LP, JRUN7, IS,JS,KR,LR, JRUN8, NTORB,NSUM)

       CALL GNDAFF(KP,LP,IP,JP, KRUN1, KR,LR,IP,JP, KRUN2,
     $             KP,LP,IQ,JQ, KRUN3, KR,LR,IQ,JQ, KRUN4,
     $             KP,LP,IR,JR, KRUN5, KR,LR,IR,JR, KRUN6,
     $             KP,LP,IS,JS, KRUN7, KR,LR,IS,JS, KRUN8, NTORB,NSUM)

      ELSE
       IF(BOOL2)THEN
	Z2 = (-1)**(M1+M2+M3   ) * OUIZA(IY, JY, KY, LX)
	Z3 = (-1)**(M1   +M3+M4) * OUIZA(IY, JX, KY, LY)
	Z4 = (-1)**(M1   +M3   ) * OUIZA(IY, JX, KY, LX)

	Z5 = (-1)**(   M2+M3+M4) * OUIZA(IX, JY, KY, LY)
	Z6 = (-1)**(   M2+M3   ) * OUIZA(IX, JY, KY, LX)
	Z7 = (-1)**(      M3+M4) * OUIZA(IX, JX, KY, LY)

	CALL COEFFS(MC2,MC3,MC4,MC5,MC6,MC7,MC8, 1,1,-1,-1,1,1,1)
	CALL GNDAFF(IP,JP,KP,LP, JRUN1, IP,JP,KQ,LQ, JRUN2,
     $              IQ,JQ,KP,LP, JRUN3, IQ,JQ,KQ,LQ, JRUN4,
     $              IR,JR,KP,LP, JRUN5, IR,JR,KQ,LQ, JRUN6,
     $              IS,JS,KP,LP, JRUN7, IS,JS,KQ,LQ, JRUN8, NTORB,NSUM)

	CALL GNDAFF(KP,LP,IP,JP, KRUN1, KQ,LQ,IP,JP, KRUN2,
     $              KP,LP,IQ,JQ, KRUN3, KQ,LQ,IQ,JQ, KRUN4,
     $              KP,LP,IR,JR, KRUN5, KQ,LQ,IR,JR, KRUN6,
     $              KP,LP,IS,JS, KRUN7, KQ,LQ,IS,JS, KRUN8, NTORB,NSUM)


       ELSE
	IF(BOOL3)THEN
	 Z2 = (-1)**(M1+M2+M3   ) * OUIZA(IY, JY, KY, LX)
	 Z3 = (-1)**(M1+M2   +M4) * OUIZA(IY, JY, KX, LY)
	 Z4 = (-1)**(M1+M2      ) * OUIZA(IY, JY, KX, LX)

	 Z5 = (-1)**(   M2+M3+M4) * OUIZA(IX, JY, KY, LY)
	 Z6 = (-1)**(   M2+M3   ) * OUIZA(IX, JY, KY, LX)
	 Z7 = (-1)**(   M2   +M4) * OUIZA(IX, JY, KX, LY)

	 CALL COEFFS(MC2,MC3,MC4,MC5,MC6,MC7,MC8, 1,-1,1,-1,1,-1,-1)
	 CALL GNDAFF(IP,JP,KP,LP, JRUN1, IP,JP,KQ,LQ, JRUN2,
     $               IP,JP,KR,LR, JRUN3, IP,JP,KS,LS, JRUN4,
     $               IR,JR,KP,LP, JRUN5, IR,JR,KQ,LQ, JRUN6,
     $               IR,JR,KR,LR, JRUN7, IR,JR,KS,LS, JRUN8, NTORB,NSUM)

	 CALL GNDAFF(KP,LP,IP,JP, KRUN1, KQ,LQ,IP,JP, KRUN2,
     $               KR,LR,IP,JP, KRUN3, KS,LS,IP,JP, KRUN4,
     $               KP,LP,IR,JR, KRUN5, KQ,LQ,IR,JR, KRUN6,
     $               KR,LR,IR,JR, KRUN7, KS,LS,IR,JR, KRUN8, NTORB,NSUM)

	ELSE
	 Z2 = (-1)**(M1+M2+M3   ) * OUIZA(IY, JY, KY, LX)
	 Z3 = (-1)**(M1+M2   +M4) * OUIZA(IY, JY, KX, LY)
	 Z4 = (-1)**(M1+M2      ) * OUIZA(IY, JY, KX, LX)

	 Z5 = (-1)**(M1   +M3+M4) * OUIZA(IY, JX, KY, LY)
	 Z6 = (-1)**(M1   +M3   ) * OUIZA(IY, JX, KY, LX)
	 Z7 = (-1)**(M1      +M4) * OUIZA(IY, JX, KX, LY)

	 CALL COEFFS(MC2,MC3,MC4,MC5,MC6,MC7,MC8, 1,-1,1, 1,-1,1,1)
	 CALL GNDAFF(IP,JP,KP,LP, JRUN1, IP,JP,KQ,LQ, JRUN2,
     $               IP,JP,KR,LR, JRUN3, IP,JP,KS,LS, JRUN4,
     $               IQ,JQ,KP,LP, JRUN5, IQ,JQ,KQ,LQ, JRUN6,
     $               IQ,JQ,KR,LR, JRUN7, IQ,JQ,KS,LS, JRUN8, NTORB,NSUM)

	 CALL GNDAFF(KP,LP,IP,JP, KRUN1, KQ,LQ,IP,JP, KRUN2,
     $               KR,LR,IP,JP, KRUN3, KS,LS,IP,JP, KRUN4,
     $               KP,LP,IQ,JQ, KRUN5, KQ,LQ,IQ,JQ, KRUN6,
     $               KR,LR,IQ,JQ, KRUN7, KS,LS,IQ,JQ, KRUN8, NTORB,NSUM)

	ENDIF
       ENDIF
      ENDIF

      ZRE1 =       0.5D0*SQRT2 * (Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8)
      ZIM1 = MC2 * 0.5D0*SQRT2 * (Z1 - Z2 + Z3 - Z4 + Z5 - Z6 + Z7 - Z8)
      ZIM2 = MC3 * 0.5D0*SQRT2 * (Z1 + Z2 - Z3 - Z4 + Z5 + Z6 - Z7 - Z8)
      ZRE2 = MC4 * 0.5D0*SQRT2 * (Z1 - Z2 - Z3 + Z4 + Z5 - Z6 - Z7 + Z8)
      ZIM3 = MC5 * 0.5D0*SQRT2 * (Z1 + Z2 + Z3 + Z4 - Z5 - Z6 - Z7 - Z8)
      ZRE3 = MC6 * 0.5D0*SQRT2 * (Z1 - Z2 + Z3 - Z4 - Z5 + Z6 - Z7 + Z8)
      ZRE4 = MC7 * 0.5D0*SQRT2 * (Z1 + Z2 - Z3 - Z4 - Z5 - Z6 + Z7 + Z8)
      ZIM4 = MC8 * 0.5D0*SQRT2 * (Z1 - Z2 - Z3 + Z4 - Z5 + Z6 + Z7 - Z8)

      FATIMA(JRUN1) = DREAL(ZRE1)
      FATIMA(JRUN2) = DIMAG(ZIM1)

      FATIMA(JRUN3) = DIMAG(ZIM2)
      FATIMA(JRUN4) = DREAL(ZRE2)

      FATIMA(JRUN5) = DIMAG(ZIM3)
      FATIMA(JRUN6) = DREAL(ZRE3)

      FATIMA(JRUN7) = DREAL(ZRE4)
      FATIMA(JRUN8) = DIMAG(ZIM4)

      FATIMA(KRUN1) = FATIMA(JRUN1)
      FATIMA(KRUN2) = FATIMA(JRUN2)

      FATIMA(KRUN3) = FATIMA(JRUN3)
      FATIMA(KRUN4) = FATIMA(JRUN4)

      FATIMA(KRUN5) = FATIMA(JRUN5)
      FATIMA(KRUN6) = FATIMA(JRUN6)

      FATIMA(KRUN7) = FATIMA(JRUN7)
      FATIMA(KRUN8) = FATIMA(JRUN8)

      RETURN
      END



      SUBROUTINE GNDAFF(I1,J1,K1,L1,JRUN1, I2,J2,K2,L2,JRUN2,
     $                  I3,J3,K3,L3,JRUN3, I4,J4,K4,L4,JRUN4,
     $                  I5,J5,K5,L5,JRUN5, I6,J6,K6,L6,JRUN6,
     $                  I7,J7,K7,L7,JRUN7, I8,J8,K8,L8,JRUN8,NTORB,NSUM)
      IMPLICIT REAL*8 (A-H, O-Z)

      JRUN1 = (L1-1)*NTORB*NSUM + (K1-1)*NSUM + NSUM -
     $                           ((NTORB-J1+1)*(NTORB-J1+2))/2 + I1-J1+1
      JRUN2 = (L2-1)*NTORB*NSUM + (K2-1)*NSUM + NSUM -
     $                           ((NTORB-J2+1)*(NTORB-J2+2))/2 + I2-J2+1
      JRUN3 = (L3-1)*NTORB*NSUM + (K3-1)*NSUM + NSUM -
     $                           ((NTORB-J3+1)*(NTORB-J3+2))/2 + I3-J3+1
      JRUN4 = (L4-1)*NTORB*NSUM + (K4-1)*NSUM + NSUM -
     $                           ((NTORB-J4+1)*(NTORB-J4+2))/2 + I4-J4+1
      JRUN5 = (L5-1)*NTORB*NSUM + (K5-1)*NSUM + NSUM -
     $                           ((NTORB-J5+1)*(NTORB-J5+2))/2 + I5-J5+1
      JRUN6 = (L6-1)*NTORB*NSUM + (K6-1)*NSUM + NSUM -
     $                           ((NTORB-J6+1)*(NTORB-J6+2))/2 + I6-J6+1
      JRUN7 = (L7-1)*NTORB*NSUM + (K7-1)*NSUM + NSUM -
     $                           ((NTORB-J7+1)*(NTORB-J7+2))/2 + I7-J7+1
      JRUN8 = (L8-1)*NTORB*NSUM + (K8-1)*NSUM + NSUM -
     $                           ((NTORB-J8+1)*(NTORB-J8+2))/2 + I8-J8+1

      RETURN
      END


