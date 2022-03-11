C.... SECOND CASE : TWO ORBITALS ARE REAL

      SUBROUTINE CMINT2(M1,M2,M3,M4,BOOL1,BOOL2,BOOL3,BOOL4,BOOL5,BOOL6,
     $             IY,JY,KY,LY, IX,JX,KX,LX,
     $             IP,IQ,IR,IS, JP,JQ,JR,JS, KP,KQ,KR,KS, LP,LQ,LR,LS,
     $                                       NTORB, NSUM, OUIZA, FATIMA)

      IMPLICIT REAL*8 (A-H, O-Z)
      COMPLEX*16 OUIZA
      COMPLEX*16 Z1,Z2,Z3,Z4
      COMPLEX*16 ZRE1,ZIM1, ZRE2,ZIM2
      LOGICAL*1 BOOL1, BOOL2, BOOL3, BOOL4, BOOL5, BOOL6

      DIMENSION OUIZA(20,20,20,20), FATIMA(*)
      DATA SQRT2/0.707106781186547524D0/


      Z1 = (-1)**(M1+M2+M3+M4) * OUIZA(IY, JY, KY, LY)
      Z4 =                       OUIZA(IX, JX, KX, LX)

      IF(BOOL1)THEN
       Z2  = (-1)**(M1   +M3+M4) * OUIZA(IY, JX, KY, LY)
       Z3  = (-1)**(   M2+M3+M4) * OUIZA(IX, JY, KY, LY)
       MC2 =  1
       MC3 = -1
       MC4 =  1
       CALL HNDAFF(IP,JP,KP,LP, JRUN1, IQ,JQ,KP,LP, JRUN2,
     $             IR,JR,KP,LP, JRUN3, IS,JS,KP,LP, JRUN4, NTORB,NSUM)

       CALL HNDAFF(KP,LP,IP,JP, KRUN1, KP,LP,IQ,JQ, KRUN2,
     $             KP,LP,IR,JR, KRUN3, KP,LP,IS,JS, KRUN4, NTORB,NSUM)

      ELSE
       IF(BOOL2)THEN
	Z2  = (-1)**(M1+M2   +M4) * OUIZA(IY, JY, KX, LY)
	Z3  = (-1)**(   M2+M3+M4) * OUIZA(IX, JY, KY, LY)
	MC2 = -1
	MC3 = -1
	MC4 = -1
	CALL HNDAFF(IP,JP,KP,LP, JRUN1, IP,JP,KR,LR, JRUN2,
     $              IR,JR,KP,LP, JRUN3, IR,JR,KR,LR, JRUN4, NTORB,NSUM)

	CALL HNDAFF(KP,LP,IP,JP, KRUN1, KR,LR,IP,JP, KRUN2,
     $              KP,LP,IR,JR, KRUN3, KR,LR,IR,JR, KRUN4, NTORB,NSUM)
       ELSE
	IF(BOOL3)THEN
	 Z2  = (-1)**(M1+M2+M3   ) * OUIZA(IY, JY, KY, LX)
	 Z3  = (-1)**(   M2+M3+M4) * OUIZA(IX, JY, KY, LY)
	 MC2 =  1
	 MC3 = -1
	 MC4 =  1
	 CALL HNDAFF(IP,JP,KP,LP, JRUN1, IP,JP,KQ,LQ, JRUN2,
     $               IR,JR,KP,LP, JRUN3, IR,JR,KQ,LQ, JRUN4, NTORB,NSUM)

	 CALL HNDAFF(KP,LP,IP,JP, KRUN1, KQ,LQ,IP,JP, KRUN2,
     $               KP,LP,IR,JR, KRUN3, KQ,LQ,IR,JR, KRUN4, NTORB,NSUM)
	ELSE
	 IF(BOOL4)THEN
	  Z2  = (-1)**(M1+M2   +M4) * OUIZA(IY, JY, KX, LY)
	  Z3  = (-1)**(M1   +M3+M4) * OUIZA(IY, JX, KY, LY)
	  MC2 = -1
	  MC3 =  1
	  MC4 =  1
	  CALL HNDAFF(IP,JP,KP,LP, JRUN1, IP,JP,KR,LR, JRUN2,
     $               IQ,JQ,KP,LP, JRUN3, IQ,JQ,KR,LR, JRUN4, NTORB,NSUM)

	  CALL HNDAFF(KP,LP,IP,JP, KRUN1, KR,LR,IP,JP, KRUN2,
     $               KP,LP,IQ,JQ, KRUN3, KR,LR,IQ,JQ, KRUN4, NTORB,NSUM)
	 ELSE
	  IF(BOOL5)THEN
	   Z2  = (-1)**(M1+M2+M3   ) * OUIZA(IY, JY, KY, LX)
	   Z3  = (-1)**(M1   +M3+M4) * OUIZA(IY, JX, KY, LY)
	   MC2 =  1
	   MC3 =  1
	   MC4 = -1
	   CALL HNDAFF(IP,JP,KP,LP, JRUN1, IP,JP,KQ,LQ, JRUN2,
     $               IQ,JQ,KP,LP, JRUN3, IQ,JQ,KQ,LQ, JRUN4, NTORB,NSUM)

	   CALL HNDAFF(KP,LP,IP,JP, KRUN1, KQ,LQ,IP,JP, KRUN2,
     $               KP,LP,IQ,JQ, KRUN3, KQ,LQ,IQ,JQ, KRUN4, NTORB,NSUM)
	  ELSE
	   Z2 = (-1)**(M1+M2+M3   ) * OUIZA(IY, JY, KY, LX)
	   Z3 = (-1)**(M1+M2   +M4) * OUIZA(IY, JY, KX, LY)
	   MC2 =  1
	   MC3 = -1
	   MC4 =  1
	   CALL HNDAFF(IP,JP,KP,LP, JRUN1, IP,JP,KQ,LQ, JRUN2,
     $              IP,JP,KR,LR, JRUN3, IP,JP,KS,LS, JRUN4, NTORB,NSUM)

	   CALL HNDAFF(KP,LP,IP,JP, KRUN1, KQ,LQ,IP,JP, KRUN2,
     $              KR,LR,IP,JP, KRUN3, KS,LS,IP,JP, KRUN4, NTORB,NSUM)
	  ENDIF
	 ENDIF
	ENDIF
       ENDIF
      ENDIF

      ZRE1 =       0.5D0 * (Z1 + Z2 + Z3 + Z4)
      ZIM1 = MC2 * 0.5D0 * (Z1 - Z2 + Z3 - Z4)
      ZIM2 = MC3 * 0.5D0 * (Z1 + Z2 - Z3 - Z4)
      ZRE2 = MC4 * 0.5D0 * (Z1 - Z2 - Z3 + Z4)

      FATIMA(JRUN1) = DREAL(ZRE1)
      FATIMA(JRUN2) = DIMAG(ZIM1)
      FATIMA(JRUN3) = DIMAG(ZIM2)
      FATIMA(JRUN4) = DREAL(ZRE2)

      FATIMA(KRUN1) = FATIMA(JRUN1)
      FATIMA(KRUN2) = FATIMA(JRUN2)
      FATIMA(KRUN3) = FATIMA(JRUN3)
      FATIMA(KRUN4) = FATIMA(JRUN4)

      RETURN
      END



      SUBROUTINE HNDAFF(I1,J1,K1,L1, JRUN1, I2,J2,K2,L2, JRUN2,
     $           I3,J3,K3,L3, JRUN3, I4,J4,K4,L4, JRUN4, NTORB,NSUM)
      IMPLICIT REAL*8 (A-H, O-Z)
c
      JRUN1 = (L1-1)*NTORB*NSUM + (K1-1)*NSUM + NSUM -
     $                  ((NTORB-J1+1)*(NTORB-J1+2))/2 + I1-J1+1
      JRUN2 = (L2-1)*NTORB*NSUM + (K2-1)*NSUM + NSUM -
     $                  ((NTORB-J2+1)*(NTORB-J2+2))/2 + I2-J2+1
      JRUN3 = (L3-1)*NTORB*NSUM + (K3-1)*NSUM + NSUM -
     $                  ((NTORB-J3+1)*(NTORB-J3+2))/2 + I3-J3+1
      JRUN4 = (L4-1)*NTORB*NSUM + (K4-1)*NSUM + NSUM -
     $                  ((NTORB-J4+1)*(NTORB-J4+2))/2 + I4-J4+1
c
      RETURN
      END






