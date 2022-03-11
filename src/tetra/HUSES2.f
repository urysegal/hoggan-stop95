C.... THESE ARE SOME GAUNT COEFFICIENTS

      SUBROUTINE HUSES2(L3,M3,YLMAC, L4,M4,YLMAD, ARGNT3,ARGNT4,ARGNT5)
      INCLUDE "../lib95/SIZE.INCL"
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION ARGNT3(*), ARGNT4(*), ARGNT5(*)
      DIMENSION YLMAC(0:*), YLMAD(0:*)
      COMMON/XGONE/GONE((LMDMAX+1)**2 * (LMDMAX+6+1)**2, LMDMAX+3+1)
      COMMON/POWUN/YIWAN(0:100), LMAXX2

      CPIN = MCLOCK()
cccc      WRITE(*, *)'ENTERING THE JOURNEY '

       K4 = 0
       K6 = 0

      DO 5 LP3=0, L3
      DO 5 MP3=MAX0(-LP3, -L3+LP3+M3), MIN0(LP3, L3-LP3+M3)
       K4   = K4 + 1
       AUX4 = YIWAN(IABS(MP3)) * ARGNT3(K4)

       K5   = 0

      DO 5 LP4=0, L4
      DO 5 MP4=MAX0(-LP4, -L4+LP4+M4), MIN0(LP4, L4-LP4+M4)
       K5   = K5 + 1
       AUX5 = YIWAN(IABS(MP4)) * ARGNT4(K5) * AUX4

      DO 5 LP34=IABS(LP4-LP3), LP4+LP3, 2
       K6      = K6 + 1
       AUX6    = ARGNT5(K6) * AUX5
C       PR1(K6) = AUX6

 5    CONTINUE

      INDEX = 0

      DO 10 L=0, LMDMAX
       JJ1 = L * l * LMAXX2
       JJ2 = -LMAXX2

      DO 10 M=-L, L

       JJ2  = JJ2 + LMAXX2
       KPR1 = 0

      DO 10 LP3=0, L3
      DO 10 MP3=MAX0(-LP3, -L3+LP3+M3), MIN0(LP3, L3-LP3+M3)

      DO 10 LP4=0, L4
      DO 10 MP4=MAX0(-LP4, -L4+LP4+M4), MIN0(LP4, L4-LP4+M4)

      DO 10 LP34=IABS(LP4-LP3), LP4+LP3, 2

       KPR1 = KPR1 + 1
C      AUX6 = PR1(KPR1)

       JJ3  = LP34 * LP34
       JJ4  = MP3-MP4+LP34+1
       JTT  = JJ1 + JJ2 + JJ3 + JJ4
       JAL  = 1

      DO 10 LMD=IABS(L-LP34), L+LP34, 2

       AUX7 = GONE(JTT, JAL) * AUX6
       JAL  = JAL + 1

       KK3  = LMD*LMD
       KK4  = M-(MP3-MP4) + LMD + 1

      DO 10 LMD4=0, LMDMAX

       NYL4L = (LMD4*(LMD4+1))/2

       KK1   = LMD4 * LMD4 * LMAXX2
       KK2   = -LMAXX2

      DO 10 MU4=-LMD4, LMD4

       KK2 = KK2 + LMAXX2

       IMU4 = IABS(MU4)
       NYL4 = NYL4L + IMU4
       MUPW = (MU4 - IMU4)/2
       AUX8 = YIWAN(IABS(MU4)) * YLMAD(NYL4) * AUX7

       KTT  = KK1 + KK2 + KK3 + KK4
       KAL  = 1

       MU3  = MU4 - M + MP3 - MP4
       IMU3 = IABS(MU3)
       MUPW = (MU3 - IMU3)/2

       LMD3X = LMD+LMD4
       LMD3M = 2*((LMD3X + MAX0(IABS(LMD-LMD4), IMU3)) /2 ) - LMD3X

      DO 10 LMD3=IABS(LMD-LMD4), LMD+LMD4, 2

       IF(LMD3 .GE. LMD3M .AND. AUX8 .NE. 0.D0)THEN
	NYL3 = (LMD3*(LMD3+1))/2 + IMU3
	AUX9 = YIWAN(IABS(MUPW)) * YLMAC(NYL3) * GONE(KTT, KAL) * AUX8
C       PR2(INDEX) = AUX9
       ELSE
C       PR2(INDEX) = 0.D0
       ENDIF

       KAL   = KAL + 1
       INDEX = INDEX + 1

 10   CONTINUE

      CPOUT = MCLOCK()
cccc      WRITE(*, *)'CA C EST L INDICE ', INDEX, 1.D-2*(CPOUT-CPIN)

      STOP
      RETURN
      END
