cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C.... FOUR CENTER INTEGRALS
      SUBROUTINE FCINTS(NEPSIL, N1,L1,M1, N2,L2,M2, N3,L3,M3, N4,L4,M4,
     $                  ZETA1,ZETA2,ZETA3,ZETA4,
     $                  AB,PHIAB,YLMAB, AC,PHIAC,YLMAC, AD,PHIAD,YLMAD,
     $                  ARGNT1, ARGNT2, ARGNT3, ARGNT4, ARGNT5,
     $                  RADAR, TETRAR, TETRAI)

      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER YIWAN
      REAL*4 GONE
      INCLUDE "../lib95/SIZE.INCL"
      INCLUDE "../lib95/CONST.INCL"
c      PARAMETER IGN1=(LMDMAX+1)**2 * (LMDMAX+6+1)**2
c      PARAMETER IGN2=(LMDMAX+1)**2* (LMDMAX+3+1)
      DIMENSION YLMAB(0:*), YLMAC(0:*), YLMAD(0:*)
      DIMENSION ARGNT1(*),ARGNT2(*),ARGNT3(*),ARGNT4(*),ARGNT5(*)
      DIMENSION RADAR(*)
      DIMENSION ARGNTA(100), ARGNTB(100), ARGNTC(100), INX(10),
     $                                ARVALI(0:100), ARVALR(0:100)
      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:20)
      COMMON/XGONE/GONE((LMDMAX+1)**2 * (LMDMAX+6+1)**2, LMDMAX+3+1)
      COMMON/POWUN/YIWAN(0:100), LMAXX2
c      IGN1=(LMDMAX+1)**2 * (LMDMAX+6+1)**2
c      IGN2=(LMDMAX+1)**2* (LMDMAX+3+1)
c      COMMON/XGONE/GONE(IGN1,IGN2)
C     WRITE(*, 3)N1,L1,M1, N2,L2,M2, N3,L3,M3, N4,L4,M4
C3    FORMAT(4(3(I3, 1X), 3X))
      IGN1=(LMDMAX+1)**2 * (LMDMAX+6+1)**2
      IGN2= (LMDMAX+3+1)
c
      A1 = (2.d0*ZETA1)**N1 * DSQRT((2.D0*ZETA1)/FACT(2*N1))
      A2 = (2.d0*ZETA2)**N2 * DSQRT((2.D0*ZETA2)/FACT(2*N2))
      A3 = (2.d0*ZETA3)**N3 * DSQRT((2.D0*ZETA3)/FACT(2*N3))
      A4 = (2.d0*ZETA4)**N4 * DSQRT((2.D0*ZETA4)/FACT(2*N4))
      PI7 = PI5*PI2

c      CSTE = 16384.D0*PI7 * DFAC(L2)*DFAC(L3)*DFAC(L4)*A1*A2*A3*A4*
      CSTE = PI7 * DFAC(L2)*DFAC(L3)*DFAC(L4)*A1*A2*A3*A4*
     $       (-AB)**L2 * (-AC)**L3 * (-AD)**L4
      CALL RAZ0(ARVALR, 0, LMDMAX)
      CALL RAZ0(ARVALI, 0, LMDMAX)

      CALL RAZ1(INX, 1, 10)

      DO 5 L=0, LMDMAX

       AUX0  = 1.D0/DBLE(L+L+1)

       INX(1) = INX(1)+INX(2)+INX(3)+INX(4)+INX(5) +
     $                      INX(6)+INX(7)+INX(8)+INX(9)+INX(10)

       II1 = L*L  * LMAXX2
       JJ1 = II1

       II2 = -LMAXX2
       JJ2 = II2

      DO 5 M=-L, L

       II2 = II2 + LMAXX2
       JJ2 = II2

       ABANG = (M2-M1-M) * PHIAB

       CALL RAZ1(INX, 2, 10)

       K1 = 0
       K2 = 0

      DO 5 LP2=0, L2
       INX(2) = INX(2)+INX(3)+INX(4)+INX(5) +
     $                      INX(6)+INX(7)+INX(8)+INX(9)+INX(10)

      DO 5 MP2=MAX(-LP2, -L2+LP2+M2), MIN(LP2, L2-LP2+M2)

       K1   = K1 + 1
       AUX1 = ARGNT1(K1) * AUX0
       CALL RAZ1(INX, 3, 10)

      DO 5 LP12=IABS(LP2-L1), L1+LP2, 2
       K2     = K2 + 1
       AUX2   = ARGNT2(K2) * AUX1
       INX(3)=INX(3)+INX(4)+INX(5)+INX(6)+INX(7)+INX(8)+INX(9)+INX(10)

       CALL RAZ1(INX, 4, 10)

       II3  = LP12*LP12
       II4  = MP2-M1+LP12+1
       ITT  = II1 + II2 + II3 +II4

       IAL  = 1
      DO 5 LMD2=IABS(L-LP12), L+LP12, 2

       NYL2 = (LMD2*(LMD2+1))/2 + IABS(M-(MP2-M1))
       MUPW = (M-(MP2-M1) - IABS(M-(MP2-M1)))/2
       AUX3 = YIWAN(IABS(MUPW)) * YLMAB(NYL2) * GONE(ITT, IAL)*AUX2
c       write(*,*)"FG",gone(itt,ial),itt,ial
       IAL  = IAL + 1

       INX(4) = INX(4)+INX(5)+INX(6)+INX(7)+INX(8)+INX(9)+INX(10)
       CALL RAZ1(INX, 5, 10)

       K4 = 0
       K6 = 0

      DO 5 LP3=0, L3
       INX(5) = INX(5)+INX(6)+INX(7)+INX(8)+INX(9)+INX(10)

      DO 5 MP3=MAX(-LP3, -L3+LP3+M3), MIN(LP3, L3-LP3+M3)
       K4   = K4 + 1
       AUX4 = YIWAN(IABS(MP3)) * ARGNT3(K4) * AUX3

       CALL RAZ1(INX, 6, 10)
       K5 = 0

      DO 5 LP4=0, L4
       INX(6) = INX(6)+INX(7)+INX(8)+INX(9)+INX(10)

      DO 5 MP4=MAX(-LP4, -L4+LP4+M4), MIN(LP4, L4-LP4+M4)
       K5   = K5 + 1
       AUX5 = YIWAN(IABS(MP4)) * ARGNT4(K5) * AUX4

       CALL RAZ1(INX, 7, 10)

      DO 5 LP34=IABS(LP4-LP3), LP4+LP3, 2
       K6   = K6 + 1
       AUX6 = ARGNT5(K6) * AUX5

       INX(7) = INX(7)+INX(8)+INX(9)+INX(10)
       CALL RAZ1(INX, 8, 10)

       JJ3 = LP34*LP34
       JJ4 = MP3-MP4+LP34+1
       JTT = JJ1 + JJ2 + JJ3 + JJ4
       JAL = 1
      DO 5 LMD=IABS(L-LP34), L+LP34, 2

       KK3 = LMD*LMD
c       KK4 = M-(MP3-MP4) + LMD + 1
c Index redone to match JJ4
c Modif Hoggan 24.09.2003 to be verified
c
c controle underflow
c
       KK4 = MP3-MP4 + LMD +1
       AUX7 = GONE(JTT, JAL) * AUX6

       JAL  = JAL + 1

       INX(8) = INX(8)+INX(9)+INX(10)
       CALL RAZ1(INX, 9, 10)

      DO 5 LMD4=0, LMDMAX
       INX(9) = INX(9)+INX(10)
       KK1    = LMD4 * LMD4 * LMAXX2
       KK2    = -LMAXX2

       NYL4L  = (LMD4*(LMD4+1))/2

      DO 5 MU4=-LMD4, LMD4

       KK2   = KK2 + LMAXX2

       ADANG = (M4 - MP4 + MU4) * PHIAD

       NYL4 = NYL4L + IABS(MU4)
       MUPW = (MU4 - IABS(MU4))/2
cc       If(DABS(AUX7).LT.1.0d-15)AUX7=1.0d0
       AUX8 = YIWAN(IABS(MUPW)) * YLMAD(NYL4) * AUX7

       ACANG  = (M - MU4 + MP4 - M3) * PHIAC
       ANGLE  = ABANG + ACANG + ADANG
c        write(*,*)"FNZ",AUX8,AUX7,AUX6

       IF(AUX8 .NE. 0.D0)THEN
        COSINE = DCOS(ANGLE)
	SINE   = DSIN(ANGLE)
c        write(*,*)"FNZ",cosine,sine	
	AUX8R  = COSINE
	AUX8I  = SINE
       ELSE
	AUX8R = 0.D0
	AUX8I = 0.D0
       ENDIF
c       write(*,*)"FNZ",AUX7,AUX8,AUX8R,AUX8I
       CALL RAZ1(INX, 10, 10)

       KTT = KK1 + KK2 + KK3 + KK4
c       write(*,*)"FK",KTT
cb       KTT = 2
cc       KAL = 1
       MU3    = MU4 - M + MP3 - MP4
       MUPW   = (MU3 - IABS(MU3))/2
       IMU3   = IABS(MU3)

       LMD3X = LMD+LMD4
       LMD3M = 2*((LMD3X + MAX0(IABS(LMD-LMD4), IMU3)+1)/2)-LMD3X
      DO 5 LMD3=IABS(LMD-LMD4), LMD3X, 2
       INX(10) = INX(10) + 1
       INDIC   = INX(1)+INX(2)+INX(3)+INX(4)+INX(5) +
     $                     INX(6)+INX(7)+INX(8)+INX(9)+INX(10)
       IF(LMD3 .GE. LMD3M)THEN
	IF(AUX8R*AUX8I .NE. 0.D0)THEN
	 NYL3 = (LMD3*(LMD3+1))/2 + IMU3
	 AUX9 = YIWAN(IABS(MUPW)) * YLMAC(NYL3) * GONE(KTT, KAL)
c         write(*,*)"FIN",GONE(KTT,KAL),KTT,KAL
c         If(KTT.EQ.2)STOP

         ARVALR(L) = ARVALR(L) + AUX8R*RADAR(INDIC)
         ARVALI(L) = ARVALI(L) + AUX8I*RADAR(INDIC)
c         ARVALR(L) = ARVALR(L) +  AUX9*RADAR(INDIC)
c  	 ARVALI(L) = ARVALI(L) +  AUX9*RADAR(INDIC)
	ELSE
	 IF(AUX8R .NE. 0.D0)THEN
	  NYL3 = (LMD3*(LMD3+1))/2 + IMU3
	  AUX9 = YIWAN(IABS(MUPW)) * YLMAC(NYL3) * GONE(KTT, KAL)
	  ARVALR(L) = ARVALR(L) + AUX8R * RADAR(INDIC)
c	  ARVALR(L) = ARVALR(L) + AUX9*RADAR(INDIC)

	 ELSE
	  NYL3 = (LMD3*(LMD3+1))/2 + IMU3
	  AUX9 = YIWAN(IABS(MUPW)) * YLMAC(NYL3) * GONE(KTT, KAL)
	  ARVALI(L) = ARVALI(L) + AUX8I *RADAR(INDIC)
c	  ARVALI(L) = ARVALI(L) + AUX9*RADAR(INDIC)
	 ENDIF
	ENDIF
       ENDIF

       KAL    = KAL + 1

 5    CONTINUE

      CSTX = CSTE/DSQRT(AB*AC*AD)

c      tetrar = 0.d0
c      do 25 i=0, lmdmax
c      tetrar = tetrar + cstx * arvalr(i)
c25    continue

c      tetrai = 0.d0

c      do 30 i=0, lmdmax
c      tetrai = tetrai + cstx * arvali(i)
c30    continue
c
c      write(*,*)"FF",tetrar,tetrai,cstx,gone(2,1),arvalr(lmdmax)
      CALL ACCEL(NEPSIL, ARVALR, LMDMAX, CSTX, TETRAR)
      CALL ACCEL(NEPSIL, ARVALI, LMDMAX, CSTX, TETRAI)
c
      RETURN
      END
