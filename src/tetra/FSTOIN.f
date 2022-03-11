C.... FOUR CENTER INTEGRALS
      SUBROUTINE FCINTS(N1,L1,M1, N2,L2,M2, N3,L3,M3, N4,L4,M4,
     $                  ZETA1,ZETA2,ZETA3,ZETA4,
     $                  AB,PHIAB,YLMAB, AC,PHIAC,YLMAC, AD,PHIAD,YLMAD,
     $                  ARGNT1, ARGNT2, ARGNT3, ARGNT4, ARGNT5,
     $                  RADAR, VALUE)

      IMPLICIT REAL*8 (A-H, O-Z)
c     REAL*4 GONE
      INCLUDE "../lib95/SIZE.INCL"
      INCLUDE "../lib95/CONST.INCL"
      DIMENSION YLMAB(0:*), YLMAC(0:*), YLMAD(0:*)
      DIMENSION ARGNT1(*), ARGNT2(*), ARGNT3(*), ARGNT4(*), ARGNT5(*)
      DIMENSION RADAR(*)
      DIMENSION ARGNTA(100), ARGNTB(100), ARGNTC(100), INX(10),
     $                                                    ARVALR(0:100)
      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:20)
      COMMON/XGONE/GONE((LMDMAX+1)**2 * (LMDMAX+6+1)**2, LMDMAX+3+1)

      WRITE(*, 3)N1,L1,M1, N2,L2,M2, N3,L3,M3, N4,L4,M4
 3    FORMAT(4(3(I3, 1X), 3X))

      A1 = (2.d0*ZETA1)**N1 * DSQRT((2.D0*ZETA1)/FACT(2*N1))
      A2 = (2.d0*ZETA2)**N2 * DSQRT((2.D0*ZETA2)/FACT(2*N2))
      A3 = (2.d0*ZETA3)**N3 * DSQRT((2.D0*ZETA3)/FACT(2*N3))
      A4 = (2.d0*ZETA4)**N4 * DSQRT((2.D0*ZETA4)/FACT(2*N4))

      CSTE = 16384.D0*PI7 * DFAC(L2)*DFAC(L3)*DFAC(L4) * A1*A2*A3*A4 *
     $       (-AB)**L2 * (-AC)**L3 * (-AD)**L4

      CALL RAZ0(ARVALR, 0, LMDMAX)
      CALL RAZ1(INX, 1, 10)

      DO 5 L=0, LMDMAX

       AUX0  = 1.D0/DBLE(L+L+1)

       INX(1) = INX(1)+INX(2)+INX(3)+INX(4)+INX(5) +
     $                               INX(6)+INX(7)+INX(8)+INX(9)+INX(10)

      DO 5 M=-L, L

       ABANG = (M2-M1-M) * PHIAB

       CALL RAZ1(INX, 2, 10)

       K1 = 0
       K2 = 0

      DO 5 LP2=0, L2
       INX(2) = INX(2)+INX(3)+INX(4)+INX(5) +
     $                               INX(6)+INX(7)+INX(8)+INX(9)+INX(10)

      DO 5 MP2=MAX(-LP2, -L2+LP2+M2), MIN(LP2, L2-LP2+M2)

       K1   = K1 + 1
       AUX1 = ARGNT1(K1)
       CALL RAZ1(INX, 3, 10)

      DO 5 LP12=IABS(LP2-L1), L1+LP2, 2
       K2     = K2 + 1
       AUX2   = ARGNT2(K2)
       INX(3) = INX(3)+INX(4)+INX(5)+INX(6)+INX(7)+INX(8)+INX(9)+INX(10)

       IF(AUX2 .NE. 0.D0)THEN
	CALL GAUNT(LP12,MP2-M1, L,M, LMD2M,LMD2X,MOX, NGAUNT, ARGNTA)
       ENDIF

       K3 = 0
       CALL RAZ1(INX, 4, 10)
       IAL = 1
      DO 5 LMD2=IABS(L-LP12), L+LP12, 2

       IF(LMD2 .GE. LMD2M)THEN
	K3   = K3 + 1
	NYL2 = (LMD2*(LMD2+1))/2 + IABS(M-(MP2-M1))
	MUPW = (M-(MP2-M1) - IABS(M-(MP2-M1)))/2
	AUX3 = (-1)**MUPW * ARGNTA(K3) * YLMAB(NYL2)
	II1  = L*L  *(LMDMAX+6+1)**2
	II2  = (M+L)*(LMDMAX+6+1)**2
	II3  = LP12**2
	II4  = MP2-M1+LP12+1
	IIT  = II1 + II2 + II3 +II4
	if(argnta(k3) .ne. gone(iit, ial))then
	 wRITE(*, 1)L,M, LP12,MP2-M1, ARGNTA(K3), GONE(IIT,IAL),IIT,IAL
 1       FORMAT(4(I3, 2X), 2(D16.10, 2X), 2(I4, 2X))
	endif
       ELSE
	AUX3 = 0.D0
       ENDIF

       IAL = IAL + 1

       INX(4) = INX(4)+INX(5)+INX(6)+INX(7)+INX(8)+INX(9)+INX(10)
       CALL RAZ1(INX, 5, 10)

       K4 = 0
       K6 = 0

      DO 5 LP3=0, L3
       INX(5) = INX(5)+INX(6)+INX(7)+INX(8)+INX(9)+INX(10)

      DO 5 MP3=MAX(-LP3, -L3+LP3+M3), MIN(LP3, L3-LP3+M3)
       K4   = K4 + 1
       AUX4 = (-1)**MP3 * ARGNT3(K4)

       CALL RAZ1(INX, 6, 10)
       K5 = 0

      DO 5 LP4=0, L4
       INX(6) = INX(6)+INX(7)+INX(8)+INX(9)+INX(10)

      DO 5 MP4=MAX(-LP4, -L4+LP4+M4), MIN(LP4, L4-LP4+M4)
       K5   = K5 + 1
       AUX5 = (-1)**MP4 * ARGNT4(K5)

       CALL RAZ1(INX, 7, 10)

      DO 5 LP34=IABS(LP4-LP3), LP4+LP3, 2
       K6   = K6 + 1
       AUX6 = ARGNT5(K6)

       IF(AUX6 .NE. 0.D0)THEN
	CALL GAUNT(LP34, MP3-MP4, L, M, LMDM,LMDX,MOX, NGAUNT, ARGNTB)
       ENDIF
       K7 = 0

       INX(7) = INX(7)+INX(8)+INX(9)+INX(10)
       CALL RAZ1(INX, 8, 10)

       JAL = 1
      DO 5 LMD=IABS(L-LP34), L+LP34, 2
       IF(LMD .GE. LMDM)THEN
	K7   = K7 + 1
	AUX7 = ARGNTB(K7)
	JJ1 = L*L*(LMDMAX+6+1)**2
	JJ2 = (M+L)*(LMDMAX+6+1)**2
	JJ3 = LP34**2
	JJ4 = MP3-MP4+LP34+1
	JTT = JJ1 + JJ2 + JJ3 + JJ4
	IF(ARGNTB(K7) .NE. GONE(JTT, JAL))THEN
	 WRITE(*, 2)L,M, LP34,MP3-MP4, JTT,JAL, ARGNTB(K7),GONE(JTT,JAL)
 2       FORMAT('SECOND ', 6(I5, 2X), D16.10, 2X, D16.10)
	ENDIF
       ELSE
	AUX7 = 0.D0
       ENDIF
       JAL = JAL + 1

       INX(8) = INX(8)+INX(9)+INX(10)
       CALL RAZ1(INX, 9, 10)

      DO 5 LMD4=0, LMDMAX
       INX(9) = INX(9)+INX(10)
      DO 5 MU4=-LMD4, LMD4

       ADANG = (M4 - MP4 + MU4) * PHIAD

       IF(AUX7 .NE. 0.D0)THEN
	NYL4 = (LMD4*(LMD4+1))/2 + IABS(MU4)
	MUPW = (MU4 - IABS(MU4))/2
	AUX8 = (-1)**MUPW * YLMAD(NYL4)
	CALL GAUNT(LMD,M-(MP3-MP4), LMD4,MU4, LMD3M,LMD3X,MOX, NGAUNT,
     $                                                           ARGNTC)
       ENDIF

       CALL RAZ1(INX, 10, 10)
       K9 = 0

       KAL = 1

      DO 5 LMD3=IABS(LMD-LMD4), LMD+LMD4, 2
       INX(10) = INX(10) + 1
       INDIC   = INX(1)+INX(2)+INX(3)+INX(4)+INX(5) +
     $                               INX(6)+INX(7)+INX(8)+INX(9)+INX(10)

       ACANG = (M - MU4 + MP4 - M3) * PHIAC

       IF(LMD3 .GE. LMD3M)THEN
	K9    = K9 + 1
	MU3   = MU4 - M + MP3 - MP4
	NYL3  = (LMD3*(LMD3+1))/2 + IABS(MU3)
	MUPW  = (MU3 - IABS(MU3))/2
	AUX9  = (-1)**MUPW * YLMAC(NYL3) * ARGNTC(K9)
	ANGLE  = ABANG + ACANG + ADANG
	COSINE = DCOS(ANGLE)
	SINE   = DSIN(ANGLE)

	KK1 = LMD4*LMD4 *(LMDMAX+6+1)**2
	KK2 = (MU4+LMD4)*(LMDMAX+6+1)**2
	KK3 = LMD*LMD
	KK4 = M-(MP3-MP4) + LMD + 1
	KTT = KK1 + KK2 + KK3 + KK4
	IF(ARGNTC(K9) .NE. GONE(KTT, KAL))THEN
	 WRITE(*, 4)LMD4,MU4, LMD,M-(MP3-MP4), KTT,JAL,
     $                                          ARGNTC(K9),GONE(KTT,KAL)
 4       FORMAT('THIRD ', 6(I5, 2X), D16.10, 2X, D16.10)
	ENDIF

       ELSE
	AUX9 = 0.D0
       ENDIF

       KAL = KAL + 1

       ARVALR(L) = ARVALR(L) +
     $               AUX0*AUX1*AUX2*AUX3*AUX4*AUX5*AUX6*AUX7*AUX8*AUX9 *
     $               RADAR(INDIC) * COSINE

 5    CONTINUE

      SUM = 0.D0
      DO 25 I=0, LMDMAX
       SUM = SUM + CSTE/DSQRT(AB*AC*AD) * ARVALR(I)
       WRITE(*, *)I, SUM
 25   CONTINUE

      RETURN
      END
