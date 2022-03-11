CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C* THIS ROUTINE COMPUTES A THREE-CENTER TWO-ELECTRON COULOMB         *C
C* INTEGRAL, GIVEN THE RADIAL INTEGRALS STORED IN THE 'RADAR' ARRAY  *C
C* AND THE COEFFICIENTS INVOLVED IN THE SOLID SPHERICAL HARMONICS    *C
C* ADDITION AND MULTIPLICATION THEOREMS.                             *C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TCTEI(NEPSIL, N1,L1,M1, N2,L2,M2, N3,L3,M3, N4,L4,M4,
     $       ZETA1,ZETA2,ZETA3,ZETA4, AB,YLMAB,PHIAB, AC,YLMAC,PHIAC,
     $                   ARGNT1, ARGNT2, ARGNT3, RADAR, TRI32R,TRI32I)

      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      INCLUDE "../lib95/CONST.INCL"

      DIMENSION YLMAB(0:*),YLMAC(0:*), ARGNT1(*),ARGNT2(*),ARGNT3(*),
     $          ARGNT4(L1MAX+1), ARGNT5(L1MAX+L2MAX+L3MAX+L4MAX+1),
     $          ARGNT6(LDEV+1), ARVALR(0:LDEV), ARVALI(0:LDEV),
     $          RADAR(*), INX(7)

      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:100)
c      COMMON/NORMC/A12,A34

C.... THE NORMALIZATION CONSTANTS

c      A12 = ANP(ZETA1,ZETA2,N1,N2)
c      A34 = ANP(ZETA3,ZETA4,N3,N4)
       A1 = (2.d0*ZETA1)**N1 * DSQRT((2.D0*ZETA1)/FACT(2*N1))
       A2 = (2.d0*ZETA2)**N2 * DSQRT((2.D0*ZETA2)/FACT(2*N2))
       A3 = (2.d0*ZETA3)**N3 * DSQRT((2.D0*ZETA3)/FACT(2*N3))
       A4 = (2.d0*ZETA4)**N4 * DSQRT((2.D0*ZETA4)/FACT(2*N4))
      AP1 = A12*A34

c       CSTE =  DFAC(L3)*DFAC(L4) * AP1*
       CSTE =  DFAC(L3)*DFAC(L4) * A1*A2*A3*A4*
     $                                     (-AB)**L3 * (-AC)**L4

C.... INITIALIZATION OF THE RESULTING AND THE INDEXING ARRAY

      CALL RAZ0(ARVALR, 0, LMDMAX)
      CALL RAZ0(ARVALI, 0, LMDMAX)

      CALL RAZ1(INX, 1, 7)

      CALL GAUNT(L2,M2, L1,M1, LMIN1,LMAX1, MG1, NGAUNT, ARGNT4)
      K1 = 1

      DO 5 L=IABS(L1-L2), L1+L2, 2
       INX(1) = INX(1)+INX(2)+INX(3)+INX(4)+INX(5)+INX(6)+INX(7)
       CALL RAZ1(INX, 2, 7)
       IF(L .GE. LMIN1)THEN
	AUX1 = ARGNT4(K1)/DBLE(L+L+1)
	K1   = K1 + 1
       ELSE
	AUX1 = 0.D0
       ENDIF

       K2 = 1
       K3 = 1

      DO 5 LP3=0, L3
       INX(2) = INX(2)+INX(3)+INX(4)+INX(5)+INX(6)+INX(7)
      DO 5 MP3=MAX(-LP3, -L3+LP3+M3), MIN(LP3, L3-LP3+M3)

       CALL RAZ1(INX, 3, 7)
       AUX2 = ARGNT1(K3)
       K3   = K3 + 1
       K4   = 1

      DO 5 LP4=0, L4
       INX(3) = INX(3)+INX(4)+INX(5)+INX(6)+INX(7)
      DO 5 MP4=MAX(-LP4, -L4+LP4+M4), MIN(LP4, L4-LP4+M4)

       CALL RAZ1(INX, 4, 7)
       AUX3 =  ARGNT2(K4)
       K4   = K4 + 1

      DO 5 LP34=IABS(LP3-LP4), LP3+LP4, 2
       INX(4) = INX(4)+INX(5)+INX(6)+INX(7)
       CALL RAZ1(INX, 5, 7)
       AUX4   = ARGNT3(K2)
       K2     = K2 + 1

       MP34 = MP4-MP3
       IF((AUX1 .NE. 0.D0).AND.(LP34 .GE. IABS(MP34)))THEN
	CALL GAUNT(L,MG1,LP34,MP34,LMIN2,LMAX2,MG2,NGAUNT,ARGNT5)
       ENDIF

       K5 = 1

      DO 5 LP=IABS(L-LP34), L+LP34, 2
       INX(5) = INX(5) + INX(6) + INX(7)
       CALL RAZ1(INX, 6, 7)

       IF(LP .GE. LMIN2)THEN
	AUX5 = ARGNT5(K5)
	K5   = K5 + 1
       ELSE
	AUX5 = 0.D0
       ENDIF

      DO 5 LMD1=0, LMDMAX
       INX(6) = INX(6) + INX(7)
      DO 5 MU1=-LMD1, LMD1
       CALL RAZ1(INX, 7, 7)
       IF(AUX5 .NE. 0.D0)THEN
	CALL GAUNT(LP,MG2,LMD1,MU1,LMIN3,LMAX3,MG3,NGAUNT,ARGNT6)
       ENDIF

       K6 = 1

      DO 5 LMD2=IABS(LP-LMD1), LP+LMD1, 2
       INX(7) = INX(7) + 1
       INDIC  = INX(1)+INX(2)+INX(3)+INX(4)+INX(5)+INX(6)+INX(7)
       IF(LMD2 .GE. LMIN3)THEN
	ANGLE  = (MU1-(M3-MP3))*PHIAB + ((M4-MP4)-MG3)*PHIAC
	COSINE = DCOS(ANGLE)
	SINE   = DSIN(ANGLE)
	NYL1   = (LMD1*(LMD1+1))/2 + IABS(MU1)
	NYL2   = (LMD2*(LMD2+1))/2 + IABS(MG3)
	MUPW   = (MU1 - IABS(MU1))/2 + (MG3 - IABS(MG3))/2
	RAD67  = (-1)**MUPW*ARGNT6(K6)*YLMAB(NYL1)*YLMAC(NYL2)
	AUX6   = RAD67 * COSINE
	AUX7   = RAD67 * SINE
	K6     = K6 + 1
       ELSE
	AUX6 = 0.D0
	AUX7 = 0.D0
       ENDIF

       ARVALR(LMD1) = ARVALR(LMD1) + AUX1 * AUX2 * AUX3 * AUX4 * AUX5 *
     $                                            AUX6 * RADAR(INDIC)

       ARVALI(LMD1) = ARVALI(LMD1) + AUX1 * AUX2 * AUX3 * AUX4 * AUX5 *
     $                                            AUX7 * RADAR(INDIC)
 5    CONTINUE
c
c wat      C32 = C3/17.8d0
c hcn      C32 = C3/3.0d0
c      C32 = C3/1.0d0
      CALL ACCEL(NEPSIL,ARVALR,LMDMAX,C3*CSTE/DSQRT(AB*AC),TRI32R)
c
      CALL ACCEL(NEPSIL, ARVALI, LMDMAX,C3*CSTE/DSQRT(AB*AC), TRI32I)
c
c      TRI32R = DABS(TRI32R)
c      TRI32I = DABS(TRI32R)
      RETURN
      END

      FUNCTION ANP1(ZNA,ZNB,NA,NB)
c
c     calculates the normalisation constant
c

      IMPLICIT REAL*8 (A-H, O-Z)
      COMMON/FACT0/FACT(0:30)
c
c    calculates norm product
c
      IF(NA.EQ.NB)THEN
		  ZP = ZNA*ZNB

		  IF(NA.EQ.1)THEN

		  ANP = 4.0d0*ZP*DSQRT(ZP)

		  ELSE
		  ANP = 2.0d0*(4.0d0*ZP)**NA*DSQRT(ZP)/FACT(2*NA)
		  ENDIF
      ELSE
      ZNA2 = 2.0d0*ZNA
      ZNB2 = 2.0d0*ZNB
      ZP = ZNA*ZNB
      ANP = 2.0d0*ZNA2**NA*DSQRT(ZP/(FACT(2*NA)*FACT(2*NB)))*ZNB2**NB
      ENDIF
      RETURN
      END
      FUNCTION ANP(ZNA,ZNB,NA,NB)
c
c     calculates the normalisation constant
c

      IMPLICIT REAL*8 (A-H, O-Z)
      COMMON/FACT0/FACT(0:30)
c
c    calculates norm product
c

      IF(NA.EQ.NB)THEN
		  ZP = ZNA*ZNB

		  IF(NA.EQ.1)THEN

		  ANP = 4.0d0*ZP*DSQRT(ZP)

		  ELSE
c                  ANP = 2.0d0*(4.0d0*ZP)**NA*DSQRT(ZP)/FACT(2*NA)
		  ANP = 2.0d0*4**NA*ZP**NA*DSQRT(ZP)/FACT(2*NA)
		  ENDIF
      ELSE
      ZNA2 = 2.0d0*ZNA
      ZNB2 = 2.0d0*ZNB
		      IF(NA.EQ.1)ANP = ZNA2*ZNB2**NB*DSQRT(ZNB2/FACT(2*NB))

c
		      IF(NB.EQ.1)ANP = ZNB2*ZNA2**NA*DSQRT(ZNA2/FACT(2*NA))
c
      NAB = NA+NB
      ZP = ZNA*ZNB
c
      ANP=2.0d0*2**NAB*ZNA**NA*DSQRT(ZP/(FACT(2*NA)*FACT(2*NB)))*ZNB**NB
c      ANP=2.0d0*ZNA2**NA*DSQRT(ZP/(FACT(2*NA)*FACT(2*NB)))*ZNB2**NB
      ENDIF
      RETURN
      END
      FUNCTION ANP2(ZETA1,ZETA2,N1,N2)
c
c     calculates the normalisation constant
c

      IMPLICIT REAL*8 (A-H, O-Z)
      COMMON/FACT0/FACT(0:30)
c
c    calculates norm product
c

      ZNA=ZETA1
      ZNB=ZETA2
      NA=N1
      NB=N2
      IF(NA.EQ.NB)THEN
		  ZP = ZNA*ZNB

		  IF(NA.EQ.1)THEN

		  ANP = 4.0d0*ZP*DSQRT(ZP)

		  ELSE
		  ANP = 2.0d0*4**NA*ZP**NA*DSQRT(ZP)/FACT(2*NA)
		  ENDIF
c      ENDIF
      ELSE
      ZNA2 = 2.0d0*ZNA
      ZNB2 = 2.0d0*ZNB
c                      IF(NA.EQ.1)ANP = ZNA2*ZNB2**NB*DSQRT(ZNB2/FACT(2*NB))

c
c                      IF(NB.EQ.1)ANP = ZNB2*ZNA2**NA*DSQRT(ZNA2/FACT(2*NA))
c
      NAB = NA+NB
      ZP = ZNA*ZNB
c
c      ANP=2.0d0*2**NAB*ZNA**NA*DSQRT(ZP/(FACT(2*NA)*FACT(2*NB)))*ZNB**NB
       A1 = (2.d0*ZETA1)**N1 * DSQRT((2.D0*ZETA1)/FACT(2*N1))
       A2 = (2.d0*ZETA2)**N2 * DSQRT((2.D0*ZETA2)/FACT(2*N2))
c       A3 = (2.d0*ZETA3)**N3 * DSQRT((2.D0*ZETA3)/FACT(2*N3))
c       A4 = (2.d0*ZETA4)**N4 * DSQRT((2.D0*ZETA4)/FACT(2*N4))
      ANP = A1*A2
      ENDIF
      RETURN
      END




