      SUBROUTINE BASIS(NBASE, NATOM, NUMAT, NORBAT, NLMAT, ZETAAT,
     $                                                    ZETAMY, LM123)
       IMPLICIT REAL*8 (A-H, O-Z)
       DIMENSION NUMAT(*), NORBAT(*), NLMAT(*), ZETAAT(*), ZETAMY(*)


      IF(NBASE .NE. 0)THEN
       DO 5 I=1, NATOM
	CALL SBASE(I, NUMAT(I), NLMAT, ZETAAT, ZETAMY, LM123)
 5     CONTINUE
      ELSE

       DO 15 I=1, N_ORB
	READ(*, *, END=110)NAT, N, L, M, ZETA

	NORBAT(NAT)   = NORBAT(NAT) + 1
	NDEX          = 5*(I-1)
	NLMAT(NDEX+1) = NUMAT(NAT)
	NLMAT(NDEX+2) = I
	NLMAT(NDEX+3) = N
	NLMAT(NDEX+4) = L
	NLMAT(NDEX+5) = M
	ZETAAT(I)     = ZETA
	ZETAMY(NAT)   = ZETAMY(NAT) + ZETA
	LM123         = MAX(LM123, L)
 15    CONTINUE
      ENDIF

 110  RETURN
      END




      SUBROUTINE SBASE(I, NUMAT, NLMAT, ZETAAT, ZETAMY, LM123)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER CELC, RELC
      DIMENSION NLMAT(*), ZETAAT(*), ZETAMY(*)
      DIMENSION LEVEL(36), NL(2, 100), NL1(2, 100)
      DATA LEVEL/1,0,2,
     $           2,0,2, 2,1,6,
     $           3,0,2, 3,1,6,
     $           4,0,2, 3,2,10,
     $           4,1,6, 5,0,2,
     $           4,2,10, 5,1,6, 6,0,2/

      RELC = NUMAT
      CELC = 0
      LVL  = 0

      DO 5 WHILE(RELC .GT. 0)
       N = LEVEL(3*LVL + 1)
       L = LEVEL(3*LVL + 2)

       LM123 = MAX0(L, LM123)

       NL(1, LVL+1) = N
       NL(2, LVL+1) = L

       CELC = CELC + LEVEL(3*LVL + 3)

       RELC = RELC - LEVEL(3*LVL + 3)
       LVL  = LVL + 1

 5    CONTINUE

      IF(NUMAT .GT. CELC)then
       N = LEVEL(3*LVL + 1)
       L = LEVEL(3*LVL + 2)

       LM123 = MAX0(L, LM123)

       NL(1, LVL+1) = N
       NL(2, LVL+1) = L

      ENDIF

c.... ORDONNER SELON LA SYMMETRIE

      NL1(1, 1) = NL(1, 1)
      NL1(2, 1) = NL(2, 1)

      NL1(1, 2) = NL(1, 2)
      NL1(2, 2) = NL(2, 2)

      K = 3

      DO 10 L=0, LM123
      DO 10 LV=3, LVL

       LL = NL(2, LV)

       IF(LL .EQ. L)THEN
	NL1(1, K) = NL(1, LV)
	NL1(2, K) = NL(2, LV)
	K = K + 1
       ENDIF

 10   CONTINUE

      DO 15 K=1, LVL
       NDEX = 5*(K-1)
       NLMAT(NDEX + 1) = NUMAT
       NLMAT(NDEX + 2) = K
       NLMAT(NDEX + 3) = NL1(1, K)
       NLMAT(NDEX + 4) = NL1(2, K)
       NLMAT(NDEX + 5) = 0

c       ZETAAT(

 15   CONTINUE

      ORBIT = LVL

      DO 20 L=1, LM123
      DO 20 K=1, LVL
       LL = NL1(2, K)
       IF(LL .GE. L)THEN
	WRITE(*, *)NL1(1, K), NL1(2, K), L
	WRITE(*, *)NL1(1, K), NL1(2, K), -L
	ORBIT = ORBIT + 2
       ENDIF
 20   CONTINUE

      RETURN
      END


