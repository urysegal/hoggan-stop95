
      SUBROUTINE JACOBI(NDIM,H,C,EV,N,NN,DH,JT,NS)
CR8   1
      IMPLICIT REAL*8    (A-H,O-Z)
C*****************************************************
C     JACOBI DIAGONALISIERUNG
C     H=UNTERES LINKES DREIECK DER MATRIX,DIE DIAGONALISIERT WERDEN SOLL
C      H(1,1),H(2,1),H(2,2),H(3,1),...H(N,N)
C     C=(INPUT)UNITAERE MATRIX//(OUTPUT)EIGENVEKTORMATRIX(NN,N)
C     DH=SCHWELLE, ALLE H(I,J)'S .GE. TH*NORM(H) WERDEN WEGROTIERT
C      WENN DH=0 : DEFAULT VALUE 1.0D-10
C     JT=MAXIMALZAHL DER ITERATIONEN; DEFAULT VALUE : N*(N*(N+1)/2+1)
C     EV=EIGENWERTE, EV(1),EV(2),...EV(N)
C     NS=ZAHL DER EIGENWERTE UND -VEKTOREN, DIE GEORDNET AM ANFANG
C     STEHEN SOLLEN:   NS.GT.0 : GROESSTER EV AM ANFANG
C     NS.LT.0 : KLEINSTER EV AM ANFANG // NS.EQ.0 : KEINE ORDNUNG
c      DIMENSION H(2  ),C(NDIM  ,2 ),EV(2   )
c      DIMENSION H(2  ),C(NDIM  ,2 ),EV(2   )
      DIMENSION H(NDIM),C(NDIM  ,NDIM ),EV(NDIM )
      DATA SMALL/1.D-30/,SMALL1/1.D-10/,ZERO/0.D0/
CR8   2
C      SQRT(XXX)=DSQRT(XXX)
C      ABS(XXX)=DABS(XXX)
C     FIND SUM OF SQUARES
cc      write(*,*) 'NDIM,N,NN,JT,NS', NDIM,N,NN,JT,NS
c
      RTWO=SQRT(2.D0)/2.D0
      IT=JT
      SUMSQ = SMALL
      IJ=1
      DO 10 I=1,N
      JM=I-1
      IF(JM.EQ.0) GO TO 20
      DO 30 J=1,JM
      SUMSQ = SUMSQ + H(IJ)*H(IJ)
      IJ = IJ + 1
30    CONTINUE

   20 EV(I) = H(IJ)
      H(IJ) = 0.D0
      IJ = IJ + 1
10    CONTINUE
      IF(N.LE.1) GO TO 200
      IF(IT.EQ.0) IT=IJ*N
      ITO = IT
      AN = 1.9D0/DBLE(N*(N-1))
      Q = 0.5D0
   40 THF =   SQRT(AN*SUMSQ)
      TH = DH*THF
      IF(TH.LE.0.0D0) TH=THF*SMALL1
C     FIND PIVOT
   50 IJ = 1
      IC = 0
      DO 1761 I=1,N
      DO 1760 J=1,I
      IF(THF.LE.  ABS(H(IJ))) GO TO 70
      GO TO 60
C     FIND SIN AND COS
   70 A = EV(J) - EV(I)
      SUMSQ = SUMSQ - H(IJ)**2
      IC = 1
      IF(  ABS(A).GT.(TH*  ABS(H(IJ)))) GO TO 80
      XCOS = RTWO
      XSIN = RTWO
  100 COSS = XCOS*XCOS
      SICO = XCOS*XSIN
      GO TO 110
   80 TL = 2.0D0*H(IJ)/A
C     COMPUTE SIN AND COS
      XSIN = 0.5D0/  SQRT(TL*TL + 1.0D0)
      COSS = XSIN + 0.5D0
      XCOS =   SQRT(COSS)
      SICO = XSIN*TL
      XSIN = SICO/XCOS
  110 SINS = XSIN*XSIN
      B = EV(J)
      EV(J) = COSS*B + 2.0D0*SICO*H(IJ) + SINS*EV(I)
      EV(I) = COSS*EV(I) - 2.0D0*SICO*H(IJ) + SINS*B
      H(IJ) =(COSS-SINS)*H(IJ) - SICO*A
C     TRANSFORM H(IJ) AND C(IJ)
      JL = J*(J-1)/2
      IL = I*(I-1)/2
      DO 120 L=1,N
      JL = JL + 1
      IL = IL + 1
      IF(J-L) 130,120,150
  130 JL = JL + L - 2
      IF(I-L) 140,120,150
  140 IL = IL + L - 2
  150 A = H(IL)
      H(IL) = XCOS*A - XSIN*H(JL)
      H(JL) = XSIN*A + XCOS*H(JL)
  120 CONTINUE
      DO 180 L=1,NN
      A = C(L,I)
      C(L,I) = XCOS*A - XSIN*C(L,J)
      C(L,J) = XCOS*C(L,J) + XSIN*A
180   CONTINUE
      IT = IT - 1
      IF (IT.EQ.0) GO TO 200
   60 IJ = IJ + 1
1760  CONTINUE
1761  CONTINUE
c 
      IF (IC.GT.0) GO TO 65
c      write(*,*)"jacobi",DH,THF,TH
      IF (THF.LE.TH) GO TO 200
c      GO TO 200
      SUMSQ = 0.D0
   65 SUMSQ = MAX(SUMSQ,ZERO)
      AN = AN*Q
      THF = MAX( SQRT(SUMSQ*AN),TH)
      GO TO 50
  200 IF(NS.EQ.0) GO TO 310
      NA = IABS(NS)
      IF(NS.GT.N) NA = NS-N
      A = DBLE(  NS/NA)
      DO 300 I=1,NA
      J=I
      XCOS = EV(I)
      DO 250 K=I,N
      COSS = XCOS**2
      XSIN =(EV(K))**2
      IF((NS.GT.N).AND.((COSS-XSIN).GE.0.0D0)) GO TO 250
      IF((NS.LE.N).AND.(A*(XCOS-EV(K)).GE.0.0D0)) GO TO 250
      XCOS = EV(K)
      J = K
  250 CONTINUE
      IF(I.EQ.J) GO TO 300
      EV(J) = EV(I)
      EV(I) = XCOS
      DO 280 K=1,NN
      XSIN = C(K,I)
      C(K,I) = C(K,J)
      C(K,J) = XSIN
  280 CONTINUE
  300 CONTINUE
  310 ITA = ITO - IT
 
      RETURN
      END
