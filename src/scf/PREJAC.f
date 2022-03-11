      SUBROUTINE PREJAC(NATO,EIGV,EIGF,HCORP,SOVLP,HT)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION EIGF(N_ORB,N_ORB), HT((N_ORB*(N_ORB+1))/2),
     $ HCORP(N_ORB,N_ORB),SOVLP(N_ORB,N_ORB),EIGV(N_ORB)
c       
      DO 104 I=1,NATO
      DO 104 J=1,NATO
      EIGF(I,J)=0.0d0
      EIGF(I,I)=1.0d0
c 
 104  CONTINUE
      DO 101 I=1,NATO
      DO 103 J=1,NATO
 
c      DO 103 K=1,NATO
c      HCORP(I,J) = HCORP(I,J) - SOVLP(I,K)*EIGF(K,J)
      HCORP(I,J) = HCORP(I,J) - SOVLP(I,J)
 103  CONTINUE
      HCORP(I,I) = HCORP(I,I) - SOVLP(I,I)*EIGV(I)

 101  CONTINUE
      IJ=1
      DO 102 I=1,NATO
      DO 102 J=1,I
      HT(IJ)=0.0d0
      HT(IJ)=HCORP(I,J)
c      write(*,*)HT(IJ),IJ,SOVLP(I,J),I,J
 
      IJ=IJ+1
 102  CONTINUE
      RETURN
      END