      SUBROUTINE FUSED0(DFAC, L1,M1, L2,M2, AB,YLMAB, ARGNT1, ARGNT2)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION DFAC(0:*), YLMAB(0:*), ARGNT1(*), ARGNT2(*), ARGNT(100)

c      lp = l1+l2
      K1   = 1
      K3   = 1

      if(m1.eq.m2)go to 4
c      if(lp.eq.1)then
c      if(iabs(m1-m2).eq.1)then
      if(m1.eq.0.or.m2.eq.0)then
c
      do 3 k11 = 0,l1
      do 3 k12 = 0,l2
      argnt1(k1) = 0.0d0
      argnt2(k3) = 0.0d0
      write(*,*)"argnt", argnt1(k1), argnt2(k2)
      k1 = k1 + 1
      k3 = k3 + 1
  3   continue
c      if(lp.eq.1)go to 10
      go to 10
c      else
c      go to 4
      endif
  4   continue
c
c
      ABLP = -AB

      DO 5 LP2=0, L2
       ABLP = -ABLP/AB
      DO 5 MP2=MAX(-LP2, -L2+LP2+M2), MIN(LP2, L2-LP2+M2)
       CALL GAUNT(LP2,MP2, L2,M2, LMIN,LMAX,M12, NGAUNT, ARGNT)
       NYL        = ((L2-LP2)*(L2-LP2+1))/2 + IABS(M2-MP2)
       MUPW       = (M2-MP2 - IABS(M2-MP2))/2
       ARGNT1(K1) = ABLP*YLMAB(NYL)*ARGNT(1)/(DFAC(LP2)*DFAC(L2-LP2)) *
     $              (-1)**MUPW
       K1         = K1 + 1
       CALL GAUNT(L1,M1, LP2,MP2, LMIN,LMAX,M12, NGAUNT, ARGNT)
       K2 = 1
      DO 5 LP12=IABS(L1-LP2), L1+LP2, 2
       IF(LP12 .GE. LMIN)THEN
	ARGNT2(K3) = ARGNT(K2)
	K2         = K2 + 1
       ELSE
	ARGNT2(K3) = 0.D0
       ENDIF
       K3 = K3 + 1
 5    CONTINUE
 10   continue
      RETURN
      END


