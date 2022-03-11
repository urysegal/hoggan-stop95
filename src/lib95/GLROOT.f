      SUBROUTINE GLROOT(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
c
c      write(*,*)"in glroot", lmherm,isp,N, N+isp
      DO 5 I=1, N
       RLEAG0(I+ISP) =0.0d0
       RLEAG0(I+ISP) = T(I)
       Y(I)          = 0.D0
 5    CONTINUE
c
      RETURN
      END





