      SUBROUTINE DUMMY0(AB, LEGP)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "SIZE.INCL"
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      EXTERNAL GLROOT

      ISP  = 0
      XXXX = DGLNQ(GLROOT, 0.D0, AB, LEGP)

      RETURN
      END


      SUBROUTINE DUMMY1(RMIN, RMAX, LEGP, LEGQ)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "SIZE.INCL"
      COMMON/GLAE0/RLEAG0(LMHERM), ISP
      EXTERNAL GLROOT

C.... COMPUTATION AND STORAGE OF THE ROOTS

      ISP  = 0
      XXXX = DGLNQ(GLROOT, 0.D0, RMIN, LEGP)
      ISP  = ISP + LEGP
      XXXX = DGLNQ(GLROOT, RMIN, RMAX, LEGQ)

      RETURN
      END

