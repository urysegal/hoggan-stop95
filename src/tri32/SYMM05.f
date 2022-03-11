      SUBROUTINE SYMM05(NORB1,NORB2,NORB3, NLMA,NLMB,NLMC, OUIZA)

      IMPLICIT REAL*8 (A-H, O-Z)
      COMPLEX*16 OUIZA

      DIMENSION NLMA(*), NLMB(*), NLMC(*)
      DIMENSION OUIZA(20,20,20,20)

C.... DEFINITION OF SOME CONSTANTS

      NSUM = NORB1 * (NORB1 + 1)
      NSX  = NORB2 * NSUM
      NSY  = 2 * (NORB1 + 1)

C.... INITIALISATION OF THE RUNNING INDEX USED FOR THE TEMPORARY STORAGE

      IRUN = 1

      DO 5 NAT4=1, NORB3
       LST = 5*(NAT4-1)
       L   = NLMC(LST+2)
       M4  = NLMC(LST+5)

       IF(M4.NE.0)THEN
	IS4 = ISIGN(1, M4)
       ELSE
	IS4 = 0
       ENDIF

       LX   = NAT4 + IS4
       NSU4 = (LX-1) * NSX


      DO 5 NAT3=1, NORB2
       KST = 5*(NAT3-1)
       K   = NLMB(KST+2)
       M3  = NLMB(KST+5)

       IF(M3.NE.0)THEN
	IS3 = ISIGN(1, M3)
       ELSE
	IS3 = 0
       ENDIF

       KX    = NAT3 + IS3
       NSU3  = (KX-1) * NSUM
       NSU43 = NSU4 + NSU3

       NSU2  = -NSY

      DO 5 NAT2=1, NORB1
       JST    = 5*(NAT2-1)
       J      = NLMA(JST+2)

       NSU2   = NSU2 + (NSY + 2 - 2*NAT2)
       NSU432 = NSU43 + NSU2

       NSU1 = 0

      DO 5 NAT1=NAT2, NORB1
       IST  = 5*(NAT1-1)
       I    = NLMA(IST+2)

       NSU1 = NSU1 + 2

       OUIZA(NAT2,NAT1,NAT3,NAT4) = (-1)**(M3+M4) *
     $                DCONJG(OUIZA(NAT1, NAT2, NAT3+IS3,NAT4+IS4))

 5    CONTINUE

      RETURN
      END
