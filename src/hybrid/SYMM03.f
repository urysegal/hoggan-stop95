      SUBROUTINE SYMM03(NORB1,NORB2, NLMA,NLMB, OUIZA)
      IMPLICIT REAL*8 (A-H, O-Z)
      COMPLEX*16 OUIZA

      DIMENSION NLMA(*), NLMB(*)
      DIMENSION OUIZA(20,20,20,20)

C.... COMBINATION OF THE ORBITALS

      DO 5 NAT4=1, NORB2
       LST = 5*(NAT4-1)
       M4  = NLMB(LST+5)

       IF(M4.NE.0)THEN
	IS4 = ISIGN(1, M4)
       ELSE
	IS4 = 0
       ENDIF

      DO 5 NAT3=1, NORB1
       KST = 5*(NAT3-1)
       M3  = NLMA(KST+5)

       IF(M3.NE.0)THEN
	IS3 = ISIGN(1, M3)
       ELSE
	IS3 = 0
       ENDIF

      DO 5 NAT2=1, NORB1
       JST    = 5*(NAT2-1)

      DO 5 NAT1=NAT2, NORB1
       IST  = 5*(NAT1-1)

       OUIZA(NAT2,NAT1,NAT3,NAT4) = (-1)**(M3+M4) *
     $                DCONJG(OUIZA(NAT1, NAT2, NAT3+IS3,NAT4+IS4))

 5    CONTINUE

      RETURN
      END
