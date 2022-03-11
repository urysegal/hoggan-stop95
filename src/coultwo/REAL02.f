CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE TRANSFORMS THE IMAGINARY TWO-CENTER COULOMB INTEGRALS  * C
C * INTO REAL INTEGRALS.                                                * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   NORB1, NORB2 : NUMBER OF ORBITALS ATTACHED TO THE CENTERS (1), (2)* C
C *                                                                     * C
C *   NLMA, NLMB : ONE-DIMENSIONAL ARRAYS CONTAINING THE DESRIPTION OF  * C
C *    THE SLATER ORBITALS                                              * C
C *                                                                     * C
C *   NTORB : TOTAL NUMBER OF THE ORBITALS DESCRIBING THE SYSTEM        * C
C *                                                                     * C
C *   NSUM : (NTORB * (NTORB + 1))/2                                    * C
C *                                                                     * C
C *   OUIZA : FOUR-DIMENSIONAL COMPLEX ARRAY CONTAINING THE ORIGINAL    * C
C *    IMAGINARY INTEGRALS                                              * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   FATIMA : ONE-DIMENSIONAL ARRAY, FOR REAL VALUES                   * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE REAL02(NORB1,NORB2, NLMA,NLMB,
     $                                       NTORB, NSUM, OUIZA, FATIMA)
      IMPLICIT REAL*8 (A-H, O-Z)
      COMPLEX*16 OUIZA
      LOGICAL BOOLI0,BOOLI1,BOOLJ0,BOOLJ1, BOOLK0,BOOLK1,BOOLL0,BOOLL1

      DIMENSION NLMA(*), NLMB(*), FATIMA(*), OUIZA(20,20,20,20)



      DO 5 NAT4=1, NORB2
       LST = 5*(NAT4-1)
       L   = NLMB(LST+2)
       M4  = NLMB(LST+5)

       CALL ISBOOL(M4, IS4, BOOLL0, BOOLL1)

       LY = NAT4
       LX = NAT4 + IS4
       L1 = L + IS4

      DO 5 NAT3=NAT4, NORB2
       KST = 5*(NAT3-1)
       K   = NLMB(KST+2)
       M3  = NLMB(KST+5)

       CALL ISBOOL(M3, IS3, BOOLK0, BOOLK1)

       KY = NAT3
       KX = NAT3 + IS3
       K1 = K + IS3

       CALL MINMAX(K,L, K1,L1, KP,LP, KQ,LQ, KR,LR, KS,LS)

      DO 5 NAT2=1, NORB1
       JST = 5*(NAT2-1)
       J   = NLMA(JST+2)
       M2  = NLMA(JST+5)

       CALL ISBOOL(M2, IS2, BOOLJ0, BOOLJ1)

       JY = NAT2
       JX = NAT2 + IS2
       J1 = J + IS2

      DO 5 NAT1=NAT2, NORB1
       IST = 5*(NAT1-1)
       I   = NLMA(IST+2)
       M1  = NLMA(IST+5)

       CALL ISBOOL(M1, IS1, BOOLI0, BOOLI1)

       IY = NAT1
       IX = NAT1 + IS1
       I1 = I + IS1

       CALL MINMAX(I,J, I1,J1, IP,JP, IQ,JQ, IR,JR, IS,JS)

C.... GENERAL CONDITION : IN ORDER TO BE COMBINED PROPERLY, THE MAGNETIC
C.... NUMBERS OF THE ORBITALS SHOULD BE POSITIVE INTEGERS

       IF(M1.GE.0 .AND. M2.GE.0 .AND. M3.GE.0 .AND. M4.GE.0)THEN
	CALL CMBINT(M1,M2,M3,M4, BOOLI0,BOOLI1,BOOLJ0,BOOLJ1,
     $          BOOLK0,BOOLK1,BOOLL0,BOOLL1, IP,IQ,IR,IS, JP,JQ,JR,JS,
     $             KP,KQ,KR,KS, LP,LQ,LR,LS, IY,JY,KY,LY, IX,JX,KX,LX,
     $                                       NTORB, NSUM, OUIZA, FATIMA)
       ENDIF


 5    CONTINUE

      RETURN
      END









