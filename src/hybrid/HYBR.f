CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            NOV 16 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE COMPUTES THE TWO-CENTER HYBRID INTEGRALS USED IN    * C
C * THE SCF PROCEDURE. IT IS CALLED FROM THE ROUTINE TOELC.             * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   NEPSIL : ORDER OF THE EPSILON ALGORITHM USED TO ACCELERATION      * C
C *                                                                     * C
C *   NTORB : TOTAL NUMBER OF ATOMIC ORBITALS INVOLVED IN THE SYSTEM    * C
C *                                                                     * C
C *   NSUM : (NTORB*(NTORB+1))/2                                        * C
C *                                                                     * C
C *   LM123 : IS THE HIGHEST VALUE OF THE QUANTUM NUMBER L              * C
C *                                                                     * C
C *   IC1, IC2 : ARE THE ATOMS THE ORBITALS OF WHICH ARE INTERACTING    * C
C *                                                                     * C
C *   NORBAT : IS AN ARRAY CONTAINING THE NUMBER OF ORBITALS ATTACHED TO* C
C *    EACH ATOM                                                        * C
C *                                                                     * C
C *   NLMAT : IS AN ARRAY CONTAINING THE FOLLOWING INFORMATIONS :       * C
C *    NLMAT(1) -> ATOMIC NUMBER                                        * C
C *    NLMAT(2) -> ORDER OF THE ORBITAL AS OCCURING IN THE DATA         * C
C *    NLMAT(3), NLMAT(4), NLMAT(5) -> N, L, M                          * C
C *                                                                     * C
C *   ZETAAT : IS AN ARRAY CONTAINING THE SLATER ZETA EXPONENTS         * C
C *                                                                     * C
C *   XYZAT : IS AN ARRAY CONTAINING THE CARTESIAN COORDINATES OF THE   * C
C *    ATOMS                                                            * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   FATIMA : IS A ONE-DIMENSIONAL ARRAY CONTAINING THE TWO-ELECTRONS  * C
C *    INTEGRALS <I J | K L>                                            * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE HYBR(NEPSIL, NTORB,NSUM, LM123, IC1,IC2, NORBAT,
     $                                      NLMAT, ZETAAT,XYZAT, FATIMA)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER ZETST1, ZETST2
      DIMENSION NORBAT(*), NLMAT(*), ZETAAT(*), XYZAT(*), FATIMA(*)
      DIMENSION XYZV(3)


ccc   WRITE(*, 1)IC1, IC1, IC1, IC2
c1    FORMAT('HYBR : ', 4(I2, 2X))

C.... COMPUTATION OF THE COORDINATE OF THE CENTER LABELED IC2 TO WHICH THE
C.... ORBITAL X4(2) IS ATTACHED

      NIC1 = 3*(IC1-1)+1
      NIC2 = 3*(IC2-1)+1

      CALL CARTCO(XYZAT(NIC1),XYZAT(NIC2), XYZV, RVMOD,THETRV,PHIRV)


C.... COMPUTATION OF THE STARTING INDICES FOR THE ORBITALS AND THE SLATER
C.... EXPONENTS

      CALL NDXSTR(NORBAT, IC1, NLMST1, ZETST1)
      CALL NDXSTR(NORBAT, IC2, NLMST2, ZETST2)

C.... COMPUTATION OF THE TWO-CENTER HYBRID INTEGRALS

      CALL HYBRID(NEPSIL, NTORB,NSUM, LM123, NORBAT(IC1),NORBAT(IC2),
     $     NLMAT(NLMST1),NLMAT(NLMST2), ZETAAT(ZETST1),ZETAAT(ZETST2),
     $                                       RVMOD,THETRV,PHIRV, FATIMA)

      RETURN
      END
