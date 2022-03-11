CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUROUTINE RETURNS THE TWO-CENTER OVERLAP INTEGRALS NEEDED IN   * C
C * IN THE SCF ROUTINE. IT IS CALLED FROM THE ONELEC ROUTINE.           * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LM123 : IS THE HIGHEST VALUE OF THE QUANTUM NUMBER L              * C
C *                                                                     * C
C *   IC1, IC2 : ARE THE ATOMS THE ORBITALS OF WHICH ARE OVERLAPING     * C
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
C *   SOVL : IS A TWO-DIMENSIONAL ARRAY CONTAINING THE REAL OVELAPS     * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE OVERLP(LM123,IC1,IC2,NORBAT,NLMAT,ZETAAT,XYZAT,XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      INTEGER ZETST1, ZETST2
      DIMENSION XCORE(N_ORB, N_ORB)
      DIMENSION NORBAT(*), NLMAT(*), ZETAAT(*), XYZAT(*)
      DIMENSION XYZV(3)

C     WRITE(*, 1)IC1, IC2
C1    FORMAT('OVERLP : ', 2(I2, 2X))

C.... COMPUTATION OF THE COORDINATES OF THE CENTER LABELED IC2 TO WHICH THE
C.... ORBITAL X2(1) IS ATTACHED

      NIC1 = 3*(IC1-1)+1
      NIC2 = 3*(IC2-1)+1

      CALL CARTCO(XYZAT(NIC1),XYZAT(NIC2), XYZV, RV,THETRV,PHIRV)

C.... DETERMINATION OF THE STARTING INDICES OF THE ORBITALS AND THE SLATER
C.... EXPONENTS

      CALL NDXSTR(NORBAT, IC1, NLMST1, ZETST1)
      CALL NDXSTR(NORBAT, IC2, NLMST2, ZETST2)

      CALL OVELP1(LM123, NORBAT(IC1),NORBAT(IC2),
     $  NLMAT(NLMST1),NLMAT(NLMST2), ZETAAT(ZETST1),ZETAAT(ZETST2),
     $                                    RV,THETRV,PHIRV, XCORE)

      RETURN
      END
