CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 23 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS THE ONE-CENTER ONE-ELECTRON INTEGRALS NEEDED   * C
C * IN THE SCF ROUTINE. IT IS CALLED FROM THE ONELEC ROUTINE.           * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   IC1, IC2, IOP : ARE THE RANKS OF THE ATOMS THE ORBITALS OF WHICH  * C
C *    ARE INTERACTING OF THE NUCLEUS IOP                               * C
C *                                                                     * C
C *   NUMAT : IS AN ARRAY CONTAINING THE ATOMIC NUMBERS                 * C
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
C * OUTPUT :                                                            * C
C *   XCORE : IS A TWO-DIMENSIONAL ARRAY CONTAINING THE CORE INTEGRALS  * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE OCOEIS(IC1,IOP,IC2,NUMAT,NORBAT,NLMAT,ZETAAT,XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER ZETSTR
      INCLUDE "../lib95/SIZE.INCL"

      DIMENSION NUMAT(*), NORBAT(*), NLMAT(*), ZETAAT(*)
      DIMENSION XCORE(N_ORB, N_ORB)

ccc   WRITE(*, 1)IC1, IOP, IC2
c1    FORMAT('OCOEIS : ', 3(I2, 2X))

C.... DETERMINING THE STARTING INDICES FOR THE ORBITALS AND THE SLATER
C.... EXPONENTS

      CALL NDXSTR(NORBAT, IC1, NLMSTR, ZETSTR)

      CALL OCOEI(NUMAT(IOP),NORBAT(IC1),NLMAT(NLMSTR),ZETAAT(ZETSTR),
     $                                                            XCORE)
      RETURN
      END
