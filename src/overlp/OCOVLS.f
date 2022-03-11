CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 20 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUROUTINE RETURNS THE TWO-CENTER OVERLAP INTEGRALS NEEDED IN   * C
C * IN THE SCF ROUTINE. IT IS CALLED FROM THE ONELEC ROUTINE.           * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   IC1 : IS THE OF THE CENTER THE ORBITALS OF WHICH ARE OVERLAPING   * C
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

      SUBROUTINE OCOVLS(IC1, NORBAT,NLMAT,ZETAAT, XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"

      INTEGER ZETSTR
      DIMENSION XCORE(N_ORB, N_ORB)
      DIMENSION NORBAT(*), NLMAT(*), ZETAAT(*)


cccc  WRITE(*, 1)IC1, IC1
c1    FORMAT('OCOVLS : ', 2(I2, 2X))

C.... DETERMINATION OF THE STARTING INDICES OF THE ORBITALS AND THE SLATER
C.... EXPONENTS

      CALL NDXSTR(NORBAT, IC1, NLMSTR, ZETSTR)

      CALL OCOVL(NORBAT(IC1),NLMAT(NLMSTR),ZETAAT(ZETSTR),XCORE)

      RETURN
      END
