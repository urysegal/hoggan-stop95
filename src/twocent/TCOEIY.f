CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 22 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE COMPUTES THE TWO-CENTER ONE-ELECTRON INTEGRALS OF   * C
C * THE FIRST KIND INVOLVED IN THE SCF PROCEDURE. IT IS CALLED FROM THE * C
C * THE ONELEC ROUTINE.                                                 * C
C *                       < XA | 1/RB | XA >                            * C
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
C *   XYZAT : IS AN ARRAY CONTAINING THE CARTESIAN COORDINATES OF THE   * C
C *    ATOMS INVOLVED IN THE MOLECULE                                   * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   XCORE : IS A TWO-DIMENSIONAL ARRAY CONTAINING THE IMAGINARY CORE  * C
C *    INTEGRALS                                                        * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TCOEIY(IC1,IOP,IC2,NUMAT,NORBAT,NLMAT,ZETAAT,XYZAT,
     $                                                            XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER ZETSTR
      INCLUDE "../lib95/SIZE.INCL"

      DIMENSION XCORE(N_ORB, N_ORB)
      DIMENSION NUMAT(*),NORBAT(*),NLMAT(*),ZETAAT(*),XYZAT(*)
      DIMENSION XYZV(3)


ccc   WRITE(*, 1)IC1, IOP, IC2
c1    FORMAT('TCOEIY :', 3(I2, 2x))

C.... COMPUTATION OF THE COORDINATES OF THE COULOMB OPERATOR WITH RESPECT
C.... TO THE COORDINATES OF THE CENTER LABELED IC1

      NIC1 = 3*(IC1-1) + 1
      NIOP = 3*(IOP-1) + 1
      CALL CARTCO(XYZAT(NIC1), XYZAT(NIOP), XYZV, RVMOD,THETRV,PHIRV)


C.... DETERMINING THE STARTING INDEX FOR THE ORBITALS AND THE SLATER EXPONENTS

      CALL NDXSTR(NORBAT, IC1, NLMSTR, ZETSTR)

      CALL TCOEI2(NUMAT(IOP),NORBAT(IC1),NLMAT(NLMSTR),ZETAAT(ZETSTR),
     $                                        RVMOD,THETRV,PHIRV, XCORE)

      RETURN
      END
