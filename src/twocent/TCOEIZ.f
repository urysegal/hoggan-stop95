CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 22 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE COMPUTES THE TWO-CENTER ONE-ELECTRON INTEGRALS      * C
C * OF THE SECOND KIND INVOLVED IN THE SCF PROCEDURE. ITS IS CALLED     * C
C * CALLED FROM THE INTGRL ROUTINE.                                     * C
C *             < XA | 1/RA | XB > <==> < XB | 1/RA | XA >              * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LM123 : IS THE HIGHEST VALUE OF THE QUANTUM NUMBER L              * C
C *                                                                     * C
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

      SUBROUTINE TCOEIZ(LM123,IC1,IOP,IC2,NUMAT,NORBAT,NLMAT,ZETAAT,
     $                                                      XYZAT,XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER ZETST1, ZETST2
      INCLUDE "../lib95/SIZE.INCL"

      DIMENSION XCORE(N_ORB, N_ORB)
      DIMENSION NUMAT(*),NORBAT(*),NLMAT(*),ZETAAT(*),XYZAT(*)
      DIMENSION XYZV(3)


ccc   WRITE(*, 1)IC1, IOP, IC2
c1    FORMAT('TCOEIZ :', 3(I2, 2x))

C.... COMPUTING THE COORDINATES OF THE CENTER LABELED IC2 TO WHICH THE
C.... ATOMIC ORBITAL XB IS ATTACHED

      NIC1 = 3*(IC1-1) + 1
      NIC2 = 3*(IC2-1) + 1
c
      CALL CARTCO(XYZAT(NIC1), XYZAT(NIC2), XYZV, RVMOD,THETRV,PHIRV)
c
C.... DETERMINING THE STARTING INDEX FOR THE ORBITALS AND THE SLATER EXPONENTS
c
      CALL NDXSTR(NORBAT, IC1, NLMST1, ZETST1)
c
      CALL NDXSTR(NORBAT, IC2, NLMST2, ZETST2) 
c
      CALL TCOEI(LM123, NUMAT(IOP),NORBAT(IC1),NORBAT(IC2),
     $     NLMAT(NLMST1),NLMAT(NLMST2), ZETAAT(ZETST1),ZETAAT(ZETST2),
     $                                        RVMOD,THETRV,PHIRV,XCORE)
c
c
      RETURN
      END
