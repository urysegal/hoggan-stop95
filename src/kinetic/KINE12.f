CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 23 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE COMPUTES THE ONE-ELECTRON INTEGRALS NEEDED IN THE   * C
C * COMPUTATION OF THE ONE-CENTER KINETIC ENERGY, I.E. :                * C
C *                     < XA | -1/2 D | XA >                            * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   IC1, IC2 : ARE THE RANKS OF THE ATOMS THE ORBITALS OF WHICH  ARE  * C
C *    INTERACTING OF THE NUCLEUS IOP                                   * C
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

      SUBROUTINE KINET1(IC1, IC2, NORBAT, NLMAT, ZETAAT, XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER ZETSTR
      INCLUDE "../lib95/SIZE.INCL"

      DIMENSION NORBAT(*), NLMAT(*), ZETAAT(*)
      DIMENSION XCORE(N_ORB, N_ORB)


ccc   WRITE(*, 1)IC1, IC2
c1    FORMAT('KINETIC 1 : ', 2(I2, 2X))

C.... DETERMINING THE STARTING INDEX FOR THE ORBITALS AND THE SLATER EXPONENTS

      CALL NDXSTR(NORBAT, IC1, NLMSTR, ZETSTR)

      CALL EKINE1(NORBAT(IC1),NLMAT(NLMSTR),ZETAAT(ZETSTR),XCORE)

      RETURN
      END




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C * THIS SUBROUTINE COMPUTES THE ONE-ELECTRON INTEGRALS NEEDED IN THE   * C
C * COMPUTATION OF THE TWO-CENTER KINETIC ENERGY, I.E. :                * C
C *            < XA | -1/2 D | XB > = < -1/2 D. XB | XA >               * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE KINET2(IC1,IC2, NORBAT,NLMAT,ZETAAT,XYZAT,XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER ZETST1, ZETST2
      INCLUDE "../lib95/SIZE.INCL"

      DIMENSION NORBAT(*), NLMAT(*), ZETAAT(*), XYZAT(*)
      DIMENSION XYZV(3)
      DIMENSION XCORE(N_ORB, N_ORB)


ccc   WRITE(*, 1)IC1, IC2
c1    FORMAT('KINETIC 2 : ', 2(I2, 2X))

C.... COMPUTATION OF THE COORDINATES OF THE CENTER LABELED IC1 TO WHICH THE
C.... ATOMIC ORBITAL XA IS ATTACHED (WITH RESPECT TO XB)

      NIC1 = 3*(IC1-1)+1
      NIC2 = 3*(IC2-1)+1
      CALL CARTCO(XYZAT(NIC2), XYZAT(NIC1), XYZV, RVMOD,THETRV,PHIRV)

C.... DETERMINING THE STARTING INDEX FOR THE ORBITALS AND THE SLATER EXPONENTS

      CALL NDXSTR(NORBAT, IC1, NLMST1, ZETST1)
      CALL NDXSTR(NORBAT, IC2, NLMST2, ZETST2)

      CALL EKINE2(NORBAT(IC1),NORBAT(IC2), NLMAT(NLMST1),NLMAT(NLMST2),
     $   ZETAAT(ZETST1),ZETAAT(ZETST2), RVMOD,THETRV,PHIRV, XCORE)

      RETURN
      END
