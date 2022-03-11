CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE RETURNS THE ONE-CENTER TWO-ELECTRON INTEGRALS NEEDED   * C
C * IN THE SCF ROUTINE. IT IS CALLED FROM THE TOELEC ROUTINE.           * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   LM123 : IS THE HIGHEST VALUE OF THE QUANTUM NUMBER L              * C
C *                                                                     * C
C *   IC1, IC2, IC3, IC4 : ARE THE CENTERS TO WHICH THE FOUR ORBITALS   * C
C *    ARE ATTACHED. IN THIS CASE IC1=IC2=IC3=IC4                       * C
C *                                                                     * C
C *   NORBAT : IS AN ARRAY CONTAINING THE NUMBER OF ORBITALS DESCRIBING * C
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
C *   TWOEL : IS A ONE-DIMENSIONAL ARRAY CONTAINING THE TWO-ELECTRON    * C
C *    INTEGRALS OF THE FORM <I J | K L>.                               * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE OCTEIS(NTORB,NSUM, LM123,IC1,IC2,IC3,IC4,NORBAT,NLMAT,
     $                                                   ZETAAT, FATIMA)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER ZETSTR
      COMPLEX*16 TWOEL
      DIMENSION NORBAT(*), NLMAT(*), ZETAAT(*), FATIMA(*)

cccc  WRITE(*, 1)IC1, IC2, IC3, IC4
c1    FORMAT('OCTEIS : ', 4(I2, 2X))

C.... DETERMINING THE STARTING INDICES FOR THE ORBITALS AND THE SLATER
C.... EXPONENTS

      CALL NDXSTR(NORBAT, IC1, NLMSTR, ZETSTR)

      CALL OCTEI(NTORB,NSUM, LM123,NORBAT(IC1),NLMAT(NLMSTR),
     $                                           ZETAAT(ZETSTR), FATIMA)

      RETURN
      END
