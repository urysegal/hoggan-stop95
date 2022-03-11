CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 19 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE COMPUTES ALL THE INTEGRALS NEEDED IN THE SCF        * C
C * PROCEDURE.                                                          * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   NATOM : IS THE NUMBER OF ATOMS INVOLVED IN THE MOLECULE           * C
C *                                                                     * C
C *   NCHARG : IS THE ELECTRONIC CHARGE OF THE MOLECULE                 * C
C *                                                                     * C
C *   NTORB : IS THE TOTAL NUMBER OF ATOMIC ORBITALS USED               * C
C *                                                                     * C
C *   NEPSIL : IS THE ORDER OF THE EPSILON-ALGORITHM ACCELERATOR        * C
C *                                                                     * C
C *   LM123 : IS THE MAXIMUM OF ALL L's                                 * C
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
C *   ZETAMY : IS AN ARRAY CONTAINING THE ARITHMETIC MEAN OF THE ZETAs  * C
C *                                                                     * C
C *   NUMAT : IS AN ARRAY CONTAINING THE ATOMIC NUMBERS                 * C
C *                                                                     * C
C *   XYZAT : IS AN ARRAY CONTAINING THE CARTESIAN COORDINATES OF THE   * C
C *    ATOMS INVOLVED IN THE MOLECULE                                   * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE INTGRL(NATOM,NTORB,NCHARG,NEPSIL,LM123,NORBAT,NLMAT,
     $                                        ZETAAT,ZETAMY,NUMAT,XYZAT)
      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*4 GONE
      INTEGER YIWAN

      INCLUDE "./lib95/SIZE.INCL"

      DIMENSION SOVL(N_ORB, N_ORB), HCORE(N_ORB, N_ORB), TEINT(NTBIEL)
      DIMENSION HZCORE(N_ORB, N_ORB)

      DIMENSION NORBAT(*),NLMAT(*),ZETAAT(*),ZETAMY(*),NUMAT(*),XYZAT(*)
      DIMENSION FATIMA(NTBIEL)

      COMMON/OPTION/NWRITE
      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:100)
      COMMON/XGONE/GONE((LMDMAX+1)**2 * (LMDMAX+6+1)**2, LMDMAX+3+1)
      COMMON/POWUN/YIWAN(0:100), LMAXX2


      CALL YGAUNT(GONE, YIWAN, LMAXX2)


C.... COMPUTATION OF THE NUCLEAR REPULSION ENERGY WITH ENUC FUNCTION

      WRITE(*, 1)
 1    FORMAT(/, 27HNUCLEAR REPULSION ENERGY : , /, 43(1H-), /)

      ENC = ENUC(NATOM, NUMAT, XYZAT)
c
      WRITE(*, 2)ENC
 2    FORMAT(/, 33HTOTAL NUCLEAR REPULSION ENERGY : , D16.10, /,
     $          49(1H-), //)

C.... COMPUTING AND COLLECTING ALL SORTS OF ONE-ELCTRON INTEGRALS
C....                       + OVERLAP INTEGRALS

cc      call cpu_time(cpin)

      CALL ONELEC(NATOM,NTORB,NEPSIL,LM123,NORBAT,NLMAT,ZETAAT,
     $     ZETAMY,NUMAT,XYZAT,SOVL,HCORE,HZCORE)
      
cc      call cpu_time(cpout)
      WRITE(*, 3) CPOUT-CPIN
 3    FORMAT(45HELAPSED TIME FOR THE ONE-ELECTRON INTEGRALS :,D20.12, /)

C.... OUTPUTING SOVL AND HCORE
      write(*,*) 'OVERLAP MATRIX'
      call MATOUT(SOVL, NTORB, NTORB)
      
      write(*,*) 'HCORE MATRIX'
      call MATOUT(HCORE, NTORB, NTORB)
      
      
C.... COMPUTING AND COLLECTING ALL SORTS OF TWO-ELCTRON INTEGRALS
cc      call cpu_time(cpin)
      CALL TOELEC(NATOM,NTORB,NEPSIL,LM123,NORBAT,NLMAT,ZETAAT,
     $     ZETAMY, XYZAT, FATIMA)
cc      call cpu_time(cpout)


      WRITE(*, 4) CPOUT-CPIN
 4    FORMAT(45HELAPSED TIME FOR THE TWO-ELECTRON INTEGRALS :,D20.12, /)

c      IF(NWRITE .NE. 0)THEN
c      CALL STPOUT(SOVL, HCORE, FATIMA)
c      ENDIF

cc      call cpu_time(cpin)
      CALL SCFCLO(NATOM,NCHARG,NUMAT,NORBAT,ENC,SOVL,HCORE,
     &  HZCORE,FATIMA)
cc      call cpu_time(cpout)

c      WRITE(*, 5) CPOUT-CPIN
c 5    FORMAT('ELAPSED TIME FOR THE SCF :',D20.12, /)

      RETURN
      END

