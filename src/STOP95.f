CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *                                                                     * C
C *              STOP !                                                 * C
C *                                                                     * C
C *         SLATER TYPE ORBITAL PACKAGE                                 * C
C *                                                                     * C
C *         PROJECT DIRECTED BY PHILIP HOGGAN                           * C
C *                                                                     * C
C *        BASED ON THE ORIGINAL METHOD IN Int. J. Quantum Chem         * C
C *        8 ICQC PRAGUE 19-23/6/94 Vol 51 By A.Bouferguene, M.Fares    * C
C *        and P. E. Hoggan.                                            * C
C *                                                                     * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 19 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C *     Any correspondence shold be addressed to P. E. Hoggan           * C
C *     URA 414, ISMRa, 6 Bd MARECHAL JUIN                              * C
C *     14050 CAEN CEDEX FRANCE                                         * C
C *                                                                     * C
C * THIS WILL BE SLATER TYPE ORBITAL PACKAGE.                           * C
C *   THIS ROUTINE READS THE INPUT DATA AND CALLS THE INTGRL            * C
C * ROUTINE.                                                            * C
C *                                                                     * C
C * DESCRIPTION :                                                       * C
C *   RTPAT : IS AN ARRAY CONTAINING THE SPHERICAL COORDINATES OF THE   * C
C *    ATOMS INVOLVED IN THE MOLECULE                                   * C
C *                                                                     * C
C *   XYZAT : IS AN ARRAY CONTAINING THE CARTESIAN COORDINATES OF THE   * C
C *    ATOMS INVOLVED IN THE MOLECULE                                   * C
C *                                                                     * C
C *   NUMAT : IS AN ARRAY CONTAINING THE ATOMIC NUMBERS                 * C
C *                                                                     * C
C *   NLMAT : IS AN ARRAY CONTAINING THE FOLLOWING INFORMATIONS :       * C
C *    NLMAT(1) -> ATOMIC NUMBER                                        * C
C *    NLMAT(2) -> ORDER OF THE ORBITAL AS OCCURING IN THE DATA         * C
C *    NLMAT(3), NLMAT(4), NLMAT(5) -> N, L, M                          * C
C *                                                                     * C
C *   ZETAAT : IS AN ARRAY CONTAINING THE SLATER ZETA EXPONENTS         * C
C *                                                                     * C
C *   NORBAT : IS AN ARRAY CONTAINING THE NUMBER OF ORBITALS ATTACHED TO* C
C *    EACH ATOM                                                        * C
C *                                                                     * C
C *   ZETAMY : IS AN ARRAY CONTAINING THE ARITHMETIC MEAN OF THE ZETAs  * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PROGRAM STOP95

      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "./lib95/SIZE.INCL"
      INTEGER DAYS, HOURS, SECDS
      CHARACTER*60 MOLEC
c      REAL*4 ZETA,ZETAAT
      DIMENSION RTPAT(3*MAXAT), XYZAT(3*MAXAT), NUMAT(MAXAT),
     $       NLMAT(5*N_ORB), ZETAAT(N_ORB), NORBAT(MAXAT), ZETAMY(MAXAT)
      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:100)

      DATA RADEG/1.74532925199432958D-2/

C.... INITIALIZING THE TIMER
      call cpu_time(cpin)

C.... HANDLING THE UNDEFLOW ERROR VIA THE IBM E.S.S.L SOFT

C.... CALL ERRSET(208, 300, -1)

C.... SOME FREQUENTLY USED CONSTANTS

      CALL SCONST()

C.... READING THE INPUT PARAMETER FOR THE EPSILON ALGORITHM. THIS PARAMETER
C.... IS USED FOR THE EVALUATION OF THE INTEGRALS :
C.... TRI31, TRI32, EXCH2 AND EXCH3

      READ(*, '(A60)')MOLEC

      READ(*, *, END=100)NEPSIL
C.... READING THE NUMBER OF ATOMS INVOLVED IN THE MOLECULE AND ITS CHARGE

      READ(*, *, END=100)NATOM, NCHARG

C.... READING THE INPUT QUANTUM NUMBERS (BASIS SET)
C.... THE GEOMETRY IS ASSUMED TO BE GIVEN IN SPHERICAL COORDINATES

      IF(NCHARG.EQ.5)then
      DO 95 I=1, NATOM
       NDEX  = 3*(I-1)
       READ(*, *, END=100)XYZAT(NDEX+1),XYZAT(NDEX+2),XYZAT(NDEX+3)
   95 CONTINUE
      goto 1234
      endif
      DO 5 I=1, NATOM

       NSTR = 3*(I-1) + 1
       READ(*, *, END=100)RTPAT(NSTR), RTPAT(NSTR+1), RTPAT(NSTR+2)

C.... COMPUTATION OF THE COORDINATES OF THE CENTERS CONSTITUTING THE THREE-
C.... ATOMIC MOLECULE
       NDEX  = 3*(I-1)
       R     = RTPAT(NDEX+1)
       THETA = RADEG * RTPAT(NDEX+2)
       PHI   = RADEG * RTPAT(NDEX+3)

       XYZAT(NDEX+1) = R * DSIN(THETA) * DCOS(PHI)
       XYZAT(NDEX+2) = R * DSIN(THETA) * DSIN(PHI)
       XYZAT(NDEX+3) = R * DCOS(THETA)

 5    CONTINUE
 1234 continue
C.... READING THE ATOMIC NUMBERS THAT CHARACTERIZE THE ATOMS

      READ(*, *, END=100)(NUMAT(I), I=1, NATOM)

C.... THE FOLLOWING PARAMETERS ARE USED TO DETERMINE THE MAXIMUMS OF THE
C.... SECOND QUANTUM NUMBERS

      LM123  = -1

C.... INITIALIZATION OF SOME DIMENSIONAL ARRAYS

      CALL RAZ0(ZETAMY, 1, NATOM)
      CALL RAZ1(NORBAT, 1, NATOM)

C.... READING THE LOGICAL VARIABLE THAT DETERMINES IF A STANDARD BASIS IS
C.... REQUESTED. IF NBASE = 0, THEN A NON-STANDARD OTHERWISE A STANDARD.

c     READ(*, *, END=100)NBASE
c
      CALL BASIS(1, NATOM,NUMAT,NORBAT,NLMAT,ZETAAT,ZETAMY,LM123)
c
C.... READING THE ATOMIC ORBITALS AND THE SLATER EXPONENTS

      DO 15 I=1, N_ORB
       READ(*, *, END=110)NAT, N, L, M, ZETA

       NORBAT(NAT)   = NORBAT(NAT) + 1
       NDEX          = 5*(I-1)
       NLMAT(NDEX+1) = NUMAT(NAT)
       NLMAT(NDEX+2) = I
       NLMAT(NDEX+3) = N
       NLMAT(NDEX+4) = L
       NLMAT(NDEX+5) = M
       ZETAAT(I)     = ZETA
       ZETAMY(NAT)   = ZETAMY(NAT) + ZETA
       LM123         = MAX(LM123, L)
 15   CONTINUE

 110  CONTINUE

C.... NUMBER OF ATOMIC ORBITALS INVOLVED IN THE DESCRIPTION OF THE MOLECULE

      NTORB = I-1
C.... THIS STEP COMPUTES THE ARITHMETIC MEAN OF THE SLATER EXPONENTS

      DO 20 I=1, NATOM
       ZETAMY(I) = ZETAMY(I)/DBLE(NORBAT(I))
 20   CONTINUE

C.... OUTPUTING SOME INFORMATIONS OF THE MOLECULAR SYSTEM UNDER STUDY

      CALL OUT(MOLEC, NTORB, NLMAT, ZETAAT, NATOM, XYZAT, NUMAT)

      
C.... AND NOW, LADIES AND GENTLEMEN THE EVALUATION OF THE INTEGRALS
C.... I DO NOT KNOW, HOW I WILL DO THIS, I WILL TRY TO DO MY BEST ...

      CALL INTGRL(NATOM,NTORB,NCHARG,NEPSIL,LM123,NORBAT,NLMAT,
     $     ZETAAT,ZETAMY,NUMAT,XYZAT)
      call cpu_time(cpout)
      
      DAYS  = (CPOUT-CPIN)/86400
      HOURS = (CPOUT-CPIN - 86400*DAYS)/3600
      MINUT = (CPOUT-CPIN - 86400*DAYS - 3600*HOURS)/60
      SECDS = IDINT(CPOUT-CPIN - 86400*DAYS - 3600*HOURS - 60*MINUT)
      
      WRITE(*, 2)DAYS, HOURS, MINUT, SECDS
      
      write(*,*) 'TIME =', (CPOUT-CPIN)

 2    FORMAT(29HTOTAL TIME FOR THE PROCESS : ,
     $     I2, 4Hdays, 2X, I2, 1Hh, 2X, I2, 3Hmin, 2X, I2, 3Hsec)

 100  END


