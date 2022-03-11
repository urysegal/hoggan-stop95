CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 17 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS ROUTINE COMPUTES A STREAM OF TWO-CENTER ONE-ELECTRON INTEGRALS * C
C * OVER A SLATER TYPE ORBITAL BASIS SET :                              * C
C *                    < XA | 1/RB | XA >                               * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   NUMAT : IS AN ARRAY CONTAINING THE ATOMIC NUMBERS OF THE ATOMS    * C
C *    INVOLVED IN THE MOLECULE                                         * C
C *                                                                     * C
C *   NORB1 : TOTAL NUMBER OF STOs ATTACHED TO THE CENTER (1)           * C
C *                                                                     * C
C *   NLMA, ZETAA : ARRAYS CONTAINING THE DESCRIPTION OF THE ORBITALS,  * C
C *    i.e. QUANTUM NUMBERS AND SLATER EXPONENTS                        * C
C *                                                                     * C
C *   AB, THETAB, PHIAB : SPHERICAL COORDINATES OF THE VECTOR SEPARATING* C
C *    THE CENTERS                                                      * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   XCORE : IS A TWO-DIMENSIONAL ARRAY CONTAINING THE CORE INTEGRALS  * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TCOEI2(NUMAT,NORB1,NLMA,ZETAA,AB,THETAB,PHIAB,XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      INCLUDE "../lib95/CONST.INCL"

      DIMENSION XCORE(N_ORB,N_ORB)
      DIMENSION NLMA(*),ZETAA(*),YLMAB(0:((LDEV+1)*(LDEV+2))/2),
     $                                                    ARGNT(L1MAX+1)

      COMMON/STOP09/ZETAS, ABP, N12, L
      COMMON/FACT0/FACT(0:30)
      COMMON/DFACT/DFAC(0:100)
      COMMON/GAUSOR/NORDER(0:99)
      EXTERNAL RAD11, RAD12


C.... PASSAGE TO COMMON

      ABP = AB

C.... COMPUTATION OF THE REQUIRED SPHERICAL HARMONICS

      CALL YLMV0(THETAB, YLMAB)

C.... COMBINIG THE GIVEN SET OF ATOMIC ORBITALS

      DO 5 NAT1=1, NORB1
       IST   = 5*(NAT1-1)
       LINE  = NLMA(IST+2)
       N1    = NLMA(IST+3)
       L1    = NLMA(IST+4)
       M1    = NLMA(IST+5)
       ZETA1 = ZETAA(NAT1)

      DO 5 NAT2=NAT1, NORB1
       JST   = 5*(NAT2-1)
       ICOL  = NLMA(JST+2)
       N2    = NLMA(JST+3)
       L2    = NLMA(JST+4)
       M2    = NLMA(JST+5)
       ZETA2 = ZETAA(NAT2)

       N12   = N1 + N2
       ZETAS = ZETA1 + ZETA2

C.... NORMALIZATION CONSTANTS

       A1 = (2.D0*ZETA1)**N1 * DSQRT((2.D0*ZETA1)/FACT(2*N1))
       A2 = (2.D0*ZETA2)**N2 * DSQRT((2.D0*ZETA2)/FACT(2*N2))

       CSTE   = 4.D0 * PI * A1 * A2
       ANGLE  = (M2-M1)*PHIAB
       COSINE = DCOS(ANGLE)
       SINE   = DSIN(ANGLE)

       IF(DABS(SINE) .LE. 1.D-10)THEN
	SINE = 0.D0
       ELSE
	IF(DABS(COSINE) .LE. 1.D-10)THEN
	 COSINE = 0.D0
	ENDIF
       ENDIF

       TCOEIV = 0.D0

       CALL GAUNT(L2,M2, L1,M1, LMIN,LMAX,M12, NGAUNT,ARGNT)

       KGN = 1

       DO 10 L=LMIN, LMAX, 2
	NYL    = (L*(L+1))/2 + IABS(M12)
	MUPW   = (M12 - IABS(M12))/2
	AUX    = (-1)**MUPW * 1.D0/DBLE(L+L+1)*YLMAB(NYL)*ARGNT(KGN)
        KGN    = KGN + 1
	XINT1  = DGLNQ(RAD11, 0.D0, AB, LEG15)
	XINT2  = DGLGQ(RAD12, AB, ZETAS, NORDER((N1+N2-L-1)/2+1))
	TCOEIV = TCOEIV + AUX * (XINT1 + XINT2)
 10    CONTINUE

       TCOEIR = COSINE * CSTE * TCOEIV
       TCOEII = SINE   * CSTE * TCOEIV

cccc   WRITE(*, 1)LINE,ICOL, NUMAT*TCOEIR, NUMAT*TCOEII

       XCORE(LINE, ICOL) = XCORE(LINE, ICOL) - NUMAT*TCOEIR
       IF(LINE.NE.ICOL)THEN
	XCORE(ICOL, LINE) = XCORE(ICOL, LINE) - NUMAT*TCOEII
       ENDIF

 5    CONTINUE

c1    FORMAT(2Z2, 1X, 2(D16.10, 2X))

      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C * THIS ROUTINE COMPUTES THE RADIAL INTEGRAL INVOLVED IN THE ONE-CENTER* C
C * ONE-ELCTRON INTEGRALS                                               * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RAD11(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION T(*), Y(*)
      COMMON/STOP09/ZETAS, AB, N12, L


      DO 5 I=1, N
       R    = T(I)
       Y(I) = R**(N12-1) * (R/AB)**(L+1) * DEXP(-ZETAS*R)
 5    CONTINUE

      RETURN
      END



      SUBROUTINE RAD12(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION T(*), Y(*)
      COMMON/STOP09/ZETAS, AB, N12, L


      DO 5 I=1, N
       R = T(I)
       Y(I) = R**(N12-1) * (AB/R)**L * DEXP(-ZETAS*R)
 5    CONTINUE

      RETURN
      END
