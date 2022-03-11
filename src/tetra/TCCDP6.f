CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 19 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE FIRST PART OF THE RADIAL INTEGRAL       * C
C * INVOLVED IN THE T.C.C.D.P., I.E. [R1=0, R1=R2]                      * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TCCDP6(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "../lib95/SIZE.INCL"
      DIMENSION T(*), Y(*)
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))
      COMMON/EXCH04/ZETA1, R2, N1, LP2, L, LMD2, LMDL12, IO, ISTR

      INDO = (IO-1)*LEA12*LMDL12 + LMD2*LEA12 + ISTR

      DO 5 I=1, N
       R1 = T(I)
       Y(I) = 1.D0/DSQRT(R1) * R1**(N1+LP2) * (R1/R2)**(L+1) *
     $                                  DEXP(-ZETA1*R1) * GEGIN0(INDO+I)
 5    CONTINUE

      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 19 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE RETURNS THE SECOND PART OF THE RADIAL INTEGRAL      * C
C * INVOLVED IN THE T.C.C.D.P., I.E. [R1=R2, R1=+\INFTY[                * C
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TCCDP7(T, Y, N)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION T(*), Y(*)
      INCLUDE "../lib95/SIZE.INCL"
      COMMON/COMG00/GEGIN0(LMHERM*LMHERM*(LMDMAX+L1MAX+L2MAX+1))
      COMMON/EXCH04/ZETA1, R2, N1, LP2, L, LMD2, LMDL12, IO, ISTR


      INDO = (IO-1)*LEA12*LMDL12 + LMD2*LEA12 + ISTR

      DO 5 I=1, N
       R1   = T(I)
       Y(I) = 1.D0/DSQRT(R1) * R1**(N1+LP2) * (R2/R1)**L *
     $                                  DEXP(-ZETA1*R1) * GEGIN0(INDO+I)
 5    CONTINUE

      RETURN
      END
