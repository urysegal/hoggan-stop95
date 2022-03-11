      SUBROUTINE TOELEC(NATOM,NTORB,NEPSIL,LM123,NORBAT,NLMAT,ZETAAT,
     $                                            ZETAMY, XYZAT, FATIMA)
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "./lib95/SIZE.INCL"

      DIMENSION NORBAT(*), NLMAT(*), ZETAAT(*), ZETAMY(*), XYZAT(*),
     $                                                         FATIMA(*)

      COMMON/INTEGN/NINT

      NSUM = (NTORB * (NTORB+1))/2

      write(*,*) ' ONE-CENTER TWO-ELECTRONS INTEGRALS '
      write(*,*) '------------------------------------'
      DO 5 IC1=1, NATOM
       IC2 = IC1
       IC3 = IC1
       IC4 = IC1
       CALL OCTEIS(NTORB,NSUM, LM123, IC1,IC2,IC3,IC4,NORBAT,NLMAT,
     $      ZETAAT, FATIMA)
 5    CONTINUE
      write(*,*) '------------------------------------'
      write(*,*)

      DO 10 IC1=1, NATOM
      DO 10 IC2=IC1+1, NATOM
            
       write(*,*) ' TWO-CENTER TWO-ELECTRONS INTEGRALS '
       write(*,*) '------------------------------------'
       CALL TCCIS(NEPSIL, NTORB,NSUM, LM123, IC2,IC1, NORBAT,
     $      NLMAT, ZETAAT,XYZAT, FATIMA)

       write(*,*)
       write(*,*) '   HYBRID TWO-ELECTRONS INTEGRALS   '
       write(*,*) '------------------------------------' 
       
       CALL HYBR(NEPSIL, NTORB,NSUM, LM123, IC1,IC2, NORBAT,
     $      NLMAT, ZETAAT,XYZAT, FATIMA)
       
       CALL HYBR(NEPSIL, NTORB,NSUM, LM123, IC2,IC1, NORBAT,
     $      NLMAT, ZETAAT,XYZAT, FATIMA)
       
       write(*,*)
       write(*,*) '   EXCH2 TWO-ELECTRONS INTEGRALS    '
       write(*,*) '------------------------------------' 
       CALL EXCH2(NEPSIL, NTORB,NSUM, LM123, IC1,IC2, NORBAT,
     $      NLMAT, ZETAAT, ZETAMY, XYZAT, FATIMA)
       
 10   CONTINUE

      DO 15 IC1=1, NATOM
      DO 15 IC2=IC1+1, NATOM
      DO 15 IC3=IC2+1, NATOM
      CALL TCTEIS(NEPSIL, NTORB,NSUM, LM123, IC1,IC2,IC3, NORBAT,
     $                           NLMAT, ZETAAT, ZETAMY, XYZAT, FATIMA)
      CALL TCTEIS(NEPSIL, NTORB,NSUM, LM123, IC2,IC1,IC3, NORBAT,
     $                           NLMAT, ZETAAT, ZETAMY, XYZAT, FATIMA)
      CALL TCTEIS(NEPSIL, NTORB,NSUM, LM123, IC3,IC1,IC2, NORBAT,
     $                           NLMAT, ZETAAT, ZETAMY, XYZAT, FATIMA)

      CALL EXCH3(NEPSIL, NTORB,NSUM, LM123, IC1,IC2,IC3, NORBAT,
     $                             NLMAT, ZETAAT, ZETAMY, XYZAT, FATIMA)

      CALL EXCH3(NEPSIL, NTORB,NSUM, LM123, IC2,IC3,IC1, NORBAT,
     $                             NLMAT, ZETAAT, ZETAMY, XYZAT, FATIMA)

      CALL EXCH3(NEPSIL, NTORB,NSUM, LM123, IC3,IC2,IC1, NORBAT,
     $                             NLMAT, ZETAAT, ZETAMY, XYZAT, FATIMA)

 15   CONTINUE

      DO 20 IC1=1, NATOM
      DO 20 IC2=IC1+1, NATOM
      DO 20 IC3=IC2+1, NATOM
      DO 20 IC4=IC3+1, NATOM

       CALL TETRA(NEPSIL, NTORB,NSUM, LM123, IC1,IC2, IC3,IC4,
     $                     NORBAT, NLMAT, ZETAAT, ZETAMY, XYZAT, FATIMA)

      CALL TETRA(NEPSIL, NTORB,NSUM, LM123, IC1,IC3, IC2,IC4,
     $                     NORBAT, NLMAT, ZETAAT, ZETAMY, XYZAT, FATIMA)

       CALL TETRA(NEPSIL, NTORB,NSUM, LM123, IC1,IC4, IC2,IC3,
     $                     NORBAT, NLMAT, ZETAAT, ZETAMY, XYZAT, FATIMA)

 20   CONTINUE


 100  RETURN
      END
