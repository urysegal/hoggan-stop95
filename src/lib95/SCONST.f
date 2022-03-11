      SUBROUTINE SCONST()                                                       
      IMPLICIT REAL*8 (A-H, O-Z)                                                
      COMMON/FACT0/FACT(0:30)                                                   
      COMMON/DFACT/DFAC(0:100)
      COMMON/CNP0/CNP(0:30, 0:30)
      COMMON/GAMMA0/GAMA(0:30)

      DATA PI, SPI/3.14159265358979324D0, 1.77245385090551603D0/

      FACT(0) = 1.0D0                                                           
      DO 5 I=1, 30                                                              
       FACT(I) = FACT(I-1) * DBLE(I)                                            
 5    CONTINUE                                                                  
                                                                                
      DFAC(0) = 1.0D0                                                           
      DO 10 I=1, 100
       DFAC(I) = DFAC(I-1) * DBLE(2*I+1)                                        
 10   CONTINUE                                                                  
                                                                                
      CNP(0,0) = 1.D0                                                           
      DO 15 I=1, 30                                                             
       CNP(I, 0) = 1.D0                                                         
       CNP(I, I) = 1.D0                                                         
       DO 15 J=1, I                                                             
        CNP(I, J) = CNP(I-1, J-1) + CNP(I-1, J)                                 
 15   CONTINUE                                                                  
                                                                                
      GAMA(0) = SPI
      DO 2 I=1, 30
       GAMA(I) = GAMA(I-1) * (DBLE(I) - 0.5D0)
 2    CONTINUE

      RETURN                                                                    
      END                                                                       
