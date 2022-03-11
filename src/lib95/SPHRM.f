      SUBROUTINE SPHRM(L, M, COSINE, YLM)
      IMPLICIT REAL*8 (A-H, O-Z)                                                
      DIMENSION YLM(0:*)                                                        
      DATA SPI/1.77245385090551603D0/                                           
                                                                                
                                                                                
      DM = DBLE(M)                                                              
      MA = IABS(M)                                                              
                                                                                
      DO 5 I=0, MA-1                                                            
       YLM(I) = 0.D0                                                            
 5    CONTINUE                                                                  
                                                                                
      SINE2 = 1.D0 - COSINE*COSINE                                              
      PROD = 1.D0                                                               
                                                                                
      DO 10 I=1, MA                                                             
       DII  = DBLE(I)                                                           
       PROD = (1.D0 - 0.5D0/DII) * SINE2 * PROD                                 
 10   CONTINUE                                                                  
                                                                                
      YLM(MA) = 0.5D0/SPI * DSQRT(2.D0*MA + 1.D0) * DSQRT(PROD)                 
                                                                                
      IF(M .GE. 0)THEN
       IF(2*(MA/2) .EQ. MA)THEN
	DMA = 1.D0
       ELSE
	DMA = -1.D0
       ENDIF
       YLM(MA) = DMA * YLM(MA)
      ENDIF
                                                                                
      DO 20 LL=MA+1, L                                                          
       DLL = DBLE(LL)                                                           
       CO1 = (2.D0*DLL + 1.D0)*(2.D0*DLL - 1.D0)                                
       CO2 = (2.D0*DLL + 1.D0)/(2.D0*DLL - 3.D0) *                              
     $                         (DLL + DM - 1.D0)*(DLL - DM - 1.D0)              
       YLM(LL) = 1.D0/DSQRT((DLL + DM)*(DLL - DM)) *                            
     $       (DSQRT(CO1) * COSINE * YLM(LL-1) - DSQRT(CO2) * YLM(LL-2))         
 20   CONTINUE                                                                  
                                                                                
      RETURN                                                                    
      END                                                                       
