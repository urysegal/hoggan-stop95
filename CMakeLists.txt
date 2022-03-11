cmake_minimum_required(VERSION 3.5)

project(STOP)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#    set(dialect "-std=legacy -g -fsanitize=address")
    set(dialect "-std=legacy -g")
    set(bounds "-fbounds-check")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

add_compile_options(-fdefault-integer-8 -m64)




file (GLOB Base_OBJ 
	src/CPUTIM.f
	src/BASIS.f
	src/ONELEC.f 
	src/INTGRL.f
	src/NDXSTR.f
	src/OUT.f
	src/TOELEC.f
	src/STOP95.f
	src/YGAUNT.f
	src/dglgq.f
	src/dglnq.f
 	src/SCFCLO.f     
)

file (GLOB Overlap_OBJ 
	src/overlp/OVERLP.f
        src/overlp/GEGEN0.f
	src/overlp/OVELP1.f
	src/overlp/ROVELP.f
        src/overlp/OCOVLS.f
        src/overlp/OCOVL.f 
)

file (GLOB One_Center_OBJ
	 src/monocent/OCOEIS.f 
	 src/monocent/OCOEI.f
)

file (GLOB Two_Center_One_e_OBJ
	 src/twocent/TCOEIZ.f  
	 src/twocent/GEGEN8.f  
	 src/twocent/TCOEI.f   
	 src/twocent/RTCOEI.f
)

file (GLOB Two_Center_Two_e_OBJ
	 src/twocent/TCOEIY.f  
	 src/twocent/TCOEI2.f
)

file (GLOB Three_Center_One_e_OBJ
	 src/tri31/TRI31.f  
	 src/tri31/GEGEN1.f 
	 src/tri31/RTCNAI.f 
	 src/tri31/TRIMON.f
)


file (GLOB Kinetic_OBJ
	 src/kinetic/KINE12.f 
	 src/kinetic/EKINE1.f 
	 src/kinetic/DUMBID.f 
	 src/kinetic/EKINE2.f 
	 src/kinetic/OVLBID.f
)


file (GLOB One_Center_2_OBJ 
	 src/monocent/OCTEIS.f  
	 src/monocent/OCTEI.f  
	 src/monocent/REAL01.f
	 )


file (GLOB Coulomb_Two_e_OBJ
	 src/coultwo/TCCIS.f   
	 src/coultwo/GEGEN3.f  
	 src/coultwo/GUSED.f   
	 src/coultwo/RTCCI.f   
	 src/coultwo/TCCINS.f  
	 src/coultwo/TCCI.f    
	 src/coultwo/REAL02.f
	 )

file (GLOB Hybrid_OBJ 
	 src/hybrid/HYBR.f    
	 src/hybrid/GEGEN2.f  
	 src/hybrid/HYBRD.f   
	 src/hybrid/HYBRID.f  
	 src/hybrid/RHYB.f    
	 src/hybrid/SYMM03.f  
	 src/hybrid/REAL03.f
	 )

file (GLOB Exchange_OBJ
	 src/exch2/EXCH2.f   
	 src/exch2/DUMMY2.f  
	 src/exch2/EXCH21.f  
	 src/exch2/FUN05.f   
	 src/exch2/FUN06.f   
	 src/exch2/GEGEN5.f  
	 src/exch2/RADP.f    
	 src/exch2/TCCDP.f   
	 src/exch2/REAL04.f
	 )



file (GLOB TRI_32_OBJ
	 src/tri32/TCTEIS.f   
	 src/tri32/GEGEN4.f   
	 src/tri32/HUSED.f    
	 src/tri32/RTCTEI.f   
	 src/tri32/TCTEI.f    
	 src/tri32/TRI32.f    
	 src/tri32/SYMM05.f   
	 src/tri32/REAL05.f
	 )

 file (GLOB Exchange_3_Center_OBJ
	 src/exch3/EXCH3.f    
	 src/exch3/DUMMY3.f   
	 src/exch3/EXCH31.f   
	 src/exch3/FUN07.f    
	 src/exch3/FUN08.f    
	 src/exch3/GEGEN6.f   
	 src/exch3/RADP1.f    
	 src/exch3/TCCDP2.f   
	 src/exch3/TIZIRI.f   
	 src/exch3/REAL06.f
	 )

file (GLOB Exchange_4_Center_OBJ
	 src/tetra/BESPRO.f   
	 src/tetra/DUMMY4.f   
	 src/tetra/DUMMY5.f   
	 src/tetra/FCINTS.f   
	 src/tetra/FUN09.f    
	 src/tetra/FUN10.f    
	 src/tetra/GEGDER.f   
	 src/tetra/GEG.f      
	 src/tetra/GEGEN9.f   
	 src/tetra/HUSES2.f   
	 src/tetra/OUADHIAS.f 
	 src/tetra/RADFC.f    
	 src/tetra/REAL07.f   
	 src/tetra/STOPAT.f   
	 src/tetra/TCCDP5.f   
	 src/tetra/TCCDP6.f   
	 src/tetra/TETRA.f    
	 src/tetra/TETRAC.f   
	 src/tetra/YLMV2.f
	 )

 file (GLOB SCF_OBJ
	 src/scf/CHDENS.f     
	 src/scf/DENSM.f      
	 src/scf/ECORE.f      
	 src/scf/ENUC.f       
	 src/scf/eign.f       
	 src/scf/FOCKER.f     
	 src/scf/jacobi.f     
	 src/scf/MATOUT.f     
	 src/scf/PREJAC.f     
	 src/scf/SAVEM.f      
 #	 src/scf/SCFCLO.f     
	 src/scf/TEENER.f     
	 src/scf/VECTOUT.f
	 )

file (GLOB Lib_OBJ
	 src/lib95/ABIN03.f
	 src/lib95/ACCEL.f
         src/lib95/BAD01.f
         src/lib95/BA01.f
	 src/lib95/BESSI.f
	 src/lib95/BESSK.f
         src/lib95/BLKDT1.f
	 src/lib95/BLKDT2.f
	 src/lib95/BLKDT3.f
	 src/lib95/BLKDT4.f
	 src/lib95/CARTCO.f
	 src/lib95/CM2REA.f
	 src/lib95/CMBINT.f
	 src/lib95/CMINT0.f
	 src/lib95/CMINT1.f
	 src/lib95/CMINT2.f
	 src/lib95/CMINT3.f
	 src/lib95/CMINT4.f
	 src/lib95/COEFFS.f
	 src/lib95/DUMMY.f
	 src/lib95/EPS94.f
	 src/lib95/FUSED.f
	 src/lib95/GAUNT.f
	 src/lib95/GLROOT.f
	 src/lib95/HERMITE.f
	 src/lib95/ISBOOL.f
	 src/lib95/OCCDP.f
	 src/lib95/OCCDP1.f
	 src/lib95/OCCDP2.f
	 src/lib95/RAZ.f
	 src/lib95/SCONST.f
	 src/lib95/SPHRM.f
	 src/lib95/YLMV.f
	 )

add_executable(stop ${Base_OBJ} ${Overlap_OBJ}  ${One_Center_OBJ} ${Two_Center_One_e_OBJ} ${Two_Center_Two_e_OBJ} ${Three_Center_One_e_OBJ} ${Kinetic_OBJ} ${One_Center_2_OBJ} ${Coulomb_Two_e_OBJ} ${Hybrid_OBJ} ${Exchange_OBJ} ${TRI_32_OBJ} ${Exchange_3_Center_OBJ} ${Exchange_4_Center_OBJ} ${SCF_OBJ}  ${Lib_OBJ} )

set(BLA_VENDER OpenBLAS)
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
    message("OpenBLAS found.")
    target_link_libraries(stop ${BLAS_LIBRARIES})
endif(BLAS_FOUND)

target_link_libraries(stop ${BLAS_LIBRARIES})

