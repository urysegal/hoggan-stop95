      BLOCK DATA BLKDT1
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "SIZE.INCL"
      DIMENSION ABSC14(L001), ABSC24(L002), ABSC20(L003)
      COMMON/GLE01/RLEG01(L001+L002+L003)

      DATA ABSC14/ .228603188386461D-02, .119275193894044D-01,
     A.287997808217058D-01, .521178491980524D-01, .807918939403076D-01,
     B.113481271845352D-00, .148657508548776D-00, .184675824784557D-00,
     C.219852061487982D-00, .252541439393026D-00, .281215484135281D-00,
     D.304533552511628D-00, .321405813943929D-00, .331047301449469D-00/

      DATA ABSC24/ .334135463333830D+00, .337545240671448D+00,
     A.343620907999545D+00, .352264078832600D+00, .363333002337683D+00,
     B.376645968070241D+00, .391984391343837D+00, .409096421435193D+00,
     C.427701082062326D+00, .447492886717306D+00, .468146855421064D+00,
     D.489323851189566D+00, .510676148810434D+00, .531853144578936D+00,
     E.552507113282694D+00, .572298917937674D+00, .590903578564807D+00,
     F.608015608656163D+00, .623354031929759D+00, .636666997662317D+00,
     G.647735921167400D+00, .656379092000455D+00, .662454759328552D+00,
     H.665864536666170D+00/

      DATA ABSC20/ .667811900135817D+00, .672671345453681D+00,
     A.681294261958112D+00, .693480504696297D+00, .708944682256641D+00,
     B.727324386545581D+00, .748188833008195D+00, .771048985214097D+00,
     C.795369024809726D+00, .820578913144417D+00, .846087753522250D+00,
     D.871297641856941D+00, .895617681452570D+00, .918477833658471D+00,
     E.939342280121086D+00, .957721984410025D+00, .973186161970370D+00,
     F.985372404708554D+00, .993995321212986D+00, .998854766530849D+00/

      EQUIVALENCE (ABSC14(1),RLEG01(1)),(ABSC24(1),RLEG01(L001+1)),
     $            (ABSC20(1),RLEG01(L001+L002+1))


      DIMENSION ABSD08(L004), ABSD14(L005), ABSD24(L006)
      COMMON/GLE02/RLEG02(L004+L005+L006)

      DATA ABSD08/ .661835725041062D-02, .338889204310622D-01,
     A.790779316806118D-01, .136094226250725D+00, .197239107082608D+00,
     B.254255401652721D+00, .299444412902271D+00, .326714976082923D+00/

      DATA ABSD14/ .335619365217198D+00, .345260852722738D+00,
     A.362133114155039D+00, .385451182531386D+00, .414125227273641D+00,
     B.446814605178685D+00, .481990841882109D+00, .518009158117891D+00,
     C.553185394821315D+00, .585874772726359D+00, .614548817468614D+00,
     D.637866885844961D+00, .654739147277262D+00, .664380634782802D+00/

      DATA ABSD24/ .667468796667163D+00, .670878574004782D+00,
     A.676954241332878D+00, .685597412165933D+00, .696666335671016D+00,
     B.709979301403574D+00, .725317724677171D+00, .742429754768527D+00,
     C.761034415395659D+00, .780826220050639D+00, .801480188754397D+00,
     D.822657184522899D+00, .844009482143768D+00, .865186477912269D+00,
     E.885840446616027D+00, .905632251271007D+00, .924236911898140D+00,
     F.941348941989496D+00, .956687365263092D+00, .970000330995650D+00,
     G.981069254500733D+00, .989712425333789D+00, .995788092661885D+00,
     H.999197869999503D+00/


      EQUIVALENCE (ABSD08(1),RLEG02(1)), (ABSD14(1),RLEG02(L004+1)),
     $            (ABSD24(1),RLEG02(L004+L005+1))
      END


