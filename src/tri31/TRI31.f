CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C *********************************************************************** C
C *             USAGE STRICTEMENT RESERVE POUR STOP-95                  * C
C *                                                                     * C
C * PROGRAMMER : AHMED BOUFERGUENE            OCT 20 1994               * C
C *     LABORATOIRE DE CATALYSE ET SPECTROCHIMIE, CAEN, FRANCE          * C
C *                                                                     * C
C * THIS SUBROUTINE COMPUTES THE THREE-CENTER ONE-ELECTRON INTEGRALS    * C
C * INVOLVED IN THE SCF PROCEDURE. ITS IS CALLED  FROM THE INTGR        * C
C * ROUTINE.                                                            * C
C *              < XA | 1/RB | XC > <==> < XC | 1/RB | XA >             * C
C *                                                                     * C
C * INPUT :                                                             * C
C *   NEPSIL : IS THE ORDER OF THE EPSILON ALGORITHM                    * C
C *                                                                     * C
C *   LM123 : IS THE HIGHEST QUANTUM NUMBER L                           * C
C *                                                                     * C
C *   IC1, IC2, IOP : ARE THE RANKS OF THE ATOMS THE ORBITALS OF WHICH  * C
C *    ARE INTERACTING OF THE NUCLEUS IOP                               * C
C *                                                                     * C
C *   NUMAT : IS AN ARRAY CONTAINING THE ATOMIC NUMBERS                 * C
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
C *   XYZAT : IS AN ARRAY CONTAINING THE CARTESIAN COORDINATES OF THE   * C
C *    ATOMS INVOLVED IN THE MOLECULE                                   * C
C *                                                                     * C
C * OUTPUT :                                                            * C
C *   XCORE : IS A TWO-DIMENSIONAL ARRAY CONTAINING THE IMAGINARY CORE  * C
C *    INTEGRALS
C *********************************************************************** C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TRI31(NEPSIL,LM123,IC1,IOP,IC2,NUMAT,NORBAT,NLMAT,
     $                                             ZETAAT, XYZAT, XCORE)
      IMPLICIT REAL*8 (A-H, O-Z)
      INTEGER ZETST1, ZETST2
      INCLUDE "../lib95/SIZE.INCL"
      INCLUDE "../lib95/CONST.INCL"

      DIMENSION XCORE(N_ORB, N_ORB)
      DIMENSION NUMAT(*),NORBAT(*),NLMAT(*),ZETAAT(*),XYZAT(*)
      DIMENSION XYZV(3)
      COMMON/FACT0/FACT(0:30)
cc      COMMON/CNP0/CNP(0:30,0:30),CNPB(0:901)
cc      COMMON/DFACT/DFACT0(0:20)

cc      dimension FACT0(0:30), Dfact0(0:30), Cnp0(0:(30+1)*(30+2)/2)
      dimension Dfact0(0:100),Cnp0(0:(30+1)*(30+2)/2)
cc      dimension Cnp0(0:(30+1)*(30+2)/2)
      dimension xh(256), xbar(256)

c     WRITE(*, 1)IC1, IOP, IC2
c1    FORMAT('TRI31 :', 3(I2, 2X))

C.... COMPUTATION OF THE COORDINATES OF THE CENTER LABELED IOP TO WHICH THE
C.... COULOMB OPERATOR 1/RB IS ATTACHED.

      NIC1 = 3*(IC1-1)+1
      NIOP = 3*(IOP-1)+1
      NIC2 = 3*(IC2-1)+1

      CALL CARTCO(XYZAT(NIC1), XYZAT(NIOP), XYZV, ac, thetac, phiac)

C.... COMPUTATION OF THE COORDINATES OF THE CENTER LABELED IC2 TO WHICH THE
C.... ATOMIC ORBITAL XC IS ATTACHED

      CALL CARTCO(XYZAT(NIC1), XYZAT(NIC2), XYZV, ab, thetab, phiab)

C.... DETERMINING THE STARTING INDEX FOR THE ORBITALS AND THE SLATER EXPONENTS

      CALL NDXSTR(NORBAT, IC1, NLMST1, ZETST1)
      CALL NDXSTR(NORBAT, IC2, NLMST2, ZETST2)

C.... COMPUTATION OF THE THREE-CENTER ONE-ELECTRON INTEGRALS
c      CALL TRIMON(NEPSIL,LM123, NUMAT(IOP),NORBAT(IC1),NORBAT(IC2),
c     $     NLMAT(NLMST1),NLMAT(NLMST2), ZETAAT(ZETST1),ZETAAT(ZETST2),
c     $                     RVMOD,THETRV,PHIRV, RVCL,THETCL,PHICL, XCORE)

      
      norb1 = NORBAT(IC1)
      norb2 = NORBAT(IC2)

      nrac = 20
      nsd = 6
      jc = 6

c      call factoriel(30, Fact0)
      call dfactoriel(101, Dfact0)
      call combinatoire(30, Cnp0)
      call glncoefs(nrac, xbar, xh, ifail)
      call trimon(lm123, NUMAT(IOP), norb1, norb2, 
     $     NLMAT(NLMST1),NLMAT(NLMST2),ZETAAT(ZETST1),ZETAAT(ZETST2),
     $     ab, thetab, phiab, ac, thetac, phiac, 
     $     Fact,Dfact0,Cnp0, xh,xbar, nrac,nsd,jc, Xcore)
   
 
      RETURN
      END

c**********************************************************************c

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction factoriel
C
C     Programmeur Lilian Berlu 17/07/2002
C     
C  cette fonction calcul n! pour n=0 à n=nmax.
C  le résultat est rangé dans le tableau fact tq: fact[n] = n!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine factoriel(nmax, Fact)
      implicit double precision (a-h, o-z)
      implicit integer (i-n)

      dimension Fact(0:*) 
      
      if(nmax.LT.0) then
         write(6,*) 'erreur dans la fonction factoriel: nmax < 0'
         stop
      endif

      Fact(0) = 1.0d0
      do i = 1, nmax
         Fact(i) = dble(i) * Fact(i-1)
      end do

      return
      end

c**********************************************************************c

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction dfactoriel
C
C     Programmeur Lilian Berlu 17/07/2002
C     
C  fonction double_factoriel
C  cette fonction calcul n!! pour n=0 à n=nmax.
C  le résultat est rangé dans le tableau dfact tq: dfact[n] = n!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dfactoriel(imax, Dfact)
      implicit double precision (a-h, o-z)
      implicit integer (i-n)

      dimension Dfact(0:*)
      
      if(imax.LT.0) then 
         write(6,*) 'erreur dans la fonction dfactoriel: nmax<0'
         stop
      endif

      Dfact(0) = 1.0d0

      if(imax.GT.0) then
         Dfact(1) = 1.0d0
        
         do i=2, imax
            Dfact(i) = dble(i) * Dfact(i-2)
         end do
      endif

      return
      end

c**********************************************************************c

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction combinatoire
C
C     Programmeur Lilian Berlu 17/07/2002
C     
C     calcul des cnp basé sur la relation (n+1; p+1) = (n; p) + (n; p+1)
C     le résultats est rangé dans le tableau 
C     cnp (0 : nu * (nu + 1)/2 + nu)
C     tq :      (n ; p) = cnp [ n * (n + 1)/2 + p]
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine combinatoire(nu, Cnp)
      implicit double precision (a-h, o-z)
      implicit integer (i-n)

      dimension Cnp(0:*)

      do n=0, nu 
         do np=0, n
            if((np.EQ.0).OR.(np.EQ.n)) then
               Cnp(n*(n+1)/2+np) = 1.0d0
            else
               Cnp(n*(n+1)/2+np)=Cnp(n*(n-1)/2+np-1) + Cnp(n*(n-1)/2+np)
            endif
         end do
      end do
      return
      end

c**********************************************************************c

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C                         Fonction glncoefs
c
c     Programmeur Lilian Berlu 02/12/2002
c
c     n : nombre de points de quadrature
c     RacLeg, Weight :  racines et poids pour la quadrature de
c     Gauss_Legendre.
c
c     Les valeurs possible pour n sont:
c     1,2,3,4,5,6,8,10,12,14,16,20,24, 32, 40, 48, 64, 96, 256
c
c     Si n prend une valeur non permise alors ifail est 
c     initialisé à 1, sinon ifail == 0.
c
c     -----------------------------------------
c     Pour évaluer \int_{a}{b} f(x), 
c     Le changement de variables s'effectue comme suit:
c     
c      do i = 1, n
c         RacLeg(i) = 0.5 * ( (b-a)*RacLeg(i) + (b+a) )
c         Weight(i) = 0.5 * (b-a) * Weight(i)
c      end do
c
c     L'intégrale vaut:
c     \sum_{i=1}^{n} Weight(i) * f(RacLag(i))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine glncoefs(n, RacLeg, Weight, ifail)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      dimension RacLeg(256), Weight(256)

      ifail = 0
 
      if(n.EQ.1) then
c     Racines et Poids de Gauss-Legendre pour N = 1
c     i       ième racine            ième poids
c     Racines et Poids de Gauss-Legendre pour N = 1
c     i       ième racine            ième poids
      RacLeg(  1) =  .00000000000000D+00
      Weight(  1) =  .20000000000000D+01
      
      else if(n.EQ.2) then
      
c     Racines et Poids de Gauss-Legendre pour N = 2
c     i       ième racine            ième poids
      RacLeg(  1) = -.57735026918963D+00
      Weight(  1) =  .10000000000000D+01
      RacLeg(  2) =  .57735026918963D+00
      Weight(  2) =  .10000000000000D+01
      
      else if(n.EQ.3) then
      
c     Racines et Poids de Gauss-Legendre pour N = 3
c     i       ième racine            ième poids
      RacLeg(  1) = -.77459666924148D+00
      Weight(  1) =  .55555555555556D+00
      RacLeg(  2) =  .00000000000000D+00
      Weight(  2) =  .88888888888889D+00
      RacLeg(  3) =  .77459666924148D+00
      Weight(  3) =  .55555555555556D+00
      
      else if(n.EQ.4) then
      
c     Racines et Poids de Gauss-Legendre pour N = 4
c     i       ième racine            ième poids
      RacLeg(  1) = -.86113631159405D+00
      Weight(  1) =  .34785484513745D+00
      RacLeg(  2) = -.33998104358486D+00
      Weight(  2) =  .65214515486255D+00
      RacLeg(  3) =  .33998104358486D+00
      Weight(  3) =  .65214515486255D+00
      RacLeg(  4) =  .86113631159405D+00
      Weight(  4) =  .34785484513745D+00
      
      else if(n.EQ.5) then
      
c     Racines et Poids de Gauss-Legendre pour N = 5
c     i       ième racine            ième poids
      RacLeg(  1) = -.90617984593866D+00
      Weight(  1) =  .23692688505619D+00
      RacLeg(  2) = -.53846931010568D+00
      Weight(  2) =  .47862867049937D+00
      RacLeg(  3) =  .00000000000000D+00
      Weight(  3) =  .56888888888889D+00
      RacLeg(  4) =  .53846931010568D+00
      Weight(  4) =  .47862867049937D+00
      RacLeg(  5) =  .90617984593866D+00
      Weight(  5) =  .23692688505619D+00
      
      else if(n.EQ.6) then
      
c     Racines et Poids de Gauss-Legendre pour N = 6
c     i       ième racine            ième poids
      RacLeg(  1) = -.93246951420315D+00
      Weight(  1) =  .17132449237917D+00
      RacLeg(  2) = -.66120938646626D+00
      Weight(  2) =  .36076157304814D+00
      RacLeg(  3) = -.23861918608320D+00
      Weight(  3) =  .46791393457269D+00
      RacLeg(  4) =  .23861918608320D+00
      Weight(  4) =  .46791393457269D+00
      RacLeg(  5) =  .66120938646626D+00
      Weight(  5) =  .36076157304814D+00
      RacLeg(  6) =  .93246951420315D+00
      Weight(  6) =  .17132449237917D+00
      
      else if(n.EQ.8) then
      
c     Racines et Poids de Gauss-Legendre pour N = 8
c     i       ième racine            ième poids
      RacLeg(  1) = -.96028985649754D+00
      Weight(  1) =  .10122853629038D+00
      RacLeg(  2) = -.79666647741363D+00
      Weight(  2) =  .22238103445337D+00
      RacLeg(  3) = -.52553240991633D+00
      Weight(  3) =  .31370664587789D+00
      RacLeg(  4) = -.18343464249565D+00
      Weight(  4) =  .36268378337836D+00
      RacLeg(  5) =  .18343464249565D+00
      Weight(  5) =  .36268378337836D+00
      RacLeg(  6) =  .52553240991633D+00
      Weight(  6) =  .31370664587789D+00
      RacLeg(  7) =  .79666647741363D+00
      Weight(  7) =  .22238103445337D+00
      RacLeg(  8) =  .96028985649754D+00
      Weight(  8) =  .10122853629038D+00
      
      else if(n.EQ.10) then
      
c     Racines et Poids de Gauss-Legendre pour N = 10
c     i       ième racine            ième poids
      RacLeg(  1) = -.97390652851717D+00
      Weight(  1) =  .66671344308688D-01
      RacLeg(  2) = -.86506336668898D+00
      Weight(  2) =  .14945134915058D+00
      RacLeg(  3) = -.67940956829902D+00
      Weight(  3) =  .21908636251598D+00
      RacLeg(  4) = -.43339539412925D+00
      Weight(  4) =  .26926671931000D+00
      RacLeg(  5) = -.14887433898163D+00
      Weight(  5) =  .29552422471475D+00
      RacLeg(  6) =  .14887433898163D+00
      Weight(  6) =  .29552422471475D+00
      RacLeg(  7) =  .43339539412925D+00
      Weight(  7) =  .26926671931000D+00
      RacLeg(  8) =  .67940956829902D+00
      Weight(  8) =  .21908636251598D+00
      RacLeg(  9) =  .86506336668898D+00
      Weight(  9) =  .14945134915058D+00
      RacLeg( 10) =  .97390652851717D+00
      Weight( 10) =  .66671344308688D-01
      
      else if(n.EQ.12) then
      
c     Racines et Poids de Gauss-Legendre pour N = 12
c     i       ième racine            ième poids
      RacLeg(  1) = -.98156063424672D+00
      Weight(  1) =  .47175336386512D-01
      RacLeg(  2) = -.90411725637047D+00
      Weight(  2) =  .10693932599532D+00
      RacLeg(  3) = -.76990267419430D+00
      Weight(  3) =  .16007832854335D+00
      RacLeg(  4) = -.58731795428662D+00
      Weight(  4) =  .20316742672307D+00
      RacLeg(  5) = -.36783149899818D+00
      Weight(  5) =  .23349253653835D+00
      RacLeg(  6) = -.12523340851147D+00
      Weight(  6) =  .24914704581340D+00
      RacLeg(  7) =  .12523340851147D+00
      Weight(  7) =  .24914704581340D+00
      RacLeg(  8) =  .36783149899818D+00
      Weight(  8) =  .23349253653835D+00
      RacLeg(  9) =  .58731795428662D+00
      Weight(  9) =  .20316742672307D+00
      RacLeg( 10) =  .76990267419430D+00
      Weight( 10) =  .16007832854335D+00
      RacLeg( 11) =  .90411725637047D+00
      Weight( 11) =  .10693932599532D+00
      RacLeg( 12) =  .98156063424672D+00
      Weight( 12) =  .47175336386512D-01
      
      else if(n.EQ.14) then
      
c     Racines et Poids de Gauss-Legendre pour N = 14
c     i       ième racine            ième poids
      RacLeg(  1) = -.98628380869681D+00
      Weight(  1) =  .35119460331752D-01
      RacLeg(  2) = -.92843488366357D+00
      Weight(  2) =  .80158087159760D-01
      RacLeg(  3) = -.82720131506977D+00
      Weight(  3) =  .12151857068790D+00
      RacLeg(  4) = -.68729290481169D+00
      Weight(  4) =  .15720316715819D+00
      RacLeg(  5) = -.51524863635815D+00
      Weight(  5) =  .18553839747794D+00
      RacLeg(  6) = -.31911236892789D+00
      Weight(  6) =  .20519846372130D+00
      RacLeg(  7) = -.10805494870734D+00
      Weight(  7) =  .21526385346316D+00
      RacLeg(  8) =  .10805494870734D+00
      Weight(  8) =  .21526385346316D+00
      RacLeg(  9) =  .31911236892789D+00
      Weight(  9) =  .20519846372130D+00
      RacLeg( 10) =  .51524863635815D+00
      Weight( 10) =  .18553839747794D+00
      RacLeg( 11) =  .68729290481169D+00
      Weight( 11) =  .15720316715819D+00
      RacLeg( 12) =  .82720131506977D+00
      Weight( 12) =  .12151857068790D+00
      RacLeg( 13) =  .92843488366357D+00
      Weight( 13) =  .80158087159760D-01
      RacLeg( 14) =  .98628380869681D+00
      Weight( 14) =  .35119460331752D-01
      
      else if(n.EQ.16) then
      
c     Racines et Poids de Gauss-Legendre pour N = 16
c     i       ième racine            ième poids
      RacLeg(  1) = -.98940093499165D+00
      Weight(  1) =  .27152459411754D-01
      RacLeg(  2) = -.94457502307323D+00
      Weight(  2) =  .62253523938648D-01
      RacLeg(  3) = -.86563120238783D+00
      Weight(  3) =  .95158511682493D-01
      RacLeg(  4) = -.75540440835500D+00
      Weight(  4) =  .12462897125553D+00
      RacLeg(  5) = -.61787624440264D+00
      Weight(  5) =  .14959598881658D+00
      RacLeg(  6) = -.45801677765723D+00
      Weight(  6) =  .16915651939500D+00
      RacLeg(  7) = -.28160355077926D+00
      Weight(  7) =  .18260341504492D+00
      RacLeg(  8) = -.95012509837637D-01
      Weight(  8) =  .18945061045507D+00
      RacLeg(  9) =  .95012509837637D-01
      Weight(  9) =  .18945061045507D+00
      RacLeg( 10) =  .28160355077926D+00
      Weight( 10) =  .18260341504492D+00
      RacLeg( 11) =  .45801677765723D+00
      Weight( 11) =  .16915651939500D+00
      RacLeg( 12) =  .61787624440264D+00
      Weight( 12) =  .14959598881658D+00
      RacLeg( 13) =  .75540440835500D+00
      Weight( 13) =  .12462897125553D+00
      RacLeg( 14) =  .86563120238783D+00
      Weight( 14) =  .95158511682493D-01
      RacLeg( 15) =  .94457502307323D+00
      Weight( 15) =  .62253523938648D-01
      RacLeg( 16) =  .98940093499165D+00
      Weight( 16) =  .27152459411754D-01
      
      else if(n.EQ.20) then
      
c     Racines et Poids de Gauss-Legendre pour N = 20
c     i       ième racine            ième poids
      RacLeg(  1) = -.99312859918509D+00
      Weight(  1) =  .17614007139152D-01
      RacLeg(  2) = -.96397192727791D+00
      Weight(  2) =  .40601429800387D-01
      RacLeg(  3) = -.91223442825133D+00
      Weight(  3) =  .62672048334109D-01
      RacLeg(  4) = -.83911697182222D+00
      Weight(  4) =  .83276741576705D-01
      RacLeg(  5) = -.74633190646015D+00
      Weight(  5) =  .10193011981724D+00
      RacLeg(  6) = -.63605368072652D+00
      Weight(  6) =  .11819453196152D+00
      RacLeg(  7) = -.51086700195083D+00
      Weight(  7) =  .13168863844918D+00
      RacLeg(  8) = -.37370608871542D+00
      Weight(  8) =  .14209610931838D+00
      RacLeg(  9) = -.22778585114165D+00
      Weight(  9) =  .14917298647260D+00
      RacLeg( 10) = -.76526521133497D-01
      Weight( 10) =  .15275338713073D+00
      RacLeg( 11) =  .76526521133497D-01
      Weight( 11) =  .15275338713073D+00
      RacLeg( 12) =  .22778585114165D+00
      Weight( 12) =  .14917298647260D+00
      RacLeg( 13) =  .37370608871542D+00
      Weight( 13) =  .14209610931838D+00
      RacLeg( 14) =  .51086700195083D+00
      Weight( 14) =  .13168863844918D+00
      RacLeg( 15) =  .63605368072652D+00
      Weight( 15) =  .11819453196152D+00
      RacLeg( 16) =  .74633190646015D+00
      Weight( 16) =  .10193011981724D+00
      RacLeg( 17) =  .83911697182222D+00
      Weight( 17) =  .83276741576705D-01
      RacLeg( 18) =  .91223442825133D+00
      Weight( 18) =  .62672048334109D-01
      RacLeg( 19) =  .96397192727791D+00
      Weight( 19) =  .40601429800387D-01
      RacLeg( 20) =  .99312859918509D+00
      Weight( 20) =  .17614007139152D-01
      
      else if(n.EQ.24) then
      
c     Racines et Poids de Gauss-Legendre pour N = 24
c     i       ième racine            ième poids
      RacLeg(  1) = -.99518721999702D+00
      Weight(  1) =  .12341229799987D-01
      RacLeg(  2) = -.97472855597131D+00
      Weight(  2) =  .28531388628934D-01
      RacLeg(  3) = -.93827455200273D+00
      Weight(  3) =  .44277438817420D-01
      RacLeg(  4) = -.88641552700440D+00
      Weight(  4) =  .59298584915437D-01
      RacLeg(  5) = -.82000198597390D+00
      Weight(  5) =  .73346481411080D-01
      RacLeg(  6) = -.74012419157855D+00
      Weight(  6) =  .86190161531953D-01
      RacLeg(  7) = -.64809365193698D+00
      Weight(  7) =  .97618652104114D-01
      RacLeg(  8) = -.54542147138884D+00
      Weight(  8) =  .10744427011597D+00
      RacLeg(  9) = -.43379350762605D+00
      Weight(  9) =  .11550566805373D+00
      RacLeg( 10) = -.31504267969616D+00
      Weight( 10) =  .12167047292780D+00
      RacLeg( 11) = -.19111886747362D+00
      Weight( 11) =  .12583745634683D+00
      RacLeg( 12) = -.64056892862606D-01
      Weight( 12) =  .12793819534675D+00
      RacLeg( 13) =  .64056892862606D-01
      Weight( 13) =  .12793819534675D+00
      RacLeg( 14) =  .19111886747362D+00
      Weight( 14) =  .12583745634683D+00
      RacLeg( 15) =  .31504267969616D+00
      Weight( 15) =  .12167047292780D+00
      RacLeg( 16) =  .43379350762605D+00
      Weight( 16) =  .11550566805373D+00
      RacLeg( 17) =  .54542147138884D+00
      Weight( 17) =  .10744427011597D+00
      RacLeg( 18) =  .64809365193698D+00
      Weight( 18) =  .97618652104114D-01
      RacLeg( 19) =  .74012419157855D+00
      Weight( 19) =  .86190161531953D-01
      RacLeg( 20) =  .82000198597390D+00
      Weight( 20) =  .73346481411080D-01
      RacLeg( 21) =  .88641552700440D+00
      Weight( 21) =  .59298584915437D-01
      RacLeg( 22) =  .93827455200273D+00
      Weight( 22) =  .44277438817420D-01
      RacLeg( 23) =  .97472855597131D+00
      Weight( 23) =  .28531388628934D-01
      RacLeg( 24) =  .99518721999702D+00
      Weight( 24) =  .12341229799987D-01
      
      else if(n.EQ.32) then
      
c     Racines et Poids de Gauss-Legendre pour N = 32
c     i       ième racine            ième poids
      RacLeg(  1) = -.99726386184948D+00
      Weight(  1) =  .70186100094701D-02
      RacLeg(  2) = -.98561151154527D+00
      Weight(  2) =  .16274394730906D-01
      RacLeg(  3) = -.96476225558751D+00
      Weight(  3) =  .25392065309262D-01
      RacLeg(  4) = -.93490607593774D+00
      Weight(  4) =  .34273862913021D-01
      RacLeg(  5) = -.89632115576605D+00
      Weight(  5) =  .42835898022227D-01
      RacLeg(  6) = -.84936761373257D+00
      Weight(  6) =  .50998059262376D-01
      RacLeg(  7) = -.79448379596794D+00
      Weight(  7) =  .58684093478536D-01
      RacLeg(  8) = -.73218211874029D+00
      Weight(  8) =  .65822222776362D-01
      RacLeg(  9) = -.66304426693022D+00
      Weight(  9) =  .72345794108849D-01
      RacLeg( 10) = -.58771575724076D+00
      Weight( 10) =  .78193895787070D-01
      RacLeg( 11) = -.50689990893223D+00
      Weight( 11) =  .83311924226947D-01
      RacLeg( 12) = -.42135127613064D+00
      Weight( 12) =  .87652093004404D-01
      RacLeg( 13) = -.33186860228213D+00
      Weight( 13) =  .91173878695764D-01
      RacLeg( 14) = -.23928736225214D+00
      Weight( 14) =  .93844399080805D-01
      RacLeg( 15) = -.14447196158280D+00
      Weight( 15) =  .95638720079275D-01
      RacLeg( 16) = -.48307665687738D-01
      Weight( 16) =  .96540088514728D-01
      RacLeg( 17) =  .48307665687738D-01
      Weight( 17) =  .96540088514728D-01
      RacLeg( 18) =  .14447196158280D+00
      Weight( 18) =  .95638720079275D-01
      RacLeg( 19) =  .23928736225214D+00
      Weight( 19) =  .93844399080805D-01
      RacLeg( 20) =  .33186860228213D+00
      Weight( 20) =  .91173878695764D-01
      RacLeg( 21) =  .42135127613064D+00
      Weight( 21) =  .87652093004404D-01
      RacLeg( 22) =  .50689990893223D+00
      Weight( 22) =  .83311924226947D-01
      RacLeg( 23) =  .58771575724076D+00
      Weight( 23) =  .78193895787070D-01
      RacLeg( 24) =  .66304426693022D+00
      Weight( 24) =  .72345794108849D-01
      RacLeg( 25) =  .73218211874029D+00
      Weight( 25) =  .65822222776362D-01
      RacLeg( 26) =  .79448379596794D+00
      Weight( 26) =  .58684093478536D-01
      RacLeg( 27) =  .84936761373257D+00
      Weight( 27) =  .50998059262376D-01
      RacLeg( 28) =  .89632115576605D+00
      Weight( 28) =  .42835898022227D-01
      RacLeg( 29) =  .93490607593774D+00
      Weight( 29) =  .34273862913021D-01
      RacLeg( 30) =  .96476225558751D+00
      Weight( 30) =  .25392065309262D-01
      RacLeg( 31) =  .98561151154527D+00
      Weight( 31) =  .16274394730906D-01
      RacLeg( 32) =  .99726386184948D+00
      Weight( 32) =  .70186100094701D-02
      
      else if(n.EQ.40) then
 
c     Racines et Poids de Gauss-Legendre pour N = 40
c     i       ième racine            ième poids
      RacLeg(  1) = -.99823770971056D+00
      Weight(  1) =  .45212770985332D-02
      RacLeg(  2) = -.99072623869946D+00
      Weight(  2) =  .10498284531153D-01
      RacLeg(  3) = -.97725994998377D+00
      Weight(  3) =  .16421058381908D-01
      RacLeg(  4) = -.95791681921379D+00
      Weight(  4) =  .22245849194167D-01
      RacLeg(  5) = -.93281280827868D+00
      Weight(  5) =  .27937006980023D-01
      RacLeg(  6) = -.90209880696887D+00
      Weight(  6) =  .33460195282548D-01
      RacLeg(  7) = -.86595950321226D+00
      Weight(  7) =  .38782167974472D-01
      RacLeg(  8) = -.82461223083331D+00
      Weight(  8) =  .43870908185673D-01
      RacLeg(  9) = -.77830565142652D+00
      Weight(  9) =  .48695807635072D-01
      RacLeg( 10) = -.72731825518993D+00
      Weight( 10) =  .53227846983937D-01
      RacLeg( 11) = -.67195668461418D+00
      Weight( 11) =  .57439769099392D-01
      RacLeg( 12) = -.61255388966798D+00
      Weight( 12) =  .61306242492929D-01
      RacLeg( 13) = -.54946712509513D+00
      Weight( 13) =  .64804013456601D-01
      RacLeg( 14) = -.48307580168618D+00
      Weight( 14) =  .67912045815234D-01
      RacLeg( 15) = -.41377920437160D+00
      Weight( 15) =  .70611647391287D-01
      RacLeg( 16) = -.34199409082576D+00
      Weight( 16) =  .72886582395804D-01
      RacLeg( 17) = -.26815218500725D+00
      Weight( 17) =  .74723169057968D-01
      RacLeg( 18) = -.19269758070137D+00
      Weight( 18) =  .76110361900626D-01
      RacLeg( 19) = -.11608407067526D+00
      Weight( 19) =  .77039818164248D-01
      RacLeg( 20) = -.38772417506051D-01
      Weight( 20) =  .77505947978425D-01
      RacLeg( 21) =  .38772417506051D-01
      Weight( 21) =  .77505947978425D-01
      RacLeg( 22) =  .11608407067526D+00
      Weight( 22) =  .77039818164248D-01
      RacLeg( 23) =  .19269758070137D+00
      Weight( 23) =  .76110361900626D-01
      RacLeg( 24) =  .26815218500725D+00
      Weight( 24) =  .74723169057968D-01
      RacLeg( 25) =  .34199409082576D+00
      Weight( 25) =  .72886582395804D-01
      RacLeg( 26) =  .41377920437160D+00
      Weight( 26) =  .70611647391287D-01
      RacLeg( 27) =  .48307580168618D+00
      Weight( 27) =  .67912045815234D-01
      RacLeg( 28) =  .54946712509513D+00
      Weight( 28) =  .64804013456601D-01
      RacLeg( 29) =  .61255388966798D+00
      Weight( 29) =  .61306242492929D-01
      RacLeg( 30) =  .67195668461418D+00
      Weight( 30) =  .57439769099392D-01
      RacLeg( 31) =  .72731825518993D+00
      Weight( 31) =  .53227846983937D-01
      RacLeg( 32) =  .77830565142652D+00
      Weight( 32) =  .48695807635072D-01
      RacLeg( 33) =  .82461223083331D+00
      Weight( 33) =  .43870908185673D-01
      RacLeg( 34) =  .86595950321226D+00
      Weight( 34) =  .38782167974472D-01
      RacLeg( 35) =  .90209880696887D+00
      Weight( 35) =  .33460195282548D-01
      RacLeg( 36) =  .93281280827868D+00
      Weight( 36) =  .27937006980023D-01
      RacLeg( 37) =  .95791681921379D+00
      Weight( 37) =  .22245849194167D-01
      RacLeg( 38) =  .97725994998377D+00
      Weight( 38) =  .16421058381908D-01
      RacLeg( 39) =  .99072623869946D+00
      Weight( 39) =  .10498284531153D-01
      RacLeg( 40) =  .99823770971056D+00
      Weight( 40) =  .45212770985332D-02
 
      else if(n.EQ.48) then
 
c     Racines et Poids de Gauss-Legendre pour N = 48
c     i       ième racine            ième poids
      RacLeg(  1) = -.99877100725243D+00
      Weight(  1) =  .31533460523058D-02
      RacLeg(  2) = -.99353017226635D+00
      Weight(  2) =  .73275539012763D-02
      RacLeg(  3) = -.98412458372283D+00
      Weight(  3) =  .11477234579235D-01
      RacLeg(  4) = -.97059159254625D+00
      Weight(  4) =  .15579315722944D-01
      RacLeg(  5) = -.95298770316043D+00
      Weight(  5) =  .19616160457356D-01
      RacLeg(  6) = -.93138669070655D+00
      Weight(  6) =  .23570760839324D-01
      RacLeg(  7) = -.90587913671557D+00
      Weight(  7) =  .27426509708357D-01
      RacLeg(  8) = -.87657202027425D+00
      Weight(  8) =  .31167227832798D-01
      RacLeg(  9) = -.84358826162439D+00
      Weight(  9) =  .34777222564770D-01
      RacLeg( 10) = -.80706620402944D+00
      Weight( 10) =  .38241351065831D-01
      RacLeg( 11) = -.76715903251574D+00
      Weight( 11) =  .41545082943465D-01
      RacLeg( 12) = -.72403413092381D+00
      Weight( 12) =  .44674560856694D-01
      RacLeg( 13) = -.67787237963266D+00
      Weight( 13) =  .47616658492490D-01
      RacLeg( 14) = -.62886739677651D+00
      Weight( 14) =  .50359035553854D-01
      RacLeg( 15) = -.57722472608397D+00
      Weight( 15) =  .52890189485194D-01
      RacLeg( 16) = -.52316097472223D+00
      Weight( 16) =  .55199503699984D-01
      RacLeg( 17) = -.46690290475096D+00
      Weight( 17) =  .57277292100403D-01
      RacLeg( 18) = -.40868648199072D+00
      Weight( 18) =  .59114839698396D-01
      RacLeg( 19) = -.34875588629216D+00
      Weight( 19) =  .60704439165894D-01
      RacLeg( 20) = -.28736248735546D+00
      Weight( 20) =  .62039423159893D-01
      RacLeg( 21) = -.22476379039469D+00
      Weight( 21) =  .63114192286254D-01
      RacLeg( 22) = -.16122235606889D+00
      Weight( 22) =  .63924238584648D-01
      RacLeg( 23) = -.97004699209463D-01
      Weight( 23) =  .64466164435950D-01
      RacLeg( 24) = -.32380170962869D-01
      Weight( 24) =  .64737696812684D-01
      RacLeg( 25) =  .32380170962869D-01
      Weight( 25) =  .64737696812684D-01
      RacLeg( 26) =  .97004699209463D-01
      Weight( 26) =  .64466164435950D-01
      RacLeg( 27) =  .16122235606889D+00
      Weight( 27) =  .63924238584648D-01
      RacLeg( 28) =  .22476379039469D+00
      Weight( 28) =  .63114192286254D-01
      RacLeg( 29) =  .28736248735546D+00
      Weight( 29) =  .62039423159893D-01
      RacLeg( 30) =  .34875588629216D+00
      Weight( 30) =  .60704439165894D-01
      RacLeg( 31) =  .40868648199072D+00
      Weight( 31) =  .59114839698396D-01
      RacLeg( 32) =  .46690290475096D+00
      Weight( 32) =  .57277292100403D-01
      RacLeg( 33) =  .52316097472223D+00
      Weight( 33) =  .55199503699984D-01
      RacLeg( 34) =  .57722472608397D+00
      Weight( 34) =  .52890189485194D-01
      RacLeg( 35) =  .62886739677651D+00
      Weight( 35) =  .50359035553854D-01
      RacLeg( 36) =  .67787237963266D+00
      Weight( 36) =  .47616658492490D-01
      RacLeg( 37) =  .72403413092381D+00
      Weight( 37) =  .44674560856694D-01
      RacLeg( 38) =  .76715903251574D+00
      Weight( 38) =  .41545082943465D-01
      RacLeg( 39) =  .80706620402944D+00
      Weight( 39) =  .38241351065831D-01
      RacLeg( 40) =  .84358826162439D+00
      Weight( 40) =  .34777222564770D-01
      RacLeg( 41) =  .87657202027425D+00
      Weight( 41) =  .31167227832798D-01
      RacLeg( 42) =  .90587913671557D+00
      Weight( 42) =  .27426509708357D-01
      RacLeg( 43) =  .93138669070655D+00
      Weight( 43) =  .23570760839324D-01
      RacLeg( 44) =  .95298770316043D+00
      Weight( 44) =  .19616160457356D-01
      RacLeg( 45) =  .97059159254625D+00
      Weight( 45) =  .15579315722944D-01
      RacLeg( 46) =  .98412458372283D+00
      Weight( 46) =  .11477234579235D-01
      RacLeg( 47) =  .99353017226635D+00
      Weight( 47) =  .73275539012763D-02
      RacLeg( 48) =  .99877100725243D+00
      Weight( 48) =  .31533460523058D-02
 
      else if(n.EQ.64) then
 
c     Racines et Poids de Gauss-Legendre pour N = 64
c     i       ième racine            ième poids
      RacLeg(  1) = -.99930504173577D+00
      Weight(  1) =  .17832807216964D-02
      RacLeg(  2) = -.99634011677196D+00
      Weight(  2) =  .41470332605625D-02
      RacLeg(  3) = -.99101337147674D+00
      Weight(  3) =  .65044579689784D-02
      RacLeg(  4) = -.98333625388463D+00
      Weight(  4) =  .88467598263639D-02
      RacLeg(  5) = -.97332682778991D+00
      Weight(  5) =  .11168139460131D-01
      RacLeg(  6) = -.96100879965205D+00
      Weight(  6) =  .13463047896719D-01
      RacLeg(  7) = -.94641137485840D+00
      Weight(  7) =  .15726030476025D-01
      RacLeg(  8) = -.92956917213194D+00
      Weight(  8) =  .17951715775697D-01
      RacLeg(  9) = -.91052213707850D+00
      Weight(  9) =  .20134823153530D-01
      RacLeg( 10) = -.88931544599511D+00
      Weight( 10) =  .22270173808383D-01
      RacLeg( 11) = -.86599939815409D+00
      Weight( 11) =  .24352702568711D-01
      RacLeg( 12) = -.84062929625258D+00
      Weight( 12) =  .26377469715055D-01
      RacLeg( 13) = -.81326531512280D+00
      Weight( 13) =  .28339672614259D-01
      RacLeg( 14) = -.78397235894334D+00
      Weight( 14) =  .30234657072402D-01
      RacLeg( 15) = -.75281990726053D+00
      Weight( 15) =  .32057928354852D-01
      RacLeg( 16) = -.71988185017161D+00
      Weight( 16) =  .33805161837142D-01
      RacLeg( 17) = -.68523631305423D+00
      Weight( 17) =  .35472213256882D-01
      RacLeg( 18) = -.64896547125466D+00
      Weight( 18) =  .37055128540240D-01
      RacLeg( 19) = -.61115535517239D+00
      Weight( 19) =  .38550153178616D-01
      RacLeg( 20) = -.57189564620263D+00
      Weight( 20) =  .39953741132720D-01
      RacLeg( 21) = -.53127946401989D+00
      Weight( 21) =  .41262563242624D-01
      RacLeg( 22) = -.48940314570705D+00
      Weight( 22) =  .42473515123654D-01
      RacLeg( 23) = -.44636601725346D+00
      Weight( 23) =  .43583724529323D-01
      RacLeg( 24) = -.40227015796399D+00
      Weight( 24) =  .44590558163757D-01
      RacLeg( 25) = -.35722015833767D+00
      Weight( 25) =  .45491627927418D-01
      RacLeg( 26) = -.31132287199021D+00
      Weight( 26) =  .46284796581314D-01
      RacLeg( 27) = -.26468716220877D+00
      Weight( 27) =  .46968182816210D-01
      RacLeg( 28) = -.21742364374001D+00
      Weight( 28) =  .47540165714830D-01
      RacLeg( 29) = -.16964442042399D+00
      Weight( 29) =  .47999388596458D-01
      RacLeg( 30) = -.12146281929612D+00
      Weight( 30) =  .48344762234803D-01
      RacLeg( 31) = -.72993121787799D-01
      Weight( 31) =  .48575467441503D-01
      RacLeg( 32) = -.24350292663424D-01
      Weight( 32) =  .48690957009140D-01
      RacLeg( 33) =  .24350292663424D-01
      Weight( 33) =  .48690957009140D-01
      RacLeg( 34) =  .72993121787799D-01
      Weight( 34) =  .48575467441503D-01
      RacLeg( 35) =  .12146281929612D+00
      Weight( 35) =  .48344762234803D-01
      RacLeg( 36) =  .16964442042399D+00
      Weight( 36) =  .47999388596458D-01
      RacLeg( 37) =  .21742364374001D+00
      Weight( 37) =  .47540165714830D-01
      RacLeg( 38) =  .26468716220877D+00
      Weight( 38) =  .46968182816210D-01
      RacLeg( 39) =  .31132287199021D+00
      Weight( 39) =  .46284796581314D-01
      RacLeg( 40) =  .35722015833767D+00
      Weight( 40) =  .45491627927418D-01
      RacLeg( 41) =  .40227015796399D+00
      Weight( 41) =  .44590558163757D-01
      RacLeg( 42) =  .44636601725346D+00
      Weight( 42) =  .43583724529323D-01
      RacLeg( 43) =  .48940314570705D+00
      Weight( 43) =  .42473515123654D-01
      RacLeg( 44) =  .53127946401989D+00
      Weight( 44) =  .41262563242624D-01
      RacLeg( 45) =  .57189564620263D+00
      Weight( 45) =  .39953741132720D-01
      RacLeg( 46) =  .61115535517239D+00
      Weight( 46) =  .38550153178616D-01
      RacLeg( 47) =  .64896547125466D+00
      Weight( 47) =  .37055128540240D-01
      RacLeg( 48) =  .68523631305423D+00
      Weight( 48) =  .35472213256882D-01
      RacLeg( 49) =  .71988185017161D+00
      Weight( 49) =  .33805161837142D-01
      RacLeg( 50) =  .75281990726053D+00
      Weight( 50) =  .32057928354852D-01
      RacLeg( 51) =  .78397235894334D+00
      Weight( 51) =  .30234657072402D-01
      RacLeg( 52) =  .81326531512280D+00
      Weight( 52) =  .28339672614259D-01
      RacLeg( 53) =  .84062929625258D+00
      Weight( 53) =  .26377469715055D-01
      RacLeg( 54) =  .86599939815409D+00
      Weight( 54) =  .24352702568711D-01
      RacLeg( 55) =  .88931544599511D+00
      Weight( 55) =  .22270173808383D-01
      RacLeg( 56) =  .91052213707850D+00
      Weight( 56) =  .20134823153530D-01
      RacLeg( 57) =  .92956917213194D+00
      Weight( 57) =  .17951715775697D-01
      RacLeg( 58) =  .94641137485840D+00
      Weight( 58) =  .15726030476025D-01
      RacLeg( 59) =  .96100879965205D+00
      Weight( 59) =  .13463047896719D-01
      RacLeg( 60) =  .97332682778991D+00
      Weight( 60) =  .11168139460131D-01
      RacLeg( 61) =  .98333625388463D+00
      Weight( 61) =  .88467598263639D-02
      RacLeg( 62) =  .99101337147674D+00
      Weight( 62) =  .65044579689784D-02
      RacLeg( 63) =  .99634011677196D+00
      Weight( 63) =  .41470332605625D-02
      RacLeg( 64) =  .99930504173577D+00
      Weight( 64) =  .17832807216964D-02
 
      else if(n.EQ.96) then
 
c     Racines et Poids de Gauss-Legendre pour N = 96
c     i       ième racine            ième poids
      RacLeg(  1) = -.99968950388323D+00
      Weight(  1) =  .79679206555201D-03
      RacLeg(  2) = -.99836437586318D+00
      Weight(  2) =  .18539607889469D-02
      RacLeg(  3) = -.99598184298721D+00
      Weight(  3) =  .29107318179349D-02
      RacLeg(  4) = -.99254390032376D+00
      Weight(  4) =  .39645543384447D-02
      RacLeg(  5) = -.98805412632962D+00
      Weight(  5) =  .50142027429275D-02
      RacLeg(  6) = -.98251726356301D+00
      Weight(  6) =  .60585455042360D-02
      RacLeg(  7) = -.97593917458514D+00
      Weight(  7) =  .70964707911539D-02
      RacLeg(  8) = -.96832682846326D+00
      Weight(  8) =  .81268769256988D-02
      RacLeg(  9) = -.95968829144874D+00
      Weight(  9) =  .91486712307834D-02
      RacLeg( 10) = -.95003271778444D+00
      Weight( 10) =  .10160770535008D-01
      RacLeg( 11) = -.93937033975276D+00
      Weight( 11) =  .11162102099838D-01
      RacLeg( 12) = -.92771245672231D+00
      Weight( 12) =  .12151604671088D-01
      RacLeg( 13) = -.91507142312090D+00
      Weight( 13) =  .13128229566962D-01
      RacLeg( 14) = -.90146063531585D+00
      Weight( 14) =  .14090941772315D-01
      RacLeg( 15) = -.88689451740242D+00
      Weight( 15) =  .15038721026995D-01
      RacLeg( 16) = -.87138850590930D+00
      Weight( 16) =  .15970562902562D-01
      RacLeg( 17) = -.85495903343460D+00
      Weight( 17) =  .16885479864245D-01
      RacLeg( 18) = -.83762351122819D+00
      Weight( 18) =  .17782502316045D-01
      RacLeg( 19) = -.81940031073793D+00
      Weight( 19) =  .18660679627411D-01
      RacLeg( 20) = -.80030874413914D+00
      Weight( 20) =  .19519081140145D-01
      RacLeg( 21) = -.78036904386743D+00
      Weight( 21) =  .20356797154333D-01
      RacLeg( 22) = -.75960234117665D+00
      Weight( 22) =  .21172939892191D-01
      RacLeg( 23) = -.73803064374440D+00
      Weight( 23) =  .21966644438744D-01
      RacLeg( 24) = -.71567681234897D+00
      Weight( 24) =  .22737069658329D-01
      RacLeg( 25) = -.69256453664217D+00
      Weight( 25) =  .23483399085926D-01
      RacLeg( 26) = -.66871831004392D+00
      Weight( 26) =  .24204841792365D-01
      RacLeg( 27) = -.64416340378497D+00
      Weight( 27) =  .24900633222484D-01
      RacLeg( 28) = -.61892584012547D+00
      Weight( 28) =  .25570036005349D-01
      RacLeg( 29) = -.59303236477757D+00
      Weight( 29) =  .26212340735672D-01
      RacLeg( 30) = -.56651041856140D+00
      Weight( 30) =  .26826866725592D-01
      RacLeg( 31) = -.53938810832436D+00
      Weight( 31) =  .27412962726029D-01
      RacLeg( 32) = -.51169417715467D+00
      Weight( 32) =  .27970007616848D-01
      RacLeg( 33) = -.48345797392060D+00
      Weight( 33) =  .28497411065085D-01
      RacLeg( 34) = -.45470942216774D+00
      Weight( 34) =  .28994614150555D-01
      RacLeg( 35) = -.42547898840730D+00
      Weight( 35) =  .29461089958168D-01
      RacLeg( 36) = -.39579764982891D+00
      Weight( 36) =  .29896344136328D-01
      RacLeg( 37) = -.36569686147231D+00
      Weight( 37) =  .30299915420828D-01
      RacLeg( 38) = -.33520852289263D+00
      Weight( 38) =  .30671376123669D-01
      RacLeg( 39) = -.30436494435450D+00
      Weight( 39) =  .31010332586314D-01
      RacLeg( 40) = -.27319881259105D+00
      Weight( 40) =  .31316425596861D-01
      RacLeg( 41) = -.24174315616384D+00
      Weight( 41) =  .31589330770727D-01
      RacLeg( 42) = -.21003131046057D+00
      Weight( 42) =  .31828758894411D-01
      RacLeg( 43) = -.17809688236762D+00
      Weight( 43) =  .32034456231993D-01
      RacLeg( 44) = -.14597371465490D+00
      Weight( 44) =  .32206204794030D-01
      RacLeg( 45) = -.11369585011067D+00
      Weight( 45) =  .32343822568576D-01
      RacLeg( 46) = -.81297495464426D-01
      Weight( 46) =  .32447163714064D-01
      RacLeg( 47) = -.48812985136050D-01
      Weight( 47) =  .32516118713869D-01
      RacLeg( 48) = -.16276744849603D-01
      Weight( 48) =  .32550614492363D-01
      RacLeg( 49) =  .16276744849603D-01
      Weight( 49) =  .32550614492363D-01
      RacLeg( 50) =  .48812985136050D-01
      Weight( 50) =  .32516118713869D-01
      RacLeg( 51) =  .81297495464426D-01
      Weight( 51) =  .32447163714064D-01
      RacLeg( 52) =  .11369585011067D+00
      Weight( 52) =  .32343822568576D-01
      RacLeg( 53) =  .14597371465490D+00
      Weight( 53) =  .32206204794030D-01
      RacLeg( 54) =  .17809688236762D+00
      Weight( 54) =  .32034456231993D-01
      RacLeg( 55) =  .21003131046057D+00
      Weight( 55) =  .31828758894411D-01
      RacLeg( 56) =  .24174315616384D+00
      Weight( 56) =  .31589330770727D-01
      RacLeg( 57) =  .27319881259105D+00
      Weight( 57) =  .31316425596861D-01
      RacLeg( 58) =  .30436494435450D+00
      Weight( 58) =  .31010332586314D-01
      RacLeg( 59) =  .33520852289263D+00
      Weight( 59) =  .30671376123669D-01
      RacLeg( 60) =  .36569686147231D+00
      Weight( 60) =  .30299915420828D-01
      RacLeg( 61) =  .39579764982891D+00
      Weight( 61) =  .29896344136328D-01
      RacLeg( 62) =  .42547898840730D+00
      Weight( 62) =  .29461089958168D-01
      RacLeg( 63) =  .45470942216774D+00
      Weight( 63) =  .28994614150555D-01
      RacLeg( 64) =  .48345797392060D+00
      Weight( 64) =  .28497411065085D-01
      RacLeg( 65) =  .51169417715467D+00
      Weight( 65) =  .27970007616848D-01
      RacLeg( 66) =  .53938810832436D+00
      Weight( 66) =  .27412962726029D-01
      RacLeg( 67) =  .56651041856140D+00
      Weight( 67) =  .26826866725592D-01
      RacLeg( 68) =  .59303236477757D+00
      Weight( 68) =  .26212340735672D-01
      RacLeg( 69) =  .61892584012547D+00
      Weight( 69) =  .25570036005349D-01
      RacLeg( 70) =  .64416340378497D+00
      Weight( 70) =  .24900633222484D-01
      RacLeg( 71) =  .66871831004392D+00
      Weight( 71) =  .24204841792365D-01
      RacLeg( 72) =  .69256453664217D+00
      Weight( 72) =  .23483399085926D-01
      RacLeg( 73) =  .71567681234897D+00
      Weight( 73) =  .22737069658329D-01
      RacLeg( 74) =  .73803064374440D+00
      Weight( 74) =  .21966644438744D-01
      RacLeg( 75) =  .75960234117665D+00
      Weight( 75) =  .21172939892191D-01
      RacLeg( 76) =  .78036904386743D+00
      Weight( 76) =  .20356797154333D-01
      RacLeg( 77) =  .80030874413914D+00
      Weight( 77) =  .19519081140145D-01
      RacLeg( 78) =  .81940031073793D+00
      Weight( 78) =  .18660679627411D-01
      RacLeg( 79) =  .83762351122819D+00
      Weight( 79) =  .17782502316045D-01
      RacLeg( 80) =  .85495903343460D+00
      Weight( 80) =  .16885479864245D-01
      RacLeg( 81) =  .87138850590930D+00
      Weight( 81) =  .15970562902562D-01
      RacLeg( 82) =  .88689451740242D+00
      Weight( 82) =  .15038721026995D-01
      RacLeg( 83) =  .90146063531585D+00
      Weight( 83) =  .14090941772315D-01
      RacLeg( 84) =  .91507142312090D+00
      Weight( 84) =  .13128229566962D-01
      RacLeg( 85) =  .92771245672231D+00
      Weight( 85) =  .12151604671088D-01
      RacLeg( 86) =  .93937033975276D+00
      Weight( 86) =  .11162102099838D-01
      RacLeg( 87) =  .95003271778444D+00
      Weight( 87) =  .10160770535008D-01
      RacLeg( 88) =  .95968829144874D+00
      Weight( 88) =  .91486712307834D-02
      RacLeg( 89) =  .96832682846326D+00
      Weight( 89) =  .81268769256988D-02
      RacLeg( 90) =  .97593917458514D+00
      Weight( 90) =  .70964707911539D-02
      RacLeg( 91) =  .98251726356301D+00
      Weight( 91) =  .60585455042360D-02
      RacLeg( 92) =  .98805412632962D+00
      Weight( 92) =  .50142027429275D-02
      RacLeg( 93) =  .99254390032376D+00
      Weight( 93) =  .39645543384447D-02
      RacLeg( 94) =  .99598184298721D+00
      Weight( 94) =  .29107318179349D-02
      RacLeg( 95) =  .99836437586318D+00
      Weight( 95) =  .18539607889469D-02
      RacLeg( 96) =  .99968950388323D+00
      Weight( 96) =  .79679206555201D-03
 
      else if(n.EQ.256) then
 
c     Racines et Poids de Gauss-Legendre pour N = 256
c     i       ième racine            ième poids
      RacLeg(  1) = -.99995605001899D+00
      Weight(  1) =  .11278901782227D-03
      RacLeg(  2) = -.99976843740926D+00
      Weight(  2) =  .26253494429645D-03
      RacLeg(  3) = -.99943093746626D+00
      Weight(  3) =  .41246325442618D-03
      RacLeg(  4) = -.99894352584341D+00
      Weight(  4) =  .56234895403141D-03
      RacLeg(  5) = -.99830626647301D+00
      Weight(  5) =  .71215416347332D-03
      RacLeg(  6) = -.99751925275672D+00
      Weight(  6) =  .86185370142009D-03
      RacLeg(  7) = -.99658260202338D+00
      Weight(  7) =  .10114243932084D-02
      RacLeg(  8) = -.99549645448110D+00
      Weight(  8) =  .11608435575677D-02
      RacLeg(  9) = -.99426097292241D+00
      Weight(  9) =  .13100886819025D-02
      RacLeg( 10) = -.99287634260882D+00
      Weight( 10) =  .14591373333107D-02
      RacLeg( 11) = -.99134277120758D+00
      Weight( 11) =  .16079671307493D-02
      RacLeg( 12) = -.98966048874507D+00
      Weight( 12) =  .17565557363307D-02
      RacLeg( 13) = -.98782974756486D+00
      Weight( 13) =  .19048808534997D-02
      RacLeg( 14) = -.98585082228613D+00
      Weight( 14) =  .20529202279661D-02
      RacLeg( 15) = -.98372400976032D+00
      Weight( 15) =  .22006516498399D-02
      RacLeg( 16) = -.98144962902546D+00
      Weight( 16) =  .23480529563273D-02
      RacLeg( 17) = -.97902802125762D+00
      Weight( 17) =  .24951020347037D-02
      RacLeg( 18) = -.97645954971923D+00
      Weight( 18) =  .26417768254275D-02
      RacLeg( 19) = -.97374459970437D+00
      Weight( 19) =  .27880553253277D-02
      RacLeg( 20) = -.97088357848074D+00
      Weight( 20) =  .29339155908297D-02
      RacLeg( 21) = -.96787691522849D+00
      Weight( 21) =  .30793357411993D-02
      RacLeg( 22) = -.96472506097571D+00
      Weight( 22) =  .32242939617942D-02
      RacLeg( 23) = -.96142848853073D+00
      Weight( 23) =  .33687685073156D-02
      RacLeg( 24) = -.95798769241118D+00
      Weight( 24) =  .35127377050563D-02
      RacLeg( 25) = -.95440318876972D+00
      Weight( 25) =  .36561799581425D-02
      RacLeg( 26) = -.95067551531663D+00
      Weight( 26) =  .37990737487663D-02
      RacLeg( 27) = -.94680523123913D+00
      Weight( 27) =  .39413976414088D-02
      RacLeg( 28) = -.94279291711746D+00
      Weight( 28) =  .40831302860527D-02
      RacLeg( 29) = -.93863917483781D+00
      Weight( 29) =  .42242504213815D-02
      RacLeg( 30) = -.93434462750200D+00
      Weight( 30) =  .43647368779681D-02
      RacLeg( 31) = -.92990991933401D+00
      Weight( 31) =  .45045685814479D-02
      RacLeg( 32) = -.92533571558332D+00
      Weight( 32) =  .46437245556801D-02
      RacLeg( 33) = -.92062270242515D+00
      Weight( 33) =  .47821839258927D-02
      RacLeg( 34) = -.91577158685749D+00
      Weight( 34) =  .49199259218139D-02
      RacLeg( 35) = -.91078309659507D+00
      Weight( 35) =  .50569298807868D-02
      RacLeg( 36) = -.90565797996014D+00
      Weight( 36) =  .51931752508693D-02
      RacLeg( 37) = -.90039700577030D+00
      Weight( 37) =  .53286415939159D-02
      RacLeg( 38) = -.89500096322308D+00
      Weight( 38) =  .54633085886443D-02
      RacLeg( 39) = -.88947066177761D+00
      Weight( 39) =  .55971560336829D-02
      RacLeg( 40) = -.88380693103316D+00
      Weight( 40) =  .57301638506014D-02
      RacLeg( 41) = -.87801062060471D+00
      Weight( 41) =  .58623120869227D-02
      RacLeg( 42) = -.87208259999549D+00
      Weight( 42) =  .59935809191153D-02
      RacLeg( 43) = -.86602375846655D+00
      Weight( 43) =  .61239506555679D-02
      RacLeg( 44) = -.85983500490338D+00
      Weight( 44) =  .62534017395424D-02
      RacLeg( 45) = -.85351726767950D+00
      Weight( 45) =  .63819147521079D-02
      RacLeg( 46) = -.84707149451730D+00
      Weight( 46) =  .65094704150537D-02
      RacLeg( 47) = -.84049865234576D+00
      Weight( 47) =  .66360495937811D-02
      RacLeg( 48) = -.83379972715550D+00
      Weight( 48) =  .67616333001738D-02
      RacLeg( 49) = -.82697572385081D+00
      Weight( 49) =  .68862026954463D-02
      RacLeg( 50) = -.82002766609892D+00
      Weight( 50) =  .70097390929698D-02
      RacLeg( 51) = -.81295659617643D+00
      Weight( 51) =  .71322239610754D-02
      RacLeg( 52) = -.80576357481300D+00
      Weight( 52) =  .72536389258339D-02
      RacLeg( 53) = -.79844968103217D+00
      Weight( 53) =  .73739657738123D-02
      RacLeg( 54) = -.79101601198955D+00
      Weight( 54) =  .74931864548059D-02
      RacLeg( 55) = -.78346368280818D+00
      Weight( 55) =  .76112830845457D-02
      RacLeg( 56) = -.77579382641133D+00
      Weight( 56) =  .77282379473816D-02
      RacLeg( 57) = -.76800759335245D+00
      Weight( 57) =  .78440334989397D-02
      RacLeg( 58) = -.76010615164266D+00
      Weight( 58) =  .79586523687543D-02
      RacLeg( 59) = -.75209068657549D+00
      Weight( 59) =  .80720773628735D-02
      RacLeg( 60) = -.74396240054911D+00
      Weight( 60) =  .81842914664383D-02
      RacLeg( 61) = -.73572251288592D+00
      Weight( 61) =  .82952778462352D-02
      RacLeg( 62) = -.72737225964965D+00
      Weight( 62) =  .84050198532215D-02
      RacLeg( 63) = -.71891289345997D+00
      Weight( 63) =  .85135010250225D-02
      RacLeg( 64) = -.71034568330454D+00
      Weight( 64) =  .86207050884010D-02
      RacLeg( 65) = -.70167191434869D+00
      Weight( 65) =  .87266159616988D-02
      RacLeg( 66) = -.69289288774258D+00
      Weight( 66) =  .88312177572488D-02
      RacLeg( 67) = -.68400992042608D+00
      Weight( 67) =  .89344947837582D-02
      RacLeg( 68) = -.67502434493116D+00
      Weight( 68) =  .90364315486629D-02
      RacLeg( 69) = -.66593750918205D+00
      Weight( 69) =  .91370127604508D-02
      RacLeg( 70) = -.65675077629297D+00
      Weight( 70) =  .92362233309563D-02
      RacLeg( 71) = -.64746552436372D+00
      Weight( 71) =  .93340483776233D-02
      RacLeg( 72) = -.63808314627291D+00
      Weight( 72) =  .94304732257378D-02
      RacLeg( 73) = -.62860504946901D+00
      Weight( 73) =  .95254834106293D-02
      RacLeg( 74) = -.61903265575926D+00
      Weight( 74) =  .96190646798407D-02
      RacLeg( 75) = -.60936740109633D+00
      Weight( 75) =  .97112029952663D-02
      RacLeg( 76) = -.59961073536297D+00
      Weight( 76) =  .98018845352573D-02
      RacLeg( 77) = -.58976412215445D+00
      Weight( 77) =  .98910956966958D-02
      RacLeg( 78) = -.57982903855908D+00
      Weight( 78) =  .99788230970349D-02
      RacLeg( 79) = -.56980697493657D+00
      Weight( 79) =  .10065053576306D-01
      RacLeg( 80) = -.55969943469448D+00
      Weight( 80) =  .10149774199095D-01
      RacLeg( 81) = -.54950793406272D+00
      Weight( 81) =  .10232972256478D-01
      RacLeg( 82) = -.53923400186606D+00
      Weight( 82) =  .10314635267934D-01
      RacLeg( 83) = -.52887917929482D+00
      Weight( 83) =  .10394750983212D-01
      RacLeg( 84) = -.51844501967367D+00
      Weight( 84) =  .10473307384170D-01
      RacLeg( 85) = -.50793308822862D+00
      Weight( 85) =  .10550292686581D-01
      RacLeg( 86) = -.49734496185218D+00
      Weight( 86) =  .10625695341897D-01
      RacLeg( 87) = -.48668222886689D+00
      Weight( 87) =  .10699504038980D-01
      RacLeg( 88) = -.47594648878698D+00
      Weight( 88) =  .10771707705805D-01
      RacLeg( 89) = -.46513935207848D+00
      Weight( 89) =  .10842295511115D-01
      RacLeg( 90) = -.45426243991759D+00
      Weight( 90) =  .10911256866049D-01
      RacLeg( 91) = -.44331738394753D+00
      Weight( 91) =  .10978581425730D-01
      RacLeg( 92) = -.43230582603374D+00
      Weight( 92) =  .11044259090814D-01
      RacLeg( 93) = -.42122941801762D+00
      Weight( 93) =  .11108280009010D-01
      RacLeg( 94) = -.41008982146872D+00
      Weight( 94) =  .11170634576553D-01
      RacLeg( 95) = -.39888870743546D+00
      Weight( 95) =  .11231313439650D-01
      RacLeg( 96) = -.38762775619452D+00
      Weight( 96) =  .11290307495876D-01
      RacLeg( 97) = -.37630865699872D+00
      Weight( 97) =  .11347607895545D-01
      RacLeg( 98) = -.36493310782365D+00
      Weight( 98) =  .11403206043039D-01
      RacLeg( 99) = -.35350281511297D+00
      Weight( 99) =  .11457093598091D-01
      RacLeg(100) = -.34201949352237D+00
      Weight(100) =  .11509262477039D-01
      RacLeg(101) = -.33048486566242D+00
      Weight(101) =  .11559704854044D-01
      RacLeg(102) = -.31890066184011D+00
      Weight(102) =  .11608413162253D-01
      RacLeg(103) = -.30726861979932D+00
      Weight(103) =  .11655380094945D-01
      RacLeg(104) = -.29559048446014D+00
      Weight(104) =  .11700598606621D-01
      RacLeg(105) = -.28386800765708D+00
      Weight(105) =  .11744061914061D-01
      RacLeg(106) = -.27210294787634D+00
      Weight(106) =  .11785763497343D-01
      RacLeg(107) = -.26029706999194D+00
      Weight(107) =  .11825697100824D-01
      RacLeg(108) = -.24845214500106D+00
      Weight(108) =  .11863856734071D-01
      RacLeg(109) = -.23656994975828D+00
      Weight(109) =  .11900236672766D-01
      RacLeg(110) = -.22465226670913D+00
      Weight(110) =  .11934831459564D-01
      RacLeg(111) = -.21270088362263D+00
      Weight(111) =  .11967635904906D-01
      RacLeg(112) = -.20071759332313D+00
      Weight(112) =  .11998645087806D-01
      RacLeg(113) = -.18870419342139D+00
      Weight(113) =  .12027854356583D-01
      RacLeg(114) = -.17666248604490D+00
      Weight(114) =  .12055259329560D-01
      RacLeg(115) = -.16459427756755D+00
      Weight(115) =  .12080855895725D-01
      RacLeg(116) = -.15250137833866D+00
      Weight(116) =  .12104640215340D-01
      RacLeg(117) = -.14038560241138D+00
      Weight(117) =  .12126608720527D-01
      RacLeg(118) = -.12824876727061D+00
      Weight(118) =  .12146758115794D-01
      RacLeg(119) = -.11609269356033D+00
      Weight(119) =  .12165085378536D-01
      RacLeg(120) = -.10391920481051D+00
      Weight(120) =  .12181587759482D-01
      RacLeg(121) = -.91730127163520D-01
      Weight(121) =  .12196262783115D-01
      RacLeg(122) = -.79527289100233D-01
      Weight(122) =  .12209108248037D-01
      RacLeg(123) = -.67312521165716D-01
      Weight(123) =  .12220122227304D-01
      RacLeg(124) = -.55087655694634D-01
      Weight(124) =  .12229303068710D-01
      RacLeg(125) = -.42854526536379D-01
      Weight(125) =  .12236649395040D-01
      RacLeg(126) = -.30614968779979D-01
      Weight(126) =  .12242160104273D-01
      RacLeg(127) = -.18370818478814D-01
      Weight(127) =  .12245834369748D-01
      RacLeg(128) = -.61239123751895D-02
      Weight(128) =  .12247671640290D-01
      RacLeg(129) =  .61239123751895D-02
      Weight(129) =  .12247671640290D-01
      RacLeg(130) =  .18370818478814D-01
      Weight(130) =  .12245834369748D-01
      RacLeg(131) =  .30614968779979D-01
      Weight(131) =  .12242160104273D-01
      RacLeg(132) =  .42854526536379D-01
      Weight(132) =  .12236649395040D-01
      RacLeg(133) =  .55087655694634D-01
      Weight(133) =  .12229303068710D-01
      RacLeg(134) =  .67312521165716D-01
      Weight(134) =  .12220122227304D-01
      RacLeg(135) =  .79527289100233D-01
      Weight(135) =  .12209108248037D-01
      RacLeg(136) =  .91730127163520D-01
      Weight(136) =  .12196262783115D-01
      RacLeg(137) =  .10391920481051D+00
      Weight(137) =  .12181587759482D-01
      RacLeg(138) =  .11609269356033D+00
      Weight(138) =  .12165085378536D-01
      RacLeg(139) =  .12824876727061D+00
      Weight(139) =  .12146758115794D-01
      RacLeg(140) =  .14038560241138D+00
      Weight(140) =  .12126608720527D-01
      RacLeg(141) =  .15250137833866D+00
      Weight(141) =  .12104640215340D-01
      RacLeg(142) =  .16459427756755D+00
      Weight(142) =  .12080855895725D-01
      RacLeg(143) =  .17666248604490D+00
      Weight(143) =  .12055259329560D-01
      RacLeg(144) =  .18870419342139D+00
      Weight(144) =  .12027854356583D-01
      RacLeg(145) =  .20071759332313D+00
      Weight(145) =  .11998645087806D-01
      RacLeg(146) =  .21270088362263D+00
      Weight(146) =  .11967635904906D-01
      RacLeg(147) =  .22465226670913D+00
      Weight(147) =  .11934831459564D-01
      RacLeg(148) =  .23656994975828D+00
      Weight(148) =  .11900236672766D-01
      RacLeg(149) =  .24845214500106D+00
      Weight(149) =  .11863856734071D-01
      RacLeg(150) =  .26029706999194D+00
      Weight(150) =  .11825697100824D-01
      RacLeg(151) =  .27210294787634D+00
      Weight(151) =  .11785763497343D-01
      RacLeg(152) =  .28386800765708D+00
      Weight(152) =  .11744061914061D-01
      RacLeg(153) =  .29559048446014D+00
      Weight(153) =  .11700598606621D-01
      RacLeg(154) =  .30726861979932D+00
      Weight(154) =  .11655380094945D-01
      RacLeg(155) =  .31890066184011D+00
      Weight(155) =  .11608413162253D-01
      RacLeg(156) =  .33048486566242D+00
      Weight(156) =  .11559704854044D-01
      RacLeg(157) =  .34201949352237D+00
      Weight(157) =  .11509262477039D-01
      RacLeg(158) =  .35350281511297D+00
      Weight(158) =  .11457093598091D-01
      RacLeg(159) =  .36493310782365D+00
      Weight(159) =  .11403206043039D-01
      RacLeg(160) =  .37630865699872D+00
      Weight(160) =  .11347607895545D-01
      RacLeg(161) =  .38762775619452D+00
      Weight(161) =  .11290307495876D-01
      RacLeg(162) =  .39888870743546D+00
      Weight(162) =  .11231313439650D-01
      RacLeg(163) =  .41008982146872D+00
      Weight(163) =  .11170634576553D-01
      RacLeg(164) =  .42122941801762D+00
      Weight(164) =  .11108280009010D-01
      RacLeg(165) =  .43230582603374D+00
      Weight(165) =  .11044259090814D-01
      RacLeg(166) =  .44331738394753D+00
      Weight(166) =  .10978581425730D-01
      RacLeg(167) =  .45426243991759D+00
      Weight(167) =  .10911256866049D-01
      RacLeg(168) =  .46513935207848D+00
      Weight(168) =  .10842295511115D-01
      RacLeg(169) =  .47594648878698D+00
      Weight(169) =  .10771707705805D-01
      RacLeg(170) =  .48668222886689D+00
      Weight(170) =  .10699504038980D-01
      RacLeg(171) =  .49734496185218D+00
      Weight(171) =  .10625695341897D-01
      RacLeg(172) =  .50793308822862D+00
      Weight(172) =  .10550292686581D-01
      RacLeg(173) =  .51844501967367D+00
      Weight(173) =  .10473307384170D-01
      RacLeg(174) =  .52887917929482D+00
      Weight(174) =  .10394750983212D-01
      RacLeg(175) =  .53923400186606D+00
      Weight(175) =  .10314635267934D-01
      RacLeg(176) =  .54950793406272D+00
      Weight(176) =  .10232972256478D-01
      RacLeg(177) =  .55969943469448D+00
      Weight(177) =  .10149774199095D-01
      RacLeg(178) =  .56980697493657D+00
      Weight(178) =  .10065053576306D-01
      RacLeg(179) =  .57982903855908D+00
      Weight(179) =  .99788230970349D-02
      RacLeg(180) =  .58976412215445D+00
      Weight(180) =  .98910956966958D-02
      RacLeg(181) =  .59961073536297D+00
      Weight(181) =  .98018845352573D-02
      RacLeg(182) =  .60936740109633D+00
      Weight(182) =  .97112029952663D-02
      RacLeg(183) =  .61903265575926D+00
      Weight(183) =  .96190646798407D-02
      RacLeg(184) =  .62860504946901D+00
      Weight(184) =  .95254834106293D-02
      RacLeg(185) =  .63808314627291D+00
      Weight(185) =  .94304732257378D-02
      RacLeg(186) =  .64746552436372D+00
      Weight(186) =  .93340483776233D-02
      RacLeg(187) =  .65675077629297D+00
      Weight(187) =  .92362233309563D-02
      RacLeg(188) =  .66593750918205D+00
      Weight(188) =  .91370127604508D-02
      RacLeg(189) =  .67502434493116D+00
      Weight(189) =  .90364315486629D-02
      RacLeg(190) =  .68400992042608D+00
      Weight(190) =  .89344947837582D-02
      RacLeg(191) =  .69289288774258D+00
      Weight(191) =  .88312177572488D-02
      RacLeg(192) =  .70167191434869D+00
      Weight(192) =  .87266159616988D-02
      RacLeg(193) =  .71034568330454D+00
      Weight(193) =  .86207050884010D-02
      RacLeg(194) =  .71891289345997D+00
      Weight(194) =  .85135010250225D-02
      RacLeg(195) =  .72737225964965D+00
      Weight(195) =  .84050198532215D-02
      RacLeg(196) =  .73572251288592D+00
      Weight(196) =  .82952778462352D-02
      RacLeg(197) =  .74396240054911D+00
      Weight(197) =  .81842914664383D-02
      RacLeg(198) =  .75209068657549D+00
      Weight(198) =  .80720773628735D-02
      RacLeg(199) =  .76010615164266D+00
      Weight(199) =  .79586523687543D-02
      RacLeg(200) =  .76800759335245D+00
      Weight(200) =  .78440334989397D-02
      RacLeg(201) =  .77579382641133D+00
      Weight(201) =  .77282379473816D-02
      RacLeg(202) =  .78346368280818D+00
      Weight(202) =  .76112830845457D-02
      RacLeg(203) =  .79101601198955D+00
      Weight(203) =  .74931864548059D-02
      RacLeg(204) =  .79844968103217D+00
      Weight(204) =  .73739657738123D-02
      RacLeg(205) =  .80576357481300D+00
      Weight(205) =  .72536389258339D-02
      RacLeg(206) =  .81295659617643D+00
      Weight(206) =  .71322239610754D-02
      RacLeg(207) =  .82002766609892D+00
      Weight(207) =  .70097390929698D-02
      RacLeg(208) =  .82697572385081D+00
      Weight(208) =  .68862026954463D-02
      RacLeg(209) =  .83379972715550D+00
      Weight(209) =  .67616333001738D-02
      RacLeg(210) =  .84049865234576D+00
      Weight(210) =  .66360495937811D-02
      RacLeg(211) =  .84707149451730D+00
      Weight(211) =  .65094704150537D-02
      RacLeg(212) =  .85351726767950D+00
      Weight(212) =  .63819147521079D-02
      RacLeg(213) =  .85983500490338D+00
      Weight(213) =  .62534017395424D-02
      RacLeg(214) =  .86602375846655D+00
      Weight(214) =  .61239506555679D-02
      RacLeg(215) =  .87208259999549D+00
      Weight(215) =  .59935809191153D-02
      RacLeg(216) =  .87801062060471D+00
      Weight(216) =  .58623120869227D-02
      RacLeg(217) =  .88380693103316D+00
      Weight(217) =  .57301638506014D-02
      RacLeg(218) =  .88947066177761D+00
      Weight(218) =  .55971560336829D-02
      RacLeg(219) =  .89500096322308D+00
      Weight(219) =  .54633085886443D-02
      RacLeg(220) =  .90039700577030D+00
      Weight(220) =  .53286415939159D-02
      RacLeg(221) =  .90565797996014D+00
      Weight(221) =  .51931752508693D-02
      RacLeg(222) =  .91078309659507D+00
      Weight(222) =  .50569298807868D-02
      RacLeg(223) =  .91577158685749D+00
      Weight(223) =  .49199259218139D-02
      RacLeg(224) =  .92062270242515D+00
      Weight(224) =  .47821839258927D-02
      RacLeg(225) =  .92533571558332D+00
      Weight(225) =  .46437245556801D-02
      RacLeg(226) =  .92990991933401D+00
      Weight(226) =  .45045685814479D-02
      RacLeg(227) =  .93434462750200D+00
      Weight(227) =  .43647368779681D-02
      RacLeg(228) =  .93863917483781D+00
      Weight(228) =  .42242504213815D-02
      RacLeg(229) =  .94279291711746D+00
      Weight(229) =  .40831302860527D-02
      RacLeg(230) =  .94680523123913D+00
      Weight(230) =  .39413976414088D-02
      RacLeg(231) =  .95067551531663D+00
      Weight(231) =  .37990737487663D-02
      RacLeg(232) =  .95440318876972D+00
      Weight(232) =  .36561799581425D-02
      RacLeg(233) =  .95798769241118D+00
      Weight(233) =  .35127377050563D-02
      RacLeg(234) =  .96142848853073D+00
      Weight(234) =  .33687685073156D-02
      RacLeg(235) =  .96472506097571D+00
      Weight(235) =  .32242939617942D-02
      RacLeg(236) =  .96787691522849D+00
      Weight(236) =  .30793357411993D-02
      RacLeg(237) =  .97088357848074D+00
      Weight(237) =  .29339155908297D-02
      RacLeg(238) =  .97374459970437D+00
      Weight(238) =  .27880553253277D-02
      RacLeg(239) =  .97645954971923D+00
      Weight(239) =  .26417768254275D-02
      RacLeg(240) =  .97902802125762D+00
      Weight(240) =  .24951020347037D-02
      RacLeg(241) =  .98144962902546D+00
      Weight(241) =  .23480529563273D-02
      RacLeg(242) =  .98372400976032D+00
      Weight(242) =  .22006516498399D-02
      RacLeg(243) =  .98585082228613D+00
      Weight(243) =  .20529202279661D-02
      RacLeg(244) =  .98782974756486D+00
      Weight(244) =  .19048808534997D-02
      RacLeg(245) =  .98966048874507D+00
      Weight(245) =  .17565557363307D-02
      RacLeg(246) =  .99134277120758D+00
      Weight(246) =  .16079671307493D-02
      RacLeg(247) =  .99287634260882D+00
      Weight(247) =  .14591373333107D-02
      RacLeg(248) =  .99426097292241D+00
      Weight(248) =  .13100886819025D-02
      RacLeg(249) =  .99549645448110D+00
      Weight(249) =  .11608435575677D-02
      RacLeg(250) =  .99658260202338D+00
      Weight(250) =  .10114243932084D-02
      RacLeg(251) =  .99751925275672D+00
      Weight(251) =  .86185370142009D-03
      RacLeg(252) =  .99830626647301D+00
      Weight(252) =  .71215416347332D-03
      RacLeg(253) =  .99894352584341D+00
      Weight(253) =  .56234895403141D-03
      RacLeg(254) =  .99943093746626D+00
      Weight(254) =  .41246325442618D-03
      RacLeg(255) =  .99976843740926D+00
      Weight(255) =  .26253494429645D-03
      RacLeg(256) =  .99995605001899D+00
      Weight(256) =  .11278901782227D-03

      else
         ifail = 1  
         
c.... n ne possède pas une des valeurs permises
         write(*,*) 'fatal error in GLNCOEFS'
         write(*,*)
         write(*,*) 'L''erreur est humaine, mais'
         write(*,*) '      pour que ce soit vraiment le bordel,'
         write(*,*) '               il faut ajouter un ordinateur !'
      end if

      return
      end
