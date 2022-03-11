cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Programme de quadrature de Gauss_Laguerre
c
c     Func : fonction à intégrer
c     a : borne de l'intervalle 
c     b : la constante de l'exponentiel
c       - b < 0 => on intègre sur ]-\infty, a]
c       - b > 0 => on intègre sur [a, +\infty[
c     n : nombre de points de quadrature
c
c     Les valeurs possible pour n sont:
c     1,2,3,4,5,6,8,10,12,14,16,20,24, 32, 40, 48, 64
c
c     -----------------------------------------
c
c     Func doit être définie comme suit:
c
c     subroutine Func( X(), FX(), N)
c     implicit double precision (a-h,o-z)
c     implicit integer (i-n)
c     dimension X(*), FX(*)
c     
c     do i = 1, N
c     FX(i) = Func(X(i))
c     // où Func est la fonction elle même
c     // et X(i) la variable : ex Func(x) == cos(x) * exp(-b.x)
c     end do
c     -----------------------------------------
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function dglgq(Func, a, b, n)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      dimension RacLag(64), Weight(64), Func_Lag(64)
      external Func

      if(n.EQ.1) then
c     Racines et Poids de Gauss-Laguerre pour N = 1
c     i       ième racine            ième poids
         RacLag(  1) =  .10000000000000D+01
         Weight(  1) =  .10000000000000D+01
 
      else if(n.EQ.2) then 
c     Racines et Poids de Gauss-Laguerre pour N = 2
c     i       ième racine            ième poids
         RacLag(  1) =  .58578643762690D+00
         Weight(  1) =  .85355339059327D+00
         RacLag(  2) =  .34142135623731D+01
         Weight(  2) =  .14644660940673D+00
 
      else if(n.EQ.3) then 
c     Racines et Poids de Gauss-Laguerre pour N = 3
c     i       ième racine            ième poids
         RacLag(  1) =  .41577455678348D+00
         Weight(  1) =  .71109300992917D+00
         RacLag(  2) =  .22942803602790D+01
         Weight(  2) =  .27851773356924D+00
         RacLag(  3) =  .62899450829375D+01
         Weight(  3) =  .10389256501586D-01
 
      else if(n.EQ.4) then 
c     Racines et Poids de Gauss-Laguerre pour N = 4
c     i       ième racine            ième poids
         RacLag(  1) =  .32254768961939D+00
         Weight(  1) =  .60315410434163D+00
         RacLag(  2) =  .17457611011583D+01
         Weight(  2) =  .35741869243780D+00
         RacLag(  3) =  .45366202969211D+01
         Weight(  3) =  .38887908515005D-01
         RacLag(  4) =  .93950709123011D+01
         Weight(  4) =  .53929470556133D-03
 
      else if(n.EQ.5) then 
c     Racines et Poids de Gauss-Laguerre pour N = 5
c     i       ième racine            ième poids
         RacLag(  1) =  .26356031971814D+00
         Weight(  1) =  .52175561058281D+00
         RacLag(  2) =  .14134030591065D+01
         Weight(  2) =  .39866681108318D+00
         RacLag(  3) =  .35964257710407D+01
         Weight(  3) =  .75942449681708D-01
         RacLag(  4) =  .70858100058588D+01
         Weight(  4) =  .36117586799220D-02
         RacLag(  5) =  .12640800844276D+02
         Weight(  5) =  .23369972385776D-04
 
      else if(n.EQ.6) then 
c     Racines et Poids de Gauss-Laguerre pour N = 6
c     i       ième racine            ième poids
         RacLag(  1) =  .22284660417926D+00
         Weight(  1) =  .45896467394996D+00
         RacLag(  2) =  .11889321016726D+01
         Weight(  2) =  .41700083077212D+00
         RacLag(  3) =  .29927363260593D+01
         Weight(  3) =  .11337338207404D+00
         RacLag(  4) =  .57751435691045D+01
         Weight(  4) =  .10399197453149D-01
         RacLag(  5) =  .98374674183826D+01
         Weight(  5) =  .26101720281493D-03
         RacLag(  6) =  .15982873980602D+02
         Weight(  6) =  .89854790642962D-06
 
      else if(n.EQ.8) then 
c     Racines et Poids de Gauss-Laguerre pour N = 8
c     i       ième racine            ième poids
         RacLag(  1) =  .17027963230510D+00
         Weight(  1) =  .36918858934164D+00
         RacLag(  2) =  .90370177679938D+00
         Weight(  2) =  .41878678081434D+00
         RacLag(  3) =  .22510866298661D+01
         Weight(  3) =  .17579498663717D+00
         RacLag(  4) =  .42667001702877D+01
         Weight(  4) =  .33343492261216D-01
         RacLag(  5) =  .70459054023935D+01
         Weight(  5) =  .27945362352257D-02
         RacLag(  6) =  .10758516010181D+02
         Weight(  6) =  .90765087733582D-04
         RacLag(  7) =  .15740678641278D+02
         Weight(  7) =  .84857467162725D-06
         RacLag(  8) =  .22863131736889D+02
         Weight(  8) =  .10480011748715D-08
 
      else if(n.EQ.10) then 
c     Racines et Poids de Gauss-Laguerre pour N = 10
c     i       ième racine            ième poids
         RacLag(  1) =  .13779347054049D+00
         Weight(  1) =  .30844111576502D+00
         RacLag(  2) =  .72945454950317D+00
         Weight(  2) =  .40111992915527D+00
         RacLag(  3) =  .18083429017403D+01
         Weight(  3) =  .21806828761181D+00
         RacLag(  4) =  .34014336978549D+01
         Weight(  4) =  .62087456098678D-01
         RacLag(  5) =  .55524961400638D+01
         Weight(  5) =  .95015169751811D-02
         RacLag(  6) =  .83301527467645D+01
         Weight(  6) =  .75300838858754D-03
         RacLag(  7) =  .11843785837900D+02
         Weight(  7) =  .28259233495996D-04
         RacLag(  8) =  .16279257831378D+02
         Weight(  8) =  .42493139849627D-06
         RacLag(  9) =  .21996585811981D+02
         Weight(  9) =  .18395648239796D-08
         RacLag( 10) =  .29920697012274D+02
         Weight( 10) =  .99118272196090D-12
 
      else if(n.EQ.12) then 
c     Racines et Poids de Gauss-Laguerre pour N = 12
c     i       ième racine            ième poids
         RacLag(  1) =  .11572211735802D+00
         Weight(  1) =  .26473137105544D+00
         RacLag(  2) =  .61175748451513D+00
         Weight(  2) =  .37775927587314D+00
         RacLag(  3) =  .15126102697764D+01
         Weight(  3) =  .24408201131988D+00
         RacLag(  4) =  .28337513377435D+01
         Weight(  4) =  .90449222211681D-01
         RacLag(  5) =  .45992276394183D+01
         Weight(  5) =  .20102381154634D-01
         RacLag(  6) =  .68445254531152D+01
         Weight(  6) =  .26639735418653D-02
         RacLag(  7) =  .96213168424569D+01
         Weight(  7) =  .20323159266300D-03
         RacLag(  8) =  .13006054993306D+02
         Weight(  8) =  .83650558568198D-05
         RacLag(  9) =  .17116855187462D+02
         Weight(  9) =  .16684938765409D-06
         RacLag( 10) =  .22151090379397D+02
         Weight( 10) =  .13423910305150D-08
         RacLag( 11) =  .28487967250984D+02
         Weight( 11) =  .30616016350350D-11
         RacLag( 12) =  .37099121044467D+02
         Weight( 12) =  .81480774674263D-15
 
      else if(n.EQ.14) then 
c     Racines et Poids de Gauss-Laguerre pour N = 14
c     i       ième racine            ième poids
         RacLag(  1) =  .99747507032598D-01
         Weight(  1) =  .23181557714486D+00
         RacLag(  2) =  .52685764885190D+00
         Weight(  2) =  .35378469159754D+00
         RacLag(  3) =  .13006291212515D+01
         Weight(  3) =  .25873461024543D+00
         RacLag(  4) =  .24308010787308D+01
         Weight(  4) =  .11548289355692D+00
         RacLag(  5) =  .39321028222932D+01
         Weight(  5) =  .33192092159337D-01
         RacLag(  6) =  .58255362183017D+01
         Weight(  6) =  .61928694370066D-02
         RacLag(  7) =  .81402401415651D+01
         Weight(  7) =  .73989037786739D-03
         RacLag(  8) =  .10916499507366D+02
         Weight(  8) =  .54907194668417D-04
         RacLag(  9) =  .14210805011161D+02
         Weight(  9) =  .24095857640854D-05
         RacLag( 10) =  .18104892220218D+02
         Weight( 10) =  .58015439816765D-07
         RacLag( 11) =  .22723381628270D+02
         Weight( 11) =  .68193146924850D-09
         RacLag( 12) =  .28272981723248D+02
         Weight( 12) =  .32212077518949D-11
         RacLag( 13) =  .35149443660592D+02
         Weight( 13) =  .42213524405166D-14
         RacLag( 14) =  .44366081711117D+02
         Weight( 14) =  .60523750222892D-18
 
      else if(n.EQ.16) then 
c     Racines et Poids de Gauss-Laguerre pour N = 16
c     i       ième racine            ième poids
         RacLag(  1) =  .87649410478928D-01
         Weight(  1) =  .20615171495780D+00
         RacLag(  2) =  .46269632891508D+00
         Weight(  2) =  .33105785495088D+00
         RacLag(  3) =  .11410577748312D+01
         Weight(  3) =  .26579577764421D+00
         RacLag(  4) =  .21292836450984D+01
         Weight(  4) =  .13629693429638D+00
         RacLag(  5) =  .34370866338932D+01
         Weight(  5) =  .47328928694125D-01
         RacLag(  6) =  .50780186145498D+01
         Weight(  6) =  .11299900080339D-01
         RacLag(  7) =  .70703385350482D+01
         Weight(  7) =  .18490709435263D-02
         RacLag(  8) =  .94383143363919D+01
         Weight(  8) =  .20427191530828D-03
         RacLag(  9) =  .12214223368866D+02
         Weight(  9) =  .14844586873981D-04
         RacLag( 10) =  .15441527368782D+02
         Weight( 10) =  .68283193308712D-06
         RacLag( 11) =  .19180156856753D+02
         Weight( 11) =  .18810248410797D-07
         RacLag( 12) =  .23515905693992D+02
         Weight( 12) =  .28623502429739D-09
         RacLag( 13) =  .28578729742882D+02
         Weight( 13) =  .21270790332241D-11
         RacLag( 14) =  .34583398702287D+02
         Weight( 14) =  .62979670025179D-14
         RacLag( 15) =  .41940452647688D+02
         Weight( 15) =  .50504737000355D-17
         RacLag( 16) =  .51701160339543D+02
         Weight( 16) =  .41614623703728D-21
 
      else if(n.EQ.20) then 
c     Racines et Poids de Gauss-Laguerre pour N = 20
c     i       ième racine            ième poids
         RacLag(  1) =  .70539889691989D-01
         Weight(  1) =  .16874680185111D+00
         RacLag(  2) =  .37212681800161D+00
         Weight(  2) =  .29125436200607D+00
         RacLag(  3) =  .91658210248327D+00
         Weight(  3) =  .26668610286700D+00
         RacLag(  4) =  .17073065310283D+01
         Weight(  4) =  .16600245326951D+00
         RacLag(  5) =  .27491992553094D+01
         Weight(  5) =  .74826064668792D-01
         RacLag(  6) =  .40489253138509D+01
         Weight(  6) =  .24964417309283D-01
         RacLag(  7) =  .56151749708616D+01
         Weight(  7) =  .62025508445722D-02
         RacLag(  8) =  .74590174536711D+01
         Weight(  8) =  .11449623864769D-02
         RacLag(  9) =  .95943928695811D+01
         Weight(  9) =  .15574177302781D-03
         RacLag( 10) =  .12038802546964D+02
         Weight( 10) =  .15401440865225D-04
         RacLag( 11) =  .14814293442631D+02
         Weight( 11) =  .10864863665180D-05
         RacLag( 12) =  .17948895520519D+02
         Weight( 12) =  .53301209095567D-07
         RacLag( 13) =  .21478788240285D+02
         Weight( 13) =  .17579811790506D-08
         RacLag( 14) =  .25451702793187D+02
         Weight( 14) =  .37255024025123D-10
         RacLag( 15) =  .29932554631701D+02
         Weight( 15) =  .47675292515782D-12
         RacLag( 16) =  .35013434240479D+02
         Weight( 16) =  .33728442433624D-14
         RacLag( 17) =  .40833057056729D+02
         Weight( 17) =  .11550143395004D-16
         RacLag( 18) =  .47619994047347D+02
         Weight( 18) =  .15395221405823D-19
         RacLag( 19) =  .55810795750064D+02
         Weight( 19) =  .52864427255691D-23
         RacLag( 20) =  .66524416525616D+02
         Weight( 20) =  .16564566124990D-27
 
      else if(n.EQ.24) then 
c     Racines et Poids de Gauss-Laguerre pour N = 24
c     i       ième racine            ième poids
         RacLag(  1) =  .59019852181508D-01
         Weight(  1) =  .14281197333478D+00
         RacLag(  2) =  .31123914619848D+00
         Weight(  2) =  .25877410751742D+00
         RacLag(  3) =  .76609690554594D+00
         Weight(  3) =  .25880670727287D+00
         RacLag(  4) =  .14255975908036D+01
         Weight(  4) =  .18332268897778D+00
         RacLag(  5) =  .22925620586322D+01
         Weight(  5) =  .98166272629919D-01
         RacLag(  6) =  .33707742642090D+01
         Weight(  6) =  .40732478151409D-01
         RacLag(  7) =  .46650837034672D+01
         Weight(  7) =  .13226019405120D-01
         RacLag(  8) =  .61815351187368D+01
         Weight(  8) =  .33693490584783D-02
         RacLag(  9) =  .79275392471722D+01
         Weight(  9) =  .67216256409355D-03
         RacLag( 10) =  .99120980150777D+01
         Weight( 10) =  .10446121465928D-03
         RacLag( 11) =  .12146102711730D+02
         Weight( 11) =  .12544721977993D-04
         RacLag( 12) =  .14642732289597D+02
         Weight( 12) =  .11513158127373D-05
         RacLag( 13) =  .17417992646509D+02
         Weight( 13) =  .79608129591336D-07
         RacLag( 14) =  .20491460082616D+02
         Weight( 14) =  .40728589875500D-08
         RacLag( 15) =  .23887329848170D+02
         Weight( 15) =  .15070082262926D-09
         RacLag( 16) =  .27635937174333D+02
         Weight( 16) =  .39177365150585D-11
         RacLag( 17) =  .31776041352375D+02
         Weight( 17) =  .68941810529581D-13
         RacLag( 18) =  .36358405801652D+02
         Weight( 18) =  .78198003824595D-15
         RacLag( 19) =  .41451720484871D+02
         Weight( 19) =  .53501888130100D-17
         RacLag( 20) =  .47153106445156D+02
         Weight( 20) =  .20105174645555D-19
         RacLag( 21) =  .53608574544695D+02
         Weight( 21) =  .36057658645530D-22
         RacLag( 22) =  .61058531447219D+02
         Weight( 22) =  .24518188458784D-25
         RacLag( 23) =  .69962240035105D+02
         Weight( 23) =  .40883015936807D-29
         RacLag( 24) =  .81498279233949D+02
         Weight( 24) =  .55753457883283D-34
 
      else if(n.EQ.32) then 
c     Racines et Poids de Gauss-Laguerre pour N = 32
c     i       ième racine            ième poids
         RacLag(  1) =  .44489365833267D-01
         Weight(  1) =  .10921834195238D+00
         RacLag(  2) =  .23452610951962D+00
         Weight(  2) =  .21044310793881D+00
         RacLag(  3) =  .57688462930189D+00
         Weight(  3) =  .23521322966985D+00
         RacLag(  4) =  .10724487538178D+01
         Weight(  4) =  .19590333597288D+00
         RacLag(  5) =  .17224087764446D+01
         Weight(  5) =  .12998378628607D+00
         RacLag(  6) =  .25283367064258D+01
         Weight(  6) =  .70578623865717D-01
         RacLag(  7) =  .34922132730220D+01
         Weight(  7) =  .31760912509175D-01
         RacLag(  8) =  .46164567697498D+01
         Weight(  8) =  .11918214834839D-01
         RacLag(  9) =  .59039585041742D+01
         Weight(  9) =  .37388162946115D-02
         RacLag( 10) =  .73581267331862D+01
         Weight( 10) =  .98080330661496D-03
         RacLag( 11) =  .89829409242126D+01
         Weight( 11) =  .21486491880136D-03
         RacLag( 12) =  .10783018632540D+02
         Weight( 12) =  .39203419679879D-04
         RacLag( 13) =  .12763697986743D+02
         Weight( 13) =  .59345416128686D-05
         RacLag( 14) =  .14931139755523D+02
         Weight( 14) =  .74164045786676D-06
         RacLag( 15) =  .17292454336715D+02
         Weight( 15) =  .76045678791208D-07
         RacLag( 16) =  .19855860940336D+02
         Weight( 16) =  .63506022266258D-08
         RacLag( 17) =  .22630889013197D+02
         Weight( 17) =  .42813829710409D-09
         RacLag( 18) =  .25628636022459D+02
         Weight( 18) =  .23058994918913D-10
         RacLag( 19) =  .28862101816323D+02
         Weight( 19) =  .97993792887271D-12
         RacLag( 20) =  .32346629153965D+02
         Weight( 20) =  .32378016577293D-13
         RacLag( 21) =  .36100494805752D+02
         Weight( 21) =  .81718234434207D-15
         RacLag( 22) =  .40145719771539D+02
         Weight( 22) =  .15421338333938D-16
         RacLag( 23) =  .44509207995755D+02
         Weight( 23) =  .21197922901636D-18
         RacLag( 24) =  .49224394987309D+02
         Weight( 24) =  .20544296737881D-20
         RacLag( 25) =  .54333721333397D+02
         Weight( 25) =  .13469825866374D-22
         RacLag( 26) =  .59892509162134D+02
         Weight( 26) =  .56612941303974D-25
         RacLag( 27) =  .65975377287935D+02
         Weight( 27) =  .14185605454630D-27
         RacLag( 28) =  .72687628090663D+02
         Weight( 28) =  .19133754944542D-30
         RacLag( 29) =  .80187446977914D+02
         Weight( 29) =  .11922487600982D-33
         RacLag( 30) =  .88735340417892D+02
         Weight( 30) =  .26715112192401D-37
         RacLag( 31) =  .98829542868284D+02
         Weight( 31) =  .13386169421063D-41
         RacLag( 32) =  .11175139809794D+03
         Weight( 32) =  .45105361938990D-47
 
      else if(n.EQ.40) then 
c     Racines et Poids de Gauss-Laguerre pour N = 40
c     i       ième racine            ième poids
         RacLag(  1) =  .35700394308888D-01
         Weight(  1) =  .88412106190342D-01
         RacLag(  2) =  .18816228315870D+00
         Weight(  2) =  .17681473909572D+00
         RacLag(  3) =  .46269428131458D+00
         Weight(  3) =  .21136311701596D+00
         RacLag(  4) =  .85977296397293D+00
         Weight(  4) =  .19408119531860D+00
         RacLag(  5) =  .13800108205273D+01
         Weight(  5) =  .14643428242413D+00
         RacLag(  6) =  .20242091359228D+01
         Weight(  6) =  .93326798435771D-01
         RacLag(  7) =  .27933693535068D+01
         Weight(  7) =  .50932204361044D-01
         RacLag(  8) =  .36887026779083D+01
         Weight(  8) =  .23976193015685D-01
         RacLag(  9) =  .47116411465550D+01
         Weight(  9) =  .97746252467145D-02
         RacLag( 10) =  .58638508783437D+01
         Weight( 10) =  .34579399930185D-02
         RacLag( 11) =  .71472479081023D+01
         Weight( 11) =  .10622468938969D-02
         RacLag( 12) =  .85640170175862D+01
         Weight( 12) =  .28327168532432D-03
         RacLag( 13) =  .10116634048452D+02
         Weight( 13) =  .65509405003246D-04
         RacLag( 14) =  .11807892294005D+02
         Weight( 14) =  .13116069073268D-04
         RacLag( 15) =  .13640933712537D+02
         Weight( 15) =  .22684528787794D-05
         RacLag( 16) =  .15619285893339D+02
         Weight( 16) =  .33796264822007D-06
         RacLag( 17) =  .17746905950096D+02
         Weight( 17) =  .43228213222821D-07
         RacLag( 18) =  .20028232834575D+02
         Weight( 18) =  .47284937709908D-08
         RacLag( 19) =  .22468249983498D+02
         Weight( 19) =  .44031741042329D-09
         RacLag( 20) =  .25072560772426D+02
         Weight( 20) =  .34724414848038D-10
         RacLag( 21) =  .27847480009169D+02
         Weight( 21) =  .23053815449168D-11
         RacLag( 22) =  .30800145739445D+02
         Weight( 22) =  .12797725976766D-12
         RacLag( 23) =  .33938657084914D+02
         Weight( 23) =  .58941771723512D-14
         RacLag( 24) =  .37272245880476D+02
         Weight( 24) =  .22322175799046D-15
         RacLag( 25) =  .40811492823887D+02
         Weight( 25) =  .68803364842843D-17
         RacLag( 26) =  .44568603175334D+02
         Weight( 26) =  .17056037368181D-18
         RacLag( 27) =  .48557763533060D+02
         Weight( 27) =  .33537119406662D-20
         RacLag( 28) =  .52795611187217D+02
         Weight( 28) =  .51461995601367D-22
         RacLag( 29) =  .57301863323394D+02
         Weight( 29) =  .60447625115877D-24
         RacLag( 30) =  .62100179072775D+02
         Weight( 30) =  .53105847773213D-26
         RacLag( 31) =  .67219370927127D+02
         Weight( 31) =  .33925280532805D-28
         RacLag( 32) =  .72695158847612D+02
         Weight( 32) =  .15217354931815D-30
         RacLag( 33) =  .78572802911571D+02
         Weight( 33) =  .45852916145027D-33
         RacLag( 34) =  .84911231135705D+02
         Weight( 34) =  .87621586574863D-36
         RacLag( 35) =  .91789874671236D+02
         Weight( 35) =  .98274157251479D-39
         RacLag( 36) =  .99320808717447D+02
         Weight( 36) =  .58011520191698D-42
         RacLag( 37) =  .10767244063939D+03
         Weight( 37) =  .15309086846067D-45
         RacLag( 38) =  .11712230951269D+03
         Weight( 38) =  .13819863056493D-49
         RacLag( 39) =  .12820184198826D+03
         Weight( 39) =  .25666336050124D-54
         RacLag( 40) =  .14228004446916D+03
         Weight( 40) =  .27003609402170D-60
 
      else if(n.EQ.48) then 
c     Racines et Poids de Gauss-Laguerre pour N = 48
c     i       ième racine            ième poids
         RacLag(  1) =  .29811235829960D-01
         Weight(  1) =  .74262005828026D-01
         RacLag(  2) =  .15710799061788D+00
         Weight(  2) =  .15227194980935D+00
         RacLag(  3) =  .38626503757646D+00
         Weight(  3) =  .19040908826391D+00
         RacLag(  4) =  .71757469411697D+00
         Weight(  4) =  .18663305948481D+00
         RacLag(  5) =  .11513938340264D+01
         Weight(  5) =  .15342420015758D+00
         RacLag(  6) =  .16881858234190D+01
         Weight(  6) =  .10877969280749D+00
         RacLag(  7) =  .23285270066532D+01
         Weight(  7) =  .67460738609219D-01
         RacLag(  8) =  .30731108616526D+01
         Weight(  8) =  .36881194115821D-01
         RacLag(  9) =  .39227524130465D+01
         Weight(  9) =  .17856844269157D-01
         RacLag( 10) =  .48783933559213D+01
         Weight( 10) =  .76776165144976D-02
         RacLag( 11) =  .59411080546246D+01
         Weight( 11) =  .29357859037395D-02
         RacLag( 12) =  .71121105358907D+01
         Weight( 12) =  .99906553781589D-03
         RacLag( 13) =  .83927625990912D+01
         Weight( 13) =  .30259801699226D-03
         RacLag( 14) =  .97845831846873D+01
         Weight( 14) =  .81538711803554D-04
         RacLag( 15) =  .11289259168010D+02
         Weight( 15) =  .19531587157281D-04
         RacLag( 16) =  .12908657778286D+02
         Weight( 16) =  .41541829450522D-05
         RacLag( 17) =  .14644840883210D+02
         Weight( 17) =  .78337003802776D-06
         RacLag( 18) =  .16500081428965D+02
         Weight( 18) =  .13073947749206D-06
         RacLag( 19) =  .18476882386874D+02
         Weight( 19) =  .19270714080170D-07
         RacLag( 20) =  .20577998634022D+02
         Weight( 20) =  .25026389371263D-08
         RacLag( 21) =  .22806462290521D+02
         Weight( 21) =  .28557855087716D-09
         RacLag( 22) =  .25165612156439D+02
         Weight( 22) =  .28546224120592D-10
         RacLag( 23) =  .27659128044481D+02
         Weight( 23) =  .24910106849372D-11
         RacLag( 24) =  .30291071001009D+02
         Weight( 24) =  .18903366069715D-12
         RacLag( 25) =  .33065930662499D+02
         Weight( 25) =  .12421626859492D-13
         RacLag( 26) =  .35988681327479D+02
         Weight( 26) =  .70342315202126D-15
         RacLag( 27) =  .39064848764198D+02
         Weight( 27) =  .34145491485919D-16
         RacLag( 28) =  .42300590362903D+02
         Weight( 28) =  .14123154148957D-17
         RacLag( 29) =  .45702792038511D+02
         Weight( 29) =  .49442180080975D-19
         RacLag( 30) =  .49279186382837D+02
         Weight( 30) =  .14539524813679D-20
         RacLag( 31) =  .53038498087817D+02
         Weight( 31) =  .35610683650041D-22
         RacLag( 32) =  .56990624814804D+02
         Weight( 32) =  .71940559964947D-24
         RacLag( 33) =  .61146864786140D+02
         Weight( 33) =  .11855372283506D-25
         RacLag( 34) =  .65520206929019D+02
         Weight( 34) =  .15734913570756D-27
         RacLag( 35) =  .70125706236113D+02
         Weight( 35) =  .16572854409195D-29
         RacLag( 36) =  .74980977518911D+02
         Weight( 36) =  .13614341627163D-31
         RacLag( 37) =  .80106857350324D+02
         Weight( 37) =  .85461558139632D-34
         RacLag( 38) =  .85528311116034D+02
         Weight( 38) =  .40000905324813D-36
         RacLag( 39) =  .91275707993668D+02
         Weight( 39) =  .13550199911030D-38
         RacLag( 40) =  .97386667713582D+02
         Weight( 40) =  .32016367953549D-41
         RacLag( 41) =  .10390883335718D+03
         Weight( 41) =  .50358691660611D-44
         RacLag( 42) =  .11090422088498D+03
         Weight( 42) =  .49624875407027D-47
         RacLag( 43) =  .11845642504628D+03
         Weight( 43) =  .28235107161201D-50
         RacLag( 44) =  .12668342576889D+03
         Weight( 44) =  .82684460695051D-54
         RacLag( 45) =  .13576258957786D+03
         Weight( 45) =  .10490648478213D-57
         RacLag( 46) =  .14598643270946D+03
         Weight( 46) =  .43465744227388D-62
         RacLag( 47) =  .15791561202298D+03
         Weight( 47) =  .34347364383966D-67
         RacLag( 48) =  .17299632814856D+03
         Weight( 48) =  .13190660883980D-73
 
      else if(n.EQ.64) then 
c     Racines et Poids de Gauss-Laguerre pour N = 64
c     i       ième racine            ième poids
         RacLag(  1) =  .22415874146705D-01
         Weight(  1) =  .56252842339030D-01
         RacLag(  2) =  .11812251209677D+00
         Weight(  2) =  .11902398731243D+00
         RacLag(  3) =  .29036574401804D+00
         Weight(  3) =  .15749640386214D+00
         RacLag(  4) =  .53928622122798D+00
         Weight(  4) =  .16754705041577D+00
         RacLag(  5) =  .86503700464811D+00
         Weight(  5) =  .15335285577924D+00
         RacLag(  6) =  .12678140407752D+01
         Weight(  6) =  .12422105360933D+00
         RacLag(  7) =  .17478596260594D+01
         Weight(  7) =  .90342300986485D-01
         RacLag(  8) =  .23054637393075D+01
         Weight(  8) =  .59477755768355D-01
         RacLag(  9) =  .29409651567253D+01
         Weight(  9) =  .35627518904036D-01
         RacLag( 10) =  .36547526502073D+01
         Weight( 10) =  .19480410431166D-01
         RacLag( 11) =  .44472663433131D+01
         Weight( 11) =  .97435948993820D-02
         RacLag( 12) =  .53189992544964D+01
         Weight( 12) =  .44643103641663D-02
         RacLag( 13) =  .62704990469237D+01
         Weight( 13) =  .18753595813231D-02
         RacLag( 14) =  .73023700025874D+01
         Weight( 14) =  .72264698157501D-03
         RacLag( 15) =  .84152752394830D+01
         Weight( 15) =  .25548753283350D-03
         RacLag( 16) =  .96099391927961D+01
         Weight( 16) =  .82871435343969D-04
         RacLag( 17) =  .10887150383886D+02
         Weight( 17) =  .24656863967886D-04
         RacLag( 18) =  .12247764504244D+02
         Weight( 18) =  .67267138788297D-05
         RacLag( 19) =  .13692707845548D+02
         Weight( 19) =  .16817853699641D-05
         RacLag( 20) =  .15222981111525D+02
         Weight( 20) =  .38508129815467D-06
         RacLag( 21) =  .16839663652649D+02
         Weight( 21) =  .80687280409905D-07
         RacLag( 22) =  .18543918170859D+02
         Weight( 22) =  .15457237067577D-07
         RacLag( 23) =  .20336995948730D+02
         Weight( 23) =  .27044801476175D-08
         RacLag( 24) =  .22220242665951D+02
         Weight( 24) =  .43167754754272D-09
         RacLag( 25) =  .24195104875933D+02
         Weight( 25) =  .62777525417615D-10
         RacLag( 26) =  .26263137227118D+02
         Weight( 26) =  .83063173762890D-11
         RacLag( 27) =  .28426010527501D+02
         Weight( 27) =  .99840317872202D-12
         RacLag( 28) =  .30685520767526D+02
         Weight( 28) =  .10883538871167D-12
         RacLag( 29) =  .33043599236438D+02
         Weight( 29) =  .10740174034416D-13
         RacLag( 30) =  .35502323891141D+02
         Weight( 30) =  .95757372315744D-15
         RacLag( 31) =  .38063932165646D+02
         Weight( 31) =  .76970280236486D-16
         RacLag( 32) =  .40730835444459D+02
         Weight( 32) =  .55648811374540D-17
         RacLag( 33) =  .43505635466422D+02
         Weight( 33) =  .36097564090104D-18
         RacLag( 34) =  .46391142978616D+02
         Weight( 34) =  .20950953695489D-19
         RacLag( 35) =  .49390399025625D+02
         Weight( 35) =  .10847933010975D-20
         RacLag( 36) =  .52506699341346D+02
         Weight( 36) =  .49946994863638D-22
         RacLag( 37) =  .55743622413278D+02
         Weight( 37) =  .20378369745988D-23
         RacLag( 38) =  .59105061919017D+02
         Weight( 38) =  .73395375642789D-25
         RacLag( 39) =  .62595264400151D+02
         Weight( 39) =  .23237830821987D-26
         RacLag( 40) =  .66218873251248D+02
         Weight( 40) =  .64382347069087D-28
         RacLag( 41) =  .69980980377147D+02
         Weight( 41) =  .15531210957883D-29
         RacLag( 42) =  .73887187232483D+02
         Weight( 42) =  .32442500920196D-31
         RacLag( 43) =  .77943677434463D+02
         Weight( 43) =  .58323862678363D-33
         RacLag( 44) =  .82157303778319D+02
         Weight( 44) =  .89632548331029D-35
         RacLag( 45) =  .86535693349457D+02
         Weight( 45) =  .11687039895507D-36
         RacLag( 46) =  .91087375613133D+02
         Weight( 46) =  .12820559843600D-38
         RacLag( 47) =  .95821940015521D+02
         Weight( 47) =  .11720949374050D-40
         RacLag( 48) =  .10075023196951D+03
         Weight( 48) =  .88353396723286D-43
         RacLag( 49) =  .10588459946880D+03
         Weight( 49) =  .54249555903062D-45
         RacLag( 50) =  .11123920752444D+03
         Weight( 50) =  .26755426666789D-47
         RacLag( 51) =  .11683044505131D+03
         Weight( 51) =  .10429170314114D-49
         RacLag( 52) =  .12267746026854D+03
         Weight( 52) =  .31529023519578D-52
         RacLag( 53) =  .12880287876924D+03
         Weight( 53) =  .72295419106475D-55
         RacLag( 54) =  .13523378794953D+03
         Weight( 54) =  .12242353012301D-57
         RacLag( 55) =  .14200312148993D+03
         Weight( 55) =  .14821685049019D-60
         RacLag( 56) =  .14915166590005D+03
         Weight( 56) =  .12325193488145D-63
         RacLag( 57) =  .15673107513267D+03
         Weight( 57) =  .66914990045713D-67
         RacLag( 58) =  .16480860265515D+03
         Weight( 58) =  .22204659418504D-70
         RacLag( 59) =  .17347494683642D+03
         Weight( 59) =  .41209460947389D-74
         RacLag( 60) =  .18285820469143D+03
         Weight( 60) =  .37743990618965D-78
         RacLag( 61) =  .19315113603707D+03
         Weight( 61) =  .14141150529176D-82
         RacLag( 62) =  .20467202848506D+03
         Weight( 62) =  .15918330640414D-87
         RacLag( 63) =  .21803185193533D+03
         Weight( 63) =  .29894843488606D-93
         RacLag( 64) =  .23480957917133D+03
         Weight( 64) =  .20890635084370-100
 
      else 
c.... n ne possède pas une des valeurs permises
         write(*,*) 'fatal error in DGLGQ'
         write(*,*) 'L''erreur est humaine, mais'
         write(*,*) '      pour que ce soit vraiment le bordel,'
         write(*,*) '             il faut ajouter un ordinateur !'
      end if

      do i = 1, n
         RacLag(i) = RacLag(i)/b + a
         Weight(i) = Weight(i) 
      end do

c.... evaluation de la fonction aux points de quadrature
      call Func(RacLag, Func_Lag, n)

      dglgq = 0.0d0
      do i =1, n
         dglgq = dglgq + Weight(i) * Func_Lag(i) 
     $         * dexp(b*RacLag(i)) * dexp(-a*b)/dabs(b)
      end do

      return
      end

      subroutine Func(x, fx, n)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      dimension x(*), fx(*)

      do i = 1, n
         fx(i) = dsin(3.0d0 * x(i)) * dexp(-1.50d0 * x(i))
      end do
      return
      end

      subroutine Func1(x, fx, n)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      dimension x(*), fx(*)

      do i = 1, n
         fx(i) = dsin(3.0d0 * x(i)) * dexp(1.50d0 * x(i))
      end do
      return
      end

      subroutine Func2(x, fx, n)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      dimension x(*), fx(*)

      do i = 1, n
         fx(i) = dsin(x(i)) * dexp(-x(i))
      end do
      return
      end


      subroutine Func3(x, fx, n)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      dimension x(*), fx(*)

      do i = 1, n
         fx(i) = dsin(x(i)) * dexp(-2.0d0*x(i))
      end do
      return
      end


c      program test
c      implicit double precision (a-h,o-z)
c      implicit integer (i-n)
c      external Func, Func1, Func2, Func3

c      xint = dglgq(Func, -2.0d0,  1.50d0, 20)
c      write(*,*) xint

c     xint = dglgq(Func1,-2.0d0, -1.50d0, 20)
c      write(*,*) xint

c      xint = dglgq(Func2, 0.0d0,  1.0d0, 20)
c      write(*,*) xint

c      xint = dglgq(Func3, 0.0d0,  2.0d0, 20)
c      write(*,*) xint
c      return
c      end
