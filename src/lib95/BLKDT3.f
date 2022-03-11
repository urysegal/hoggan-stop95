      BLOCK DATA BLKDT3
      IMPLICIT REAL*8 (A-H, O-Z)
      INCLUDE "SIZE.INCL"
      DIMENSION POW1(LTOT6), POW2(LTOT6), POW3(LTOT6), POW4(LTOT6),
     $          POW5(LTOT6), POW6(LTOT6), POW7(LTOT6), POW8(LTOT6),
     $          POW9(LTOT6), POW10(LTOT6)

      COMMON/X4POW/X4NL(10*LTOT6)

      DATA POW1/ .273104674257146D-10, .202395347413040D-07,
     A.687949771026924D-06, .737812382673878D-05, .426060385798173D-04,
     B.165842859780055D-03, .488368225171792D-03, .116316195635049D-02,
     C.233626535568468D-02, .406752851154131D-02, .625398621413452D-02,
     D.860083458106947D-02, .106712406520364D-01, .120104756817440D-01,
     E.124649427167331D-01, .129815920015549D-01, .139417831801759D-01,
     F.153983234430963D-01, .174268699533264D-01, .201248796077552D-01,
     G.236088642905854D-01, .280093301134574D-01, .334627320616066D-01,
     H.401000159068799D-01, .480317535446184D-01, .573304750782334D-01,
     I.680114913649849D-01, .800141738330133D-01, .931861769932777D-01,
     J.107273314393253D+00, .121917621236246D+00, .136665505974894D+00,
     K.150986838503578D+00, .164304460321500D+00, .176032123797871D+00,
     L.185617515420140D+00, .192586082183072D+00, .196581905995263D+00,
     M.198891679438442D+00, .204744248965846D+00, .215446242840244D+00,
     N.231279439588142D+00, .252609335831536D+00, .279841807533795D+00,
     O.313360958988255D+00, .353449913479522D+00, .400197797836263D+00,
     P.453399889457478D+00, .512461963261454D+00, .576323268976899D+00,
     Q.643414191423295D+00, .711663542973871D+00, .778566085946715D+00,
     R.841313440086147D+00, .896981941617831D+00, .942760944615913D+00,
     S.976196757133190D+00, .995426929475141D+00, .191867238303669D-08,
     T.131895790926355D-05, .391040020036561D-04, .343051089716158D-03,
     U.151346405588412D-02, .417908065611975D-02, .804016307333669D-02,
     V.113939987935775D-01, .126878503199337D-01, .142098454767706D-01,
     W.171978025103350D-01, .220738216649183D-01, .294121484255923D-01,
     X.398574438014183D-01, .539703389825902D-01, .720028607253628D-01,
     Y.936446177969677D-01, .117820047156796D+00, .142634817019054D+00,
     Z.165546585377529D+00, .183769415405175D+00, .194835397090169D+00,
     1.198483253934161D+00, .202570261046350D+00, .210008684403409D+00,
     2.220941184793960D+00, .235558782084044D+00, .254087178274428D+00,
     3.276766270782811D+00, .303823579503825D+00, .335441805741553D+00,
     4.371721391177556D+00, .412639850233104D+00, .458010653055330D+00,
     5.507445380378982D+00, .560323571688387D+00, .615774958270275D+00,
     6.672678459072591D+00, .729681334530121D+00, .785240254793436D+00,
     7.837683871040825D+00, .885294018363790D+00, .926400248103472D+00,
     8.959480362582479D+00, .983258513063104D+00, .996795338409250D+00/

      DATA POW2/ .745861631011020D-21, .409638766544451D-15,
     A.473274887455998D-12, .544367112026905D-10, .181527452346488D-08,
     B.275038541400269D-07, .238503523357446D-06, .135294573670109D-05,
     C.545813581217244D-05, .165447881922015D-04, .391123435665846D-04,
     D.739743554909204D-04, .113875377053675D-03, .144251526101765D-03,
     E.155374796931438D-03, .168521730894835D-03, .194373318243036D-03,
     F.237108364858209D-03, .303695796370150D-03, .405010779226640D-03,
     G.557378473091280D-03, .784522573404629D-03, .111975443702687D-02,
     H.160801127573202D-02, .230704934857096D-02, .328678337269594D-02,
     I.462556295768942D-02, .640226801417968D-02, .868366358262247D-02,
     J.115075639809138D-01, .148639063679048D-01, .186774605233738D-01,
     K.227970254013055D-01, .269959556815394D-01, .309873086087889D-01,
     L.344538620307461D-01, .370893990506251D-01, .386444457647303D-01,
     M.395579001498440D-01, .419202074845881D-01, .464170835539773D-01,
     N.534901791762053D-01, .638114765492497D-01, .783114372437814D-01,
     O.981950906180385D-01, .124926841338682D+00, .160158277392995D+00,
     P.205571459760054D+00, .262617263789784D+00, .332148510364220D+00,
     Q.413981821724893D+00, .506464998398123D+00, .606165150186388D+00,
     R.707808304469587D+00, .804576603588494D+00, .888798198693089D+00,
     S.952960108637356D+00, .990874771924308D+00, .368130371342769D-17,
     T.173964996640888D-11, .152912297270194D-08, .117684050155444D-06,
     U.229057344845320D-05, .174647151303543D-04, .646442222458469D-04,
     V.129823208508047D-03, .160981545741042D-03, .201919708473698D-03,
     W.295764411184486D-03, .487253602894616D-03, .865074475009069D-03,
     X.158861582638322D-02, .291279748989569D-02, .518441195263599D-02,
     Y.876931444234015D-02, .138815635120296D-01, .203446910260590D-01,
     Z.274056719301595D-01, .337711980383598D-01, .379608319592840D-01,
     1.393956020922926D-01, .410347106603862D-01, .441036475248506D-01,
     2.488150071381589D-01, .554879398169179D-01, .645602941634611D-01,
     3.765995686430241D-01, .923087674625170D-01, .112521205039154D+00,
     4.138176792658978D+00, .170271646000399D+00, .209773758312170D+00,
     5.257500814067970D+00, .313962504989631D+00, .379178799232759D+00,
     6.452496309300275D+00, .532434849961659D+00, .616602257748060D+00,
     7.701714267801942D+00, .783745498950706D+00, .858217419686175D+00,
     8.920602566181405D+00, .966797303511067D+00, .993600946674412D+00/

      DATA POW3/ .203698297778169D-31, .829089804686133D-23,
     A.325589350458147D-18, .401640795973868D-15, .773416563797040D-13,
     B.456131782555556D-11, .116477542399295D-09, .157369500993730D-08,
     C.127516536046003D-07, .672963976891914D-07, .244608057467913D-06,
     D.636241194818634D-06, .121519155288116D-05, .173252944629971D-05,
     E.193673794337441D-05, .218768035387259D-05, .270991065895574D-05,
     F.365107129315039D-05, .529246714871449D-05, .815079317177923D-05,
     G.131590727297057D-04, .219739517399493D-04, .374700427010254D-04,
     H.644812777352963D-04, .110811625725833D-03, .188432852235897D-03,
     I.314591435155088D-03, .512272185812114D-03, .809197411560337D-03,
     J.123445452882505D-02, .181217210665325D-02, .255256459275299D-02,
     K.344205079262888D-02, .443555592911846D-02, .545476174518515D-02,
     L.639524026677540D-02, .714290205368444D-02, .759679880456125D-02,
     M.786773719586068D-02, .858292139792441D-02, .100003862553061D-01,
     N.123711786633421D-01, .161193747095356D-01, .219148141488691D-01,
     O.307705077640071D-01, .441553812624271D-01, .640949899179259D-01,
     P.932060771308207D-01, .134581358588064D+00, .191424915278915D+00,
     Q.266361779089064D+00, .360432675152264D+00, .471939628417919D+00,
     R.595488639554851D+00, .721690684067087D+00, .837924229372819D+00,
     S.930276567729079D+00, .986343431710994D+00, .706321576852411D-26,
     T.229452508254507D-17, .597948277883734D-13, .403716416480359D-10,
     U.346670058159645D-08, .729864531659056D-07, .519750088605628D-06,
     V.147920548111905D-05, .204250975663390D-05, .286924785612582D-05,
     W.508649793313632D-05, .107555491358847D-04, .254436988581581D-04,
     X.633181660221126D-04, .157204667917308D-03, .373292491768556D-03,
     Y.821199099294372D-03, .163552646759738D-02, .290186128181112D-02,
     Z.453691540801469D-02, .620611332104178D-02, .739611376866028D-02,
     1.781936729397366D-02, .831241205043587D-02, .926214899408552D-02,
     2.107852455128304D-01, .130706715236259D-01, .164039429725609D-01,
     3.212001769569017D-01, .280455801500481D-01, .377443162025492D-01,
     4.513632695956478D-01, .702608665045485D-01, .960786160384280D-01,
     5.130667598542619D+00, .175920592172023D+00, .233488809274525D+00,
     6.304384520076143D+00, .388507771870368D+00, .484180913980294D+00,
     7.587814724216909D+00, .693845202140604D+00, .795052830523994D+00,
     8.883300083994095D+00, .950611679083710D+00, .990416791884071D+00/

      DATA POW4/ .556309572614419D-42, .167803919056059D-30,
     A.223989119096487D-24, .296335552656512D-20, .329522159554064D-17,
     B.756461992555875D-15, .568839306539160D-13, .183046216645767D-11,
     C.297912465441193D-10, .273730016324809D-09, .152977541927055D-08,
     D.547220527029706D-08, .129676014991167D-07, .208085027826881D-07,
     E.241413275214855D-07, .283995737837913D-07, .377809868448086D-07,
     F.562203766857335D-07, .922311367328994D-07, .164033731289770D-06,
     G.310670762265566D-06, .615475668181421D-06, .125384999924137D-05,
     H.258570026288132D-05, .532247669674170D-05, .108029449390305D-04,
     I.213958326755485D-04, .409890357253882D-04, .754060132161637D-04,
     J.132424028774825D-03, .220935712513842D-03, .348847531602188D-03,
     K.519704367147769D-03, .728781623159638D-03, .960213294816320D-03,
     L.118706860883369D-02, .137562352193651D-02, .149339318846318D-02,
     M.156482746426503D-02, .175730379555092D-02, .215454564565691D-02,
     N.286119926830254D-02, .407190453939544D-02, .613268120318672D-02,
     O.964227582148480D-02, .156067156868602D-01, .256506738174914D-01,
     P.422596250678793D-01, .689678272404331D-01, .110322632937170D+00,
     Q.171380948718661D+00, .256506794602411D+00, .367436189300487D+00,
     R.500992595876112D+00, .647343511041997D+00, .789962238000080D+00,
     S.908132968654122D+00, .981832813636049D+00, .135519970304965D-34,
     T.302638200562642D-23, .233821706564483D-17, .138495356609889D-13,
     U.524672672275879D-11, .305016274584427D-09, .417887546977044D-08,
     V.168540654673237D-07, .259150580691752D-07, .407715686701032D-07,
     W.874765869233057D-07, .237416073533785D-06, .748353847312217D-06,
     X.252370024383523D-05, .848438921714265D-05, .268781272946349D-04,
     Y.769008757886356D-04, .192697805538511D-03, .413906452945806D-03,
     Z.751070853943530D-03, .114049381694612D-02, .144102476304099D-02,
     1.155201346421425D-02, .168384747898162D-02, .194513172499626D-02,
     2.238290492189850D-02, .307891146512591D-02, .416803158247263D-02,
     3.586749391629736D-02, .852090855044904D-02, .126610215834633D-01,
     4.190928260295221D-01, .289924334316850D-01, .440050296764127D-01,
     5.663066692456672D-01, .985724545393643D-01, .143776561787597D+00,
     6.204752909930370D+00, .283486869453694D+00, .380198344260005D+00,
     7.492402913636816D+00, .614257007125491D+00, .736537139452796D+00,
     8.847509084859788D+00, .934697026076270D+00, .987242841232287D+00/

      DATA POW5/ .151930744614993D-52, .339627324946207D-38,
     A.154093263194951D-30, .218640040176482D-25, .140396338428652D-21,
     B.125453820160385D-18, .277803042542482D-16, .212912395456245D-14,
     C.696002572036868D-13, .111340464586583D-11, .956719438283987D-11,
     D.470655323234815D-10, .138380396276783D-09, .249920016644978D-09,
     E.300920264661209D-09, .368671679879234D-09, .526734326923400D-09,
     F.865699544299634D-09, .160730002549170D-08, .330115909381748D-08,
     G.733458386538048D-08, .172390611668942D-07, .419572465700597D-07,
     H.103686621671965D-06, .255647888944872D-06, .619337965598617D-06,
     I.145516248925973D-05, .327970382977881D-05, .702679809391886D-05,
     J.142055644719830D-04, .269359565158228D-04, .476754244145058D-04,
     K.784685193521444D-04, .119742071285471D-03, .169028385585468D-03,
     L.220340725804951D-03, .264925944648632D-03, .293574079388435D-03,
     M.311231162399070D-03, .359797845824903D-03, .464188764384589D-03,
     N.661736563323015D-03, .102860110126610D-02, .171618059292830D-02,
     O.302151279824974D-02, .551619230922023D-02, .102653431747764D-01,
     P.191605093342910D-01, .353433881495092D-01, .635815004564884D-01,
     Q.110268934545174D+00, .182546534243623D+00, .286073355738856D+00,
     R.421491804294220D+00, .580655439428154D+00, .744745545707856D+00,
     S.886516459045890D+00, .977342822935671D+00, .260018424374088D-43,
     T.399167048277386D-29, .914336448199582D-22, .475109830056505D-17,
     U.794073230594209D-14, .127468761291749D-11, .335988402401208D-10,
     V.192035201601564D-09, .328806377814085D-09, .579357690647708D-09,
     W.150440506618517D-08, .524068006756989D-08, .220106944320099D-07,
     X.100588240640288D-06, .457905362109422D-06, .193530205615417D-05,
     Y.720135312147886D-05, .227036645355585D-04, .590374711789307D-04,
     Z.124337215246936D-03, .209587882013405D-03, .280762631923859D-03,
     1.308048682526873D-03, .341097423379544D-03, .408494554557797D-03,
     2.526481836695614D-03, .725264634869657D-03, .105904338374917D-02,
     3.162392441005445D-02, .258885293642217D-02, .424703594248969D-02,
     4.709721185320500D-02, .119634333891437D-01, .201547723798130D-01,
     5.336470129970309D-01, .552324697975878D-01, .885340063350011D-01,
     6.137732871942590D+00, .206855077224738D+00, .298547044718769D+00,
     7.412477978807069D+00, .543798054146241D+00, .682328188726492D+00,
     8.813168324033214D+00, .919048808024259D+00, .984079062018247D+00/

      DATA POW6/ .414929965177234D-63, .687389904234490D-46,
     A.106008425131758D-36, .161315328990522D-30, .598173181155621D-26,
     B.208056203057309D-22, .135670178833796D-19, .247651598430155D-17,
     C.162604669651716D-15, .452880514194181D-14, .598331017782257D-13,
     D.404802857984243D-12, .147669051019372D-11, .300165828229557D-11,
     E.375095386130614D-11, .478594533072010D-11, .734361577952194D-11,
     F.133303215876668D-10, .280102085202221D-10, .664354293291230D-10,
     G.173161195105685D-09, .482854555069622D-09, .140400410001667D-08,
     H.415783517837642D-08, .122792163960021D-07, .355069398017553D-07,
     I.989677710729378D-07, .262422792356721D-06, .654800450875949D-06,
     J.152387798373667D-05, .328396774413209D-05, .651558600017626D-05,
     K.118477136590371D-04, .196741564003379D-04, .297544256967353D-04,
     L.408990980697855D-04, .510210497485295D-04, .577113520769832D-04,
     M.619012885831296D-04, .736665397229487D-04, .100007725255315D-03,
     N.153046061520330D-03, .259834241026415D-03, .480259079179475D-03,
     O.946824148054822D-03, .194969769443029D-02, .410816773257902D-02,
     P.868737281411651D-02, .181121420794091D-01, .366434981895396D-01,
     Q.709485973594915D-01, .129911713317418D+00, .222727012871244D+00,
     R.354606719838888D+00, .520837443469221D+00, .702117014170032D+00,
     S.865414492465796D+00, .972873365279421D+00, .498890169927277D-52,
     T.526484535442844D-35, .357542143024123D-26, .162986944935743D-20,
     U.120180129224412D-16, .532702234573894D-14, .270140154605558D-12,
     V.218804885537264D-11, .417184610594468D-11, .823258325988259D-11,
     W.258724612238000D-10, .115681837214430D-09, .647381811584633D-09,
     X.400919014840382D-08, .247133076149912D-07, .139347284410777D-06,
     Y.674367960681888D-06, .267494682621157D-05, .842079889887446D-05,
     Z.205836014194811D-04, .385158425536123D-04, .547024988789662D-04,
     1.611425048780650D-04, .690961940962314D-04, .857874039886395D-04,
     2.116321520772029D-03, .170842454078525D-03, .269089345047029D-03,
     3.449447503003946D-03, .786554565952773D-03, .142463340559802D-02,
     4.263818546355520D-02, .493658936196999D-02, .923110045985968D-02,
     5.170740213088949D-01, .309480547501554D-01, .545170240564356D-01,
     6.926499360619841D-01, .150938288803678D+00, .234431157462793D+00,
     7.345526150006201D+00, .481421164533536D+00, .632109003324215D+00,
     8.780219038383975D+00, .903662564410351D+00, .980925421645936D+00/

      DATA POW7/ .113319312979258D-73, .139124518475756D-53,
     A.729284717963178D-43, .119020447244318D-35, .254857896337284D-30,
     B.345046357100039D-26, .662570044458006D-23, .288058917723345D-20,
     C.379887656379856D-18, .184210440380632D-16, .374195393669931D-15,
     D.348164241946662D-14, .157581198028556D-13, .360513438044165D-13,
     E.467554250142899D-13, .621291896251552D-13, .102383098956613D-12,
     F.205264603407383D-12, .488130261247464D-12, .133700501693813D-11,
     G.408813915564571D-11, .135244326297316D-10, .469818130122550D-10,
     H.166729256791079D-09, .589792295653808D-09, .203562972740887D-08,
     I.673094570773892D-08, .209975429253754D-07, .610183507106042D-07,
     J.163471442046341D-06, .400373535581146D-06, .890455857437025D-06,
     K.178884882887367D-05, .323255164963831D-05, .523773474778225D-05,
     L.759158896663825D-05, .982594407993691D-05, .113450075888570D-04,
     M.123116512457023D-04, .150828003494878D-04, .215462886612569D-04,
     N.353964073395944D-04, .656365550519740D-04, .134396568802100D-03,
     O.296697723027696D-03, .689120481407611D-03, .164407967972012D-02,
     P.393885387359633D-02, .928178388888439D-02, .211185006633446D-01,
     Q.456493344026742D-01, .924534301732793D-01, .173407698645768D+00,
     R.298335399345319D+00, .467181781310289D+00, .661928499509844D+00,
     S.844814821121176D+00, .968424346768242D+00, .957206791207947D-61,
     T.694410942127286D-41, .139813286772068D-30, .559128490697140D-24,
     U.181888305812655D-19, .222620560397953D-16, .217197089568507D-14,
     V.249306260184046D-13, .529317589500243D-13, .116983735997580D-12,
     W.444949478583214D-12, .255354024454143D-11, .190408899303560D-10,
     X.159796071029205D-09, .133378558936210D-08, .100334031118867D-07,
     Y.631509299325759D-07, .315162361206169D-06, .120109911009521D-05,
     Z.340754492976715D-05, .707803386991510D-05, .106579830909079D-04,
     1.121357633218836D-04, .139968340753828D-04, .180160998500379D-04,
     2.257002146164074D-04, .402434404109865D-04, .683721523867135D-04,
     3.124391909319048D-03, .238973823702849D-03, .477881602093538D-03,
     4.980669970697145D-03, .203703349498563D-02, .422794235003969D-02,
     5.866413323769103D-02, .173409245744148D-01, .335702182133712D-01,
     6.623236162233495D-01, .110136852005960D+00, .184084781817604D+00,
     7.289441682883027D+00, .426199277275269D+00, .585585937507992D+00,
     8.748604845842409D+00, .888533909392914D+00, .977781887623797D+00/

      DATA POW8/ .309480340582438D-84, .281581552505725D-61,
     A.501711254736203D-49, .878147597682407D-41, .108584853637174D-34,
     B.572234746181605D-30, .323578156663952D-26, .335059174283290D-23,
     C.887518370652502D-21, .749281218371801D-19, .234021283340439D-17,
     D.299450305202669D-16, .168158688639894D-15, .432993788057137D-15,
     E.582803694499634D-15, .806535791101006D-15, .142740296696760D-14,
     F.316073075468577D-14, .850658258304279D-14, .269070650008445D-13,
     G.965163225266680D-13, .378810298123367D-12, .157213982059759D-11,
     H.668584584946454D-11, .283287581873584D-10, .116703619355725D-09,
     I.457781655880068D-09, .168010104969715D-08, .568606682915625D-08,
     J.175361233969556D-07, .488125890639990D-07, .121694600304939D-06,
     K.270092629232463D-06, .531122654255197D-06, .922009571542014D-06,
     L.140913188207834D-05, .189234007410500D-05, .223022321534823D-05,
     M.244868499291812D-05, .308811662985766D-05, .464206693921916D-05,
     N.818646125293500D-05, .165804065779492D-04, .376097787399197D-04,
     O.929734830175904D-04, .243569574530487D-03, .657957067291341D-03,
     P.178587591087773D-02, .475656119426623D-02, .121710833381896D-01,
     Q.293714295837082D-01, .657957356772034D-01, .135009353207663D+00,
     R.250993581122685D+00, .419053621288180D+00, .624040337466094D+00,
     S.824705488756548D+00, .963995673932480D+00, .183656623514585D-69,
     T.915898804397939D-47, .546725904607270D-35, .191809638025004D-27,
     U.275281413033112D-22, .930349277613623D-19, .174630001918491D-16,
     V.284059522776834D-15, .671590234728720D-15, .166232081182094D-14,
     W.765215325975065D-14, .563663919721994D-13, .560033480786997D-12,
     X.636906292073398D-11, .719848603879665D-10, .722433726866600D-09,
     Y.591374469705916D-08, .371324442593579D-07, .171318551790179D-06,
     Z.564107427643464D-06, .130072614649233D-05, .207655236769736D-05,
     1.240874579310231D-05, .283534233247274D-05, .378353742758691D-05,
     2.567823586680811D-05, .947969581008376D-05, .173724872724893D-04,
     3.344274848577865D-04, .726058825251155D-04, .160301467536923D-03,
     4.364536005793596D-03, .840561196290688D-03, .193644263682197D-02,
     5.439657438645431D-02, .971652879391503D-02, .206716997194627D-01,
     6.419237541249542D-01, .803648051526557D-01, .144550780978049D+00,
     7.242460629358025D+00, .377311670802766D+00, .542486957793308D+00,
     8.718271648919875D+00, .873658530555823D+00, .974648427564398D+00/

      DATA POW9/ .845205276037574D-95, .569907961444992D-69,
     A.345152142817402D-55, .647908171385398D-46, .462637046324925D-39,
     B.949010467722712D-34, .158025290074334D-29, .389728084652530D-26,
     C.207347842188915D-23, .304772271888971D-21, .146356587982517D-19,
     D.257552254029892D-18, .179446183420717D-17, .520046136180647D-17,
     E.726461466703836D-17, .104701185747246D-16, .199005426762021D-16,
     F.486699544771933D-16, .148243108421918D-15, .541501443740037D-15,
     G.227864076035848D-14, .106102226905146D-13, .526080935800394D-13,
     H.268102524914475D-12, .136067993148029D-11, .669067394101304D-11,
     I.311344131359357D-10, .134431897447496D-09, .529862829937360D-09,
     J.188115807840051D-08, .595111474506517D-08, .166314541250870D-07,
     K.407804321909287D-07, .872658210719228D-07, .162303303040505D-06,
     L.261559558850688D-06, .364438361029907D-06, .438421530468039D-06,
     M.487023070657195D-06, .632274120099145D-06, .100011588106768D-05,
     N.189336017078885D-05, .418836549347257D-05, .105247884635252D-04,
     O.291342597988703D-04, .860896450440446D-04, .263312969400801D-03,
     P.809715940576738D-03, .243755668798692D-02, .701447853645569D-02,
     Q.188979946165479D-01, .468244263646109D-01, .105113703693088D+00,
     R.211164273173867D+00, .375883530865055D+00, .588320858027968D+00,
     S.805074823714084D+00, .959587253729928D+00, .352376891499202D-78,
     T.120803197214569D-52, .213791708692134D-39, .658005053425396D-31,
     U.416628523878604D-25, .388800466951008D-21, .140405369292177D-18,
     V.323657385982346D-17, .852103637456715D-17, .236213218687954D-16,
     W.131600220540008D-15, .124422168428921D-14, .164717878602082D-13,
     X.253854567430851D-12, .388504731675298D-11, .520172950188806D-10,
     Y.553790361904949D-09, .437494633368464D-08, .244359902865615D-07,
     Z.933860584324769D-07, .239033683543121D-06, .404585905138846D-06,
     1.478095702915167D-06, .574356036444770D-06, .794575717558586D-06,
     2.125455615995214D-05, .223302559955054D-05, .441412627067521D-05,
     3.952836659652126D-05, .220593791218148D-04, .537718137336062D-04,
     4.135505831207905D-03, .346849046149148D-03, .886911356695014D-03,
     5.223102136189880D-02, .544440011821953D-02, .127291150321278D-01,
     6.282012063233124D-01, .586406982730430D-01, .113507092085793D+00,
     7.203105358575625D+00, .334031765220536D+00, .502560052292618D+00,
     8.689167542138357D+00, .859032187679215D+00, .971525009184098D+00/

       DATA POW10/ .230829511592663D-105, .115346719850117D-76,
     A.237447337620684D-61, .478034671683736D-51, .197111318441725D-43,
     B.157386609928342D-37, .771745304458603D-33, .453316881389166D-29,
     C.484419580281936D-26, .123966990543561D-23, .915312083590428D-22,
     D.221516433289268D-20, .191491340737194D-19, .624600147198260D-19,
     E.905530056837723D-19, .135918807544977D-18, .277449051159447D-18,
     F.749435701000595D-18, .258341337194563D-17, .108976513626939D-16,
     G.537961204782997D-16, .297185229915919D-15, .176041053974079D-14,
     H.107509155137451D-13, .653558431219696D-13, .383579515631834D-12,
     I.211749787014856D-11, .107564572110658D-10, .493758914527018D-10,
     J.201798061967667D-09, .725545753422296D-09, .227294609310326D-08,
     K.615730852931786D-08, .143381636357349D-07, .285705951336296D-07,
     L.485500354482528D-07, .701857561479699D-07, .861857400887673D-07,
     M.968648364482765D-07, .129454489860240D-06, .215471208980891D-06,
     N.437895279238555D-06, .105802022552583D-05, .294527582754373D-05,
     O.912953958998696D-05, .304283775923004D-04, .105377270495928D-03,
     P.367125117949451D-03, .124915508588687D-02, .404260720029844D-02,
     Q.121592379257279D-01, .333232371643581D-01, .818379648636902D-01,
     R.177655341087197D+00, .337160739337503D+00, .554645927851692D+00,
     S.785911432159264D+00, .955198993543866D+00, .676095810139834D-87,
     T.159334332430481D-58, .836011140506227D-44, .225729350616321D-34,
     U.630552295546324D-28, .162482851052528D-23, .112888206548116D-20,
     V.368775186541532D-19, .108113634091219D-18, .335655333712645D-18,
     W.226323460316360D-17, .274647275706243D-16, .484470669379313D-15,
     X.101179941551085D-13, .209677320648561D-12, .374539404855457D-11,
     Y.518594867802334D-10, .515456383343175D-09, .348542300320308D-08,
     Z.154597430953629D-07, .439270802868650D-07, .788276554848126D-07,
     1.948939908065422D-07, .116347452236164D-06, .166867801103373D-06,
     2.277183124370387D-06, .526008790592617D-06, .112157288866289D-05,
     3.263713048957069D-05, .670215952642171D-05, .180373142967993D-04,
     4.503704160892735D-04, .143123738456479D-03, .406214849682072D-03,
     5.113212148362237D-02, .305062571994144D-02, .783827027772602D-02,
     6.189703440135540D-01, .427890229736521D-01, .891303379103104D-01,
     7.170138083000765D+00, .295716323693239D+00, .465571757130776D+00,
     8.661242723210986D+00, .844650711530811D+00, .968411600302713D+00/

      EQUIVALENCE (POW1(1), X4NL(1)), (POW2(1), X4NL(LTOT6+1)),
     $         (POW3(1), X4NL(2*LTOT6+1)), (POW4(1) , X4NL(3*LTOT6+1)),
     $         (POW5(1), X4NL(4*LTOT6+1)), (POW6(1) , X4NL(5*LTOT6+1)),
     $         (POW7(1), X4NL(6*LTOT6+1)), (POW8(1) , X4NL(7*LTOT6+1)),
     $         (POW9(1), X4NL(8*LTOT6+1)), (POW10(1), X4NL(9*LTOT6+1))

      END