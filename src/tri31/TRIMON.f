CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction trimon
C
C     Programmeur Lilian Berlu 17/07/2002
C
C  Cette fonction évalue l'intégrale d'attraction coulombienne
C  tricentrique entre deux orbitales de Slater développées sur 
C  une base de fonctions B: 
C
C  Bn1l1m1(zeta1,0) et Bn2l2m2(zeta2, ab).
C
C                   <X1(A) | -ZC/RC | X2(B)>
C
C  L'approche utilisée est basée sur la transformée de Fourier
C  et la transformation non linéaire SDbar
C
C  input:
C  
C       ni, li, mi, zetai : paramètres quantiques des slaters
C       
C       ab, thetab, phiab : paramètres géométriques de X2(B)
C       
C       ac, thetac, phiac : paramètres géométriques de l'opérateur
C
C       numat est le numéro atomique de l'atome C : numat = ZC
C       
C       Ylmab est un tableau contenant les valeurs des harmoniques
C       sphérique en (thetab,0) pour l=0, l1+l2 et m=-l,l
C       
C       Fact[], Dfact[] : tableaux contenant n! et n!! respect.
C       
C       Cnp est un tableau contenant les coefficients binomiaux
C       jusqu'à n=l1+l2.
C       
C       xh[], xbar[] : racines et poids de Gauss_Legendre
C       
C       nrac : ordre de la quadrature de Gauss_Legendre
C       
C       nsd, jc : ordres pour SDbar.
C       
C  ouput :
C       - ovlpr : partie réelle de l'intégrale 
C       - ovlpi : partie imaginaire.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine trimon(lm123, numat, norb1, norb2, NLMA, NLMB,
     $     ZETAA, ZETAB, ab, thetab, phiab, ac, thetac, phiac, 
     $     Fact0,Dfact0,Cnp0, xh,xbar, nrac,nsd,jc, Xcore)
      
      implicit double precision (a-h, o-z)
      implicit integer (i-n)

      INCLUDE "../lib95/SIZE.INCL"
      INCLUDE "../lib95/CONST.INCL"

      DIMENSION XCORE(N_ORB,N_ORB)

      dimension xh(*), xbar(*)
      dimension Fact0(0:*), Dfact0(0:*), Cnp0(0:*)
      dimension NLMA(*), NLMB(*), ZETAA(*), ZETAB(*)

      dimension Ylmab(0:(l1max+l2max+1)*(l1max+l2max+1))
      dimension Cstlmp1(0:(l1max+1)*(l1max+1)),
     $          Cstlmp2(0:(l2max+1)*(l2max+1))

      dimension V(256), Thetav(256), Phiv(256)
      dimension Ylmv(0:256*(l1max+l2max+1)*(l1max+l2max+1))

      dimension argnt1(l1max+l2max+1), argnt2(l1max+l2max+1)


c.... évaluation des harmoniques sphériques 
      call harmonique (thetab, 2*lm123, Ylmab)

c.... tableaux des coordonnées de \vect{v}
      call vcoord(ab, thetab, phiab, ac, thetac, phiac, 
     $     nrac, xbar, V, Thetav, Phiv)

c.... tableau des harmoniques sphériques Ylmv.
c.... la valeur maxi pour lambda est l1+l2
      call vharmoniques(2*lm123, nrac, Thetav, Ylmv)


C.... on combine les orbitales sur les deux atomes
      DO 5 NAT2=1, NORB2
       JST   = 5*(NAT2-1)
       ICOL  = NLMB(JST+2)
       N2    = NLMB(JST+3)
       L2    = NLMB(JST+4)
       M2    = NLMB(JST+5)
       ZETA2 = ZETAB(NAT2)

      DO 5 NAT1=1, NORB1
       IST   = 5*(NAT1-1)
       LINE  = NLMA(IST+2)
       N1    = NLMA(IST+3)
       L1    = NLMA(IST+4)
       M1    = NLMA(IST+5)
       ZETA1 = ZETAA(NAT1)

c       write(*,*) n1,l1,m1, zeta1
c       write(*,*) n2,l2,m2, zeta2

c.... évaluation des facteurs du théorème de multiplication des
c.... harmoniques sphériques solides. Ces coefficients sont communs
c.... à toutes les fonctions B presentes dans les Slater.
      call Constlmp(l1, m1, Fact0, Dfact0, Cstlmp1)
      call Constlmp(l2, m2, Fact0, Dfact0, Cstlmp2)

c.... constante de normation
      cnorm1 = dble(2**n1) * dsqrt((2.0d0*zeta1)/Fact0(2*n1))*zeta1
      cnorm2 = dble(2**n2) * dsqrt((2.0d0*zeta2)/Fact0(2*n2))*zeta2
      cnormp = cnorm1 * cnorm2

c.... initialisation de la série
       totalnucs_r = 0.000d0 
       totalnucs_i = 0.000d0 
       
       do 10 lp1 = 0, l1
       do 10 mp1 = -lp1, lp1
       do 10 lp2 = 0, l2
       do 10 mp2 = -lp2, lp2
          cstmp1 = Cstlmp1(lp1*(lp1+1)+mp1)
          cstmp2 = Cstlmp2(lp2*(lp2+1)+mp2)
          cstmprod = cstmp1 * cstmp2
        if(cstmprod.NE.0.0d0) then
c... évaluation des : <lp2 mp2|lp1 mp1|l m>
             call gaunt(lp1,mp1, lp2,mp2, lmin,lmax,m, 
     $            ngt, argnt1)

c.... évaluation des : <l2-lp2 m2-mp2|l1-lp1 m1-mp1|lambda mu>
       if(iabs(m1-mp1).gt.iabs(l1-lp1))then
c            
	     if(ngt2.eq.0)ngt2=1
	     argnt2(ngt2)=0.0d0
	     go to 1234
	     endif
	       if(iabs(m2-mp2).gt.iabs(l2-lp2))then
	      if(ngt2.eq.0)ngt2=1
	       argnt2(ngt2)=0.0d0
	       go to 1234
	       endif

             call gaunt(l1-lp1,m1-mp1, l2-lp2,m2-mp2, 
     $            lambda_min,lambda_max,mu, ngt2, argnt2)
 1234       continue
            do 15 l = lmin, lmax, 2
                ngt = (l-lmin)/2 + 1
                nylm = l*(l+1) + m
                
                cste_l = ab**l * argnt1(ngt) * Ylmab(nylm)
c                write(*,*) 'cste_l =', cste_l
                
                do  15 lambda = lambda_min, lambda_max, 2
                   ngt2 = (lambda - lambda_min)/2 + 1
                   
                   cste_lb = dble((-1)**( lambda+l1+l2+lp2 
     $                  + (lambda+l1+lp1+l2+lp2)/2 ))
     $                  * argnt2(ngt2)
c                write(*,*) 'cste_lb =', cste_lb
                   
C.... Décomposition des Slater en fonctions B. 
                   np1min = (n1-l1+1)/2
                   np1max = n1-l1
                   
                   np2min = (n2-l2+1)/2
                   np2max = n2-l2
                   
                   do 15 np1 = np1min, np1max   
                   do 15 np2 = np2min, np2max
                      Ap1 = coeffbons(n1,l1,np1, Fact0, Dfact0)
                      Ap2 = coeffbons(n2,l2,np2, Fact0, Dfact0)
                         
                         cste_Ap = Ap1 * Ap2
c                         write(*,*) 'cste_Ap = ', cste_Ap

c.... constante d'intégration sur fonctions B
                         cste = 128.0d0 * pi*pi 
     $                        * zeta1**(2*np1+l1-1)*zeta2**(2*np2+l2-1)
     $                        * Dfact0(2*l1+1) * Dfact0(2*l2+1) 
     $                        * Fact0(np1+l1+np2+l2+1)
     $                        / ( dble(2**(np1+np2+l1+l2+1)) 
     $                        *   Fact0(np1+l1) * Fact0(np2+l2) )

c                         write(*,*) 'cste =', cste

                         ld = (lp1+lp2-l)/2
                         do 15 j = 0, ld
                            ncnp = ld*(ld+1)/2 + j                           
                            nu = np1+np2+l1+l2 - l - j
                            nx = l1-lp1 + l2-lp2
                            ng = 2*(np1+l1+np2+l2) - (lp1+lp2+l) + 1
                              
                            cste_j = dble((-2)**j) * Cnp0(ncnp) 
     $                           / Fact0(nu+l+1)
c                            write(*,*) 'cste_j =', cste_j

                              
c.... integrale externe sur s
                            do 15 i = 1, nrac 
                               s = 0.5d0 * (xbar(i)+1.d0)
                               a = (1.d0-s)*zeta1**2 + s*zeta2**2 
                               b = s*(1.d0-s) 

c                               write(*,*) i, s, a, b

                                 
c.... v = V[i]
c.... tethav = Tethav[i]
c.... phiv = Phiv[i]
                                 
c.... pour les molécules non linéaires, on distingue les 
c.... parties réelles et imaginaires
                               cosin = dcos( m*phiab + mu*Phiv(i) )
                               sine  = dsin( m*phiab + mu*Phiv(i) )
                                 
c.... Évaluation de l'intégale interne semi-infinie : SDbar
                               nk = 0
                               zetas = 0.00d0
                               call sdatt(nk, nx,nu,ng,lambda,
     $                              ab, V(i), a, b, zetas, 
     $                              nsd,jc,nrac, xbar, xh, Cnp0, sd) 
                                 
                               nylmv = nrac*(lambda*(lambda+1)+mu)+i-1
                               term = 0.5d0*xh(i) * s**(np2+l1+l2-lp1)
     $                              * (1.0d0-s)**(np1+l1+l2-lp2)
     $                                * Ylmv(nylmv) * sd                                 

c                               write(*,*) 'term(i) = ',
c     $                              ' V(i) =', V(i)

                               term = term * cstmprod * cste * cste_j 
     $                              * cste_Ap * cste_lb * cste_l 
                               
                               totalnucs_r = totalnucs_r 
     $                              + cosin * term
                               
                               totalnucs_i = totalnucs_i
     $                                + sine  * term
 15          continue 
          end if
 10    continue
       ovlp_r = numat * cnormp * totalnucs_r
       ovlp_i = numat * cnormp * totalnucs_i

c       write(*,*) Line, Icol, ovlp_r

c.... On place les valeurs dans le tableau Xcore.
       Xcore(Line, Icol) = Xcore(Line, Icol) - ovlp_r
       Xcore(Icol, Line) = Xcore(Icol, Line) - ovlp_i

 5    CONTINUE

      return
      end

c**********************************************************************c

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction harmonique
C
C     Programmeur Lilian Berlu 17/07/2002
C     
C     cette routine évalue les harmoniques sphériques au point cos(tetha)
C     pour toutes les combinaisons de l, m avec l <= lmax
C     le résultat est placé dans le tableau ylmab tel que
C     ylmab(l*(l+1) + m) = Y(l,m) m pouvant être négatif (m = -l,.,-1,0,1,.,l)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine harmonique (theta, lmax, Ylm)

      implicit doubleprecision (a-h, o-z)
      implicit integer (i-n)

      INCLUDE "../lib95/SIZE.INCL"

      dimension Ylm(0:*)
      dimension Plm(0:(8*l1max+1)*(8*l1max+2)/2)

      pi = dacos(-1.0d0)
      Eps = 1.0d-15

      costheta = dcos(theta)
      if(dabs(costheta).LT.Eps) then
         costheta = 0.0d0
      endif

      call legendre(costheta, lmax, Plm)

c.... Ylm[l*(l+1)+m] = ((double) pow(-1,m)) * Plm[l*(l+1)/2 + m] 
c....                * sqrt((2.*l+1.)/(4.*pi) * fact[l-m]/fact[l+m])
      do  l = 0, lmax
         poch = 1.0d0

         do m = 0, l
c            if(m.NE.0) then
c            poch = poch * dble((l+m)*(l+1-m))
c            endif

            Ylm(l*(l+1)+m) = (-1.0d0)**m * Plm(l*(l+1)/2 + m)
     $           * dsqrt((2.0d0*l + 1.0d0)/(4.0d0*pi*poch))

            Ylm(l*(l+1)-m) = (-1.0d0)**m * Ylm(l*(l+1)+m)

            poch = poch * dble((l+m+1)*(l-m))
         end do
      end do

      return
      end

c**********************************************************************c

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction Legendre
C
C     Programmeur Lilian Berlu 17/07/2002
C
C  fonction de legendre(x,l,Plm).
C  cette fonction évalue le polynome associé de Legendre au point x
C  le programme est basé sur les 3 relations suivantes
C  a :               P(l=m,m) =   (2m-1)!! * (1-x²)**m/2
C  b :             P(l=m+1,m) = x * (2m+1) * P(l=m,m)
C  c :          (l -m) P(l,m) = x * (2l-1) * P(l-1,m) - (l+m-1) * P(l-2,m)
C
C  Le résultats est placé dans Plm tel que :
C  P(l,m) = Plm (l*(l+1)/2 + m)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine legendre(x, lmax, Plm)
      implicit doubleprecision (a-h, o-z)
      implicit integer (i-n)

      dimension Plm(0:*)

      do 5 l = 0, lmax

         if (l.EQ.0) then
            dfact = 1.0d0
         else
            dfact = dfact * dble(2*l - 1)
         endif
c     évaluation des P(l=m,m)
         Plm (l*(l+1)/2 + l) = dfact * dsqrt((1.0d0 - x*x)**l)

c     évaluation des P(l=m+1,m)
         if (l.LT.lmax) then
            Plm ((l+1)*(l+2)/2+l) = dble(2*l+1) * x * Plm(l*(l+1)/2+l)
         endif

c     application de la troisième relation : évaluation des P(l+2, m)
         do 10 lj = l+2, lmax
            mj = l

            Plm(lj*(lj+1)/2+mj)=(dble(2*lj-1)* x * Plm(lj*(lj-1)/2+mj)
     $           - dble(lj+mj-1)*Plm((lj-1)*(lj-2)/2+mj)) / dble(lj-mj)
 10      continue
 5    continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Fonction vcoord
C
C     Cette fonction évalue les coordonnées du vecteur V
C     apparaissant dans l'intégrale externe sur s, pour 
C     l'approche transformée de Fourier de l'intégrale
C     tricentrique d'attraction coulombienne.
C
C          \vec{V} = (1-s) \vec{ab} - \vec{ac}
C input:
C 
C       - ab, thetab, phiab : les coordonnées de \vec{ab}
C       - ac, thetac, phiac : les coordonnées de \vec{ac}
C       - nrac : ordre de la quadrature de Gauss-Legendre
C       - xbar[] : le tableau des racines de Gauss-Legendre
C
C  ouput:
C       - V[], Thetav, Phiv[] : tableaux des coordonnées 
C       de \vec{v} pour chaque point de quadrature.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine vcoord(ab, thetab, phiab, ac, thetac, phiac, 
     $     nrac, xbar, V, Thetav, Phiv)
      
      implicit double precision (a-h, o-z)
      implicit integer (i-n)

      dimension xbar(*), V(*), Thetav(*), Phiv(*)
      
      pi = dacos(-1.0d0)
      Eps = 1.0d-15

c      write(*,*) '------------------------------'
      do i = 1, nrac

         s = 0.5d0 * (xbar(i)+1.0d0) 
      
c.... évaluation des coordonnées sphériques de
c.... \vec{v} = (1-s) \vec{ab} - \vec{ac}

         vx = (1.0d0 - s) * ab * dsin(thetab) * dcos(phiab) 
     $        - ac * dsin(thetac) * dcos(phiac)  
c         if(dabs(vx) .LT. Eps) vx = 0.0d0 
         
         vy = (1.0d0 - s) * ab * dsin(thetab) * dsin(phiab) 
     $        - ac * dsin(thetac) * dsin(phiac)   
c         if(dabs(vy) .LT. Eps) vy = 0.0d0 
         
         vz = (1.0d0 - s) * ab * dcos(thetab) - ac * dcos(thetac) 
c         if(dabs(vz) .LT. Eps) vz = 0.0d0 
         
c         write(*,*) ab, thetab, phiab
c         write(*,*) ac, thetac, phiac
c         write(*,*) s, vx, vy, vz

         vnorm = dsqrt(vx*vx + vy*vy + vz*vz) 
         V(i) = vnorm
c         if(dabs(vnorm) .LT. Eps) V(i) = 0.0d0 
         
         Thetav(i) = 0.0d0
         Phiv(i)   = 0.0d0
         
         if(vnorm .NE. 0.0d0) then
            Thetav(i) = dacos(vz / vnorm) 
            if(dabs(1-dcos(Thetav(i))) .LT. Eps) Thetav(i) = 0.0d0 
            if(dabs(1+dcos(Thetav(i))) .LT. Eps) Thetav(i) = pi 
            
            if(vx .GT. 0.0d0) then
               Phiv(i) = datan(vy / vx)
            endif
            
            if(vx .LT. 0.0d0) then
               Phiv(i) = datan(vy / vx) + pi
            endif

            if((vx .EQ. 0.0d0).AND.(vy.NE.0.000d0)) then
               Phiv(i) = dabs(vy)/vy * pi/2.0d0
            endif
         endif
         
      end do
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Fonction vharmoniques
C
C Cette fonction évalue les Ylmv(thetav)
C apparaissant dans l'intégrale externe sur s, pour 
C l'approche transformée de Fourier de l'intégrale
C tricentrique d'attraction coulombienne.
C
C input:
C       - lmax : la plus haute valeur pour lambda dans 
C       l'intégrale.
C       - nrac : ordre de la quadrature de Gauss-Legendre
C       - Thetav[] : tableaux des thetav pour chaque point
C       de quadrature
C
C  ouput:
C       - Ylmv[] : tableaux des Ylmv(thetav)
C       Ylmv(nrac*(l*(l+1)+m) + i-1) == Ylm(Thetav[i])
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine vharmoniques(lmax, nrac, Thetav, Ylmv)

      implicit double precision (a-h, o-z)
      implicit integer (i-n)

      INCLUDE "../lib95/SIZE.INCL"
      
      dimension Thetav(*), Ylmv(0:*)
      dimension Plm(0:(8*l1max+1)*(8*l1max+2)/2)
  
      pi = dacos(-1.0d0)
      Eps = 1.0d-15

      do i=1, nrac         
         dcosthetav = dcos(Thetav(i)) 
      if(dabs(dcosthetav).LT.Eps) dcosthetav = 0.0000d0 
      
      call legendre(dcosthetav, lmax, Plm) 
  
  
      do l = 0, lmax
	  poch = 1.0 
	  do m = 0, l
	      if(m.NE.0) poch = poch * dsqrt(dble((l+m)*(l+1-m))) 
	      
	      nylm1 = nrac*(l*(l+1) + m) + i - 1 
	      nylm2 = nrac*(l*(l+1) - m) + i - 1 

	      Ylmv(nylm1) = dble((-1)**m) * Plm(l*(l+1)/2 + m) 
     $             * dsqrt(dble(2*l+1)/(4.0d0*pi)) / poch 
	      
	      Ylmv(nylm2) = dble((-1)**m) * Ylmv(nylm1) 
	      
           end do	    
        end do
      end do

      return 
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction Constlmp
C
C     Programmeur Lilian Berlu 17/07/2002
C     
C  Cette fonction évalue les coefficients intervenant dans les 
C  sommes sur l' et m' de l'équation (43) de l'intégrale tricentrique
C  d'attraction nucléaire.
C
C  Ce jeu de coefficients est commun à toutes les fonctions B impliquées
C  dans la décomposition des Slaters de base. Par conséquent cette 
C  fonction peut être appellée avant la fonction faisant le calcul des
C  intégrales.
C
C  input :
C        - l, m : nombres quantiques secondaires et magnétiques
C
C  ouput :
C        - Cstlmp[] : tableau contenant les coefficients.
C	Cstlmp[lp*(lp+1)+mp] contient le terme telque (lp, mp).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine Constlmp(l, m, Fact0, Dfact0, Cstlmp)
      implicit none

      integer l,lp,m,mp,mmp,lmin,lmax,ngt
      dimension argntl(7)
c     dimension argntl(l1max+l2max+1)

      double precision Fact0(0:*), Dfact0(0:*), Cstlmp(0:*)
c      double precision gaunt1,mmp, lmin,lmax,ngt,argntl
      double precision gaunt1,argntl
c      integer lmin,lmax,ngt
    
      do lp = 0, l         
         do mp = -lp, lp 

            call gaunt(lp,mp,l,m, lmin,lmax,mmp, ngt, argntl)
            gaunt1 = argntl((l-lp - lmin)/2+1)

            Cstlmp(lp*(lp+1)+mp) = gaunt1 
     $           / (Dfact0(2*lp+1) * Dfact0(2*(l-lp)+1))
            
         end do
      end do
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction coeffbons
C
C     Programmeur Lilian Berlu 17/07/2002
C     
C     cette fonction évalue les coefficients qui apparaissent
C     dans le développement de l'OA de Slater |nlm> en fonction B |plm>.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function coeffbons(n, l, p, Fact0, Dfact0)
      implicit none
      integer n,l,p, pmax

      double precision coeffbons, Fact0(0:*), Dfact0(0:*)

      pmax = n - l
      coeffbons = (-1.0d0)**(pmax-p) * 2.0d0**(l+p) 
     $     * Fact0(pmax) * Fact0(l+p) 
     $     / (Fact0(2*p-n+l) * Dfact0(2*(pmax-p)))
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction sdatt
C
C     Programmeur Lilian Berlu 17/07/2002
C     
C
C    Cette routine évalue l'intégrale semi infinie sur x contenue 
C dans l'intégrale externe de eq(49) (sur s) au moyen de SDbar.
C
C    On utilise une quadrature de Gauss_Legendre pour évaluer la 
C fonction F(x_{i+j}) impliquée dans l'expression.
C
C F(x) = \int_{0}^{x} G(t) sint(vt) dt
C      = \int_{0}^{x} [d/(xdx)]^\lambda 
C                  [x^{nx+\lambda-1} \frac{k_\nu(R2.\gamma(s,x))}
C                     {\gamma(s,x)^{n_\gamma}}]
C      * sin(vx) dx
C
C -input:
C        - nx : tel qu'il est définie dans (45) 
C          nx = l1+l2 - (l'1+l'2)
C        - nu : nu == n1 + n2 + l1 + l2 - l - j 
C          (ici n=n1+n2 car les l sont pris nuls)
C        - ng : n^{\gamma} = 2(n1+l1+n2+l2)-(l'1+l'2)-l+1
C        - lambda : lambda == ordre de la bessel sphérique
C        - ab : position de la fonction B translatée
C        - v : norme de vecteur v : v = (1-s)*ab - ac
C        - a = (1-s)*zeta1*zeta1 + s * zeta2*zeta2
C        - b = s*(1-s)
C        - nsd: ordre de SDbar.
C        - jc : le j apparaissant dans l'équation (50)
C        - nrac ordre de la quadrature de Gaus_Legendre 
C          (limité à 25)
C        - xbar[], xh[]: racines et poids de la quadrature de GL
C        - Cnp[] : coefficients du binome
C
C -output:
C        - sd : valeur de l'approximation SDbar
C        - t : temps pour sdatt.
C
C  voir J. PHYS. A: Math.&Gen. 34, 2801-2818, 2001 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sdatt(nk, nx, nu, ng, lambda, r, v, a, b, zetas, 
     $     nsd, jc, nrac, xbar, xh, Cnp0, sd)
      implicit double precision (a-h,o-z) 
      implicit integer (i-n)

      dimension vcst(0:100), besk(0:100)
      dimension xbar(*), xh(*), Cnp0(0:*)
      
      pi = dacos(-1.0d0)
      Eps = 1.0d-15

      terme = 0.d0 
      do 5 i = 0, nsd+jc+1
         vcst(i) = terme 

c.... calcul de F[xi]    
         Ui = 0.000d0
         do 10 j = 1, nrac 
            xij = 0.50d0 * (xbar(j)+dble(2*i+1)) * pi/v

            z = G(a,b,xij)
            call bessel_red(nu+lambda, z*r, besk)

            Ui = Ui + 0.50d0 * pi/v * xh(j) * dsin(v*xij)
     $           * CALF(nk, nx, nu, ng, lambda, r,a,b,zetas, xij, 
     $                  besk, Cnp0) 

 10      continue 

         vcst(i) = vcst(i) + Ui
         terme = vcst(i)
         
 5    continue 
      
      xnum = 0.0d0 
      deno = 0.0d0 

      do 15 i = nsd+1, 0, -1 
         
       dijc = dble((1+i+jc)**nsd)
       xijc = dble(1+i+jc) * pi/v
       dxijc = xijc * xijc
       Gxij = CALF(nk, nx, nu, ng, lambda, r,a,b,zetas, xijc, besk,Cnp0)
       binom = Cnp0((nsd+1)*(nsd+2)/2+i)

       xnum = xnum + binom * dijc * vcst(i+jc)
       deno = deno + binom * dijc 
 15   continue 
      sd = 1.d0/v**(lambda+1) * xnum/deno  

      return 
      end 
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction bessel_red
C
C     Programmeur Lilian Berlu 17/07/2002
C     
C     cette routine évalue les bessels réduites.
C     pour les évaluer on utilise la relation de recurence suivante:
C     k(n+1/2)[z] = z**2 k(n-2+1/2)[z] + (2n-1) k[z](n-1+1/2)
C     Le résultat est rangé dans le tableau kred: kred[i] = k(i+1/2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bessel_red(nmax, z, bred)
      implicit double precision (a-h, o-z)
      implicit integer (i-n)

      dimension bred(0:*)

c.... initialisation de la série
      bred(0) = dexp(-z)
      
      if(nmax.GT.0) then
         bred(1) = (z+1) * dexp(-z)

         do i = 2, nmax
            bred(i) = (2*i - 1) * bred(i-1) + z*z * bred(i-2)
         end do
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction CALF
C
C     Programmeur Lilian Berlu 17/07/2002
C     
C Cette fonction évalue l'intégrande de l'eq(49)
C
C -input:
C        - nx : tel qu'il est définie dans (45) 
C          nx = l1+l2 - (l'1+l'2)
C        - n : nu == n1 + n2 + l1 + l2 - l - j 
C          (ici n=n1+n2 car les l sont pris nuls)
C        - r2 : position de la fonction B translatée
C        - v : norme de vecteur v : v = (1-s)r2 - r1
C        - a = (1-s)*zeta1*zeta1 + s * zeta2*zeta2
C        - b = s*(1-s)
C        - ng : n^{\gamma} = 2(n1+l1+n2+l2)-(l'1+l'2)-l+1
C        - nj : lambda == ordre de la bessel sphérique
C        - x : point auquel on évalue l'intégrande.
C        - Cnp[] : coefficents du binome.
C
C  voir J. PHYS. A: Math.&Gen. 34, 2801-2818, 2001
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION CALF(nk, nx, nu, ng, lambda, r,a,b,zetas, x, besk, Cnp0)
      implicit double precision (a-h, o-z)
      implicit integer (i-n)

      dimension Cnp0(0:*), besk(0:*)

      if(lambda.EQ.0) then
c.... évaluation de gamma(s,x)
         z = G(a, b, x)
         
         xpnk = (zetas*zetas + x*x)**nk
         xpnx = x**(nx-1)
         zpng = z**ng

         CALF = xpnx/xpnk * besk(nu) / zpng 
      else
         CALF = SinF(nk, nx, nu, ng, lambda, r,a,b,zetas, x, besk,Cnp0)
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction SinF
C
C     Programmeur Lilian Berlu 17/07/2002
C     
C 
C  Cette fonction évalue la dérivée 
C                 [d/(x dx)]^nj (x^{nj-1} g(x))
C 
C
C -input:
C        - nk : nk = l1+l2+k
C        - nx : nx = l3-l'3 + l4-l'4 + l +2r
C        - n : nu == n1 + n2 + l1 + l2 - l - j 
C          (ici n=n1+n2 car les l sont pris nuls)
C        - r2 : distance B3 B4 == R34
C        - v : norme de vecteur v : v = (1-s)R34 - R14
C        - a = (1-s)*zeta3*zeta3 + s * zeta3*zeta3
C        - b = s*(1-s)
C        - ng : n_{\gamma} = 2(n3+n4+l3+l4) - (l'3+l'4+l') + 1
C        - nj : lambda == ordre de la bessel sphérique
C        - zetas == zeta1 + zeta2
C        - x : point auquel on évalue l'intégrande.
C        - Cnp[] : coefficents du binome.
C************************************************************
      Function SinF(nk, nx, n, ng, nj, r, a, b, zetas, x, besk, Cnp0)
      implicit double precision (a-h, o-z)
      implicit integer (i-n)

      dimension Cnp0(0:*), besk(0:*)
c, besk(0:n+nj)

      zx = (zetas*zetas + x*x)
      xpnk = zx**nk

      z = G(a,b,x)

c.... évaluation des bessels réduites nécessaires
c      call bessel_red(n+nj, z*r, besk)

      
c.... initialisation de SinF
      SinF = 0.0000d00
      
      if(nk.EQ.0) then
         
         if((2*n+1).EQ.ng) then
            do l = 0, nj
               xmbp = (-b)**(nj-l)
               xpl1 = x**(nx+nj-1-2*l)
               zpl  = z**(ng+2*(nj-l))
               binom= Cnp0(nj*(nj+1)/2 + l)
               
               SinF = SinF + xmbp * binom * dobin(nx+nj-1,l) * xpl1
     $              * besk(n+nj-l) / zpl
c     $              * hatk(n+nj-l, z*r) / zpl
            end do

         else
            do l = 0, nj
               term = 0.0d0
               xpl2 = x**(nx+nj-1-2*l)
               binom= Cnp0(nj*(nj+1)/2 + l)
               
               do i = 0, nj-l
                  mup = dble((-1)**(nj-l-i))
                  zpi = z**(ng+2*(nj-l))
                  bpi = b**(nj-l)
                  bino= Cnp0((nj-l)*(nj-l+1)/2 + i)
                  
                  term = term +  mup * bino * dobin(2*n+1-ng,i) * bpi
     $                 * besk(n+nj-l-i) / zpi
c     $                 * hatk(n+nj-l-i, z*r) / zpi
               end do
               
               SinF = SinF + binom * dobin(nx+nj-1,l) * xpl2 * term
            end do
         endif         
         
C..... NK != 0
      else
         if((2*n+1) .EQ. ng) then
            do 5 j = 0, nj
               
c     // dérivée : [d/(x dx)]^{nj-j} (zetas^2 +x^2)^{-nk}
               cnp1 = Cnp0(nj*(nj+1)/2 + j)
               dxpnk = pochammer(dble(nk), nj-j) 
     $               * (-2.0d0/zx)**(nj-j) / xpnk
               
               terme = 0.0000d00
c     // dérivée : [d/(x dx)]^i ( x^{nx+nj-1} * k_n(g(x)) / g(x)^{ng} )
               do 10 l=0, j
                  
                  bp = (-b)**(j-l)
                  xpl1= x**(nx+nj-1-2*l)
                  zpl = z**(ng+2*(j-l))
                  cnp2 = Cnp0(j*(j+1)/2 + l)
                  
                  terme = terme + bp * cnp2 * dobin(nx+nj-1,l) * xpl1
     $                 * besk(n+j-l) / zpl
c     $                 * hatk(n+j-l, r*z) / zpl
 10            continue	    
               SinF = SinF + cnp1 * dxpnk * terme
 5          continue
            
         else
            do 15 j=0, nj
c     // dérivée : [d/(x dx)]^{nj-j} (zetas^2 +x^2)^{-nk}
               cnp1 = Cnp0(nj*(nj+1)/2 + j)
               dxpnk = pochammer(dble(nk), nj-j) 
     $               * (-2.0d0/zx)**(nj-j)/xpnk
               
               terme = 0.0000d00
               
c     // dérivée : [d/(x dx)]^j ( x^{nx+nj-1} * k_n(g(x)) / g(x)^{ng} )
               do 20 l=0, j
                  xpl2 = x**(nx+nj-1-2*l)
                  cnp2 = Cnp0(j*(j+1)/2 + l)
                  
                  term = 0.0000d00
                  do 25 i=0, j-l
                     up = dble((-1)**(j-l-i))
                     zpi = z**(ng+2*(j-l))
                     bpi = b**(j-l)
                     cnp3 = Cnp0((j-l)*(j-l+1)/2 + i)
                     
                     term = term + up * cnp3 * dobin(2*n+1-ng,i) * bpi 
     $                    * besk(n+j-l-i) / zpi 
c     $                    * hatk(n+j-l-i,r*z) / zpi 
 25               continue
                  terme = terme + cnp2 * dobin(nx+nj-1,l) * xpl2 * term
 20            continue
               SinF = SinF + cnp1 * dxpnk * terme
 15         continue
         end if
      end if

      return 
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction G
C
C     Programmeur Lilian Berlu 17/07/2002
C     
C Fonction G : Cette fonction évalue Gamma(s,x)
C 
C -input:
C 
C        - a = (1-s)*zeta1*zeta1 + s * zeta2*zeta2
C        - b = s*(1-s)
C        - x = point auquel on évalue l'intégrale
C 
C  voir J. PHYS. A: Math.&Gen. 34, 2801-2818, 2001 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function G(a, b, x)
      implicit double precision (a-h, o-z)
      implicit integer (i-n)

      G = dsqrt(a + b * x*x)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction dobin
C
C     Programmeur Lilian Berlu 17/07/2002
C     
C     Fonction dobin : dobin(n,m) = (n)!! / (n-2m)!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function dobin(n, m)
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
      
      dobin = 1.0d0

      do i=0, m-1
         dobin = dobin * dble(n-2*i)
      end do
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         Fonction pochammer
C
C     Programmeur Lilian Berlu 17/07/2002
C     cette fonction évalue la valeur du symbol de pochammer (a)n
C     (a)n = a(a+1)(a+2)....(a+n-1) et (a)0 = 1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function pochammer(a, n)
      implicit none
      double precision a, pochammer
      integer i, n

      pochammer = 1.0d0
      
      do i = 1, n
         pochammer = pochammer * (a + dble(i-1))
      end do

      return
      end
