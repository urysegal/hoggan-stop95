      subroutine eign(z,d,maxn,n,eps,fail)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  ******************************************************************  c
c  *                                                                *  c
c  *  this subroutine reduces the given lower triangle of a real    *  c
c  *  symmetric matrix a, stored row by row in the array z(n,n), to *  c
c  *  tridiagonal form using householders reduction. the diagonal of*  c
c  *  the result is stored in the array d, the subdiagonal in the   *  c
c  *  array e. the subroutine then finds the eigenvalues of a by    *  c
c  *  transforming the tridiagonal matrix using ql-transformations. *  c
c  *  the procedure will fail if any one eigenvalue takes more than *  c
c  *  30 iterations.                                                *  c
c  *                                                                *  c
c  *  parameters                                                    *  c
c  *                                                                *  c
c  *  a=z  original matrix(real,symmetric), destroyed in computation*  c
c  *       resultant eigen vectors are stored in the array z(n,n)   *  c
c  *       and the eigenvalues in the store d(n) in ascending order *  c
c  *  n    order of matrix a                                        *  c
c  *  eps  treshold of diagonalization                              *  c
c  *  fail=.true. if one eigenvalue takes more than 30 iterations   *  c
c  *  fail=.false. otherwise                                        *  c
c  *                                                                *  c
c  *  method                                                        *  c
c  *  householder-transformation and ql-transformation as found in  *  c
c  *  numerische mathematik, vol.17 and 18(1968)                    *  c
c  *                                                                *  c
c  ******************************************************************  c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8(a-h,o-z)
      include '../lib95/SIZE.INCL'
      dimension z(maxn,maxn),d(*),iord(N_ORB)
      dimension e(N_ORB)
      logical fail
c
c ******householders transformation
c
      incx = 1
      incy = 1
      i1 = n
      i2 = n-1
   10 i = i1
      i1 = i2
      i2 = i2-1
      if (i2) 120,120,20
   20 g = 0.0d0
      do 30 k=1,i2
      s=z(i,k)
   30 g = g+s*s
c
c *   if g is too small for orthogonality to be guaranteed,the trans=
c *   formation is skipped.
c
      if (g.gt.1d-20) go to 40
      e(i)=z(i,i1)
      d(i) = 0.0d0
      go to 10
   40 f = z(i,i1)
      h=g+f*f
      g=dsqrt(h)
      if (f) 60,50,50
   50 g = -g
   60 e(i) = g
      h=h-f*g
      z(i,i1)=f-g
      f=0.0d0
      do 90 j=1,i2
      g = 0.0d0
      s=z(i,j)/h
      z(j,i)=s
c
c *   calculate elements of a*u
c
      do 70 k=1,j
   70 g = g+z(j,k)*z(i,k)
      j1 = j+1
      do 80 k=j1,i1
   80 g = g+z(k,j)*z(i,k)
c
c *   calculate elements of p
c
      e(j)=g/h
   90 f = f+g*s
      j=i1
      g = 0.0d0
      s=z(i,j)/h
      z(j,i)=s
      do 100 k=1,j
  100 g = g+z(j,k)*z(i,k)
      e(j)=g/h
      f=f+g*s
c
c *   calculate skalar
c
      hh=f/(h+h)
c
c *   calculate transformed a
c
      do 110 j=1,i1
      f=z(i,j)
      g=e(j)-hh*f
      e(j)=g
      do 110 k=1,j
  110 z(j,k) = z(j,k)-f*e(k)-g*z(i,k)
c
      d(i)=h
      go to 10
c
c *******
c
  120 e(2) = z(2,1)
      d(2)=0.0d0
      d(1)=z(1,1)
      z(1,1)=1.0d0
      e(1)=0.0d0
      do 160 i=2,n
      l=i-1
      if (dabs(d(i)).lt.1.d-20) go to 150
      do 140 j=1,l
      g=0.0d0
      do 130 k=1,l
  130 g=g+z(i,k)*z(k,j)
      do 140 k=1,l
  140 z(k,j)=z(k,j)-g*z(k,i)
  150 d(i)=z(i,i)
      z(i,i)=1.0d0
      do 160 j=1,l
      z(i,j)=0.0d0
  160 z(j,i)=0.0d0
c
c******
c
      do 170 k=2,n
  170 e(k-1)=e(k)
c
c ********ql transformation
c
      fail=.false.
      e(n)=0.0d0
      b=0.0d0
      f=0.0d0
      do 310 l=1,n
      j=0
      h=eps*(dabs(d(l))+dabs(e(l)))
      if(b.lt.h)b=h
      do 180 m=l,n
c
c *   look for small subdiagonal element
c
      if(dabs(e(m)).le.b) go to 190
  180 continue
  190 if (m.eq.l) go to 300
      mm=m-1
  200 if (j.ne.30) go to 210
      fail=.true.
      return
c
c *   form shift
c
  210 j=j+1
      el=e(l)
      dl=d(l)
      p=(d(l+1)-dl)/(el+el)
      r=dsqrt(p*p+1.d00)
      if (p) 230,220,220
  220 q=p+r
      go to 240
  230 q=p-r
  240 h=dl-el/q
      do 250 i=l,n
  250 d(i)=d(i)-h
      f=f+h
c
c *   ql transformation
c
      p=d(m)
      c=1.0d0
      s=0.0d0
      kk=m
      do 280 i=l,mm
      kk1=kk
      kk=kk-1
      ekk=e(kk)
      dkk=d(kk)
      g=c*ekk
      h=c*p
      if (dabs(p).lt.dabs(ekk)) go to 260
      c=ekk/p
      r=dsqrt(c*c+1.d00)
      e(kk1)=s*p*r
      s=c/r
      c=1.0d0/r
      go to 270
  260 c=p/ekk
      r=dsqrt(c*c+1.0d0)
      e(kk1)=s*ekk*r
      s=1.0d0/r
      c=c/r
  270 p=c*dkk-s*g
      d(kk1)=h+s*(c*g+s*dkk)
c      do 280 ii=1,n
c      s1=z(ii,kk)
c      h=z(ii,kk1)
c      z(ii,kk1)=s1*s+c*h
c  280 z(ii,kk)=s1*c-s*h
      call drot(n,z(1,kk1),incx,z(1,kk),incy,c,s)
  280 continue
  290 el=s*p
      e(l)=el
      d(l)=c*p
      if (dabs(el).gt.b) go to 200
  300 d(l)=d(l)+f
  310 continue
c
c *   order eigenvalues and eigenvectors
c
      do 320 i=1,n
  320 iord(i)=i
      i1=1
      do 360 i=2,n
      i2=i
  330 if(d(i2).ge.d(i1)) go to 350
  340 temp=d(i2)
      d(i2)=d(i1)
      d(i1)=temp
      itemp=iord(i2)
      iord(i2)=iord(i1)
      iord(i1)=itemp
      i2=i1
      i1=i1-1
      if(i1.gt.0) go to 330
  350 i1=i
  360 continue
      do 380 i=1,n
      do 370 j=1,n
      jj=iord(j)
  370 e(j)=z(i,jj)
      do 380 j=1,n
  380 z(i,j)=e(j)
      return
      end
