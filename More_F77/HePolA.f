c djt/src/HePolA.f  Hermite (He) Polynomial Array
      subroutine HePolA(x,ndata,He,LdimHe, maxdeg, Hn,jn)
      real He(LdimHe,0:maxdeg), Hn(0:maxdeg)
      double precision  x(LdimHe), HeKp1,HeK,HeKm1
c
c x        Points where polynomials to be evaluated (Double)
c ndata    Number of x points
c He        Polynomials
c LdimHe    Leading dimension of He, LdimHe >= ndata
c maxdeg   Polynomials degree 0 to maxdeg returned
c Hn       Normalization constants
c jn       0, leave in standard form; 1, scale by 1/sqrt(Hn)
c
      call eenter('HePolA',ncseqn)
      if(maxdeg.lt.0) call djtf("maxdeg < 0 in Gbpol")
      call djtn(1,ndata,LdimHe,'number data, Lead. Dim He')
      do 100 n = 1, ndata
  100 He(n,0) = 1.
      if(maxdeg.eq.0) go to 2222
      do 200 n = 1, ndata
  200 He(n,1) = x(n)
      if(maxdeg.eq.1) go to 2222
      do 400 n = 1, ndata
      HeKm1 = 1.d+00
      HeK   = x(n)
      do 300 k = 1, maxdeg -1
      HeKp1 = x(n)*HeK - dble(k)*HeKm1
      HeKm1 = HeK
      HeK = HeKp1
  300 He(n,k+1) = HeKp1
  400 continue
 2222 continue
      call eleave(ncseqn)
      return
      end
