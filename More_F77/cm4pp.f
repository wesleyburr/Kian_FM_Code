      subroutine cm4pp(x, ndata, ave,var, ttl )
      real x(ndata)
      double precision avei, d, d2, cm2, cm3, cm4, anvt
      character ttl*(*)
c
c       Central 4 moments, & extremes
c x      Data array
c ndata  Number of data points
c ave    Upper and Lower Extreme Point trimmed Arithmetic average
c        The smallest and largest points are excluded
c        from ave, var (and sd) but not the other moments.
c var    Upper and Lower Extreme Point Trimmed Variance
c ttl    Data Title
c
      if(ndata.eq.0) then
         write(*,*) " * * * WARNING: cm4pp called with ndata = 0"
         write(*,"(a)") ttl
         var = 1.
         return
         endif
      call eenter('cm4pp',ncseqn)
      if(ndata.lt.0) call djtf('ndata negative')
      lmx = 1
      lmn = 1
      avei = 0.d+00
      do 100 n = 1, ndata
      if(x(n).gt.x(lmx)) lmx = n
      if(x(n).lt.x(lmn)) lmn = n
  100 avei = avei + dble(x(n))
      avnt = avei/dble(ndata)
      xmn = x(lmn)
      xmx = x(lmx)
      avei = (avei - xmn - xmx)/dble(ndata-2)
      ave = sngl(avei)
      if(xmn.eq.xmx) then
        write(*,*) "Min and Max equal in cm4pp", xmn,xmx
        write(*,"(a)") ttl
        call djtf("probable bad data in cm4pp")
        endif
c          Half range scaling on higher moments
      hr = (xmx-xmn)/2.
      scl = 1./hr
c          Delete extremes from variance calculation
      cm2 = -( dprod(xmn-ave,scl)**2 + dprod(xmx-ave, scl)**2 )
      cm3 = 0.d+00
      cm4 = 0.d+00
      do 200 n = 1, ndata
      d = dprod(x(n)-ave, scl)
      d2 = d**2
      cm2 = cm2 + d2
      cm3 = cm3 + d2*d
  200 cm4 = cm4 + d2**2
      var = sngl(cm2/dble(ndata-3))
      cm3 = cm3/dble(ndata-1)
      cm4 = cm4/dble(ndata-1)
      sd = sqrt(var)
      cskew = cm3/(var*sd)
      ckurt = cm4/var**2
c          Compensate for half-range scaling
      var = var*(hr**2)
      sd = sd*hr
      write(6, 12000) ndata,ttl, xmn,lmn, avnt, xmx,lmx, ave, var
     1, sd, cskew, ckurt
      call eleave(ncseqn)
      return
12000 Format(i9,'  points of ',a/
     1,' Min',1pe12.4,' at',i9,'  Average',1pe12.4
     2,'  Max',1pe12.4,' at',i9/
     3,' Single Trimmed:  Average',1pe12.4
     4,' Variance',1pe12.4,' Sdev',1pe12.4/
     5,' Skewness',0pf12.6,'    Kurtosis',0pf12.4)
      end
