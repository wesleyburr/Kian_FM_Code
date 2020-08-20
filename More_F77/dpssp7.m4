c  djt/src/ dpssp7.m4 Time Derivatives of DPSS's, Core
      subroutine dpssp7(efn,ndim,ndata,nefn,W,ev,efnp,x,y)
      real efn(ndim,nefn), ev(nefn), efnp(ndim,nefn)
      real*8 x(0:ndata), y(0:ndata), pi, tpW, ev8
c                      Kernel and derivative
      pi = 4.*ATAN8(1.d+00)
      tpW = 2.*pi*W
      x(0) = 2.*W
      y(0) = 0.
      do 100 n = 1, ndata
      x(n) = SIN8(tpW*dble(n))/(pi*dble(n))
  100 y(n) = ( 2.*W*COS8(tpW*dble(n)) - x(n) )/dble(n)
      nhalf = 1 + ndata/2
      do 1000 k = 1, nefn
      ev8 = 0.
      amx = 0.
      do 200 n = 1, nhalf
      if(abs(efn(n,k)).gt.amx) then
        amx = abs(efn(n,k))
        lmx = n
        endif
  200 continue
      do 300 n = 1, ndata
      j = iabs(n-lmx)
  300 ev8 = ev8 + efn(n,k)*x(j)
      ev8 = ev8/efn(lmx,k)
      ev(k) = ev8
      do 500 n = 1, ndata
      efnp(n,k) = 0.
      do 400 m = 1, ndata
      nmm = n - m
      if(nmm.ge.0) then
          efnp(n,k) = efnp(n,k) + y(nmm)*efn(m,k)
        else
          efnp(n,k) = efnp(n,k) - y(-nmm)*efn(m,k)
        endif
  400 continue
  500 efnp(n,k) = efnp(n,k)/ev8
 1000 continue
      return
      end
