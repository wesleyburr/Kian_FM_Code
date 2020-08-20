c djt/src/ dpssp.f  Time Derivatives of DPSS's, Driver
      subroutine dpssp(efn,ndim,ndata,nefn,W,ev,efnp)
      real efn(ndim,nefn),ev(nefn), efnp(ndim,nefn)
      real*8 R8
      common/cstak/R8(1000)
      save  /cstak/
c
c        Discrete Prolate Spheroidal Sequences
c
c efn    Eigenfunction array.
c        efn(t,k) is the k+1 sup st  dpss
c ndim   column dimension = max number of time points
c ndata  number of data (time points) to return
c        ndata = Slepian's N
c nefn   The number of eigenfunctions to compute.
c        Lowest order functions returned first.
c W      the bandwidth parameter ( as in Slepian. )
c ev     Integral Eigenvalues, lambda
c evp    1 - lambda(k)  = 1 - integral eigenvalue
c
      call senter('dpssp',ncseqn)
      call djta(numarg(),7)
      lx = istkg8(ndata+1)
      ly = istkg8(ndata+1)
      call dpssp7(efn,ndim,ndata,nefn,W,ev
     1, efnp, R8(lx), R8(ly) )
      call sleave(ncseqn)
      return
      end
