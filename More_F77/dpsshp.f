c djt/odin/ts/dpsshp.f DPSS Assoc. Poly, Mod Gram Schmidt, Full, Hermite
      subroutine dpsshp(dw,LdDW,nppb,nord, Poly,LdPoly,ndata,maxdeg
     1, HeScl,norm, A,HeNrm,np)
      real dw(LdDW,0:nord-1), Poly(LdPoly,0:maxdeg)
     1, A(0:nord-1,0:maxdeg),  Hn(0:maxdeg), rs(1000)
      double precision sum, scl, ds(500), cntr,alphaM
      equivalence (ds(1),rs(1))
      common/cstak/ds
c
c                   Orthogonal Polynomials on Slepian Sequence Projections
c
c              David J. Thomson, Queen's University, March 2009
c
c dw      Array of Slepian Sequences, orthonormal on (1,nppb)
c LdDW    Leading Dimension of dw
c nppb    Length (N) of Slepian Sequences
c nord    Number of Slepian Sequences
c
c Poly    Associated Polynomials, Based on Hermite He_n ( HeScl * x )
c LdPoly  Leading Dimension >= nppb, ndata
c ndata   number data points of polynomials, >= nppb
c maxdeg  Maximum Degree
c HeScl   Scaling for Hermite polynomials
c A       Orthonormal Projection Matrix, A(k,j) = Inner Product of
c        dw(t,k) with Poly(t,j) on [1,nppb]. The Polynomials are
c        defined by the requirement that
c        \sum_{k=0}^{K-1} A(k,j) * A(k,L) = \delta_{j,L}
c         Note that the polynomials are NOT orthogonal in the sense
c         \sum_{t=1}^{nppb} Poly(t,j) *Poly(t,L) \ne  \delta_{j,L}
c Hn      \sum_{t=1}^{nppb} Poly(t,j)^2
c np      Print control (not used)
c
c
      call senter('dpsshp',ncseqn)
      write(*,*) "dpsshp: LdDW,nppb,nord",LdDW,nppb,nord
     1," LdPoly,ndata,maxdeg",LdPoly,ndata,maxdeg," HeScl",HeScl
     2," norm,np",norm,np
      if(np.gt.0) write(*,*) "DPSS Associated Polynomials for nppb="
     1,nppb," nord=",nord," ndata", ndata," maxdeg =",maxdeg
     2," LdDW", LdDW," LdPoly", LdPoly
     3," Hermite Form, He_n, Scale = ", HeScl
      call djtn(4,nppb,LdDW,'nppb not in [4,LdDW]')
      if(nord.lt.1) call djtf('nord < 1')
      if(LdPoly.lt.LdDW) call djtf('LdPoly < LdDW')
      call djtn(nppb,ndata,LdPoly,'ndata not in [nppb,LdPoly]')
      if(maxdeg.gt.nord-1) call djtf('maxdeg > nord-1')
      if(HeScl.lt.1.0e-06) call djtf("HeScl < 1.e-06 ?")
c
c                                    Stack storage for t
      LDST = istkgt(LdPoly,4)
c                                   Storage for Weight
      LDSW = istkgt(LdPoly,4)
      Ktop = nord-1
      cntr = dble(1+nppb)/2.d+00
      scl = dble(HeScl)*2.d+00/dble(nppb-1)
      do 200 n = 1, ndata
      xx = (dble(n) - cntr) * scl
      ds(LDST +n-1) = xx
  200 continue
c                   Start with Hermite Polynomials for reasonable match
      call HePolA(ds(LDST),ndata,Poly,LdPoly,maxdeg, Hn,jn)
      if(norm.gt.0) then
        do 422  L = 0, maxdeg
        sum = 0.d+00
        do 405 n = 1, nppb
  405   sum = sum + dprod(Poly(n,L),Poly(n,L))
        scl = 1./dsqrt(sum)
        do 410  n = 1, nppb
  410   Poly(n,L) = Poly(n,L)*scl
  422   continue
        endif
c         Weighted Inner Products of Polynomials and Slepian Sequences
      do 400 L = 0, maxdeg
      do 320 k = 0, Ktop
  320 A(k,L) = 0.
      Kmin = mod(L,2)
      do 400 k = Kmin, Ktop, 2
      sum = 0.d+00
      do 300 n = 1, nppb
  300 sum = sum + dprod(dw(n,k),Poly(n,L))
CX  300 sum = sum + dprod(dw(n,k),Poly(n,L))*ds(LDSW +n-1)
  400 A(k,L) = sum
      if(np.gt.3) then
        write(*,*) "Raw Inner Products, alpha", alpha
        do 401 L = 0, maxdeg
  401   write(*,"(i3,1p11e12.4)") L,(A(k,L),k=0,nord-1)
        endif
c                         Degree 0, 1
      do 460 L = 0, min(1,maxdeg)
      sum = 0.d+00
      do 420 k = L, Ktop, 2
  420 sum = sum + A(k,L)**2
      scl = 1.d+00/dsqrt(sum)
      if(np.gt.1) write(*,*) L,sum, scl
      do 430 k = L, Ktop, 2
  430 A(k,L) = A(k,L) * scl
      do 440 n = 1, ndata
  440 Poly(n,L) = Poly(n,L)*scl
  460 continue
      if(maxdeg.lt.2) go to 9999
c                             Do Gram-Schmidt
      do 600 L = 2, maxdeg
      Kmin = mod(L,2)
      do 550 j = Kmin, L-2, 2
      sum = 0.d+00
      do 520 k = Kmin, Ktop, 2
  520 sum = sum + A(k,L) * A(k,j)
      do 530 k = Kmin, Ktop, 2
  530 A(k,L) = A(k,L) - sum * A(k,j)
      do 540 n = 1, ndata
  540 Poly(n,L) = Poly(n,L) - sum * Poly(n,j)
c                                       Orthonormalize
      sum = 0.d+00
      do 550 k = Kmin, Ktop, 2
  550 sum = sum + A(k,L)*A(k,L)
      scl = 1.d+00/dsqrt(sum)
      do 560 k = Kmin, Ktop, 2
  560 A(k,L) = A(k,L) * scl
      do 570 n = 1, ndata
  570 Poly(n,L) = Poly(n,L)*scl
  600 continue
 9999 do 700 L = 0, maxdeg
      sum = 0.d+00
      do 650 n = 1, nppb
  650 sum = sum + dprod(Poly(n,L),Poly(n,L))
  700 Hn(L) = sum
      if(np.gt.2) then
        write(*,*) "Non-zero dpsshp expansion coefficients * 100"
        do 810 L = 0, maxdeg
        Lmin = mod(L,2)
  810   write(*,"(i3,2p10f12.6)") L,(A(k,L),k=Lmin,min(21,Ktop),2)
        write(*,*) "L2 norms of polynomials", (Hn(L),L=0,maxdeg)
        endif
      call sleave(ncseqn)
      return
      end
