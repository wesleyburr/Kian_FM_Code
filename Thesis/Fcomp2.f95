subroutine Fcomp2(mxdeg, nord, FPcoef, Fcoef, Fp)
  implicit none
  real*8 FPcoef(0:mxdeg), Fcoef(nord), Fp(0:mxdeg) &
    , ssq1, ssq2
  integer mxdeg, nord, k, L
  
  
  ssq1 = 0.0d0
  do k = 1, nord
    ssq1 = ssq1 + Fcoef(k)**2
  end do
  
  
  ssq2 = 0.0d0
  do L = 0, mxdeg
    ssq2 = ssq2 + FPcoef(L)**2
    Fp(L) = dble(nord-L-1) * FPcoef(L)**2 / (ssq1 - ssq2)
  end do

end