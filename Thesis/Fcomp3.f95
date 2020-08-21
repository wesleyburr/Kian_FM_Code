subroutine Fcomp3(mxdeg, nord, FPcoef, Fcoef, Fp)
  implicit none
  real*8 FPcoef(mxdeg), Fcoef(nord), Fp(mxdeg) &
    , ssq1, ssq2
  integer mxdeg, nord, k, L
  
  
  ssq1 = 0.0d0
  do k = 1, nord
    ssq1 = ssq1 + Fcoef(k)**2
  end do
  
  
  ssq2 = 0.0d0
  do L = 1, mxdeg
    ssq2 = ssq2 + FPcoef(L)**2
    Fp(L) = dble(nord-L) * FPcoef(L)**2 / (ssq1 - ssq2)
  end do

end