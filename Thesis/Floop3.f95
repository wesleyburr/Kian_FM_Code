subroutine Floop3(mxdeg, nord, FPcoef, Fcoef, Fp, nfreqs)
  implicit none
  real*8 FPcoef(0:mxdeg,nfreqs), Fcoef(nord,nfreqs) &
    , Fp(0:mxdeg,nfreqs)
  integer mxdeg, nord, nfreqs, i
  
  do i=1, nfreqs
    call Fcomp2(mxdeg,nord,FPcoef(:,i),Fcoef(:,i),Fp(:,i))
  end do
end