sr.new <- floor(0.5/W)
ds.length <- floor(ndata/sr.new)
Zds <- Fit1[(0:(ds.length-1))*sr.new + 1,]
Zdot.ds <- Zdot[(0:(ds.length-1))*sr.new + 1,]
Uf <- Re(Zds)
Wf <- Im(Zds)
Ufdot <- Re(Zdot.ds)
Wfdot <- Im(Zdot.ds)
amp <- Uf^2 + Wf^2
IF.ds <- (Wfdot*Uf - Ufdot*Wf)/(2*pi*amp)
f.ind <- which.min(abs(freq-0.2))
plot(IF.ds[,f.ind], type = 'l')
plot(FrA[,f.ind], type = 'l')


#Vds <- dpss(ds.length, nord, ds.length*W)$v
Vds <- dw[(0:(ds.length-1))*sr.new + 1,]
Psi_f <- t(Vds) %*% IF.ds
Hp.ds <- dpsshp(Vds, ds.length, nord, mxdeg, HeSmx, norm = TRUE)
Hds <- Hp.ds$A
Gds <- Hp.ds$Poly
Psi_Pf <- t(Hds) %*% Psi_f

F1stuff <- FLoop2(nord,mxdeg,Hds,Psi_f,Psi_Pf,nFFT)
F1 <- matrix(data=F1stuff$F1, nrow = mxdeg+1, ncol = nFFT)
plot(freq, F1[3,], type = 'l')
abline(h = qf(1-1/ndata, 1, nord - 3), lty = 2, col = 2)
plot(freq, F1array[3, ,5,1000], type = 'l')
abline(h = qf(1-1/ndata, 1, nord - 3), lty = 2, col = 2)

Fit.ds <- Gds[,1:3] %*% Psi_Pf[1:3,]
plot(Fit.ds[,f.ind], type = 'l')
par(new = TRUE)
plot(FMpoly, type = 'l', col = 2)
par(new = TRUE)
Fit <- Poly[,1:3] %*% FPcoef[1:3,]
plot(Fit[,f.ind], type = 'l', col = 4)
Fit.ds2 <- Fit[(0:(ds.length-1))*sr.new + 1,]
plot(Fit.ds[,f.ind], type = 'l')
lines(Fit.ds2[,f.ind], col = 2)
plot(FMpoly, type = 'l')
lines(Fit[,f.ind], col = 4)
par(new = TRUE)
plot(Fit.ds[,f.ind], type = 'l', col = 2, axes = FALSE, xlab = '', ylab = '')

sig.ind <- which(F1array[3,,5,1000] > qf(1-1/ndata, 1, nord - 3) )
true.ind <- which(F1[3,sig.ind] > qf(1-1/ndata, 1, nord - 3))
freq[sig.ind[true.ind]]
