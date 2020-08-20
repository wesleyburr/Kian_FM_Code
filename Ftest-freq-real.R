source("utility.R")
##########################################################################
#set up frequency modulation and generate data
ndata <- 2000 #N
nppb <- ndata  #use Dave's block function to set this if I do this over blocks
nord <- 10 #This is K = 2NW - 1
mxdeg <- 4 #max degree of polynomial modulation

js1 <- 420
js2 <- 5318008
seed <- sample(js1:js2, size = 1)
set.seed(seed)

nFFT <- 2^ceiling(log2(2*ndata))
freq <- seq(0, 0.5, length.out = nFFT + 1)
freq <- freq[1:nFFT]

A <- 2
f <- freq[which.min(abs(freq-0.2))] #center frequency
W <- (nord + 1)/(2*ndata)
fbw <- 0.8
stm11 <- 2.0/(ndata-1)
Shift <- 0.0
dt <- 1
t <- 0:(ndata-1)        #time
tt <- (t + Shift) * stm11 - 1.0  # this runs from -1 to 1
FMpoly <- W * fbw * (1.0 - 2.0 * tt^2)
Fmod <- f + FMpoly  # quadratic frequency modulation
Pref <- cumsum(FMpoly)*2*pi*dt
PhMod <- 2*pi*f*dt*t + Pref

cd <- A*cos(PhMod)  #+ 1i * sin(PhMod) #data without noise

snr <- 3									#signal to noise ratio
sigmat2 <- ndata/(2*snr*nord)  
sigmat <- sqrt(sigmat2)
CDN <- cd + rnorm(ndata, sd = sigmat) #+ rnorm(ndata, sd = sigmat) * 1i

##########################################################################
#slepians, derivatives and special polynomials
DW <- dpss(n=ndata, k=nord, nw=ndata * W) #slepian sequences
dw <- DW$v
ev <- DW$eigen

dwp <- dpssp7(dw,ndata,nord,W,ev)

ap <- dpssap(dw,mxdeg,alpha=1)
H <- ap$U
Poly <- ap$R

HeSmx <- 0.9
hp <- dpsshp(dw,ndata,nord,mxdeg,HeSmx,norm = TRUE)
Poly <- hp$Poly
H <- hp$A

######################################################################################

tapered <- dw * CDN 
tapered <- rbind(tapered, matrix(0, nrow = 2*nFFT-ndata, ncol = nord))
yk <- mvfft(tapered)
yk <- t(yk)
yk <- yk[,1:nFFT]

DD <- dpssd(n=ndata,nw=ndata*W,k=nord)
dd <- DD$v
ddp <- DD$deriv

tprd <- dd * CDN 
tprd <- rbind(tprd, matrix(0, nrow = 2*nFFT-ndata, ncol = nord-2))
ykd <- mvfft(tprd)
ykd <- t(ykd)
ykd <- ykd[,1:nFFT]

Zd <- dd %*% ykd
Zddot <- ddp %*% ykd
ampd <- Mod(Zd)
Uf <- Re(Zd)
Wf <- Im(Zd)
Ufdot <- Re(Zddot)
Wfdot <- Im(Zddot)
IFd <- (Wfdot*Uf - Ufdot*Wf)/(2*pi*ampd^2)
Fdcoef <- t(dd) %*% IFd
APd <- dpssap(dd, mxdeg, alpha = 1)
Hd <- APd$U
Polyd <- APd$R

FPdcoef <- t(Hd) %*% Fdcoef

F1stuffd <- FLoop2(nord-2, mxdeg, Hd, Fdcoef, FPdcoef, nFFT)
F1d <- matrix(data= F1stuffd$F1, nrow = mxdeg+1, ncol = nFFT)

FitD <- Polyd %*% FPdcoef
#rewrite dpssp7 for the modified slepians. don't need to do this. just generate dwp first
# and then use that.



#############################################################################
Fit1 <- dw %*% yk #standard inverse
Zdot <- dwp %*% yk #time derivative of Fit1.

amp <- Mod(Fit1)
PhD <- Arg(Fit1) #phase of standard inverse
x1 <- Re(Fit1)
y1 <- Im(Fit1)
Xdot <- Re(Zdot)
Ydot <- Im(Zdot)
FrA <- ( x1*Ydot - Xdot*y1 ) / ((2*pi)*amp^2) #+ Freq

Fcoef <- t(dw) %*% FrA
FPcoef <- t(H) %*% Fcoef

Fit <- Poly %*% FPcoef
r <- Fcoef - H %*% FPcoef
FitC <- Fit + dw %*% r

SSQC <- colSums( Fcoef^2 )
F1DJTstuff <- FLoopDJT(nord,mxdeg,H,Fcoef,FPcoef,SSQC,nFFT)
F1DJT <- matrix(F1DJTstuff$F1mat, nrow = mxdeg+1, ncol = nFFT)

F1stuff <- FLoop2(nord, mxdeg, H, Fcoef, FPcoef, nFFT)
F1 <- matrix(data= F1stuff$F1, nrow = mxdeg+1, ncol = nFFT)

FPstuff <- FPLoop(nord,mxdeg,H,Fcoef,nFFT)
FP.temp <- matrix(data=FPstuff$FP, nrow = mxdeg + 1, ncol = nFFT)
plot(freq, FP.temp[1,], type = 'l')

FPk.temp <- array(data=0, dim = c(nord,mxdeg+1,nFFT))
for(k in 1:nord){
  FPstuff <- FPLoop(nord-1,mxdeg,H[-k,],Fcoef[-k,],nFFT)
  FPk.temp[k,,] <- FPstuff$FP
}

#drop stuff from frequencies less than W
notAnalytic <- which(freq < W)
FPk <- FPk.temp[,,-notAnalytic]
N <- nFFT - length(notAnalytic)

jkBar <- matrix(data=0, nrow = mxdeg+1, ncol = N)
jkVar <- matrix(data=0, nrow = mxdeg+1, ncol = N)
for(m in 1:(mxdeg+1)){
  for(n in 1:N){
    jkBar[m,n] <- mean(FPk[,m,n])
    jkVar[m,n] <- (nord-2+1/nord)*var(FPk[,m,n])
  }
}

FquantU <- matrix(rep(qf(0.975,1:(mxdeg+1),nord-1:(mxdeg+1)), N), nrow=mxdeg+1, ncol=N)
FquantL <- matrix(rep(qf(0.025,1:(mxdeg+1),nord-1:(mxdeg+1)), N), nrow=mxdeg+1, ncol=N)

FP <- FP.temp[,-notAnalytic]
FP.CI.U <- jkBar +  qnorm(0.975) * sqrt(jkVar)
FP.CI.L <- jkBar - qnorm(0.975) * sqrt(jkVar)

plot(freq[-notAnalytic], FP[1,], type = 'l', ylim = c(0,max(FP.CI.U[1,])))
lines(freq[-notAnalytic], FP.CI.U[1,], col = 2, lty = 2)
lines(freq[-notAnalytic], FP.CI.L[1,], col = 2, lty = 2)
abline(h = qf(0.999, 1, nord - 1), lty = 3, col = 4, lwd = 2)

FP2 <- matrix(0, nrow = mxdeg+1, ncol = nFFT)
for( p in 1:(mxdeg+1)){
  r <- Fcoef - H[,1:p, drop = FALSE] %*% FPcoef[1:p,, drop = FALSE]
  Z <- dw %*% r + Poly[,1:p, drop = FALSE] %*% FPcoef[1:p,, drop = FALSE]
  res2 <- colSums( (FrA - Z)^2 )
  ssq <- colSums( Z^2 )
  FP2[p,] <- ( ssq / p ) / (res2 / (nord - p))
}
plot(freq, FP2, type = 'l')
abline(h = qf(0.99, mxdeg+1, nord-mxdeg-1), lty = 2, col = 2)

plot(freq, FP2[1,], type = 'l', ylab = 'FP2', xlab = 'Frequency')
for(i in 2:(mxdeg+1)){
  lines(freq, FP2[i,], col = i)
}
legend(x=0.3, y=110, legend = 0:mxdeg, col = 1:(mxdeg+1), ncol = mxdeg+1, pch = 15,
       title = 'degree', bty = 'n')

plot(FMpoly, type = 'l')
#jkstuff <- jackknifeF(mxdeg,nord,H,Fcoef) 
#jkVar <- matrix(data = jkstuff$jkVar, nrow = mxdeg+1, ncol = nFFT)

#jackknife. something fishy happening here
#jkFarray <- array(0, dim = c(mxdeg+1,nFFT,nord))
#for(k in 1:nord){
#  Fk <- FLoop2(nord-1,mxdeg,H[-k,],Fcoef[-k],FPcoef - H[k,]*Fcoef[k], nFFT)
#  jkFarray[,,k] <- Fk$F1
#}

#jkFvar <- array(data=0, dim = c(mxdeg+1,nFFT))
#for(m in 1:(mxdeg+1)){
#  for(n in 1:nFFT){
#    jkFvar[m,n] <- var(jkFarray[m,n,])*(nord-1)^2 / nord
#  }
#}

df.ind <- 3
crit.F <- qf(1-1/ndata, 1, nord-df.ind)
plot(freq, F1[df.ind,], type = 'l')
abline(h = crit.F, lty = 2, col = 2)

plot(0:mxdeg, F1[,f.ind], type = 'h', lwd = 2)
crit.F <- qf(1-1/ndata, 1, nord - 0:mxdeg - 1)
points(0:mxdeg, crit.F, pch = 8, col = 2)

err <- abs(freq[f.ind] - f)
err

f.ind <- which.max(F1[df.ind,])
Fit <- Poly[,1:df.ind] %*% FPcoef[1:df.ind,f.ind]
r <- Fcoef[,f.ind] - H[,1:df.ind] %*% FPcoef[1:df.ind,f.ind]
FitC <- Fit + dw %*% r

plot(FMpoly, type = 'l')
lines(Fit, col = 2)
lines(FitC, col = 4)
lines(FrA[,f.ind], col = 6)

#fmax <- freq[f.ind]
#FitMod <- 2*pi*cumsum(FitC)
#FitPhMod <- 2*pi*fmax*t + FitMod
#recon <- cos(FitPhMod)

#plot(cd[1:100], type = 'l')
#lines(recon[1:100], col = 2)

#plot(Pref, type = 'l')
#lines(FitMod, col = 2)
