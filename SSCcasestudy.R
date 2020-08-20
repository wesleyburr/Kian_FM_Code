source("utility.R")
elec <- read.csv("SSC2020_hourly_demand.xlsx - Hourly Demand.csv", header = TRUE, sep = ',')
elec.series <- elec$Total.Energy.Use.from.Electricity..MW.
N <- nrow(elec)
K <- 12
NW <- (K + 1)/2
blk <- block(N, 24*365, 0.1)
M <- blk$nBlock
nppb <- blk$increment
st <- blk$startIdx
W <- NW / nppb

nFFT <- 2^ceiling(log2(2*nppb))
mxdeg <- 3

specArray <- array(0, dim = c(nFFT/2+1, M))
FArray <- array(0, dim = c(nFFT/2+1, M-1))
FParray <- array(0, dim = c(nFFT/2+1,mxdeg+1,M))
IFarray <- array(0, dim = c(nppb+1,nFFT/2+1,M))
start.time <- Sys.time()
for( m in 1:(M-1)){
  elec.spec <- spec.mtm(elec.series[st[m]:st[m+1]], nw = NW, k = K, deltat = 1
                        , dtUnits = 'hour', Ftest = TRUE, returnInternals = TRUE, nFFT = nFFT,
                        plot = FALSE)
  specArray[,m] <- elec.spec$spec
  FArray[,m] <- elec.spec$mtm$Ftest
  
  yk.raw <- elec.spec$mtm$eigenCoefs
  dk <- sqrt(elec.spec$mtm$eigenCoefWt)
  
  yk <- dk * yk.raw
  yk <- t(yk)
  
  
  dw <- elec.spec$mtm$dpss$v
  ev <- elec.spec$mtm$dpss$eigen
  
  Z <- dw %*% yk
  
  dwp <- dpssp7(dw,nppb+1,K,W,ev)
  
  Zdot <- dwp %*% yk
  
  U <- Re(Z)
  V <- Im(Z)
  Udot <- Re(Zdot)
  Vdot <- Im(Zdot)
  
  A <- Mod(Z)
  
  phi <- (Vdot*U - Udot*V)/(2*pi*A^2)
  IFarray[,,m] <- phi
  
  ap <- dpssap(dw,4,alpha=1)
  H <- ap$U
  Poly <- ap$R
  
  Fcoef <- t(dw) %*% phi
  FPcoef <- t(H) %*% Fcoef
  
  
  
  Fpstuff <- FLoop2(K,mxdeg,H, Fcoef, FPcoef, nFFT = nFFT/2+1)
  FP <- matrix(data = Fpstuff$F1, nrow = mxdeg+1, ncol = nFFT/2+1)
  FParray[,,m] <- FP
}
runtime <- Sys.time() - start.time
freq <- elec.spec$freq

aveSpec <- apply(specArray, MARGIN = 1, FUN = mean)
aveF <- apply(FArray, MARGIN = 1, FUN = mean)

plot(freq, aveSpec, type= 'l', log = 'y')
plot(freq[-1], aveF[-1], type = 'l', log = 'y', ylim = c(1, max(aveF)))
abline(h = qf(0.999,2,2*K-2), lty = 2, col = 2)

sig.ind <- which(aveF >= qf(0.999,2,2*K-2))
freq[sig.ind]
per <- 1/freq[sig.ind]
perD <- per / 24

IFweekly <- NULL
for(m in 1:M){
  IFweekly <- c(IFweekly, IFarray[,sig.ind[1],m])
}
plot(IFweekly, type = 'l')

plot(0,0, type = 'l', col = 'white', xlim = c(0,nppb), ylim = range(IFweekly), xlab = NA, ylab = NA)
for(m in 1:M){
  lines(IFarray[,sig.ind[1],m], col = m)
}

plot(freq, FParray[,3,1], type = 'l')

aveFP <- matrix(data = 0, nrow = mxdeg+1, ncol = nFFT/2+1)
for(p in 1:(mxdeg+1)){
  for(n in 1:(nFFT/2+1)){
    aveFP[p,n] <- mean(FParray[n,p,])
  }
}

for(p in 1:(mxdeg+1)){
  aveFP[p,which(aveFP[p,] < 0)] <- 0.00001
}

plot(freq,rep(1,length(freq)), type = 'l', col = 'white', xlim = c(0.15,0.18),ylim = c(1,80)
     , xlab = 'frequency'
     , ylab = NA
     , log = 'y')
for(p in 1:(mxdeg+1)){
  lines(freq, aveFP[p,], col = p)
}
abline(h = qf(0.99,1,K-1:(mxdeg+1)), lty = 2, col = 1:(mxdeg+1), lwd = 2)
lines(freq, aveF, col = 5)
abline( h = qf(0.99,2,2*K-2), lty = 2, col = 5, lwd = 2)

plot(freq[-1], aveSpec[-1], type = 'l', log = 'y', xlab = 'frequency', ylab = NA)
par(new = TRUE)
for(m in 1:M){
  plot(freq[-1], FArray[-1,m], axes = FALSE, xlab = NA, ylab = NA, type = 'l', col = m)
  par(new = TRUE)
}
axis(side = 4)#, col = 2, col.ticks = 2)
