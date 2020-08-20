source("utility.R")

N <- 10000
n <- 0:(N-1)
nFFT <- 2^ceiling(log2(2*N))


#multitaper and background spectrum stuff
NW <- 10
K <- 2*NW - 2
W <- NW / N
dt <- 1

sigmat2 <- 1 #really snr is (A^2/2) / (2W*S(f))
sigmat <- sqrt(sigmat2)

#ar <- c(0.5, -0.3)
#ma <- -0.1

###
DW <- dpss(n = N, nw = NW, k = K)
L <- 10*K
DS <- dpss(n = L, nw = NW, k = K)
V <- DS$v
Vdot <- dpssp11R(DS,NW)

nfreqs <- nFFT/2
num.sim <- 10000

scls <- inverseScales(V,Vdot)

seed <- 111
set.seed(seed)

freq <- 1:nfreqs / (nFFT*dt)
#band.ind <- which(freq < f + W & freq > f - W)
subBand <- sample(1:(nfreqs-1), size = 1)

NumArray <- AmpArray <- PhiArray <- array(data = 0, dim = c(L,num.sim))

start.time <- Sys.time()
for(i in 1:num.sim){
  
  #X <- arima.sim(model = list(ar = c(0.1, -0.2), ma = c(0.4, -0.3)), n = N)
  X <- rnorm(N)
  sp.X <- spec.mtm(X, k = K, nw = NW, deltat = dt, Ftest = FALSE, returnInternals = TRUE, 
                   dpssIN = DW, plot = FALSE, nFFT = nFFT, returnZeroFreq = FALSE)
  
  aS <- analSig(sp.X, dpssIN = DS, derivIN = Vdot, scale = TRUE, scaleIN = scls,
                subBand = subBand)
  
  NumArray[,i] <- aS$num
  AmpArray[,i] <- aS$amp2
  PhiArray[,i] <- aS$phi
  
  
  
}

runtime <- Sys.time() - start.time
runtime

###histograms and stuff down here

set.seed(666)
pdf("IFhistograms_downsampled.pdf", width = 11, height = 8.5)
  par(mfrow = c(2,1))
  par(oma = c(4, 4, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
  par(mar = c(0, 0, 1, 1)) 
  
  lMesh <- seq(-6,6, length.out = num.sim)
  cMesh <- seq(0,15, length.out = num.sim)
  
  low.ind <- sample(1:floor(0.05*L), size = 1)
  hist(NumArray[low.ind,], freq = FALSE, breaks = 60, xlab = '', ylab = '', main = '')
  lines(density(NumArray[low.ind,], bw = "SJ", kernel = "optcosine"), col = 2, lwd = 2)
  lines(lMesh, dlaplace(lMesh), col = 4, lwd = 2)
  
  hist(AmpArray[low.ind,], freq = FALSE, breaks = 60, xlab = '', ylab = '', main = '')
  lines(density(AmpArray[low.ind,], bw = "SJ", kernel = "optcosine"), col = 2, lwd = 2)
  lines(cMesh, dchisq(cMesh, 2), col = 4, lwd = 2)
  
  mid.ind <- sample((floor(0.05*L)+1):(floor(0.95*L)-1), size = 5)
  
  for(i in 1:length(mid.ind)){
    hist(NumArray[mid.ind[i],], freq = FALSE, breaks = 60, xlab = '', ylab = '', main = '')
    lines(density(NumArray[mid.ind[i],], bw = "SJ", kernel = "optcosine"), col = 2, lwd = 2)
    lines(lMesh, dlaplace(lMesh), col = 4, lwd = 2)
    
    hist(AmpArray[mid.ind[i],], freq = FALSE, breaks = 60, main = '')
    lines(density(AmpArray[mid.ind[i],], bw = "SJ",  kernel = "optcosine"), col = 2, lwd = 2)
    lines(cMesh, dchisq(cMesh, 2), col = 4, lwd = 2)
  }
  
  hi.ind <- sample( floor(0.95*L):L, size = 1)
  hist(NumArray[hi.ind,], freq = FALSE, breaks = 60, xlab = '', ylab  = '', main = '')
  lines(density(NumArray[hi.ind,], bw = "SJ", kernel = "optcosine"), col = 2)
  lines(lMesh, dlaplace(lMesh), col = 4, lwd = 2)
  
  hist(AmpArray[hi.ind,], freq = FALSE, breaks = 60, main = '')
  lines(density(AmpArray[hi.ind,], bw = "SJ", kernel = "optcosine"), col = 2)
  lines(cMesh, dchisq(cMesh, 2), col = 4, lwd = 2)
dev.off()

nMesh <- seq(-3,3, length.out = num.sim)
PhiEC <- crossprod(V,PhiArray)
##### eigencoefficients
png("IFeigenHist_downsampled.png", width = 800, height = 450)
  par(mfrow = c(4,3))
  par(oma = c(4, 4, 0, 0))
  par(mar = c(0, 0, 1, 1))
  
  for(k in 1:12){
    hist( (PhiEC[k,] - mean(PhiEC[1,]))/sd(PhiEC[1,]) , breaks = 100, freq = FALSE, xlab = '', 
          ylab = '', main = '', axes = FALSE)
    lines(density( (PhiEC[1,] - mean(PhiEC[1,]))/sd(PhiEC[1,])), col = 2, lwd = 2)
    lines(nMesh, dnorm(nMesh), col = 4, lwd = 2)
  }
dev.off()