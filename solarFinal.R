source("utility.R")
delT <- 20
ndec1 <- 9
#ndec2 <- 8
ndec <- ndec1 #* ndec2
#decimate to 16 minute sampling periods
nyquist <- 0.5/(ndec*delT)

load("solar_interp.RData")

yr <- 1576880
#solar.dt <- detrend(solar.interp, W = 6/(delT*length(solar.interp)), deltat = delT)
#solar.dm <- solar.dt - mean(solar.dt)
nflt <- 291
arOrder <- 2
solar.dm <- solar.interp[(yr/2+1):(1.5*yr + 2*(nflt-1) + arOrder)] - 
                mean(solar.interp[(yr/2+1):(1.5*yr + 2*(nflt-1) + arOrder)])


imp <- filterDesign(nflt = nflt, wflt = 0.5/ndec1 + 1/nflt)
solar.tmp1 <- filter(solar.dm, filter = imp, sides = 2)
solar.tmp1 <- solar.tmp1[-which(is.na(solar.tmp1))]
dec1 <- seq(1, length(solar.tmp1), by = ndec1)
solar.dec1 <- solar.tmp1[dec1] - mean(solar.tmp1[dec1])

#imp2 <- filterDesign(nflt = nflt, wflt = 0.5/ndec2 + 1/nflt)
#solar.tmp2 <- filter(solar.dec1, filter = imp2, sides = 2)
#solar.tmp2 <- solar.tmp2[-which(is.na(solar.tmp2))]
#dec2 <- seq(1, length(solar.tmp2), by = ndec2)
#solar.dec2 <- solar.tmp2[dec2]

solar <- solar.dec1 - mean(solar.dec1)

PW <- prewhiten(solar, deltat = ndec*delT, order = arOrder, rmLineComponents = TRUE)

solar.pw <- PW$pw


NW <- 5
K <- 2*NW - 1
M <- 3
N <- length(solar.pw)

pilotSP <- spec.mtm(solar.pw, nw = NW, k = K, deltat = ndec*delT, plot = FALSE, Ftest = TRUE,
                    returnInternals = TRUE, returnZeroFreq = FALSE)

freq <- pilotSP$freq
freq <- freq[-length(freq)]
nfreqs <- length(freq)
fBlx <- block(n = nfreqs, blockSize = floor(nfreqs/400), overlap = 0)
nblx <- fBlx$nBlock
idx <- (nfreqs - fBlx$startIdx)[nblx:1] #I want the frequency band minus the stuff near 0

DW <- dpss(n = N, nw = NW, k = K)
Vdot <- dpssp11R(DW,NW)
AP <- dpssap(DW$v, M, 1)
scls <- inverseScales(DW$v,Vdot)

F1 <- F2 <- F3 <- matrix(0, nrow = M+1, ncol = nblx*fBlx$blockSize)

for(i in 1:nblx){
  ind <- (idx[i]-fBlx$blockSize+1):idx[i]
  ModF <- ModulatedF12(pilotSP, mxdeg = M, adaptive = FALSE, subBand = ind,
                       derivIN = Vdot, dpssIN = DW, apIN = AP, scale = TRUE, scalesIN = scls)
  F1[,ind] <- ModF$F1
  F2[,ind] <- ModF$F2
  F3[,ind] <- ModF$F3
}

F1plotting <- F1
F2plotting <- F2
F3plotting <- F3

for(m in 1:(M+1)){
  F1plotting[m,which(F1plotting[m,] < 1)] <- 1
  F2plotting[m,which(F2plotting[m,] < 1)] <- 1
  F3plotting[m,which(F3plotting[m,] < 1)] <- 1
}
