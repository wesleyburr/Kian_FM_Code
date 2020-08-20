source("utility.R")
load("solar_interp.RData")

#N <- length(solar.interp)
delT <- 20
#maybe don't need to bother removing a trend
#solar.dt <- detrend(solar.interp, W = 6/(N*delT), deltat = delT)
#solar.dm <- solar.dt - mean(solar.dt)
solar.dm <- solar.interp - mean(solar.interp)

#low-pass filter stuff
## 12:1 decimation filter -> 4 minute sampling time and nyquist of 2.08 mHz

nflt1 <- 201
ndec1 <- 8


#not sure if I should just accept the drop off up near nyquist or what
imp1 <- filterDesign2(nflt = nflt1, wflt = 0.5/ndec1 + 1/nflt1)

solar.filt.tmp <- filter(solar.dm, filter = imp1, sides = 2)
solar.filt.tmp <- solar.filt.tmp[-which(is.na(solar.filt.tmp))]
dec.ind <- seq(1,length(solar.filt.tmp), by = ndec1)
solar.dec.tmp <- solar.filt.tmp[dec.ind]


nflt2 <- 201
ndec2 <- 6
imp2 <- filterDesign(nflt = nflt2, wflt = 0.5/ndec2 + 1/nflt2)

solar.filt <- filter(solar.dec.tmp, filter = imp2, sides = 2)
solar.filt <- solar.filt[-which(is.na(solar.filt))]
dec.ind <- seq(1, length(solar.filt), by = ndec2)
solar.dec <- solar.filt[dec.ind]





solar <- solar.dec - mean(solar.dec)

L <- length(solar)
ndec <- ndec1*ndec2
sr <- ndec*delT
nyquist <- 0.5/sr

###################################
PW <- prewhiten(solar, deltat = sr, order = 2, rmLineComponents = FALSE)
# I hope this is sufficiently prewhitened
#solar4 <- PW$pw
#save(solar4, file = "pw_4min_solar.RData")
solar16 <- PW$pw
save(solar16, file = "pw_16min_solar.RData")

#L <- length(solar4)
L <- length(solar16)
nFFT <- 2^ceiling(log2(5*L))
nfreqs <- nFFT/2
#think about size. L*nfreqs*8/(1024)^3 is the size of each matrix

#block in time instead of frequency, since the data is so nonstaionary.
#also, consider using my "empirical quantiles" for testing at the end
#use two week chunks, or maybe months. maybe try both. should be relatively quick.

chunk <- 24 #24 chunks - 24 months in two years. overlapped by a 0.5 gives us 2 week blocks
blx <- block(L, floor(L/chunk), 0.5)
idx <- blx$startIdx
blkLen <- blx$blockSize

#continue here. loop this over the blocks, store in arrays, make some heatmaps

#choose nw, k, mxdeg. maybe also see the affect of downsampling a bit and using adaptive
#weighting



M <- 6
NW <- 6   #think about these
# don't necessarily need W for anything but W = NW / (blkLen*sr) = 2.8e-06
K <- 2*NW - 2
DW <- dpss(n = blkLen, nw = NW, k = K)
AP <- dpssap(DW$v, M, alpha = 1)
Vdot <- dpssp11R(DW, NW)

HFarray <- array(data = NA, dim = c(nfreqs, blx$nBlock))
F1array <- F2array <- F3array <- array(data = NA, dim = c(M+1, nfreqs-2, blx$nBlock))

start.time <- Sys.time()
for(i in 1:blx$nBlock){
  solar.block <- solar16[idx[i]:(idx[i] + blkLen - 1)]
  modSP <- spec.mtm(solar.block, deltat = sr, nw = NW, k = K, dtUnits = 'second',
                    plot = FALSE, returnInternals = TRUE, Ftest = TRUE, nFFT = nFFT,
                    returnZeroFreq = FALSE)
  modF <- ModulatedF6(modSP, mxdeg = M, adaptive = FALSE, subBand = 1:(nfreqs-2), 
                      derivIN = Vdot, dpssIN = DW, apIN = AP)
  HFarray[,i] <- modSP$mtm$Ftest
  F1array[,,i] <- modF$ModF$F1
  F2array[,,i] <- modF$ModF$F2
  F3array[,,i] <- modF$ModF$F3
}
runtime <- Sys.time() - start.time
runtime

#let's see what we have.

sigF1array <- F1array
sigF2array <- F2array
sigF3array <- F3array
sigHFarray <- HFarray

conf <- 1-1/blkLen
for(b in 1:blx$nBlock){
  for(m in 1:(M+1)){
    sigF1array[m,which(sigF1array[m,,b] < quants$F1quant[m,conf*10000]),b] <- NA
    sigF2array[m,which(sigF2array[m,,b] < quants$F2quant[m,conf*10000]),b] <- NA
    sigF3array[m,which(sigF3array[m,,b] < quants$F3quant[m,conf*10000]),b] <- NA
  }
  sigHFarray[which(sigHFarray[,b] < qf(conf, 2, 2*K-2)),b] <- NA
}


pdf("Fstat_16min_blocked.pdf")
  fields::image.plot(x = modSP$freq, y = 1:blx$nBlock, z = HFarray)
    for(m in 1:(M+1)){
      fields::image.plot(x = modSP$freq[-c(1,nfreqs)], y = 1:blx$nBlock, z = sigF1array[m,,])
      title(paste("F1 degree ", m-1))
      
      fields::image.plot(x = modSP$freq[-c(1,nfreqs)], y = 1:blx$nBlock, z = sigF2array[m,,])
      title(paste("F2 degree ", m-1))
      
      fields::image.plot(x = modSP$freq[-c(1,nfreqs)], y = 1:blx$nBlock, z = sigF3array[m,,])
      title(paste("F3 degree ", m-1))
    }
dev.off()

















### let's try out my useless F statistic
#according to DJT 09 FM paper, modulation is about 20 nHz per year. So 40 nHz over two years?
#This is for g modes I think, but still. Take W at least 20 nHz so that a 2W band is 40 nHz wide.
# Using NW = 4, we have 4/(L*ndec*delT) = 6.337569e-08 > 20e-09 = 2e-08. Ok. The last eigenvalue
#is a bit small ~ 0.936, should I use k = 6?
NW <- 4
K <- 2*NW - 1
modSP <- spec.mtm(solar4, deltat = sr, nw = NW, k = K, dtUnits = 'second',
                  plot = FALSE, returnInternals = TRUE, Ftest = TRUE, nFFT = nFFT)

#want the same frequency resolution in the Harmonic F as the modulated F's. Think about this.

#set up frequency blocks
w <- NW / (L*sr)
freq <- modSP$freq
analFreq <- which(freq > w)
nfreqs <- length(analFreq)


#analytic signal only exists for f > W (apparently)
#for memory allocation purposes can only work with a small band at once.

numBands <- 12000
fblkLen <- floor(nfreqs/numBands)
fblx <- block(n = nfreqs, blockSize = fblkLen, overlap = 0)
fIdx <- fblx$startIdx
nfBlx <- fblx$nBlock

M <- 6
modF <- NULL

DW <- modSP$mtm$dpss
V <- DW$v

Vdot <- dpssp11R(DW,NW)
AP <- dpssap(V, maxdeg = M, alpha = 1)

begin.time <- Sys.time()
for(i in 1:nfBlx){
  subBand <- analFreq[fIdx[i]:(fIdx[i]+fblkLen-1)]
  modF.tmp <- ModulatedF4(spec = modSP, mxdeg = M, subBand = subBand, dpssIN = DW, 
                          derivIN = Vdot, apIN = AP, adaptive = TRUE)
  modF <- cbind(modF, modF.tmp$ModF)
  rm(modF.tmp)
}

subBand <- which(analFreq > analFreq[fIdx[nfBlx]+fblkLen-1])
modF.tmp <- ModulatedF4(spec = modSP, mxdeg = M, subBand = subBand, dpssIN = DW,
                        derivIN = Vdot, apIN = AP, adaptive = TRUE)
modF <- cbind(modF, modF.tmp$ModF)
rm(modF.tmp)
freqF <- freq[analFreq]

end.time <- Sys.time()
comp.time <- end.time - begin.time
#I'm dumb and forgot that I was using 7 tapers so with degree 6, there should be division
#by 0 in the denominator.
modF <- modF[1:M,]

modF.stuff <- list(modF = modF, freqF = freqF)
save(modF.stuff, file = "modulated_F_solar4.RData")

modF.plotting <- modF
for(m in 1:(M+1)){
  modF.plotting[m, which(modF.plotting[m,] < 1)] <- 1
}

H.Ftest.plotting <- modSP$mtm$Ftest
H.Ftest.plotting[which(H.Ftest.plotting < 1)] <- 1

png(filename = "modFsolar4.png", width = 600, height = 900)
par(mfrow = c(3,2))
plot(freq, H.Ftest.plotting, xlab = 'Frequency', ylab = '', type = 'l', log = 'y')
abline(h = qf(1-1/L,2,2*modSP$mtm$k - 2), lty = 2, col = 2)
title("Harmonic F-test Statistic")

for(m in 2:M){
  plot(freqF, modF.plotting[m,], type = 'l', log = 'y', ylab = '', xlab = 'Frequency')
  abline(h = qf(1-1/N, 1, modSP$mtm$k - m), lty = 2, col = 2)
  title(paste("Degree", m-1)) 
}
dev.off()

Hf.mxs <- findLocalFMax(modSP, cutoff = 1-1/L)
mxs <- findLocalFMaxM(modF.plotting, k = modSP$mtm$k, cutoff = 1-1/L)
sink("significantFrequencies.txt")
for(m in 1:M){
  cat(paste("degree ", m-1, "significant frequencies in Hz: "))
  print(freqF[ mxs[[m]] ])
  cat("\n")
}
cat("Harmonic F statistic significant frequencies in Hz: ")
print(freq[Hf.mxs])
cat("\n")
sink()
