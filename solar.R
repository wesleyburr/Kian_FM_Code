source("utility.R")
load("GOLF_PM2.RData")
library(tsinterp)
library(forecast)



#interpolating a 1ish year of the Appourchaux PM2 data October 1996 - October 1997
pm2 <- solar2$imDat
yr <- 1577880


solar <- pm2[(yr/2):(yr/2+yr-1 + 1800)] #to account for loss of data when filtering
solar[which(solar == 0)] <- NA
N <- length(solar)


gapz <- findGaps(solar) #data in this block look pretty clean

s.chunk1 <- solar[1:(gapz[3,1]-5000)]
s.chunk2 <- solar[(gapz[3,1]-4999):(gapz[3,2]+5000)]
s.chunk3 <- solar[(gapz[3,2]+5001):N]

s.chunk1 <- na.interp(s.chunk1)
s.chunk2 <- interpolate(as.matrix(s.chunk2), gap = which(is.na(s.chunk2)), delT = 20)
s.chunk3 <- na.interp(s.chunk3)


Solar_20s <- c(s.chunk1, s.chunk2[[1]], s.chunk3)

save(Solar_20s, file = "Solar_20s.RData")
###########################################################

rm(list = ls())
#setwd("~/Documents/University/Thesis/Rstuff")
source("utility.R")
load("solar_interp.RData")

N <- length(solar.interp)
delT <- 20
#maybe don't need to bother removing a trend
#solar.dt <- detrend(solar.interp, W = 6/(N*delT), deltat = delT)
#solar.dm <- solar.dt - mean(solar.dt)
solar.dm <- solar.interp - mean(solar.interp)

#low-pass filter stuff
## 12:1 decimation filter -> 4 minute sampling time and nyquist of 2.08 mHz
# I have filtered with two different low pass filters. imp is using Slepians in what I think is the 
# normal way, and imp2,3 are using Dave's modified Slepians instead and 5 is a convex combination of 
# 2 and 3.

nflt <- 291
ndec <- 12
nyquist <- 0.5/(ndec*delT)

imp <- filterDesign(nflt = nflt, ndec = ndec-1)

solar.filt <- filter(solar.dm, filter = imp, sides = 2)
solar.filt <- solar.filt[-which(is.na(solar.filt))]

dec.ind <- seq(1,length(solar.filt), by = ndec)

solar.dec <- solar.filt[dec.ind]
solar <- solar.dec - mean(solar.dec)

L <- length(solar)

#rm(solar.dt,solar.dm,solar.filt5,solar.dec)

nFFT <- 2^ceiling(log2(10*L))
sp <- spec.mtm(solar, deltat = ndec*delT, dtUnits = 'second', Ftest = TRUE, plot = FALSE,
                nFFT = nFFT)
conf <- 1-1/L
sigF <- findLocalFMax(sp, conf)

cmv.tmp <- sp$mtm$cmv
cmv.tmp[-sigF] <- 0
cmv <- c(cmv.tmp, Conj(rev(cmv.tmp[-1])[-1]))
p.recon <- Re(fft(cmv, inverse = TRUE))[1:L]
res <- solar - p.recon

nw <- 10
k <- 2*nw - 1

sp.pilot.tmp <- spec.mtm(res, nw = nw, k = k, deltat = ndec*delT, dtUnits = 'second',
                          plot = FALSE, nFFT = nFFT)
sp.pilot <- c(sp.pilot.tmp$spec, rev(sp.pilot.tmp$spec[-1])[-1])

acvf <- Re(fft(sp.pilot, inverse = TRUE))/length(sp.pilot)

order <- 3 #taking this to be 20 leaves me with exactly 2 years worth of data (21 actually)
#20 was also the amount that made the prewhitening work on a month block. this should be like 2 or 3
gam <- acvf[1:(order+1)]
Gam <- toeplitz(gam[-order])
QR <- qr(Gam)
phi <- qr.coef(QR, y = gam[-1])


fit <- filter(res, filter = c(0,phi), sides = 1)
pw.tmp <- fit - (res + p.recon)
pw <- pw.tmp[-which(is.na(pw.tmp))]


sp.int <- spec.mtm(pw, deltat = ndec*delT, dtUnits = 'second', plot = FALSE, 
                    nFFT = nFFT)
tf <- fft(c(-1,phi, rep(0,nFFT-order-1)))[1:sp.int$mtm$nfreqs]
sp.pw <- sp.int
sp.pw$spec <- sp.int$spec / Mod(tf)^2


#################################################################
######## check on a block
#################################################################
blkLen <- 7*24*3600/(ndec*delT)
solar.chk <- solar1[1:blkLen]
solar.chk <- solar.chk - mean(solar.chk)
sr <- ndec*delT

sp.chk <- spec.mtm(solar.chk, nw = 4, k = 7, deltat = sr, dtUnits = 'second',
                   plot = FALSE, Ftest = TRUE, nFFT = 2^ceiling(log2(10*blkLen)))
sigF.chk <- findLocalFMax(sp.chk,1-1/blkLen)
cmv.tmp.chk <- sp.chk$mtm$cmv
cmv.tmp.chk[-sigF.chk] <- 0
cmv.chk <- c(cmv.tmp.chk, Conj(rev(cmv.tmp.chk[-1])[-1]))
p.recon.chk <- Re( fft(cmv.chk, inverse = TRUE) )[1:blkLen]
res.chk <- solar.chk - p.recon.chk

nw <- 10
k <- 2*nw - 1
sp.p.tmp <- spec.mtm(res.chk, nw = nw, k = k, deltat = sr, dtUnits = 'second',
                     plot = FALSE, nFFT = 2^ceiling(log2(10*blkLen)))
sp.p.chk <- c(sp.p.tmp$spec, rev(sp.p.tmp$spec[-1])[-1])

acvf.chk <- Re( fft(sp.p.chk, inverse = TRUE)/length(sp.p.chk))
tau.chk <- 2
gam.chk <- acvf.chk[1:tau.chk]
Gam.chk <- toeplitz(gam.chk[-tau.chk])
QR.chk <- qr(Gam.chk)
phi.chk <- qr.coef(QR.chk, y = gam.chk[-1])
fit.chk <- filter(res.chk, filter = c(0,phi.chk), sides = 1)

pw.tmp.chk <- fit.chk - (res.chk + p.recon.chk)
pw.chk <- pw.tmp.chk[-which(is.na(pw.tmp.chk))]
pw.chk <- pw.chk - mean(pw.chk)
sp.int.chk <- spec.mtm(pw.chk, nw = 4, k = 7, deltat = ndec*delT, dtUnits = 'second',
                       plot = FALSE, nFFT = 2^ceiling(log2(10*blkLen)), Ftest = TRUE)
tf.chk <- fft(c(-1,phi.chk,rep(0,sp.int.chk$mtm$nFFT-tau.chk)))[1:sp.int.chk$mtm$nfreqs]
sp.pw.chk <- sp.int.chk
sp.pw.chk$spec <- sp.int.chk$spec / Mod(tf.chk)^2


###########################################################
### build all these spectrograms and F spectrograms
###########################################################
#set up time blocks

N <- length(pw)
sr <- delT*ndec
blkLen <- 14*24*(3600/sr)
blx <- block(n = N, blockSize = blkLen, overlap = 0)
idx <- blx$startIdx
nBlx <- blx$nBlock

#multitaper stuff
nyquist <- 0.5/sr
nw <- 2
k <- 2*nw - 1
DW <- dpss(n = blkLen, nw = nw, k = k)
nFFT <- 2^ceiling(log2(10*blkLen))

#set up frequency blocks
nfreqs <- nFFT/2+1
freq.tmp <- seq(0, nyquist, length.out = nfreqs)
freq <- c(freq.tmp, nyquist + freq.tmp[1])
fblkLen <- (nfreqs-1)/64
fblx <- block(n = nfreqs, blockSize = fblkLen, overlap = 0)
fIdx <- fblx$startIdx
nfBlx <- fblx$nBlock

spectrogram <- F_spectrogram <- array(data = 0, dim = c(nBlx,nfBlx,fblkLen))

for(i in 1:nBlx){
  pw.sub <- pw[idx[i]:(idx[i]+blkLen-1)]
  sp.sub <- spec.mtm(pw.sub, deltat = sr, plot = FALSE, Ftest = TRUE, dtUnits = 'second', nw = nw,
                     k = k, dpssIN = DW, nFFT = nFFT)
  for(j in 1:nfBlx){
    sp.dr <- dropFreqs(sp.sub, minFreq = freq[fIdx[j]], maxFreq = freq[fIdx[j]+fblkLen])
    spectrogram[i,j,] <- sp.dr$spec
    F_spectrogram[i,j,] <- sp.dr$mtm$Ftest
  }
  
}
library(fields)

pdf(file = "spectrograms.pdf")
  for(j in 1:nfBlx){
    image.plot(x = 1:nBlx, y = freq[fIdx[j]:(fIdx[j]+fblkLen-1)], z = spectrogram[,j,], xlab = "block",
               ylab = "frequency")
    title(paste("spectrogram in the frequency range ", format(freq[fIdx[j]], digits=3, scientific=TRUE), 
                " to ", format(freq[fIdx[j]+fblkLen-1], digits=3, scientific=TRUE), " Hz"))
    
  }
dev.off()

pdf(file = "F_spectrograms.pdf")
  for(j in 1:nfBlx){
    image.plot(x = 1:nBlx, y = freq[fIdx[j]:(fIdx[j]+fblkLen-1)], z = F_spectrogram[,j,], xlab = "block",
               ylab = "frequency")
    title(paste("F spectrogram in the frequency range ", format(freq[fIdx[j]], digits=3, scientific=TRUE),
                " to ", format(freq[fIdx[j]+fblkLen-1], digits=3, scientific=TRUE), " Hz"))
    
  }
dev.off()

Fquant1 <- qf(0.99, 2, 12)
Fquant2 <- qf(1-1/blkLen, 2, 12)

F_spectrogram.sig1 <- F_spectrogram
for(i in 1:nBlx){
  for(j in 1:nfBlx){
    F_spectrogram.sig1[i,j,which(F_spectrogram.sig1[i,j,] < Fquant1)] <- NA
  }
}

F_spectrogram.sig2 <- F_spectrogram
for(i in 1:nBlx){
  for(j in 1:nfBlx){
    F_spectrogram.sig2[i,j,which(F_spectrogram.sig2[i,j,] < Fquant2)] <- NA
  }
}

pdf(file = "signif_F_spectrograms1.pdf")
for(j in 1:nfBlx){
  image.plot(x = 1:52, y = freq[fIdx[j]:(fIdx[j]+fblkLen-1)], z = F_spectrogram.sig1[,j,], xlab = "block",
             ylab = "frequency")
  title(paste("F spectrogram in the frequency range ", format(freq[fIdx[j]], digits=3, scientific=TRUE),
              " to ", format(freq[fIdx[j]+fblkLen-1], digits=3, scientific=TRUE), " Hz"))
  
}
dev.off()

pdf(file = "signif_F_spectrograms2.pdf")
for(j in 1:nfBlx){
  image.plot(x = 1:52, y = freq[fIdx[j]:(fIdx[j]+fblkLen-1)], z = F_spectrogram.sig2[,j,], xlab = "block",
             ylab = "frequency")
  title(paste("F spectrogram in the frequency range ", format(freq[fIdx[j]], digits=3, scientific=TRUE),
              " to ", format(freq[fIdx[j]+fblkLen-1], digits=3, scientific=TRUE), " Hz"))
  
}
dev.off()

###########################################################
## I'll do it with the unwhitened data too just to compare
###########################################################
spectrogram.uw <- F_spectrogram.uw <- array(data = 0, dim = c(nBlx,nfBlx,fblkLen))

for(i in 1:nBlx){
  solar.sub <- solar[idx[i]:(idx[i]+blkLen-1)]
  sp.sub <- spec.mtm(solar.sub, deltat = sr, plot = FALSE, Ftest = TRUE, dtUnits = 'second', nw = nw,
                     k = k, dpssIN = DW, nFFT = nFFT)
  for(j in 1:nfBlx){
    sp.dr <- dropFreqs(sp.sub, minFreq = freq[fIdx[j]], maxFreq = freq[fIdx[j]+fblkLen])
    spectrogram.uw[i,j,] <- sp.dr$spec
    F_spectrogram.uw[i,j,] <- sp.dr$mtm$Ftest
  }
  
}

pdf(file = "uw_spectrograms.pdf")
  for(j in 1:nfBlx){
    image.plot(x = 1:52, y = freq[fIdx[j]:(fIdx[j]+fblkLen-1)], z = spectrogram.uw[,j,], xlab = "block",
               ylab = "frequency")
    title(paste("unwhitened spectrogram in the band ", format(freq[fIdx[j]], digits=3, scientific=TRUE), 
                " to ", format(freq[fIdx[j]+fblkLen-1], digits=3, scientific=TRUE), " Hz"))
    
  }
dev.off()

pdf(file = "uw_F_spectrograms.pdf")
  for(j in 1:nfBlx){
    image.plot(x = 1:52, y = freq[fIdx[j]:(fIdx[j]+fblkLen-1)], z = F_spectrogram.uw[,j,], xlab = "block",
               ylab = "frequency")
    title(paste("unwhitened F spectrogram in the band ", format(freq[fIdx[j]], digits=3, scientific=TRUE),
                " to ", format(freq[fIdx[j]+fblkLen-1], digits=3, scientific=TRUE), " Hz"))
    
  }
dev.off()

################################################
## plot spectrograms on a log scale?
################################################

#not sure how to get the legend in the way I want. similar to axis on mtm spectrum plots

pdf(file = "log_spectrograms.pdf")
  for(j in 1:nfBlx){
    image.plot(x = 1:52, y = freq[fIdx[j]:(fIdx[j]+fblkLen-1)], z = log10(spectrogram.uw[,j,]),
               xlab = "block", ylab = "frequency")
    title(paste("unwhitened log-spectrogram in the band ", 
                format(freq[fIdx[j]], digits=3, scientific=TRUE), " to ",
                format(freq[fIdx[j]+fblkLen-1], digits=3, scientific=TRUE), " Hz"))
    
  }
dev.off()

pdf(file = "log_F_spectrograms.pdf")
  for(j in 1:nfBlx){
    image.plot(x = 1:52, y = freq[fIdx[j]:(fIdx[j]+fblkLen-1)], z = log10(F_spectrogram.uw[,j,]),
               xlab = "block", ylab = "frequency")
    title(paste("unwhitened log - F spectrogram in the band ", 
                format(freq[fIdx[j]], digits=3, scientific=TRUE), " to ",
                format(freq[fIdx[j]+fblkLen-1], digits=3, scientific=TRUE), " Hz"))
    
  }
dev.off()


#############################################################
### low frequency stuff
#############################################################
rm(list = ls())
source("utility.R")
load("solar_pw_lf_series.RData")

N <- length(solar.pw.lf)
ndec.lf <- 48
delT <- 20
sr <- delT*ndec.lf
blkLen <- 14*24*(3600/sr)
blx <- block(n = N, blockSize = blkLen, overlap = 0)
idx <- blx$startIdx
nBlx <- blx$nBlock

#multitaper stuff
nyquist <- 0.5/sr
nw <- 4
k <- 2*nw - 1
DW <- dpssm(n = blkLen, nw = nw, k = k)
nFFT <- 2^ceiling(log2(10*blkLen))

#set up frequency blocks
nfreqs <- nFFT/2+1
freq.tmp <- seq(0, nyquist, length.out = nfreqs)
freq <- c(freq.tmp, nyquist + freq.tmp[2])
fblkLen <- (nfreqs-1)/64
fblx <- block(n = nfreqs, blockSize = fblkLen, overlap = 0)
fIdx <- fblx$startIdx
nfBlx <- fblx$nBlock

lf_spectrogram <- lf_F_spectrogram <- array(data = 0, dim = c(nBlx,nfBlx,fblkLen))

for(i in 1:nBlx){
  pw.sub <- pw[idx[i]:(idx[i]+blkLen-1)]
  sp.sub <- spec.mtm(pw.sub, deltat = sr, plot = FALSE, Ftest = TRUE, dtUnits = 'second', nw = nw,
                     k = k-2, dpssIN = DW, nFFT = nFFT)
  for(j in 1:nfBlx){
    sp.dr <- dropFreqs(sp.sub, minFreq = freq[fIdx[j]], maxFreq = freq[fIdx[j]+fblkLen])
    lf_spectrogram[i,j,] <- sp.dr$spec
    lf_F_spectrogram[i,j,] <- sp.dr$mtm$Ftest
  }
  
}


pdf(file = "lf_spectrograms.pdf")
for(j in 1:nfBlx){
  image.plot(x = 1:nBlx, y = freq[fIdx[j]:(fIdx[j]+fblkLen-1)], z = lf_spectrogram[,j,], xlab = "block",
             ylab = "frequency")
  title(paste("spectrogram in the frequency range ", format(freq[fIdx[j]], digits=3, scientific=TRUE), 
              " to ", format(freq[fIdx[j]+fblkLen-1], digits=3, scientific=TRUE), " Hz"))
  
}
dev.off()

pdf(file = "lf_F_spectrograms.pdf")
for(j in 1:nfBlx){
  image.plot(x = 1:nBlx, y = freq[fIdx[j]:(fIdx[j]+fblkLen-1)], z = lf_F_spectrogram[,j,], xlab = "block",
             ylab = "frequency")
  title(paste("F spectrogram in the frequency range ", format(freq[fIdx[j]], digits=3, scientific=TRUE),
              " to ", format(freq[fIdx[j]+fblkLen-1], digits=3, scientific=TRUE), " Hz"))
  
}
dev.off()

### let's try out my useless F statistic
modSP <- spec.mtm(solar.pw.lf, deltat = sr, nw = 5.5, k = 10, dtUnits = 'second',
                  plot = FALSE, returnInternals = TRUE, Ftest = TRUE)

#set up frequency blocks
w <- modSP$mtm$nw / (N*sr)
nfreqs <- length(which(freq > w))
freq <- modSP$freq

#analytic signal only exists for f > W (apparently)
numBands <- 50 #for memory allocation purposes. can only work with a small band at once.
fblkLen <- (nfreqs-1)/numBands
fblx <- block(n = nfreqs, blockSize = fblkLen, overlap = 0)
fIdx <- fblx$startIdx
nfBlx <- fblx$nBlock

M <- 4
modF <- NULL
freqF <- NULL

## THIS IS USING THE CMVS. CHANGE BACK BEFORE RUNNING AGAIN
for(i in 1:nfBlx){
  subBand <- which(freq > w)[fIdx[i]:(fIdx[i]+fblkLen-1)]
  modF.tmp <- ModulatedFcmv(spec = modSP, mxdeg = M, subBand = subBand)
  modF <- cbind(modF, modF.tmp$ModF)
  freqF <- c(freqF, modF.tmp$freq)
  rm(modF.tmp)
}

modF.stuff <- list(modF = modF, freqF = freqF)
save(modF.stuff, file = "modulated_F_solar.RData")

modF.plotting <- modF
for(m in 1:(M+1)){
  modF.plotting[m, which(modF.plotting[m,] < 1)] <- 1
}


png(filename = "modFsolar.png", width = 600, height = 400)
par(mfrow = c(2,2))
for(m in 2:(M+1)){
  plot(freqF, modF.plotting[m,], type = 'l', log = 'y', ylab = '', xlab = 'Frequency')
  abline(h = qf(1-1/N, 1, modSP$mtm$k - m), lty = 2, col = 2)
  title(paste("Degree", m-1)) 
}
dev.off()

#########
# need to set this up for the 4-minute data.
#########