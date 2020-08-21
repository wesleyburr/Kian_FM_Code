source("solarUtility.R")


load("Solar_20s.RData")


delT <- 20
ndec <- 5   #100 second sampling rate
sr <- ndec*delT
nyquist = 0.5/sr

yr <- 1577880
nflt <- 301

solar.dm <- Solar_20s - mean(Solar_20s)

imp <- filterDesign(nflt = nflt, wflt = 0.5/ndec + 1/nflt)
solar.tmp <- filter(solar.dm, filter = imp, sides = 2)
solar.tmp <- solar.tmp[-which(is.na(solar.tmp))]
dec <- seq(1, length(solar.tmp), by = ndec)
solar <- solar.tmp[dec] - mean(solar.tmp[dec])
hp <- solar - filter(solar, imp, sides = 2)
hp <- hp[-which(is.na(hp))]



N <- length(hp)
nFFT <- 2^ceiling(log2(10*N))
L <- 300
M <- 6
conf <- 1-1/N


sp.Solar <- spec.mtm(hp, nw = 5, k = 8, plot = FALSE, Ftest = TRUE,  dtUnits = 'second',
                     returnInternals = TRUE, deltat = sr, returnZeroFreq = FALSE, nFFT = nFFT)
sp.Solar.dr <- dropFreqs(sp.Solar, 0.0011, nyquist)
rayleigh <- 1/(nFFT*sr)
freq <- seq(0,nyquist, by = rayleigh)
subBand <- which(freq > 0.0011 & freq < nyquist)
freqF <- freq[subBand]
nfreqs <- length(subBand)
NW <- seq(5, 11, by = 0.5)
numNW <- length(NW)

mF3array <- array(0, dim = c(M, nfreqs, numNW))
HFarray <- array(0, dim = c(nfreqs, numNW))
sigFlist <- list()
sigMF3list <- list()

start.time <- Sys.time()
for(n in 1:numNW){
  K <- 2*NW[n] - 2
  DW <- dpss(n = N, nw = NW[n], k = K)
  DW$v <- DW$v * sqrt(sr)
  DS <- dpss(n = L, nw = NW[n], k = K)
  Vdot <- dpssp11R(DS,NW[n])
  AP <- dpssap(DS$v, M, 1)
  
  tapered <- DW$v * hp
  padded <- rbind(tapered, matrix(0, nrow = nFFT-N, ncol = K))
  yk <- mvfft(padded)[subBand, , drop = FALSE]
  
  HF <- HarmonicF(yk, DW$v)
  sigF <- findLocalFMax2(HF$HF, conf, K)
  HFarray[,n] <- HF$HF
  
  if(length(sigF) > 0){
    sigFlist[[n]] <- sigF
    cmv.tmp <- HF$cmv
    cmv.tmp[-sigF] <- 0
    cmv <- c(cmv.tmp, Conj(rev(cmv.tmp[-1])[-1]))
    p.recon <- Re(fft(cmv, inverse = TRUE))[1:N]
    hp.res <- hp - p.recon
    tapered <- DW$v * hp.res
    padded <- rbind(tapered, matrix(0, nrow = nFFT-N, ncol = K))
    yk <- mvfft(padded)[subBand, , drop = FALSE]
  } else {
    sigFlist[[n]] <- NA
  }
  
  
  mF3 <- ModulatedF18(yk, Vdot, AP, DS)
  
  mF3array[,,n] <- mF3
  
  sigMF3list[[n]] <- findLocalFMaxM(mF3, K, conf)
  
  
}
runtime <- Sys.time() - start.time
runtime

SolarFstats <- list(HF = HFarray, mF3 = mF3array)
save(SolarFstats, file = "SolarFstats.RData")


sink("sigFreqsHFsolar.txt", append = FALSE)
cat("------------------------\n")
for(n in 1:numNW){
  K <- 2*NW[n] - 2
  cat(paste("NW = ", NW[n], "\n"))
  cat(paste("Significant frequencies in microhertz of HF at the ", conf, "level: \n"))
  if(any(is.na(sigFlist[[n]]))){
    print("NA")
    cat("\n")
  } else {
    print(freqF[sigFlist[[n]]]*1e6)
    cat("\n")
    print(HFarray[sigFlist[[n]],n])
    cat(paste("Relevant quantile: ", qf(conf, 2, 2*K - 2), "\n"))
    cat("\n")
  }
  cat("------------------------\n")
}
sink()

sink("sigFreqsmF3solar.txt", append = FALSE)

for(n in 1:numNW){
  sigmF3s <- sigMF3list[[n]]
  K <- 2*NW[n]-2
  cat(paste("NW = ", NW[n], "\n"))
  cat(paste("Significant frequencies in microhertz of modified F3 at the ", conf, "level: \n"))
  for(m in 1:M){
    cat(paste("Degree ", m, ": \n"))
    if(any(is.na(sigmF3s[[m]]))){
      print("NA")
      cat("\n")
    } else {
      print(freqF[sigmF3s[[m]]]*1e6)
      cat("\n")
      print(mF3array[m,sigmF3s[[m]],n])
      cat(paste("Relevant quantile: ", qf(conf, 1, K - m), "\n"))
    }
    cat("------------------------\n")
  }
}
sink()








alpha <- 1/N

sigFreqHM <- array(data = 0, dim = c(M, nfreqs, numNW))
for(m in 1:M){
  for(n in 1:numNW){
    K <- 2*NW[n] - 2
    sigFreqHM[m,which(mF3array[m,,n] > qf(1-alpha, 1, K-m)),n] <- n
  }
}

cols <- rainbow(numNW)

for(m in 1:M){
  par(oma = c(4,4,1,4), mar = c(0,0,1,1), mfrow = c(1,1))
  if(length(which(sigFreqHM[m,,1]>0)) > 0){
    plot(freqF[which(sigFreqHM[m,,1] > 0)], sigFreqHM[m,which(sigFreqHM[m,,1] > 0),1], pch = 19,
         xlim = range(freqF), yaxt='n', ylim = c(1,numNW), col = cols[1], xaxt = 'n')
    grid()
  } else {
    plot(0, 0, col = 0, ylim = c(1, numNW), xlim = range(freqF), yaxt = 'n', xaxt = 'n')
    grid()
  }
  title(paste("Degree ", m))
  mtext(paste("Frequency (mHz)"), side = 1, outer = TRUE, line = 2.2)
  mtext("NW", side = 2, outer = TRUE, line = 2.3)
  #if(m == 1 | m == 4){
  axis(2, at = 1:numNW, labels = NW)
  #}
  #if(m > 3){
  axis(1, at = c(0.001, 0.002, 0.003, 0.004, 0.005), labels = 1:5)
  #}
  for(n in 2:numNW){
    points(freqF[which(sigFreqHM[m,,n] > 0)], sigFreqHM[m,which(sigFreqHM[m,,n] > 0),n], pch = 19,
           col = cols[n])
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("right", legend = NW, col = cols,
         pch = 19, horiz = FALSE, xpd = TRUE, bty = 'n', inset = c(0.01,0))
}




######################################
## tables of proportion of frequencies
######################################

alphas <- c(0.1, 0.01, 0.001, 0.0001, 0.00001, 1/N)
A <- length(alphas)

totProp <- array(0, dim = c(M, numNW, A))

for(a in 1:A){
  for(n in 1:numNW){
    K <- 2*NW[n] - 2
    maxes <- findLocalFMaxM(mF3array[,,n], K, 1-alphas[a])
    for(m in 1:M){
      totProp[m,n,a] <- length(maxes[[m]])/nfreqs
    }
  }
}

for(a in 1:A){
  plot(totProp[1,,a], type = 'b', ylim = range(totProp[,,a]), lwd = 2, log = 'y')
  for(m in 2:M){
    lines(totProp[m,,a], type = 'b', col = m, lwd = 2)
  }
}