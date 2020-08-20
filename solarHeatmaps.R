source("utility.R")
load("solar_pw_lf_series.RData")
solar <- solar.pw.lf
rm(solar.pw.lf)

N <- length(solar)
delT <- 20*12

#physical bandwidth is between about 125 and 250 nHz, depending on NW. NW will be chosen to be
#relatively high.



wk <- 7*24*60/4
mth <- 30*24*60/4

blkW <- block(N, wk, overlap = 0)
blkM <- block(N, mth, overlap = 0)

pdf("solar_weeks.pdf", width = 11, height = 8.5)
  par(oma = c(4, 4, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
  par(mar = c(0, 0, 1, 1))
  for(w in 1:blkW$nBlock){
    plot(solar[blkW$startIdx[w]:(blkW$startIdx[w] + blkW$blockSize - 1)], type = 'l', xlab = '',
         ylab = '', ylim = range(solar))
    title(paste("Week ", w, " of ", blkW$nBlock))
  }
dev.off()
pdf("solar_months.pdf", width = 11, height = 8.5)
  par(oma = c(4, 4, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
  par(mar = c(0, 0, 1, 1))
  for(m in 1:blkM$nBlock){
    plot(solar[blkM$startIdx[m]:(blkM$startIdx[m] + blkM$blockSize - 1)], type = 'l', xlab = '',
         ylab = '', ylim = range(solar))
    title(paste("Month ", m, " of ", blkM$nBlock))
  }
dev.off()



NW <- 8
K <- 2*NW - 2
DWw <- dpss(n = wk, nw = NW, k = K)
DWm <- dpss(n = mth, nw = NW, k = K)

VdotW <- dpssp11R(DWw,NW)
VdotM <- dpssp11R(DWm,NW)

M <- 6

APw <- dpssap(DWw$v, M, alpha = 1)
APm <- dpssap(DWm$v, M, alpha = 1)

nFFTw <- 2^ceiling(log2(2*wk))
nfreqsW <- nFFTw/2

nFFTm <- 2^ceiling(log2(2*mth))
nfreqsM <- nFFTm/2

F1arrayW <- F2arrayW <- F3arrayW <- array(NA, dim = c(M+1, nfreqsW-1, blkW$nBlock))
F1arrayM <- F2arrayM <- F3arrayM <- array(NA, dim = c(M+1, nfreqsM-1, blkM$nBlock))
HFarrayW <- array(NA, dim = c(nfreqsW, blkW$nBlock))
HFarrayM <- array(NA, dim = c(nfreqsM, blkM$nBlock))

start.time <- Sys.time()
for(w in 1:blkW$nBlock){
  solarW <- solar[blkW$startIdx[w]:(blkW$startIdx[w] + blkW$blockSize - 1)]
  solarW <- detrend(solarW, W = 6/(delT*length(solarW)), deltat = delT)
  sp.W <- spec.mtm(solarW, deltat = delT, nw = NW, k = K, dpssIN = DWw, plot = FALSE,
                   Ftest = TRUE, returnZeroFreq = FALSE, returnInternals = TRUE, nFFT = nFFTw)
  modFw <- ModulatedF8(sp.W, M, adaptive = FALSE, subBand = 1:(nfreqsW-1), derivIN = VdotW,
                       apIN = APw, dpssIN = DWw)
  
  HFarrayW[,w] <- sp.W$mtm$Ftest
  F1arrayW[,,w] <- modFw$F1
  F2arrayW[,,w] <- modFw$F2
  F3arrayW[,,w] <- modFw$F3
}

for(m in 1:blkM$nBlock){
  solarM <- solar[blkM$startIdx[m]:(blkM$startIdx[m] + blkM$blockSize - 1)]
  solarM <- detrend(solarM, W = 6/(delT*length(solarM)), deltat = delT)
  sp.M <- spec.mtm(solarM, deltat = delT, nw = NW, k = K, dpssIN = DWm, plot = FALSE,
                   Ftest = TRUE, returnZeroFreq = FALSE, returnInternals = TRUE, nFFT = nFFTm)
  modFm <- ModulatedF8(sp.M, M, adaptive = FALSE, subBand = 1:(nfreqsM-1), derivIN = VdotM,
                       apIN = APm, dpssIN = DWm)
  
  HFarrayM[,m] <- sp.M$mtm$Ftest
  F1arrayM[,,m] <- modFm$F1
  F2arrayM[,,m] <- modFm$F2
  F3arrayM[,,m] <- modFm$F3
}
runtime <- Sys.time() - start.time
runtime

pdf("weekly_heatmaps.pdf", height = 8.5, width = 11)
  for(p in 1:(M+1)){
    image.plot(x = sp.W$freq[-nfreqsW], y = 1:blkW$nBlock, z = F1arrayW[p,,])
    title(paste("F1 degree ", p-1))
    
    image.plot(x = sp.W$freq[-nfreqsW], y = 1:blkW$nBlock, z = F2arrayW[p,,])
    title(paste("F2 degree ", p-1))
    
    image.plot(x = sp.W$freq[-nfreqsW], y = 1:blkW$nBlock, z = F3arrayW[p,,])
    title(paste("F3 degree ", p-1))
  }

  image.plot(x = sp.W$freq, y = 1:blkW$nBlock, z = HFarrayW)
dev.off()

pdf("monthly_heatmaps.pdf", height = 8.5, width = 11)
  for(p in 1:(M+1)){
    image.plot(x = sp.M$freq[-nfreqsM], y = 1:blkM$nBlock, z = F1arrayM[p,,])
    title(paste("F1 degree ", p-1))
    
    image.plot(x = sp.M$freq[-nfreqsM], y = 1:blkM$nBlock, z = F2arrayM[p,,])
    title(paste("F2 degree ", p-1))
    
    image.plot(x = sp.M$freq[-nfreqsM], y = 1:blkM$nBlock, z = F3arrayM[p,,])
    title(paste("F3 degree ", p-1))
  }
  
  image.plot(x = sp.M$freq, y = 1:blkM$nBlock, z = HFarrayM)
dev.off()