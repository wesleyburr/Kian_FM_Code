#  A KAY AND MARPLE TYPE EXAMPLE
# SET UP A FEW MODULATED SIGNALS WITH DIFFERENT DEGREES AND DIFFERENT CARRIER FREQUENCIES
# TO CHECK THE PERFORMANCE

# THE MODULATION BANDWIDTH WILL BE THE SAME FOR ALL SIGNALS
##########################################################################################


source("utility.R")



N <- 2048
n <- 0:(N-1)
nFFT <- 2^ceiling(log2(2*N))
nfreqs <- nFFT/2-1
subBand <- 1:nfreqs
freq <- seq(1/nFFT,0.5, by = 1/nFFT)
freqF <- seq(1/nFFT, 0.5-1/nFFT, by = 1/nFFT)

NW <- 5:11
K <- 2*NW - 1
numNW <- length(NW)
W <- NW / N

sigmat2 <- 1  
sigmat <- sqrt(sigmat2)

ar <- c(0.5,0.3,-0.1)
ma <- c(0.6)#, 0.4, 0.2, -0.1)
ARMA <- list(ar = ar, ma = ma)


#signal setup
#amplitudes
A1 <- 0.5
A2 <- 0.4
A3 <- 0.2
A4 <- 0.2
A5 <- 0.3
A6 <- 0.2

#carrier frequencies
f1 <- freq[which.min(abs(freq-0.1))]
f2 <- freq[which.min(abs(freq-0.3))]
f3 <- freq[which.min(abs(freq-0.31))]
f4 <- freq[which.min(abs(freq-0.4))]
f5 <- freq[which.min(abs(freq-0.2))]
f6 <- freq[which.min(abs(freq-0.25))]

#time and band parameters
B1 <- 0.0025
B2 <- 0.0035
B3 <- 0.0015
stm11 <- 2.0/(N-1)
tt <- n * stm11 - 1.0  # this runs from -1 to 1

#modulating functions
FMpoly1 <- B2 * tt                       #linear
FMpoly2 <- B3 * (1.0 - 2.0 * tt^2)     #quadratic
FMpoly3 <- B3 * (4*tt^3 - 3*tt)        #cubic
FMpoly4 <- 3/5*B1 * (-1 + (32/3)*tt^2 - (32/3)*tt^4) # quartic
FMcos <- -B2 * cos(2*pi*(tt+1))     # cos( (d*pi/2) * (t+1) ) goes from 1 to -1 d times as 

#FMpoly5 <- B*fbw * (9.46*tt -19.8*tt^3 + 11.34*tt^5)

                                       # t goes from -1 to 1
Pref1 <- cumsum(FMpoly1)*2*pi
Pref2 <- cumsum(FMpoly2)*2*pi
Pref3 <- cumsum(FMpoly3)*2*pi
Pref4 <- cumsum(FMpoly4)*2*pi
Pref5 <- cumsum(FMcos)*2*pi
PhMod1 <- 2*pi*f1*n + Pref1
PhMod2 <- 2*pi*f2*n + Pref2
PhMod3 <- 2*pi*f3*n + Pref3
PhMod4 <- 2*pi*f4*n + Pref4
PhMod5 <- 2*pi*f5*n + Pref5

sig <- A1*cos(PhMod1) + A2*cos(PhMod2) + A3*cos(PhMod3) + 
  A4*cos(PhMod4) + A5*cos(PhMod5) + A6*cos(2*pi*f6*n)


num.sim <- 10000
L <- 100
M <- 6




HFarray <- array(0, dim = c(nfreqs, num.sim, numNW))
F3array <- array(0, dim = c(M+1, nfreqs, num.sim, numNW))
mF1array <- mF3array <- array(0, dim = c(M, nfreqs, num.sim, numNW))


  
band1 <- which(freq < f1 + B2 & freq > f1 - B2)
band2 <- which(freq < f2 + B3 & freq > f2 - B3)
band3 <- which(freq < f3 + B3 & freq > f3 - B3)
band4 <- which(freq < f4 + max(FMpoly4) & freq > f4 + min(FMpoly4))
band5 <- which(freq < f5 + B2 & freq > f5 - B2)
band6 <- which(freq < f6 + B1 & freq > f6 - B1)

band <- list(band1, band2, band3, band4, band5, band6)

numBand <- length(band)
conf <- 1-1/N
Fcounter <- array(0, dim = c(numBand, numNW))
F3counter <- array(0, dim = c(M+1, numBand, numNW))
mF1counter <- mF3counter <- array(0, dim = c(M, numBand, numNW))

sigFreqsHF <- array(NA, dim = c(num.sim, numNW))
sigFreqsF3 <- array(NA, dim = c(M+1, num.sim, numNW))
sigFreqsmF1 <- sigFreqsmF3 <- array(NA, dim = c(M, num.sim, numNW))

sigFreqsHFc <- array(NA, dim = c(num.sim, numNW, numBand))
sigFreqsF3c <- array(NA, dim = c(M+1, num.sim, numNW, numBand))
sigFreqsmF1c <- sigFreqsmF3c <- array(NA, dim = c(M, num.sim, numNW, numBand))

seed <- 666
set.seed(seed)


start.time <- Sys.time()
for(j in 1:numNW){
  DW <- dpss(n = N, nw = NW[j], k = K[j])
  DS <- dpss(n = L, nw = NW[j], k = K[j])
  Vdot <- dpssp11R(DS,NW[j])
  AP <- dpssap(DS$v, maxdeg = M, alpha = 1)
  for(i in 1:num.sim){
    noiz.tmp <- rgumbel(N, scale = sigmat * 6/pi^2) # this should have been sqrt(6)/pi
    noiz <- arima.sim(model=ARMA, n = N, innov = noiz.tmp)
    #noiz <- rnorm(N, sd = sigmat)
    X <- sig + as.numeric(noiz)
    
    tapered <- DW$v * X
    pad <- rbind(tapered, matrix(0, nrow = nFFT-N, ncol = K[j]))
    yk <- mvfft(pad)[subBand+1, , drop = FALSE]
    
    HFarray[,i,j] <- HarmonicF(yk, DW$v)
    
    modF <- ModulatedF17(yk, derivIN = Vdot, apIN = AP, dpssIN = DS)
    
    F3array[,,i,j] <- modF$F3
    mF1array[,,i,j] <- modF$mF1
    mF3array[,,i,j] <- modF$mF3
    
    
    
    
    
      for(m in 1:(M+1)){
        if(max(F3array[m,,i,j]) > qf(conf,1,K[j]-m)){
          sigFreqsF3[m,i,j] <- freq[which.max(F3array[m,,i,j])]
          for(l in 1:numBand){
            if(any(which.max(F3array[m,,i,j]) == band[[l]])){
              F3counter[m,l,j] <- F3counter[m,l,j] + 1
              sigFreqsF3c[m,i,j,l] <- freq[which.max(F3array[m,,i,j])]
            }
          }
        }
        if(m < M+1){
          if(max(mF1array[m,,i,j]) > qf(conf,m,K[j]-m)){
            sigFreqsmF1[m,i,j] <- freq[which.max(mF1array[m,,i,j])]
            for(l in 1:numBand){
              if(any(which.max(mF1array[m,,i,j]) == band[[l]])){
                mF1counter[m,l,j] <- mF1counter[m,l,j] + 1
                sigFreqsmF1c[m,i,j,l] <- freq[which.max(mF1array[m,,i,j])]
              }
            }
          }
          if(max(mF3array[m,,i,j]) > qf(conf,1,K[j]-m)){
            sigFreqsmF3[m,i,j] <- freq[which.max(mF3array[m,,i,j])]
            for(l in 1:numBand){
              if(any(which.max(mF3array[m,,i,j]) == band[[l]])){
                mF3counter[m,l,j] <- mF3counter[m,l,j] + 1
                sigFreqsmF3c[m,i,j,l] <- freq[which.max(mF3array[m,,i,j])]
              }
            }
          }
        }
      }
    
      if(max(HFarray[,i,j]) > qf(conf,2,2*K[j]-2)){
        sigFreqsHF[i,j] <- freq[which.max(HFarray[,i,j])]
        for(l in 1:numBand){
          if(any(which.max(HFarray[,i,j]) == band[[l]])){
            Fcounter[l,j] <- Fcounter[l,j] + 1
            sigFreqsHFc[i,j,l] <- freq[which.max(HFarray[,i,j])]
          }
        }
      }
  }
}


runtime <- Sys.time() - start.time
runtime

Fcounter <- Fcounter / num.sim
F3counter <- F3counter / num.sim
mF1counter <- mF1counter / num.sim
mF3counter <- mF3counter / num.sim

#kmTestStats <- list(F3 = F3array, mF1 = mF1array, mF3 = mF3array)
#save(kmTestStats, file = "kmTestStatsGWN2.RData")






par(mfrow = c(1,1))
par(oma = c(4,4,0,4))
par(mar = c(0,0,1,1))
plot(FMpoly1, type = 'l', ylab = '', lwd = 2)
lines(FMpoly2, col = 2, lwd = 2)
lines(FMpoly3, col = 3, lwd = 2)
lines(FMpoly4, col = 4, lwd = 2)
lines(FMcos, col = 5, lwd = 2)
abline(h = c(-B1,B1), col = 6, lty = 2)
abline(h = c(-B2,B2), col = 6, lty = 2)
abline(h = c(-B3,B3), col = 6, lty = 2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("right", legend = expression(phi[1], phi[2], phi[3], phi[4], phi[5]), col = 1:5,
       lty = 1, horiz = FALSE, xpd = TRUE, bty = 'n', inset = c(0.01,0), lwd = 2)
mtext("Modulating functions", side = 2, outer = TRUE, line = -1.75)


#layout(matrix(c(1,2,3,4,4,5), nrow = 2, ncol = 3, byrow = TRUE), widths = c(1,1.25), heights = c(1,1))
par(mfrow = c(2,3))
par(oma = c(4,4,0,5))
par(mar = c(0,0,1,1))
plot(NW, F3counter[2,1,], type = 'b', ylim = c(0,1), xaxt = 'n', las = 1, lwd = 2)
grid()
lines(NW, mF1counter[1,1,], type = 'b', col = 2, lwd = 2)
lines(NW, mF3counter[1,1,], type = 'b', col = 3, lwd = 2)
lines(NW, Fcounter[1,], type = 'b', col = 4, lwd = 2)
title(paste("linear FM, carrier ", round(f1,2)))

plot(NW, F3counter[3,2,], type = 'b', ylim = c(0,1), xaxt = 'n', yaxt = 'n', lwd = 2)
grid()
lines(NW, mF1counter[2,2,], type = 'b', col = 2, lwd = 2)
lines(NW, mF3counter[2,2,], type = 'b', col = 3, lwd = 2)
lines(NW, Fcounter[2,], type = 'b', col = 4, lwd = 2)
title(paste("quadratic FM, carrier ", round(f2,2)))

plot(NW, F3counter[4,3,], type = 'b', ylim = c(0,1), xaxt = 'n', yaxt = 'n', lwd = 2)
grid()
#axis(4, las = 1)
lines(NW, mF1counter[3,3,], type = 'b', col = 2, lwd = 2)
lines(NW, mF3counter[3,3,], type = 'b', col = 3, lwd = 2)
lines(NW, Fcounter[3,], type = 'b', col = 4, lwd = 2)
title(paste("cubic FM, carrier ", round(f3,2)))

plot(NW, F3counter[5,4,], type = 'b', ylim = c(0,1), yaxt = 'n', lwd = 2)
grid()
axis(2, las = 1)
lines(NW, mF1counter[4,4,], type = 'b', col = 2, lwd = 2)
lines(NW, mF3counter[4,4,], type = 'b', col = 3, lwd = 2)
lines(NW, Fcounter[4,], type = 'b', col = 4, lwd = 2)
title(paste("quartic FM, carrier ", round(f4,2)))

plot(NW, F3counter[7,5,], type = 'b', ylim = c(0,1), yaxt = 'n', lwd = 2)
grid()
lines(NW, mF1counter[6,5,], type = 'b', col = 2, lwd = 2)
lines(NW, mF3counter[6,5,], type = 'b', col = 3, lwd = 2)
lines(NW, Fcounter[5,], type = 'b', col = 4, lwd = 2)
title(paste("sinusoidal FM, carrier ", round(f5,2)))

plot(NW, F3counter[1,6,], type = 'b', ylim = c(0,1), yaxt = 'n', lwd = 2)
grid()
#axis(4, las = 1)
lines(NW, Fcounter[6,], type = 'b', col = 4, lwd = 2)
#legend(x=8,y=0.8,  horiz = FALSE, bty = 'n')
title(paste("no FM, carrier ", round(f6,2)))
mtext("NW", side = 1, outer = TRUE, line = 2.2)
mtext("Detection probability", side = 2, outer = TRUE, line = 2.3)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("right", legend = expression(F[3], tilde(F)[1], tilde(F)[3], HF), col = 1:4, lty = 1, lwd = 2, horiz = FALSE, xpd = TRUE, 
       bty = 'n', inset = c(0.01,0), cex = 1.3)
#par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
#plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")




#########################################################
#in band, wrong degree

#mF1

#layout(matrix(c(1,2,3,4,4,5), nrow = 2, ncol = 3, byrow = TRUE), widths = c(1,1.25), heights = c(1,1))
par(mfrow = c(2,3))
par(oma = c(4,4,0,5))
par(mar = c(0,0,1,1))

plot(NW, mF1counter[1,1,], type = 'b', ylim = c(0,1), xaxt = 'n', lwd = 2, las = 1)
grid()
lines(NW, mF1counter[2,1,], type = 'b', col = 2, lwd = 2)
lines(NW, mF1counter[3,1,], type = 'b', col = 3, lwd = 2)
lines(NW, mF1counter[4,1,], type = 'b', col = 4, lwd = 2)
lines(NW, mF1counter[5,1,], type = 'b', col = 5, lwd = 2)
lines(NW, mF1counter[6,1,], type = 'b', col = 6, lwd = 2)
title(paste("linear FM, carrier ", round(f1,2)))

plot(NW, mF1counter[1,2,], type = 'b', ylim = c(0,1), xaxt = 'n', yaxt = 'n', lwd = 2, las = 1)
grid()
lines(NW, mF1counter[2,2,], type = 'b', col = 2, lwd = 2)
lines(NW, mF1counter[3,2,], type = 'b', col = 3, lwd = 2)
lines(NW, mF1counter[4,2,], type = 'b', col = 4, lwd = 2)
lines(NW, mF1counter[5,2,], type = 'b', col = 5, lwd = 2)
lines(NW, mF1counter[6,2,], type = 'b', col = 6, lwd = 2)
title(paste("quadratic FM, carrier ", round(f2,2)))

plot(NW, mF1counter[1,3,], type = 'b', ylim = c(0,1), xaxt = 'n', yaxt = 'n', lwd = 2)
grid()
#axis(4, las = 1)
lines(NW, mF1counter[2,3,], type = 'b', col = 2, lwd = 2)
lines(NW, mF1counter[3,3,], type = 'b', col = 3, lwd = 2)
lines(NW, mF1counter[4,3,], type = 'b', col = 4, lwd = 2)
lines(NW, mF1counter[5,3,], type = 'b', col = 5, lwd = 2)
lines(NW, mF1counter[6,3,], type = 'b', col = 6, lwd = 2)
title(paste("cubic FM, carrier ", round(f3,2)))

plot(NW, mF1counter[1,4,], type = 'b', ylim = c(0,1), lwd = 2, las = 1)
grid()
lines(NW, mF1counter[2,4,], type = 'b', col = 2, lwd = 2)
lines(NW, mF1counter[3,4,], type = 'b', col = 3, lwd = 2)
lines(NW, mF1counter[4,4,], type = 'b', col = 4, lwd = 2)
lines(NW, mF1counter[5,4,], type = 'b', col = 5, lwd = 2)
lines(NW, mF1counter[6,4,], type = 'b', col = 6, lwd = 2)
legend(x = 4, y = 100, legend = 1:6, col = 1:6, lty = 1, horiz = TRUE)
title(paste("quartic FM, carrier ", round(f4,2)))

plot(NW, mF1counter[1,5,], type = 'b', ylim = c(0,1), yaxt = 'n', lwd = 2)
grid()
lines(NW, mF1counter[2,5,], type = 'b', col = 2, lwd = 2)
lines(NW, mF1counter[3,5,], type = 'b', col = 3, lwd = 2)
lines(NW, mF1counter[4,5,], type = 'b', col = 4, lwd = 2)
lines(NW, mF1counter[5,5,], type = 'b', col = 5, lwd = 2)
lines(NW, mF1counter[6,5,], type = 'b', col = 6, lwd = 2)
title(paste("sinusoidal FM, carrier ", round(f5,2)))

plot(NW, mF1counter[1,6,], type = 'b', ylim = c(0,1), yaxt = 'n', lwd = 2)
grid()
#axis(4, las = 1)
lines(NW, mF1counter[2,6,], type = 'b', col = 2, lwd = 2)
lines(NW, mF1counter[3,6,], type = 'b', col = 3, lwd = 2)
lines(NW, mF1counter[4,6,], type = 'b', col = 4, lwd = 2)
lines(NW, mF1counter[5,6,], type = 'b', col = 5, lwd = 2)
lines(NW, mF1counter[6,6,], type = 'b', col = 6, lwd = 2)
title(paste("no FM, carrier ", round(f6,2)))
#legend(x = 5, y = 1, legend = 1:6, lty = 1, col = 1:6, lwd = 2)

mtext("NW", side = 1, line = 2.2, outer = TRUE)
mtext("Detection probability", side = 2, line = 2.3, outer = TRUE)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("right", legend = 1:6, col = 1:6, lty = 1, lwd = 2, horiz = FALSE, xpd = TRUE, 
       bty = 'n', inset = c(0.01,0), cex = 1.3)
#################
#F3

par(mfrow = c(2,3))
par(oma = c(4,4,0,5))
par(mar = c(0,0,1,1))

plot(NW, F3counter[1,1,], type = 'b', ylim = c(0,1), xaxt = 'n', lwd = 2, las = 1)
grid()
lines(NW, F3counter[2,1,], type = 'b', col = 2, lwd = 2)
lines(NW, F3counter[3,1,], type = 'b', col = 3, lwd = 2)
lines(NW, F3counter[4,1,], type = 'b', col = 4, lwd = 2)
lines(NW, F3counter[5,1,], type = 'b', col = 5, lwd = 2)
lines(NW, F3counter[6,1,], type = 'b', col = 6, lwd = 2)
lines(NW, F3counter[7,1,], type = 'b', col = 7, lwd = 2)
title(paste("linear FM, carrier ", round(f1,2)))

plot(NW, F3counter[1,2,], type = 'b', ylim = c(0,1), xaxt = 'n', yaxt = 'n', lwd = 2)
grid()
lines(NW, F3counter[2,2,], type = 'b', col = 2, lwd = 2)
lines(NW, F3counter[3,2,], type = 'b', col = 3, lwd = 2)
lines(NW, F3counter[4,2,], type = 'b', col = 4, lwd = 2)
lines(NW, F3counter[5,2,], type = 'b', col = 5, lwd = 2)
lines(NW, F3counter[6,2,], type = 'b', col = 6, lwd = 2)
lines(NW, F3counter[7,2,], type = 'b', col = 7, lwd = 2)
title(paste("quadratic FM, carrier ", round(f2,2)))

plot(NW, F3counter[1,3,], type = 'b', ylim = c(0,1), xaxt = 'n', yaxt = 'n', lwd = 2)
grid()
#axis(4, las = 1)
lines(NW, F3counter[2,3,], type = 'b', col = 2, lwd = 2)
lines(NW, F3counter[3,3,], type = 'b', col = 3, lwd = 2)
lines(NW, F3counter[4,3,], type = 'b', col = 4, lwd = 2)
lines(NW, F3counter[5,3,], type = 'b', col = 5, lwd = 2)
lines(NW, F3counter[6,3,], type = 'b', col = 6, lwd = 2)
lines(NW, F3counter[7,3,], type = 'b', col = 7, lwd = 2)
title(paste("cubic FM, carrier ", round(f3,2)))

plot(NW, F3counter[1,4,], type = 'b', ylim = c(0,1), lwd = 2)
grid()
lines(NW, F3counter[2,4,], type = 'b', col = 2, lwd = 2)
lines(NW, F3counter[3,4,], type = 'b', col = 3, lwd = 2)
lines(NW, F3counter[4,4,], type = 'b', col = 4, lwd = 2)
lines(NW, F3counter[5,4,], type = 'b', col = 5, lwd = 2)
lines(NW, F3counter[6,4,], type = 'b', col = 6, lwd = 2)
lines(NW, F3counter[7,4,], type = 'b', col = 7, lwd = 2)
legend(x = 6, y = 60, legend = 1:6, col = 1:6, lty = 1, horiz = TRUE)
title(paste("quartic FM, carrier ", round(f4,2)))

plot(NW, F3counter[1,5,], type = 'b', ylim = c(0,1), yaxt = 'n', lwd = 2)
grid()
lines(NW, F3counter[2,5,], type = 'b', col = 2, lwd = 2)
lines(NW, F3counter[3,5,], type = 'b', col = 3, lwd = 2)
lines(NW, F3counter[4,5,], type = 'b', col = 4, lwd = 2)
lines(NW, F3counter[5,5,], type = 'b', col = 5, lwd = 2)
lines(NW, F3counter[6,5,], type = 'b', col = 6, lwd = 2)
lines(NW, F3counter[7,5,], type = 'b', col = 7, lwd = 2)
title(paste("sinusoidal FM, carrier ", round(f5,2)))

plot(NW, F3counter[1,6,], type = 'b', ylim = c(0,1), yaxt = 'n', lwd = 2)
grid()
#axis(4, las = 1)
lines(NW, F3counter[2,6,], type = 'b', col = 2, lwd = 2)
lines(NW, F3counter[3,6,], type = 'b', col = 3, lwd = 2)
lines(NW, F3counter[4,6,], type = 'b', col = 4, lwd = 2)
lines(NW, F3counter[5,6,], type = 'b', col = 5, lwd = 2)
lines(NW, F3counter[6,6,], type = 'b', col = 6, lwd = 2)
lines(NW, F3counter[7,6,], type = 'b', col = 7, lwd = 2)
title(paste("no FM, carrier ", round(f6,2)))
#legend(x = 9, y = 1, legend = 0:6, lty = 1, col = 1:7, lwd = 2)

mtext("NW", side = 1, line = 2.2, outer = TRUE)
mtext("Detection probability", side = 2, line = 2.3, outer = TRUE)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("right", legend = 0:6, col = 1:7, lty = 1, lwd = 2, horiz = FALSE, xpd = TRUE, 
       bty = 'n', inset = c(0.01,0), cex = 1.3)

#################
#mF3

par(mfrow = c(2,3))
par(oma = c(4,4,0,5))
par(mar = c(0,0,1,1))

plot(NW, mF3counter[1,1,], type = 'b', ylim = c(0,1), xaxt = 'n', lwd = 2, las = 1)
grid()
lines(NW, mF3counter[2,1,], type = 'b', col = 2, lwd = 2)
lines(NW, mF3counter[3,1,], type = 'b', col = 3, lwd = 2)
lines(NW, mF3counter[4,1,], type = 'b', col = 4, lwd = 2)
lines(NW, mF3counter[5,1,], type = 'b', col = 5, lwd = 2)
lines(NW, mF3counter[6,1,], type = 'b', col = 6, lwd = 2)
title(paste("linear FM, carrier ", round(f1,2)))

plot(NW, mF3counter[1,2,], type = 'b', ylim = c(0,1), xaxt = 'n', yaxt = 'n', lwd = 2)
grid()
lines(NW, mF3counter[2,2,], type = 'b', col = 2, lwd = 2)
lines(NW, mF3counter[3,2,], type = 'b', col = 3, lwd = 2)
lines(NW, mF3counter[4,2,], type = 'b', col = 4, lwd = 2)
lines(NW, mF3counter[5,2,], type = 'b', col = 5, lwd = 2)
lines(NW, mF3counter[6,2,], type = 'b', col = 6, lwd = 2)
title(paste("quadratic FM, carrier ", round(f2,2)))

plot(NW, mF3counter[1,3,], type = 'b', ylim = c(0,1), xaxt = 'n', yaxt = 'n', lwd = 2)
grid()
#axis(4, las = 1)
lines(NW, mF3counter[2,3,], type = 'b', col = 2, lwd = 2)
lines(NW, mF3counter[3,3,], type = 'b', col = 3, lwd = 2)
lines(NW, mF3counter[4,3,], type = 'b', col = 4, lwd = 2)
lines(NW, mF3counter[5,3,], type = 'b', col = 5, lwd = 2)
lines(NW, mF3counter[6,3,], type = 'b', col = 6, lwd = 2)
title(paste("cubic FM, carrier ", round(f3,2)))

plot(NW, mF3counter[1,4,], type = 'b', ylim = c(0,1), lwd = 2)
grid()
lines(NW, mF3counter[2,4,], type = 'b', col = 2, lwd = 2)
lines(NW, mF3counter[3,4,], type = 'b', col = 3, lwd = 2)
lines(NW, mF3counter[4,4,], type = 'b', col = 4, lwd = 2)
lines(NW, mF3counter[5,4,], type = 'b', col = 5, lwd = 2)
lines(NW, mF3counter[6,4,], type = 'b', col = 6, lwd = 2)
legend(x = 6, y = 60, legend = 1:6, col = 1:6, lty = 1, horiz = TRUE)
title(paste("quartic FM, carrier ", round(f4,2)))

plot(NW, mF3counter[1,5,], type = 'b', ylim = c(0,1), yaxt = 'n', lwd = 2)
grid()
lines(NW, mF3counter[2,5,], type = 'b', col = 2, lwd = 2)
lines(NW, mF3counter[3,5,], type = 'b', col = 3, lwd = 2)
lines(NW, mF3counter[4,5,], type = 'b', col = 4, lwd = 2)
lines(NW, mF3counter[5,5,], type = 'b', col = 5, lwd = 2)
lines(NW, mF3counter[6,5,], type = 'b', col = 6, lwd = 2)
title(paste("sinusoidal FM, carrier ", round(f5,2)))

plot(NW, mF3counter[1,6,], type = 'b', ylim = c(0,1), yaxt = 'n', lwd = 2)
grid()
#axis(4, las = 1)
lines(NW, mF3counter[2,6,], type = 'b', col = 2, lwd = 2)
lines(NW, mF3counter[3,6,], type = 'b', col = 3, lwd = 2)
lines(NW, mF3counter[4,6,], type = 'b', col = 4, lwd = 2)
lines(NW, mF3counter[5,6,], type = 'b', col = 5, lwd = 2)
lines(NW, mF3counter[6,6,], type = 'b', col = 6, lwd = 2)
title(paste("no FM, carrier ", round(f6,2)))
#legend(x = 5, y = 1, legend = 1:6, lty = 1, col = 1:6, lwd = 2)

mtext("NW", side = 1, line = 2.2, outer = TRUE)
mtext("Detection probability", side = 2, line = 2.3, outer = TRUE)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("right", legend = 1:6, col = 1:6, lty = 1, lwd = 2, horiz = FALSE, xpd = TRUE, 
       bty = 'n', inset = c(0.01,0), cex = 1.3)

#####################################################
# estimates of frequency

f <- c(f1,f2,f3,f4,f5,f6)

frMSEf3 <- array(0, dim = c(M+1,numNW, numBand))
for(l in 1:numBand){
  for(k in 1:numNW){
    frMSEf3[,k,l] <- apply(sigFreqsF3c[,,k,l], MARGIN = 1, FUN = function(x) mean( (x-f[l])^2, na.rm = TRUE))
  }
}

frMSEmf3 <- array(0, dim = c(M,numNW, numBand))
for(l in 1:numBand){
  for(k in 1:numNW){
    frMSEmf3[,k,l] <- apply(sigFreqsmF3c[,,k,l], MARGIN = 1, FUN = function(x) mean( (x-f[l])^2, na.rm = TRUE))
  }
}

frMSElist <- c( frMSEf3[1,,6], frMSEf3[2,,1], frMSEmf3[1,,1], frMSEf3[3,,2],
                frMSEmf3[2,,2], frMSEf3[4,,3], frMSEmf3[3,,3], frMSEf3[5,,4],
                frMSEmf3[4,,4], frMSEf3[7,,5], frMSEmf3[6,,5])
rnge <- range(frMSElist, na.rm = TRUE)

par(mfrow = c(2,3))
par(oma = c(4,4,0,5))
par(mar = c(0,0,1,1))

plot(NW, frMSEf3[1,,6], type = 'b', lwd = 2, log = 'y', xaxt = 'n', las = 1,
     ylim = rnge)
grid()
title(paste("no FM, carrier ", round(f6,2)))

plot(NW, frMSEf3[2,,1], type = 'b', lwd = 2,
     ylim = rnge, log = 'y', xaxt = 'n', yaxt= 'n')
grid()
lines(NW, frMSEmf3[1,,1], type = 'b', lwd = 2, col = 2)
title(paste("linear FM, carrier ", round(f1,2)))

plot(NW, frMSEf3[3,,2], type = 'b', lwd = 2, 
     ylim =rnge, log = 'y', xaxt = 'n', yaxt = 'n')
grid()
#axis(4, las = 1)
lines(NW, frMSEmf3[2,,2], type = 'b', lwd = 2, col = 2)
title(paste("quadratic FM, carrier ", round(f2,2)))

plot(NW, frMSEf3[4,,3], type = 'b', lwd = 2, 
     ylim = rnge, log = 'y', las = 1)
grid()
lines(NW, frMSEmf3[3,,3], type = 'b', lwd = 2, col = 2)
title(paste("cubic FM, carrier ", round(f3,2)))

plot(NW, frMSEf3[5,,4], type = 'b', lwd = 2, 
     ylim = rnge, log = 'y', yaxt = 'n')
grid()
lines(NW, frMSEmf3[4,,4], type = 'b', lwd = 2, col = 2)
title(paste("quartic FM, carrier ", round(f4,2)))

plot(NW, frMSEf3[7,,5], type = 'b', lwd = 2, 
     ylim = rnge, log = 'y', yaxt = 'n')
grid()
#axis(4, las = 1)
lines(NW, frMSEmf3[6,,5], type = 'b', lwd = 2, col = 2)
title(paste("sinusoidal FM, carrier ", round(f5,2)))

mtext("NW", side = 1, outer = TRUE, line = 2.2)
mtext("MSE (log scale)", side = 4, outer = TRUE, line = 2.2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("topright", legend = expression(F[3], tilde(F)[3]), col = 1:2, lty = 1, lwd = 2, horiz = FALSE, xpd = TRUE, 
       bty = 'n', inset = c(0.01,0), cex = 1.3)

#####################################################
## histograms?

#F3
par(mfrow = c(1,1))
par(oma = c(4,4,0,4))
par(mar = c(0,0,1,1))

hist(sigFreqsF3[1,,1], freq = FALSE, breaks = 1000, main = '', xaxt = 'n')
text("Degree 0", x = 0.1, y = 400, cex = 1)
axis(1, at = c(0,0.1,0.2,0.3,0.4,0.5), labels = c(0,0.1,0.2,0.3,0.4,0.5))
mtext("Frequency", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

par(mfrow = c(2,2))
hist(sigFreqsF3[2,,3], freq = FALSE, breaks = 1000, main = '')
text("Degree 1", x = 0.3, y = 500, cex = 1)
hist(sigFreqsF3[3,,3], freq = FALSE, breaks = 1000, main = '',  yaxt = 'n')
text("Degree 2", x = 0.1, y = 1000, cex = 1)
axis(4)
hist(sigFreqsF3[4,,3], freq = FALSE, breaks = 1000, main = '')
text("Degree 3", x = 0.1, y = 200, cex = 1)
hist(sigFreqsF3[5,,3], freq = FALSE, breaks = 1000, main = '', yaxt = 'n')
text("Degree 4", x = 0.1, y = 300, cex = 1)
axis(4)
mtext("Frequency", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

par(mfrow = c(1,1))
hist(sigFreqsF3[7,,5], freq = FALSE, breaks = 1000, main = '')
text("Sinusoidal, degree 6", x = 0.35, y = 400, cex = 1)
mtext("Frequency", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

### coloured noise: what happens to the degree 2 F3 when NW >= 8?
#histogram doesn't tell much
par(mfrow = c(1,1))
hist(sigFreqsF3[3,,4], freq = FALSE, breaks = 2000, main = '', xlim = c(0.29,0.31))
abline(v = c(f2 - W[4], f2 + W[4]), col = 3, lty = 3, lwd = 2)
abline(v = c(f2 - B3, f2 + B3), col = 2, lty = 2, lwd = 2)
mtext("Frequency", side = 1, outer = FALSE, line = 2.2)
mtext("Density", side = 2, outer = FALSE, line = 2.5)

###
#mF3

par(mfrow = c(2,2))
hist(sigFreqsmF3[1,,3], freq = FALSE, breaks = 1000, main = '')
text("Degree 1", x = 0.3, y = 600, cex = 1)
hist(sigFreqsmF3[2,,3], freq = FALSE, breaks = 1000, main = '', yaxt = 'n')
text("Degree 2", x = 0.1, y = 1000, cex = 1)
axis(4)
hist(sigFreqsmF3[3,,3], freq = FALSE, breaks = 1000, main = '')
text("Degree 3", x = 0.1, y = 600, cex = 1)
hist(sigFreqsmF3[4,,3], freq = FALSE, breaks = 1000, main = '', yaxt = 'n')
text("Degree 4", x = 0.1, y = 1000, cex = 1)
axis(4)
mtext("Frequency", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

par(mfrow = c(1,1))
hist(sigFreqsmF3[6,,5], freq = FALSE, breaks = 1000, main = '')
text("Sinusoidal, degree 6", x = 0.35, y = 400, cex = 1)
mtext("Frequency", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)


#######################
# plots of test statistics around carrier frequencies

ind <- 8008
F3d0 <- F3array[1,,ind,1]
F3 <- F3array[,,ind,3]
mF3 <- mF3array[,,ind,3]
F3sin <- F3array[7,,ind,5]
mF3sin <- mF3array[6,,ind,5]

F3d0[which(F3d0 < 1)] <- 1
F3sin[which(F3sin < 1)] <- 1
mF3sin[which(mF3sin < 1)] <- 1
for(m in 1:(M+1)){
  F3[m,which(F3[m,] < 1)] <- 1
  if(m < M+1){
    mF3[m, which(mF3[m,] < 1)] <- 1
  }
}

par(mfrow = c(1,1))
plot(freqF, F3d0, type = 'l', log = 'y', xlim = c(0.23,0.27), lwd = 2)
plot(freqF, F3[2,], type = 'l', log = 'y', xlim = c(0.07,0.13), lwd = 2)
lines(freqF, mF3[1,], col = 2, lwd = 2)

plot(freqF, F3[3,], type = 'l', log = 'y', xlim = c(0.27,0.33), lwd = 2)
lines(freqF, mF3[2,], col = 2, lwd = 2)

plot(freqF, F3[4,], type = 'l', log = 'y', xlim = c(0.28,0.34), lwd = 2)
lines(freqF, mF3[3,], col = 2, lwd = 2)

plot(freqF, F3[5,], type = 'l', log = 'y', xlim = c(0.37,0.43), lwd = 2)
lines(freqF, mF3[4,], col = 2, lwd = 2)

plot(freqF, F3sin, type = 'l', log = 'y', xlim = c(0.17,0.23), lwd = 2)
lines(freqF, mF3sin, col = 2, lwd = 2)

#what about F3(2;f) for NW > 7?
F3d2 <- F3array[3,,111,4]
F3d2[which(F3d2 < 1)] <- 1
plot(freqF, F3d2, type = 'l', log = 'y', lwd = 2, xlim = c(0.27,0.33))
abline(v = c(f2 - W[4], f2 + W[4]), col = 3, lty = 3, lwd = 2)
abline(v = c(f2 - B3, f2 + B3), col = 2, lty = 2, lwd = 2)
#########################
#out of band-detections?

#oobDetectF1.tmp <- oobDetectF3.tmp <- array(0, dim = c(M, numNW, num.sim))
#ooBand <- length(unlist(band))
#for(m in 1:(M)){
#  for(k in 1:numNW){
#    for(i in 1:num.sim){
#      #oobDetectF1[m,k] <- length( which(F1array[m, -unlist(band), ,k] > qf(conf, m, K[k]-m))) / ( num.sim*(nfreqs-ooBand) )
#      #oobDetectF3[m,k] <- length( which(F3array[m, -unlist(band), ,k] > qf(conf, 1, K[k]-m))) / ( num.sim*(nfreqs-ooBand) )
#      F1bandik <- F1array[,-unlist(band),i,k]
#      F3bandik <- F3array[,-unlist(band),i,k]
#      
#      oobDetectF1.tmp[m,k,i] <- length( findLocalFMaxM2(F1bandik, K[k], conf, m, type = 1) ) / (num.sim*(nfreqs-ooBand))
#      oobDetectF3.tmp[m,k,i] <- length( findLocalFMaxM2(F3bandik, K[k], conf, m, type = 3) ) / (num.sim*(nfreqs-ooBand))
#    }
#  }
#}
#oobDetectF1 <- oobDetectF3 <- array(0, dim = c(M, numNW))
#for(m in 1:M){
#  oobDetectF1[m,] <- rowSums(oobDetectF1.tmp[m,,])
#  oobDetectF3[m,] <- rowSums(oobDetectF3.tmp[m,,])
#}
#
#par(mfrow = c(2,3))
#plot(NW, oobDetectF1[1,], type = 'b', ylim = c(0,0.005), xaxt = 'n')
#grid()
#lines(NW, oobDetectF3[1,], type = 'b', col = 2)

#plot(NW, oobDetectF1[2,], type = 'b', ylim = c(0,0.005), xaxt = 'n', yaxt = 'n')
#grid()
#lines(NW, oobDetectF3[2,], type = 'b', col = 2)
#legend(x=4,y=0.1, legend = c("F1","F2"), lty = 1, col = 1:2)

#plot(NW, oobDetectF1[3,], type = 'b', ylim = c(0,0.005), xaxt = 'n', yaxt = 'n')
#grid()
#lines(NW, oobDetectF3[3,], type = 'b', col = 2)
#axis(4)

#plot(NW, oobDetectF1[4,], type = 'b', ylim = c(0,0.005))
#grid()
#lines(NW, oobDetectF3[4,], type = 'b', col = 2)

#plot(NW, oobDetectF1[5,], type = 'b', ylim = c(0,0.005), yaxt = 'n')
#grid()
#lines(NW, oobDetectF3[5,], type = 'b', col = 2)

#plot(NW, oobDetectF1[6,], type = 'b', ylim = c(0,0.005), yaxt = 'n')
#grid()
#axis(4)
#lines(NW, oobDetectF3[6,], type = 'b', col = 2)