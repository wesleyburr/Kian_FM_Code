source("utility.R")

N <- 1024
n <- 0:(N-1)
nFFT <- 2^ceiling(log2(2*N))


js1 <- 420
js2 <- 5318008


#amplitudes
A <- 1

#carrier frequencies
f <-  0.25

#time and band parameters
B <- 0.005
fbw <- 1
stm11 <- 2.0/(N-1)
Shift <- 0.0
dt <- 1
tt <- (n + Shift) * stm11 - 1.0  # this runs from -1 to 1


#modulating polynomials
FMpoly <- B * fbw * (1.0 - 2.0 * tt^2)     #quadratic
Pref <- cumsum(FMpoly)*2*pi*dt
PhMod <- 2*pi*f*dt*n + Pref


#carrier waves
sig <- A*cos(PhMod) 



#multitaper and background spectrum stuff
snr <- 10
NW <- 5.5
K <- (2*NW-3):(2*NW)
numK <- length(K)
W <- NW / N

sigmat2 <- 0.5*A^2/(snr*2*W)  #really snr is (A^2/2) / (2W*S(f))
sigmat <- sqrt(sigmat2)

#ar <- c(0.5, -0.3)
#ma <- -0.1

###

M <- 4
#L <- 10*K
#DS <- dpss(n = L, nw = NW, k = K)


nfreqs <- nFFT/2
num.sim <- 10000


HFarray <- array(data = 0, dim = c(num.sim, nfreqs, numK))
F1array <- F2array <- F3array  <- array(data = 0, dim = c(M+1,nfreqs-1,num.sim,numK))

seed <- 666
set.seed(seed)

freq <- 1:nfreqs / (nFFT*dt)
band.ind <- which(freq < f + W & freq > f - W)
F1counter <- F2counter <- F3counter <- matrix(0, nrow = M+1, ncol = numK)
Fcounter <- rep(0, numK)

conf <- 1-1/N

start.time <- Sys.time()
for(j in 1:numK){
  DW <- dpss(n = N, nw = NW, k = K[j])
  V <- DW$v
  Vdot <- dpssp11R(DW,NW)
  AP <- dpssap(V, maxdeg = M, alpha = 1)
  scls <- inverseScales(V,Vdot)
  for(i in 1:num.sim){
    
    #noiz <- arima.sim(model = list(ar = ar, ma = ma), n = N, rand.gen = rnorm, mean = 0, sd = sigmat)
    noiz <- rnorm(N, sd = sigmat)
    X <- sig + noiz
    sp.X <- spec.mtm(X, k = K[j], nw = NW, deltat = dt, Ftest = TRUE, returnInternals = TRUE, 
                     dpssIN = DW, plot = FALSE, nFFT = nFFT, returnZeroFreq = FALSE)
    
    modF.test <- ModulatedF12(sp.X, mxdeg = M, adaptive = FALSE, derivIN = Vdot, apIN = AP, 
                             dpssIN = DW, subBand = 1:(nfreqs-1), scale = TRUE, scalesIN = scls)
    
    HFarray[i,,j] <- sp.X$mtm$Ftest
    F1array[,,i,j] <- modF.test$F1
    F2array[,,i,j] <- modF.test$F2
    F3array[,,i,j] <- modF.test$F3
    
    
    
    for(m in 1:(M+1)){
      if(max(F1array[m,,i,j]) > qf(conf,m,K[j]-m)){
        if(sum(which.max(F1array[m,,i,j]) == band.ind) == 1){
          F1counter[m,j] <- F1counter[m,j] + 1
        }
      }
      if(max(F2array[m,,i,j]) > qf(conf,m,K[j]-m)){
        if(sum(which.max(F2array[m,,i,j]) == band.ind) == 1){
          F2counter[m,j] <- F2counter[m,j] + 1
        }
      }
      if(max(F3array[m,,i,j]) > qf(conf,1,K[j]-m)){
        if(sum(which.max(F3array[m,,i,j]) == band.ind) == 1){
          F3counter[m,j] <- F3counter[m,j] + 1
        }
      }
      
    }
    if(max(HFarray[i,,j]) > qf(conf,2,2*K[j]-2)){
      if(sum(which.max(HFarray[i,,j]) == band.ind) == 1){
        Fcounter[j] <- Fcounter[j] + 1
      }
    }
    
  }
  
}
runtime <- Sys.time() - start.time
runtime

seed2 <- 420
set.seed(seed2)
l <- sample(1:(nfreqs-1), size = 1)
png("KdistsD0.png", width = 1600, height = 900)
  par(mfrow = c(3,2))
  par(oma = c(0.5, 0, 0.5, 0.5)) # make room (i.e. the 4's) for the overall x and y axis titles
  par(mar = c(0, 0, 1, 1)) 
  
  f1Mesh <- seq(0, 500, length.out = num.sim)
  f3Mesh <- seq(0, 500, length.out = num.sim)
  
  for(m in 1:(M+1)){
    for(k in 3:numK){
      hist(F1array[m,l,,k], freq = FALSE, breaks = 200, xlab = '', ylab = 'F1', main = '',
           xaxt = 'n', yaxt='n', xlim = c(0,35*m/(m+1)*(K[k]+1)/K[k]))
      lines(density(F1array[m,l,,k]), col = 2, lwd = 2)
      lines(f1Mesh, df(f1Mesh, m, K[k]-m), col = 4, lwd = 2)
      legend(legend = paste("F1, degree ",m-1, "K=",K[k]), x = 3, y = 0.2, bty = 'n')
    } 
    for(k in 3:numK){    
      hist(F2array[m,l,,k], freq = FALSE, breaks = 200, xlab = '', ylab = 'F2', main = '',
           xaxt = 'n', yaxt = 'n', xlim = c(0,35*m/(m+1)*(K[k]+1)/K[k]))
      lines(density(F2array[m,l,,k]), col = 2, lwd = 2)
      lines(f1Mesh, df(f1Mesh, m, K[k]-m), col = 4, lwd = 2)
      legend(legend = paste("F2, degree ",m-1, "K=",K[k]), x = 3, y = 0.2, bty = 'n')
    }
    for(k in 3:numK){
      hist(F3array[m,l,,k], freq = FALSE, breaks = 200, xlab = '', ylab = 'F3', main = '',
           xaxt = 'n', yaxt = 'n', xlim = c(0,35*m/(m+1)*(K[k]+1)/K[k]))
      lines(density(F3array[m,l,,k]), col = 2, lwd = 2)
      lines(f3Mesh, df(f3Mesh, 1, K[k]-m), col = 4, lwd = 2)
      legend(legend = paste("F3, degree ",m-1, "K=",K[k]), x = 3, y = 0.3, bty = 'n')
    }
  }
    
  
dev.off()

#rejection probabilities
rejProbsF1 <- rejProbsF3 <- array(0, dim = c(M+1,nfreqs-1,numK))
for(m in 1:(M+1)){
  for(k in 1:numK){
    for(j in 1:(nfreqs-1)){
      rejProbsF1[m,j,k] <- length( which(F1array[m,j,,k] > qf(1-1/N, m, K[k]-m)) )/num.sim
      rejProbsF3[m,j,k] <- length( which(F3array[m,j,,k] > qf(1-1/N, 1, K[k]-m)) )/num.sim
    }
  }
}