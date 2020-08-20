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
B <- 0.002
fbw <- 1
stm11 <- 2.0/(N-1)
Shift <- 0.0
dt <- 1
tt <- (n + Shift) * stm11 - 1.0  # this runs from -1 to 1


#modulating polynomials
#FMpoly <- B * fbw * (1.0 - 2.0 * tt^2)     #quadratic
FMpoly <- B * fbw * tt
Pref <- cumsum(FMpoly)*2*pi*dt
PhMod <- 2*pi*f*dt*n + Pref


#carrier waves
sig <- A*cos(PhMod) 



#multitaper and background spectrum stuff
snr <- 10
NW <- 4
K <- 2*NW - 1
W <- NW / N
									
sigmat2 <- 0.5*A^2/(snr*2*W)  #really snr is (A^2/2) / (2W*S(f))
sigmat <- sqrt(sigmat2)

#ar <- c(0.5, -0.3)
#ma <- -0.1

###
DW <- dpss(n = N, nw = NW, k = K)
#V <- DW$v
#Vdot <- dpssp11R(DW,NW)

#or downsample
L <- 200
DS <- dpss(n = L, nw = NW, k = K)
V <- DS$v
Vdot <- dpssp11R(DS,NW)


M <- 5
AP <- dpssapD(V, maxdeg = M, alpha = 1, deriv = FALSE)
apS <- assPolySmoother(AP,V,Vdot,smooth = 2)
scls <- inverseScales(V, Vdot)

nfreqs <- nFFT/2
num.sim <- 100

HFarray <- array(data = 0, dim = c(num.sim, nfreqs))
F1array <- F2array <- F3array  <- array(data = 0, dim = c(M+1,nfreqs-1,num.sim))

seed <- 111 #sample(js1:js2, size = 1)
set.seed(seed)

freq <- 1:nfreqs / (nFFT*dt)
band.ind <- which(freq < f + W & freq > f - W)
F1counter <- F2counter <- F3counter <- rep(0,M+1)
Fcounter <- 0

conf <- 1-1/N

start.time <- Sys.time()
for(i in 1:num.sim){

  #noiz <- arima.sim(model = list(ar = ar, ma = ma), n = N, rand.gen = rnorm, mean = 0, sd = sigmat)
  noiz <- rnorm(N, sd = sigmat)
  X <- sig + noiz
  sp.X <- spec.mtm(X, k = K, nw = NW, deltat = dt, Ftest = TRUE, returnInternals = TRUE, 
                   dpssIN = DW, plot = FALSE, nFFT = nFFT, returnZeroFreq = FALSE)
  
  modF.test <- ModulatedF13(sp.X, mxdeg = M, adaptive = FALSE, derivIN = Vdot, apIN = AP, 
                           dpssIN = DS, subBand = 1:(nfreqs-1), smooth = -1, scale =  TRUE,
                           scalesIN = scls)
  
  HFarray[i,] <- sp.X$mtm$Ftest
  F1array[,,i] <- modF.test$F1
  F2array[,,i] <- modF.test$F2
  F3array[,,i] <- modF.test$F3
  
  
  for(m in 1:(M+1)){
    if(max(F1array[m,,i]) > qf(conf,m,K-m)){
      if(sum(which.max(F1array[m,,i]) == band.ind) == 1){
        F1counter[m] <- F1counter[m] + 1
      }
    }
    if(max(F2array[m,,i]) > qf(conf,m,K-m)){
      if(sum(which.max(F2array[m,,i]) == band.ind) == 1){
        F2counter[m] <- F2counter[m] + 1
      }
    }
    if(max(F3array[m,,i]) > qf(conf,1,K-m)){
      if(sum(which.max(F3array[m,,i]) == band.ind) == 1){
        F3counter[m] <- F3counter[m] + 1
      }
    }
   
  }
  if(max(HFarray[i,]) > qf(conf,2,2*K-2)){
    if(sum(which.max(HFarray[i,]) == band.ind) == 1){
      Fcounter <- Fcounter + 1
    }
  }

}

runtime <- Sys.time() - start.time
runtime

Fcounter
F1counter
F2counter
F3counter




aveF1 <- aveF2 <- aveF3 <- matrix(NA, nrow = M+1, ncol = nfreqs-1)
for(m in 1:(M+1)){
  aveF1[m,] <- apply(F1array[m,,], MARGIN = 1, FUN = mean)
  aveF2[m,] <- apply(F2array[m,,], MARGIN = 1, FUN = mean)
  aveF3[m,] <- apply(F3array[m,,], MARGIN = 1, FUN = mean)
}


qqplot(x = qf(seq(0,1, length.out = num.sim),4,K-4), y = F1array[4,69,])
abline(a=0,b=1)

Farrays <- list(F1 = F1array, F2 = F2array, F3 = F3array, Hf = HFarray)
counters <- list(F1counter=F1counter, F2counter=F2counter, F3counter=F3counter, Fcounter=Fcounter)
params <- list(A=A,B=B,dt=dt,f=f,fbw=fbw,FMpoly=FMpoly,K=K,M=M,N=N,nFFT=nFFT,NW=NW,seed=seed,
               sigmat=sigmat,snr=snr,W=W)
Fstatstuff <- list(Farrays,counters,params)
save(Fstatstuff, file = "Fstatstuff.RData")
