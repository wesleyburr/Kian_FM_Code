source("thesisUtility.R")


N <- 10000
nFFT <- 2^ceiling(log2(2*N))

#multitaper and background spectrum stuff
NW <- 4:11
numNW <- length(NW)
K <- 2*NW - 1
W <- NW / N

ar <- c(0.6, -0.3)
ma <- c(1, 0.8, 0.6)
trnd <- 10*seq(0,1, length.out = N)^2


L <- c(0.01, 0.05, 0.1) * N
numL <- length(L)
M <- 5

nfreqs <- nFFT/2
num.sim <- 10000

seed <- 666
set.seed(seed)

freq <- 1:nfreqs / nFFT
subBand <- sample(1:(nfreqs-1), size = 1)
subBand <- 12687

sigs <- c(0.9, 0.95, 0.99, 0.999)
mF1array <- mF3array <- array(data = 0, dim = c(M,num.sim,numNW,numL))
mF1counter <- mF3counter <- array(0, dim = c(M, length(sigs), numNW, numL))

F3array <- array(data = 0, dim = c(M+1,num.sim,numNW,numL))
F3counter <- array(0, dim = c(M+1, length(sigs), numNW, numL))



start.time <- Sys.time()
for(n in 1:numNW){
  DW <- dpss(n = N, nw = NW[n], k = K[n])
  for(l in 1:numL){
    DS <- dpss(n = L[l], nw = NW[n], k = K[n])
    Vdot <- dpssp11R(DS,NW[n])
    AP <- dpssapD(DS$v, maxdeg = M, alpha = 1, deriv = FALSE)
    for(i in 1:num.sim){
      X1 <- rgumbel(n = N/4, 0, 0.1)
      X2 <- rgumbel(n = N/4, 0, 0.2)
      X3 <- rgumbel(n = N/4, 0, 0.5)
      X4 <- rgumbel(n = N/4, 0, 0.8)
      X.tmp <- c(X1,X2,X3,X4)
      X <- trnd + as.numeric(arima.sim(model = list(ar=ar,ma=ma), n = N, innov = X.tmp))
      
      tapered <- X * DW$v
      pad <- rbind(tapered, matrix(0, nrow = nFFT-N, ncol = K[n]))
      yk <- mvfft(pad)[subBand+1,,drop = FALSE]
      
      modF.test <- ModulatedF17(yk, mxdeg = M, derivIN = Vdot, apIN = AP, 
                                dpssIN = DS)
      
      
      F3array[,i,n,l] <- modF.test$F3
      mF1array[,i,n,l] <- modF.test$mF1
      mF3array[,i,n,l] <- modF.test$mF3
      
      for(m in 1:(M+1)){
        for(a in 1:length(sigs)){
          
          if(F3array[m,i,n,l] > qf(sigs[a],1,K[n]-m)){
            F3counter[m,a,n,l] <- F3counter[m,a,n,l]+1
          }
          if(m < M+1){
            if(mF1array[m,i,n,l] > qf(sigs[a], m, K[n]-m)){
              mF1counter[m,a,n,l] <- mF1counter[m,a,n,l] + 1
            }
            if(mF3array[m,i,n,l] > qf(sigs[a], 1, K[n]-m)){
              mF3counter[m,a,n,l] <- mF3counter[m,a,n,l] + 1
            }
          }
        }
      }
    }
  }
}
runtime <- Sys.time() - start.time
runtime

testStatsLNS <- list(F3 = F3array, mF1 = mF1array, mF3 = mF3array)
save(testStatsLNS, file = "testStatsLNS.RData")


rejProbsmF1 <- mF1counter/num.sim
rejProbsF3 <- F3counter/num.sim
rejProbsmF3 <- mF3counter/num.sim

rejProbsF3 <- array(0, dim = c(M+1, numNW, numL, length(sigs)))
rejProbsmF1 <- rejProbsmF3 <- array(0, dim = c(M, numNW, numL, length(sigs)))

K <- 2*NW - 1
for(m in 1:(M+1)){
  for(n in 1:numNW){
    for(k in 1:length(sigs)){
      for(l in 1:numL){
        rejProbsF3[m,n,l,k] <- length( which(F3array[m,,n,l] > qf(sigs[k], 1, K[n]-m)) )/num.sim
        if(m < M+1){
          rejProbsmF1[m,n,l,k] <- length( which(mF1array[m,,n,l] > qf(sigs[k], m, K[n]-m)) )/num.sim
          rejProbsmF3[m,n,l,k] <- length( which(mF3array[m,,n,l] > qf(sigs[k], 1, K[n]-m)) )/num.sim
        }
      }
    }
  }
}


############
# change these slightly because you have an extra dimension now
par(oma = c(4,4,1,4))
par(mar = c(0,0,1,1))
alpha <- 1-sigs
par(mfrow = c(2,3))


for(l in 1:numL){
  for(m in 1:(M+1)){
    plot(alpha, rejProbsF3[m,1,l,], type = 'b', log = 'xy', xlab = '', ylab = '', 
         ylim = range(rejProbsF3), xaxt = 'n', yaxt = 'n')
    if(m == 1 | m == 4){axis(2, at = alpha, label = alpha)}
    if(m == 3 | m == 6){axis(4, at = alpha, label = alpha)}
    if(m > 3){axis(1, at = alpha, label = alpha)}
    grid()
    text(paste(m-1), x = 0.09, y = 0.001, cex = 3)
    abline(a = 0, b = 1, col = 1, lty = 2, lwd = 2)
    for(n in 2:numNW){
      lines(alpha, rejProbsF3[m,n,l,], type = 'b', col = n)
    }
  }
  mtext(paste("L = ", L[l]), side = 3, outer = TRUE, line = -1)
  mtext("Significance level (log scale)", side = 1, outer = TRUE, line = 2.2)
  mtext("Rejection probability (log scale)", side = 4, outer = TRUE, line = 2.2)
}

for(l in 1:numL){
  par(mfrow = c(2,3))
  frame()
  for(m in 1:M){
    plot(alpha, rejProbsmF1[m,1,l,], type = 'b', log = 'xy', xlab = '', ylab = '', 
         ylim = range(rejProbsmF1), xaxt = 'n', yaxt = 'n', las = 1)
    if(m == 1 | m == 3){axis(2, at = alpha, label = alpha)}
    if(m == 2 | m == 5){axis(4, at = alpha, label = alpha)}
    if(m >= 3){axis(1, at = alpha, label = alpha)}
    grid()
    text(paste(m), x = 0.065, y = 0.0006, cex = 3)
    abline(a = 0, b = 1, col = 1, lty = 2, lwd = 2)
    for(n in 2:numNW){
      lines(alpha, rejProbsmF1[m,n,l,], type = 'b', col = n)
    }
  }
  mtext(paste("L = ", L[l]), side = 3, outer = TRUE, line = -1)
  mtext("Significance level (log scale)", side = 1, outer = TRUE, line = 2.2)
  mtext("Rejection probability (log scale)", side = 4, outer = TRUE, line = 2.2)
}

for(l in 1:numL){
  par(mfrow = c(2,3))
  frame()
  for(m in 1:M){
    plot(alpha, rejProbsmF3[m,1,l,], type = 'b', log = 'xy', xlab = '', ylab = '', 
         ylim = range(rejProbsmF3), xaxt = 'n', yaxt = 'n')
    if(m == 1 | m == 3){axis(2, at = alpha, label = alpha)}
    if(m == 2 | m == 5){axis(4, at = alpha, label = alpha)}
    if(m >= 3){axis(1, at = alpha, label = alpha)}
    grid()
    text(paste(m), x = 0.065, y = 0.0006, cex = 3)
    abline(a = 0, b = 1, col = 1, lty = 2, lwd = 2)
    for(n in 2:numNW){
      lines(alpha, rejProbsmF3[m,n,l,], type = 'b', col = n)
    }
  }
  mtext(paste("L = ", L[l]), side = 3, outer = TRUE, line = -1)
  mtext("Significance level (log scale)", side = 1, outer = TRUE, line = 2.2)
  mtext("Rejection probability (log scale)", side = 4, outer = TRUE, line = 2.2)
}