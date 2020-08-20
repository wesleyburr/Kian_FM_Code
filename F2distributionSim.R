source("utility.R")

N <- 10000
nFFT <- 2^ceiling(log2(2*N))


#multitaper and background spectrum stuff
NW <- 4:11
K <- 2*NW - 1
numNW <- length(NW)
W <- NW / N
L <- N
M <- 5

sigmat2 <- 1 #really snr is (A^2/2) / (2W*S(f))
sigmat <- sqrt(sigmat2)

#ar <- c(0.5, -0.3)
#ma <- -0.1

###


nfreqs <- nFFT/2
num.sim <- 10000

#ar <- c(0.6, -0.3)
#ma <- c(1, 0.8, 0.6)
#trnd <- 10*seq(0,1, length.out = N)^2

seed <- 666
set.seed(seed)

freq <- 1:nfreqs / nFFT
subBand <- sample(1:(nfreqs-1), size = 1)

F1array <- F3array <- array(data = 0, dim = c(M+1,num.sim, numNW))
mF1array <- mF3array <- array(data = 0, dim = c(M, num.sim, numNW))

start.time <- Sys.time()
for(n in 1:numNW){
  
  DW <- dpss(n = N, nw = NW[n], k = K[n])
  #DS <- dpss(n = L, nw = NW[n], k = K[n])
  Vdot <- dpssp11R(DW,NW[n])
  AP <- dpssap(DW$v, maxdeg = M, alpha = 1)
  for(i in 1:num.sim){
    #X1 <- rgumbel(n = N/4, 0, 0.1)
    #X2 <- rgumbel(n = N/4, 0, 0.2)
    #X3 <- rgumbel(n = N/4, 0, 0.5)
    #X4 <- rgumbel(n = N/4, 0, 0.8)
    #X.tmp <- c(X1,X2,X3,X4)
    #X <- trnd + arima.sim(model = list(ar=ar,ma=ma), n = N, innov = X.tmp)
    X <- rnorm(N, sd = sigmat)
    tapered <- X * DW$v
    pad <- rbind(tapered, matrix(0, nrow = nFFT-N, ncol = K[n]))
    yk <- mvfft(pad)[subBand+1,,drop = FALSE]
    
    
    modF.test <- ModulatedF20(yk, derivIN = Vdot, apIN = AP, dpssIN = DW)
    
    F1array[,i,n] <- modF.test$F1
    F3array[,i,n] <- modF.test$F3
    mF1array[,i,n] <- modF.test$mF1
    mF3array[,i,n] <- modF.test$mF3
    
    
  }
}
runtime <- Sys.time() - start.time
runtime

sigs <- c(0.9,0.95,0.99, 0.999)
rejProbsF1 <- rejProbsF3 <- array(0, dim = c(length(sigs), M+1, numNW))
rejProbsmF1 <- rejProbsmF3 <- array(0, dim = c(length(sigs), M, numNW))

K <- 2*NW - 1
for(m in 1:(M+1)){
  for(n in 1:numNW){
    for(l in 1:length(sigs)){
      rejProbsF1[l,m,n] <- length( which(F1array[m,,n] > qf(sigs[l], m, K[n]-m)) )/num.sim
      rejProbsF3[l,m,n] <- length( which(F3array[m,,n] > qf(sigs[l], 1, K[n]-m)) )/num.sim
      if(m < M+1){
        rejProbsmF1[l,m,n] <- length( which(mF1array[m,,n] > qf(sigs[l], m, K[n]-m)) )/num.sim
        rejProbsmF3[l,m,n] <- length( which(mF3array[m,,n] > qf(sigs[l], 1, K[n]-m)) )/num.sim
      }
    }
  }
}

par(oma = c(4,4,0,4))
par(mar = c(0,0,1,1))
alpha <- 1-sigs
par(mfrow = c(2,3))
for(m in 1:(M+1)){
  plot(alpha, rejProbsF1[,m,1], type = 'b',  xlab = '', ylab = '', 
       ylim = range(rejProbsF1), xaxt = 'n', yaxt = 'n')
  text(paste(m-1), x = 0.09, y = 0.025, cex = 3)
  if(m == 1 | m == 4){axis(2, at = alpha, label = alpha)}
  if(m == 3 | m == 6){axis(4, at = alpha, label = alpha)}
  if(m > 3){axis(1, at = alpha, label = alpha)}
    
  grid()
  abline(a = 0, b = 1, col = 1, lty = 2, lwd = 2)
  for(n in 2:numNW){
    lines(alpha, rejProbsF1[,m,n], type = 'b', col = n)
  }
}
mtext("Significance level (log scale)", side = 1, outer = TRUE, line = 2.2)
mtext("Rejection probability (log scale)", side = 4, outer = TRUE, line = 2.2)

for(m in 1:(M+1)){
  plot(alpha, rejProbsF3[,m,1], type = 'b', log = 'xy', xlab = '', ylab = '', 
       ylim = range(rejProbsF3), xaxt = 'n', yaxt = 'n')
  text(paste(m-1), x = 0.05, y = 0.001, cex = 3)
  if(m == 1 | m == 4){axis(2, at = alpha, label = alpha)}
  if(m == 3 | m == 6){axis(4, at = alpha, label = alpha)}
  if(m > 3){axis(1, at = alpha, label = alpha)}
  grid()
  abline(a = 0, b = 1, col = 1, lty = 2, lwd = 2)
  for(n in 2:numNW){
    lines(alpha, rejProbsF3[,m,n], type = 'b', col = n)
  }
}
mtext("Significance level (log scale)", side = 1, outer = TRUE, line = 2.2)
mtext("Rejection probability (log scale)", side = 4, outer = TRUE, line = 2.2)

par(mfrow = c(2,3))
frame()
for(m in 1:M){
  plot(alpha, rejProbsmF1[,m,1], type = 'b', log = 'xy', xlab = '', ylab = '', 
       ylim = range(rejProbsmF1), xaxt = 'n', yaxt = 'n', las = 1)
  text(paste(m), x = 0.05, y = 0.001, cex = 3)
  if(m == 1 | m == 3){axis(2, at = alpha, label = alpha)}
  if(m == 2 | m == 5){axis(4, at = alpha, label = alpha)}
  if(m >= 3){axis(1, at = alpha, label = alpha)}
  grid()
  abline(a = 0, b = 1, col = 1, lty = 2, lwd = 2)
  for(n in 2:numNW){
    lines(alpha, rejProbsmF1[,m,n], type = 'b', col = n)
  }
}
mtext("Significance level (log scale)", side = 1, outer = TRUE, line = 2.2)
mtext("Rejection probability (log scale)", side = 4, outer = TRUE, line = 2.2)

par(mfrow = c(2,3))
rnge.tmp <- as.vector(rejProbsmF3)
rnge <- range(rnge.tmp[-which(rnge.tmp == 0)])
frame()
for(m in 1:M){
  plot(alpha, rejProbsmF3[,m,1], type = 'b', log = 'xy', xlab = '', ylab = '', 
       ylim = rnge, xaxt = 'n', yaxt = 'n')
  text(paste(m), x = 0.05, y = 0.001, cex = 3)
  if(m == 1 | m == 3){axis(2, at = alpha, label = alpha)}
  if(m == 2 | m == 5){axis(4, at = alpha, label = alpha)}
  if(m >= 3){axis(1, at = alpha, label = alpha)}
  grid()
  abline(a = 0, b = 1, col = 1, lty = 2, lwd = 2)
  for(n in 2:(numNW-1)){
    lines(alpha, rejProbsmF3[,m,n], type = 'b', col = n)
  }
}
mtext("Significance level (log scale)", side = 1, outer = TRUE, line = 2.2)
mtext("Rejection probability (log scale)", side = 4, outer = TRUE, line = 2.2)

#rejProbsNS <- list(rpF1 = rejProbsF1, rpF3 = rejProbsF3, rpmF1 = rejProbsmF1,
                   #rpmF3 = rejProbsmF3)

#save(rejProbsNS, file = "rpNS.RData")

testStatNS <- list(F1 = F1array, F3 = F3array, mF1 = mF1array, mF3 = mF3array)
save(testStatNS, file = "testStatNS.RData")


#rejProbsGaussian <- list(rpF1 = rejProbsF1, rpF3 = rejProbsF3, rpmF1 = rejProbsmF1,
 #                        rpmF3 = rejProbsmF3)
#save(rejProbsGaussian, file = "rpGaussian.RData")
######################################################
### make a handful of histograms


#F1
par(mfrow = c(2,2))
x1 <- seq(0, 10, length.out = num.sim)
hist(F1array[2,,2], breaks = 1000, freq = FALSE, main = '', xlim = c(0,10), xaxt = 'n', las = 1)
lines(x1, df(x1, 2, K[2]-2), col = 2)
text("P = 1, NW = 5", x = 6, y = 0.4)

hist(F1array[6,,2], breaks = 3000, freq = FALSE, main = '', xlim = c(0,10), yaxt = 'n', xaxt = 'n')
lines(x1, df(x1, 6, K[2]-6), col = 2)
text("P = 5, NW = 5", x = 6, y = 0.3)
axis(4)

hist(F1array[2,,7], breaks = 500, freq = FALSE, main = '', xlim = c(0,10))
lines(x1, df(x1, 2, K[7]-2), col = 2)
text("P = 1, NW = 10", x = 6, y = 0.4)

hist(F1array[6,,7], breaks = 250, freq = FALSE, main = '', xlim = c(0,10), yaxt = 'n')
lines(x1, df(x1, 6, K[7]-6), col = 2)
text("P = 5, NW = 10", x = 6, y = 0.4)
axis(4)

mtext("Density", side = 2, outer = TRUE, line = 2.2)

######
# F3
par(mfrow = c(2,2))
x1 <- seq(0, 7, length.out = num.sim)
hist(F3array[2,,2], breaks = 500, freq = FALSE, main = '', xlim = c(0,7), xaxt = 'n', las = 1)
lines(x1, df(x1, 1, K[2]-2), col = 2)
text("P = 1, NW = 5", x = 3, y = 1)

hist(F3array[6,,2], breaks = 5000, freq = FALSE, main = '', xlim = c(0,7), yaxt = 'n', xaxt = 'n')
lines(x1, df(x1, 1, K[2]-6), col = 2)
axis(4, las = 1)
text("P = 5, NW = 5", x = 3, y = 1)

hist(F3array[2,,7], breaks = 500, freq = FALSE, main = '', xlim = c(0,7))
lines(x1, df(x1, 1, K[7]-2), col = 2)
text("P = 1, NW = 10", x = 3, y = 1)

hist(F3array[5,,7], breaks = 250, freq = FALSE, main = '', xlim = c(0,7), yaxt = 'n')
lines(x1, df(x1, 1, K[7]-5), col = 2)
text("P = 5, NW = 10", x = 3, y = 1)
axis(4, las = 1)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

#######
## mF1
par(mfrow = c(2,2))
x1 <- seq(0, 10, length.out = num.sim)
hist(mF1array[1,,2], breaks = 150, freq = FALSE, main = '', xlim = c(0,10), xaxt = 'n',las=1)
lines(x1, df(x1, 1, K[2]-1), col = 2)
text("P = 1, NW = 5", x = 5, y = 1)


hist(mF1array[5,,2], breaks = 2000, freq = FALSE, main = '', xlim = c(0,10), yaxt = 'n', xaxt = 'n')
lines(x1, df(x1, 5, K[2]-5), col = 2)
axis(4,las = 1)
text("P = 5, NW = 5", x = 5, y = 0.4)


hist(mF1array[1,,7], breaks = 250, freq = FALSE, main = '', xlim = c(0,10),las=1)
lines(x1, df(x1, 1, K[7]-1), col = 2)
text("P = 1, NW = 10", x = 5, y = 1)


hist(mF1array[5,,7], breaks = 100, freq = FALSE, main = '', xlim = c(0,10), yaxt = 'n')
lines(x1, df(x1, 5, K[7]-5), col = 2)
axis(4,las = 1)
text("P = 5, NW = 10", x = 5, y = 0.4)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

######
## mF3
par(mfrow = c(2,2))
x1 <- seq(0, 7, length.out = num.sim)
hist(mF3array[1,,2], breaks = 400, freq = FALSE, main = '', xlim = c(0,7), xaxt = 'n', las = 1)
lines(x1, df(x1, 1, K[2]-1), col = 2)
text("P = 1, NW = 5", x = 3, y = 1)

hist(mF3array[5,,2], breaks = 4000, freq = FALSE, main = '', xlim = c(0,7), yaxt = 'n', 
     xaxt = 'n')
lines(x1, df(x1, 1, K[2]-5), col = 2)
axis(4, las = 1)
text("P = 5, NW = 5", x = 3, y = 1)

hist(mF3array[1,,7], breaks = 500, freq = FALSE, main = '', xlim = c(0,7), las = 1)
lines(x1, df(x1, 1, K[7]-1), col = 2)
text("P = 1, NW = 10", x = 3, y = 1)

hist(mF3array[5,,7], breaks = 500, freq = FALSE, main = '', xlim = c(0,7), yaxt = 'n')
lines(x1, df(x1, 1, K[7]-5), col = 2)
text("P = 5, NW = 10", x = 3, y = 1)
axis(4, las = 1)
mtext("Density", side = 2, outer = TRUE, line = 2.2)
####################################################
#old stuff
pdf("Fhistograms_downsampled.pdf", width = 11, height = 8.5)
  par(mfrow = c(3,1))
  par(oma = c(4, 4, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
  par(mar = c(0, 0, 1, 1)) 
  
  f1Mesh <- seq(0, max(F2array), length.out = num.sim)
  f3Mesh <- seq(0, max(F3array), length.out = num.sim)
  
  for(m in 1:(M+1)){
    hist(F1array[m,], freq = FALSE, breaks = 200, xlab = '', ylab = 'F1', main = '')
    lines(density(F1array[m,], kernel = "optcosine"), col = 2, lwd = 2)
    lines(f1Mesh, df(f1Mesh, m, K-m), col = 4, lwd = 2)
    title(paste("degree ", m-1))

    hist(F2array[m,], freq = FALSE, breaks = 200, xlab = '', ylab = 'F2', main = '')
    lines(density(F2array[m,], kernel = "optcosine"), col = 2, lwd = 2)
    lines(f1Mesh, df(f1Mesh, m, K-m), col = 4, lwd = 2)

    hist(F3array[m,], freq = FALSE, breaks = 200, xlab = '', ylab = 'F3', main = '')
    lines(density(F3array[m,], kernel = "optcosine"), col = 2, lwd = 2)
    lines(f3Mesh, df(f3Mesh, 1, K-m), col = 4, lwd = 2)
    
  }



dev.off()






#how should I estimate the quantiles?

F1sort <- F2sort <- F3sort <- array(data = 0, dim = c(M+1,num.sim)) 
for(i in 1:(M+1)){
  #for(j in 1:length(subBand)){
    F1sort[i,] <- sort(F1array[i,])
    F2sort[i,] <- sort(F2array[i,])
    F3sort[i,] <- sort(F3array[i,])
  #}
}



#compute mean of each quantile over degree and frequency? No. this is wrong
#it should be basically the same at all frequencies since I started with white noise.
F1quant.mean <- F2quant.mean <- F3quant.mean <- matrix(data = 0, nrow = M+1, ncol = num.sim)
F1quant.sd <- F2quant.sd <- F3quant.sd <- matrix(data = 0, nrow = M+1, ncol = num.sim)
for(i in 1:(M+1)){
  F1quant.mean[i,] <- apply(F1sort[i,,], MARGIN = 2, FUN = mean)
  F2quant.mean[i,] <- apply(F2sort[i,,], MARGIN = 2, FUN = mean)
  F3quant.mean[i,] <- apply(F3sort[i,,], MARGIN = 2, FUN = mean)
  F1quant.sd[i,] <- apply(F1sort[i,,], MARGIN = 2, FUN = sd)
  F2quant.sd[i,] <- apply(F2sort[i,,], MARGIN = 2, FUN = sd)
  F3quant.sd[i,] <- apply(F3sort[i,,], MARGIN = 2, FUN = sd)
}

F1CI <- F2CI <- F3CI <- array(NA, dim = c(M+1, num.sim, 2))

for(i in 1:(M+1)){
  for(j in 1:2){
    F1CI[i,,j] <- F1quant.mean[i,] +  (-1)^j * qt(0.995, length(subBand)-1)*F1quant.sd[i,]/sqrt(length(subBand))
    F2CI[i,,j] <- F2quant.mean[i,] +  (-1)^j * qt(0.995, length(subBand)-1)*F2quant.sd[i,]/sqrt(length(subBand))           
    F3CI[i,,j] <- F3quant.mean[i,] +  (-1)^j * qt(0.995, length(subBand)-1)*F3quant.sd[i,]/sqrt(length(subBand))
  }
}

quants <- list(F1quant = F1quant.mean, F2quant = F2quant.mean, F3quant = F3quant.mean,
               HFquant = HFquant.mean)
save(quants, file = "Fquants.RData")


##### let's see if I can get a better handle on the distributions

F1mean <- F2mean <- F3mean <- matrix(data = NA, nrow = M+1, ncol = nfreqs-1)


for(m in 1:(M+1)){
  F1mean[m,] <- apply(F1array[m,,], MARGIN = 1, FUN = mean)
  F2mean[m,] <- apply(F2array[m,,], MARGIN = 1, FUN = mean)
  F3mean[m,] <- apply(F3array[m,,], MARGIN = 1, FUN = mean)
}

F1var <- F2var <- F3var <- matrix(data = NA, nrow = M+1, ncol = nfreqs-1)
HFvar <- apply(HFarray, MARGIN = 1, FUN = var)

for(m in 1:(M+1)){
  F1var[m,] <- apply(F1array[m,,], MARGIN = 1, FUN = var)
  F2var[m,] <- apply(F2array[m,,], MARGIN = 1, FUN = var)
  F3var[m,] <- apply(F3array[m,,], MARGIN = 1, FUN = var)
}
