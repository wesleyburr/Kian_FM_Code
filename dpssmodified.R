library(multitaper)

nord <- 8
nw <- 4
ndata <- 1000

DW <- dpss(ndata, k = nord, nw = nw)
dw <- DW$v
eigen <- DW$eigen

D <- 2*dw[1,] - dw[2,]

dd <- matrix(data=0, nrow = ndata, ncol = nord - 2)
eigend <- numeric(nord-2)
for(k in 1:(nord-2)){
  L <- k + 2
  g <- D[k]/D[L]
  gsK <- 1/sqrt(1 + g^2)
  gsK2 <- - gsK*g
  
  dd[,k] <- gsK*dw[,k] + gsK2*dw[,L]
  eigend[k] <- gsK^2 * eigen[k] + gsK2^2 * eigen[L]
  
  plot(dw[,k], type = 'l')
  lines(dd[,k], col = 2)
}
DD <- list(v = dd, eigen = eigend)
class(DD) <- "dpss"

nFFT <- 2^ceiling(log2(10*ndata))
Vk <- mvfft(rbind(dw, matrix(data=0, nrow = nFFT-ndata, ncol = nord)))
Vdk <- mvfft(rbind(dd, matrix(data=0, nrow = nFFT-ndata, ncol = nord-2)))

freq <- seq(0,1/2, length.out = nFFT/2)

for(k in 1:(nord-2)){
  plot(freq, Mod(Vk[1:(nFFT/2),k])^2, type = 'l', xlim = range(freq[1:1000]), log = 'y')
  lines(freq, Mod(Vdk[1:(nFFT/2),k])^2, col = 2)
}


HP1 <- dpssap(dd, mxdeg, alpha = 1)
H1 <- HP1$U
G1 <- HP1$R

HP2 <- dpssap(dw, mxdeg, alpha = 1)
H2 <- HP2$U
G2 <- HP2$R

orth <- t(dd) %*% dd

Fit <- G1 %*% t(H1) %*% t(dd) %*% FrA
Fit2 <- G2 %*% t(H2) %*% t(dw) %*% FrA

plot(FMpoly, type = 'l')
lines(Fit[,f.ind], col = 2)
lines(Fit2[,f.ind], col = 4)

sum( (FMpoly - Fit[,f.ind])^2 )
sum( (FMpoly - Fit2[,f.ind])^2 )

### WRITE CODE TO GET THE EIGENVALUES (ENERGY CONCENTRATIONS) FOR THE MODIFIED SLEPIANS
plot(0:(nord-3),eigen[1:(nord-2)] - eigend, log = 'y', xlab = 'K', ylab = 'energy concentration residuals')
