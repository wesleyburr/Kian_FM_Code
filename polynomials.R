source("utility.R")
##########################################################################
#set up frequency modulation and generate data
ndata <- 2000 #N
nppb <- ndata  #use Dave's block function to set this if I do this over blocks
nord <- 10 #This is K = 2NW - 1
mxdeg <- 2 #max degree of polynomial modulation

js1 <- 420
js2 <- 5318008
seed <- sample(js1:js2, size = 1)
set.seed(seed)

nFFT <- 2^ceiling(log2(2*ndata))
freq <- seq(0, 0.5, length.out = nFFT + 1)
freq <- freq[1:nFFT]

A <- 2
f <- freq[which.min(abs(freq-0.2))] #center frequency
W <- (nord + 1)/(2*ndata)
fbw <- 0.8
stm11 <- 2.0/(ndata-1)
Shift <- 0.0
dt <- 1
t <- 0:(ndata-1)        #time
tt <- (t + Shift) * stm11 - 1.0  # this runs from -1 to 1
FMpoly <- W * fbw * (1.0 - 2.0 * tt^2)
Fmod <- f + FMpoly  # quadratic frequency modulation
Pref <- cumsum(FMpoly)*2*pi*dt
PhMod <- 2*pi*f*dt*t + Pref

cd <- A*cos(PhMod)  #+ 1i * sin(PhMod) #data without noise

snr <- 5									#signal to noise ratio
sigmat2 <- ndata/(2*snr*nord)  
sigmat <- sqrt(sigmat2)
CDN <- cd + rnorm(ndata, sd = sigmat) #+ rnorm(ndata, sd = sigmat) * 1i

##########################################################################
#slepians, derivatives and special polynomials
DW <- dpss(ndata, nord, ndata * W) #slepian sequences
dw <- DW$v
ev <- DW$eigen

dwp <- dpssp7(dw,ndata,nord,W,ev)

######################################################################################

tapered <- dw * CDN 
tapered <- rbind(tapered, matrix(0, nrow = 2*nFFT-ndata, ncol = nord))
yk <- mvfft(tapered)
yk <- t(yk)
yk <- yk[,1:nFFT]

#############################################################################
Fit1 <- dw %*% yk #standard inverse
Zdot <- dwp %*% yk #time derivative of Fit1.

amp <- Mod(Fit1)
PhD <- Arg(Fit1) #phase of standard inverse
x1 <- Re(Fit1)
y1 <- Im(Fit1)
Xdot <- Re(Zdot)
Ydot <- Im(Zdot)
FrA <- ( x1*Ydot - Xdot*y1 ) / ((2*pi)*amp^2)


#### POLYNOMIALS
#############################################################################

ap <- dpssap(dw,mxdeg,alpha=1) #CHEBYSHEV POLYNOMIALS OF THE SECOND KIND
U <- ap$U
R <- ap$R

hp <- dpsshp(dw,ndata,nord,mxdeg,1,norm=TRUE)  #HERMITE POLYNOMIALS
H <- hp$A
G <- hp$Poly

#TRY A FEW MORE. DIFFERENT GEGENBAUER POLYNOMIALS, CHEBYSHEV OF THE FIRST KIND

tp <- dpsstp(dw,mxdeg)  #CHEBYSHEV OF THE FIRST KIND
Ht <- tp$U
Tp <- tp$R

hp2 <- dpsshp2(dw,mxdeg) #HERMITE (WES' CODE)
H2 <- hp2$U
G2 <- hp2$R

lp <- dpssap(dw,mxdeg,alpha=0.5)  #LEGENDRE POLYNOMIALS
Hl <- lp$U
L <- lp$R

shp <- dpssap(dw,mxdeg,alpha = (ndata-2)/2) #SPHERICAL HARMONICS
Hs <- shp$U
S <- shp$R

Fcoef <- t(dw) %*% FrA

FPUcoef <- t(U) %*% Fcoef
FPHcoef <- t(H) %*% Fcoef
FPTcoef <- t(Ht) %*% Fcoef
FPH2coef <- t(H2) %*% Fcoef
FPLcoef <- t(Hl) %*% Fcoef
FPScoef <- t(Hs) %*% Fcoef

FitU <- R %*% FPUcoef
FitH <- G %*% FPHcoef
FitT <- Tp %*% FPTcoef
FitH2 <- G2 %*% FPH2coef
FitL <- L %*% FPLcoef
FitS <- S %*% FPScoef

f.ind <- which(freq == f)

plot(FMpoly, type = 'l')
lines(FitU[,f.ind], col = 2)
lines(FitH[,f.ind], col = 3)
lines(FitT[,f.ind], col = 4)
lines(FitH2[,f.ind], col = 5)
lines(FitL[,f.ind], col = 6)
lines(FitS[,f.ind], col = 7)

#polynomials don't seem to differ all that much - if at all. why?

plot(R[,mxdeg+1], type = 'l')
lines(G[,mxdeg+1], col = 2)
lines(Tp[,mxdeg+1], col = 3)
lines(G2[,mxdeg+1], col = 4)
lines(L[,mxdeg+1], col = 5)
lines(S[,mxdeg+1], col = 6)


#########################################
K <- 7
P <- 5
N <- 100
U <- matrix(data = 0, nrow = K, ncol = P)
raw.poly <- poly(tt, degree = P-1, raw = TRUE)
R <- cbind(rep(1,N), raw.poly)
dw <- dpss(n = N, k = K, nw = (K+1)/2)
V <- dw$v
t <- 0:(N-1)
scale <- 2/(N-1)
tt <- scale * t - 1

# Inner Products of R and V
for(L in 1:P) {
  Kmin <- ( (L-1) %% 2 ) + 1
  for(k in seq(Kmin, K, 2)) {  # loop on non-zero Slepians
    U[k, L] <- t(V[, k]) %*% R[, L]
  }
}

# Degree 0, 1 (manual) -- L = degree+1
for(L in 1:min(2,P)) {
  scl <- 1 / sqrt( sum(U[, L]^2) )
  U[, L] <- U[, L] * scl # orthonormalize
  R[, L] <- R[, L] * scl
}

# loop on higher degrees, applying Gram-Schmidt only on similar
# parity functions (as even/odd are already orthogonal in U)
if( P > 2 ) {
  for(L in 3:P) {
    if(L %% 2 == 0) {
      Kmin <- 2
    } else {
      Kmin <- 1
    }
    for(j in seq(Kmin, L-1, 2)) {
      scl <- sum( U[, L] * U[, j] )
      U[, L] <- U[, L] - scl * U[, j] # Gram-Schmidt
      R[, L] <- R[, L] - scl * R[, j]
    }
    scl <- 1 / sqrt(sum(U[, L]^2))
    U[, L] <- U[, L] * scl  # orthonormalize
    R[, L] <- R[, L] * scl
  }
}

AP <- dpssap(V,P-1,alpha=1)
H <- AP$U
G <- AP$R

plot(tt, rep(0,N), ylim = range(R), col = 'white', xlab = NA, ylab = NA)
for(p in 1:P){
  lines(tt,R[,p], col = p, lty = 1, lwd = 1)
  lines(tt,G[,p], col = p, lty = 2, lwd = 2)
}

plot(1:K, rep(0,K), ylim = range(H), col = 'white', xlab = NA, ylab = NA)
for(p in 1:P){
  lines(1:K, U[,p], col = p, lty = 1, lwd = 1)
  lines(1:K, H[,p], col = p, lty = 2, lwd = 2)
}