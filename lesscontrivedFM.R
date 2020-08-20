source("utility.R")
##########################################################################
#set up frequency modulation and generate data
N <- 2000 #N
K <- 10 #This is K = 2NW - 1
M <- 4 #max degree of polynomial modulation

js1 <- 420
js2 <- 5318008
seed <- sample(js1:js2, size = 1)
set.seed(seed)

nFFT <- 2^ceiling(log2(2*N))

A <- 1
f <- 0.3 #center frequency
B <- 0.002
rho <- 0.8
stm11 <- 2.0/(N-1)
dt <- 1
t <- 0:(N-1)        #time
tt <- t * stm11 - 1.0  # this runs from -1 to 1
phi.prime <- B * rho * (1.0 - 2.0 * tt^2)
Fmod <- f + phi.prime  # quadratic frequency modulation
phi <- cumsum(phi.prime)*2*pi*dt
Phi <- 2*pi*f*dt*t + phi

sig <- A*cos(Phi) 

sigma <- 3
noise <- rnorm(N, mean = 0, sd = sigma)

x <- sig + noise

w <- (K+1)/(2*N)
snr <- A^2 / (4*w*sigma^2)

sp.X <- spec.mtm(x, k = K, nw = N*w, deltat = 1, Ftest = TRUE, nFFT = nFFT, returnInternals = TRUE)
plot(sp.X, Ftest = TRUE)
abline(h = qf(c(0.99,0.999,1-1/N),2,2*K-2), lty = 2, col = 2)

V <- sp.X$mtm$dpss$v
lambda <- sp.X$mtm$dpss$eigen
Y.unweighted <- t(sp.X$mtm$eigenCoefs)
dk <- t(sqrt(sp.X$mtm$eigenCoefWt))
Y <- dk * Y.unweighted

Z <- V %*% Y
U <- Re(Z)
W <- Im(Z)

Hp <- dpsshp(V,N,K,M,1,norm = TRUE)
H <- Hp$A
G <- Hp$Poly

theta.temp <- atan2(W,U)
theta <- apply(theta.temp, MARGIN = 2, FUN = unwrap)
Tcoef <- t(V) %*% theta
TPcoef <- t(H) %*% Tcoef
SSQC <- colSums(Tcoef^2)

FTdjtstuff <- FLoopDJT(K,M,H,Tcoef,TPcoef,SSQC,nFFT/2)
FTdjt <- matrix(data = FTdjtstuff$F1mat, nrow = M+1, ncol = nFFT/2)

FTstuff <- FLoop2(K,M,H,Tcoef,TPcoef,nFFT/2)
FT <- matrix(data = FTstuff$F1, nrow = M+1, ncol = nFFT/2)

amp <- Mod(Z)

freq <- sp.X$freq[1:(nFFT/2)]
f.ind <- which.min(abs(freq-f))

Vdot <- dpssp7(V,N,K,w,lambda)
Zdot <- Vdot %*% Y
Udot <- Re(Zdot)
Wdot <- Im(Zdot)

psi <- ( Wdot*U - Udot*W) / (2*pi*amp^2)

plot(phi.prime, type = 'l')
lines(psi[,f.ind], col = 2)



Psi <- t(V) %*% psi
Psi <- Psi[,1:(nFFT/2)]
PsiP <- t(H) %*% Psi
PsiP <- PsiP[,1:(nFFT/2)]

F1stuff <- FLoop2(K,M,H,Psi,PsiP,nFFT/2)
F1 <- matrix(data= F1stuff$F1, nrow = M+1, ncol = nFFT/2)

Fit <- G[,1:3] %*% PsiP[1:3,]
FitC <- Fit + V %*% (Psi - H[,1:3] %*% PsiP[1:3,])
lines(Fit[,f.ind], col = 6)
lines(FitC[,f.ind], col = 5)

png('Fcurvelinear.PNG')
  plot(freq,rep(1,nFFT/2), type = 'l', col = 0, ylim = c(0.1,450), xlab = 'Frequency',
       ylab = 'Polynomial F', xlim = c(f-w,f+w))
  for(m in 1:(M+1)){
    lines(freq, F1[m,], col = m)
  }
  legend(legend = 0:M, fill = 1:(M+1), x = 0.290, y = 450, ncol = M+1, border = 'white', bty = 'n')
  title("Polynomial F at different degrees in a band around f=0.3,
        normal scale")
dev.off()

png('Fcurvelog.PNG')
  plot(freq,rep(1,nFFT/2), type = 'l', col = 0, log = 'y', ylim = c(0.1,450), xlab = 'Frequency',
       ylab = 'Polynomial F', xlim = c(f-w,f+w))
  for(m in 1:(M+1)){
    lines(freq, F1[m,], col = m)
  }
  legend(legend = 0:M, fill = 1:(M+1), x = 0.289, y = 450, ncol = M+1, border = 'white', bty = 'n')
  title("Polynomial F at different degrees in a band around f=0.3,
        log scale")
dev.off()
