source("utility.R")
##########################################################################
#set up frequency modulation and generate data
ndata <- 200 #N
nord <- 10 #This is K = 2NW - 1
mxdeg <- 4 #max degree of polynomial modulation
js1 <- 420
js2 <- 5318008
seed <- sample(js1:js2, size = 1)
set.seed(seed)

f <- 0.2  #center frequency
W <- (nord + 1)/(2*ndata)
fbw <- 0.8
stm11 <- 2.0/(ndata-1)
Shift <- 0.0
dt <- 1
t <- 0:(ndata-1)        #time
tt <- (t + Shift) * stm11 - 1.0  # this runs from -1 to 1
FMpoly <- W * fbw * (1.0 - 2.0 * tt^2)
Fmod <- f + FMpoly  # quadratic frequency modulation
PMpoly <- cumsum(FMpoly) - FMpoly[1]
PhMod <- 2*pi*f*dt*t + 2*pi*PMpoly
cd <- complex(real = cos(PhMod),
              imaginary = sin(PhMod)) #data without noise

snr <- 1 											#signal to noise ratio
sigmat2 <- ndata/(2*snr*nord)  
sigmat <- sqrt(sigmat2)
noise <- complex(real = rnorm(ndata, sd = sigmat), imaginary = rnorm(ndata, sd = sigmat))
nFFT <- 2^ceiling(log2(2*ndata))
noise <- c(noise, rep(0, nFFT - ndata))
AR <- c(0, -0.1, rep(0, nFFT - 2))
ARnoise <- fft( fft(AR) * fft(noise), inverse = TRUE) / ndata + noise
ARnoise <- ARnoise[1:ndata]
CDN <- cd + ARnoise
f.demod <- f
CDNdemod <- CDN * exp( -1i*2*pi*f.demod*dt*t )

# logical check: CDN should look (barring noise)
# like exp( i*2*pi*(Fmod - f)*dt*t )
plot(Re(CDNdemod), type = "l")
lines(Re(exp(1i*2*pi*(Fmod-f)*dt*t)), col = "red")
plot(Im(CDNdemod), type = "l")
lines(Im(exp(1i*2*pi*(Fmod-f)*dt*t)), col = "red")

##########################################################################
#slepians, derivatives and special polynomials
DW <- dpss(ndata, nord, ndata * W) #slepian sequences
dw <- DW$v
ev <- DW$eigen

dwp <- dpssp7(dw,ndata,nord,W,ev)

HeSmx <- 0.9
hp <- dpsshp(dw,ndata,nord,mxdeg,HeSmx,norm = TRUE)
Poly <- hp$Poly
H <- hp$A
#########################################################################
#inversion and fitting

yk <- t(dw) %*% CDNdemod

Fit1 <- dw %*% yk #standard inverse
Zdot <- dwp %*% yk #time derivative of Fit1.

# look at the structure of Fit1 spectrally
sp0 <- spec.mtm(as.vector(Fit1), deltat = 1)
plot(sp0$freq, sp0$spec, log = 'y', type = "l",
     xlim = c(-0.03, 0.03))

# compare Fit1 to the actual demod signal
plot(Re(exp(1i*2*pi*(Fmod-f)*dt*t)), type = "l",
     ylab = "Magnitude",
     main = "Compare FM to Projection (Real)")
lines(Re(Fit1), col = "red")

plot(Im(exp(1i*2*pi*(Fmod-f)*dt*t)), type = "l",
     ylab = "Magnitude",
     main = "Compare FM to Projection (Imag)")
lines(Im(Fit1), col = "red")

amp <- Mod(Fit1)
PhD <- Arg(Fit1) 
PhD <- unwrap(PhD)
x1 <- Re(Fit1)
y1 <- Im(Fit1)
Xdot <- Re(Zdot)
Ydot <- Im(Zdot)
FrA <- ( x1*Ydot - Xdot*y1 ) / ((2*pi)*amp^2) 

Fcoef <- t(dw) %*% FrA
FPcoef <- t(H) %*% Fcoef

F1stuff <- F1comp(mxdeg,nord,H,FPcoef,Fcoef)
F1 <- F1stuff$F1
res <- matrix(data = F1stuff$res, nrow = nord, 
              ncol = mxdeg+1)

alpha <- 1 - 1/ndata
alpha <- 0.99
sigF <- qf(alpha, 1, nord - 1:(mxdeg+1))
plot(0:mxdeg, F1, type = 'h', lwd = 3, 
     ylim = c(0, max(c(F1,sigF))))
points(0:mxdeg, sigF, pch = 8, col = 2) 
title("New one, projected thing, differenced ssqr")

Fit <- Poly %*% FPcoef
FitC <- Fit + dw %*% Fcoef - dw %*% H %*% FPcoef

# Compare FM poly to the various fits
plot(t, FMpoly, type = 'l', ylim = range(FMpoly))
lines(Fit, col = 2)
lines(FitC, col = 4)
lines(FrA, col = 5)

sum( (FMpoly - Fit)^2 )
sum( (FMpoly - FitC)^2 )

FitPh1 <- 2*pi*( cumsum(Fit) - Fit[1] )
FitPhC <- 2*pi*( cumsum(FitC) - FitC[1] )
ahat <- mean( amp * cos(PhD - FitPh) )

sigRecon <- ahat * cos(2*pi*f.demod*t + FitPh)

plot(Re(cd[1:100]), type = 'l')
lines(sigRecon[1:100], col = 2)

plot(2*pi*PMpoly, type = 'l')
lines(FitPh1, col = 2)
lines(FitPhC, col = 4)

