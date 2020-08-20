source("utility.R")
##########################################################################
#set up frequency modulation and generate data
ndata <- 2000 #N
nppb <- ndata  #use Dave's block function to set this if I do this over blocks
nord <- 15 #This is K = 2NW - 1
mxdeg <- 6 #max degree of polynomial modulation
#CONSIDER ADDING A VARIABLE P=mxdeg+1
js1 <- 69
js2 <- 5318008
seed <- sample(js1:js2, size = 1)
set.seed(seed)

W <- (nord + 1)/(2*ndata)
fbw <- 0.8
stm11 <- 2.0/(ndata-1)
Shift <- 0.0
dt <- 1
t <- 0:(ndata-1)        #time
tt <- (t + Shift) * stm11 - 1.0  # this runs from -1 to 1
FMpoly <- W * fbw * (1.0 - 2.0 * tt^2)

nFFT <- 2^ceiling(log2(2*ndata))
freq <- seq(-0.5/dt, 0.5/dt, length.out = nFFT + 1)
freq <- freq[1:nFFT]

f.ind <- which.min(abs(freq-0.2))
f <- freq[f.ind] #center frequency

Fmod <- f + FMpoly  # quadratic frequency modulation
PhMod <- 2*pi*f*dt*t + cumsum(FMpoly)*2*pi*dt

cd <- cos(PhMod)  + 1i * sin(PhMod) #data without noise

snr <- 5									#signal to noise ratio
sigmat2 <- ndata/(2*snr*nord)  
sigmat <- sqrt(sigmat2)
CDN <- cd + rnorm(ndata, sd = sigmat) + rnorm(ndata, sd = sigmat) * 1i

##########################################################################
#slepians, derivatives and special polynomials
DW <- dpss(ndata, nord, ndata * W) #slepian sequences
dw <- DW$v
ev <- DW$eigen

dwp <- dpssp7(dw,ndata,nord,W,ev)

HeSmx <- 0.9
hp <- dpsshp(dw,nppb,nord,mxdeg,HeSmx,norm = TRUE)
Poly <- hp$Poly
H <- hp$A

######################################################################################

tapered <- dw * CDN 
tapered <- rbind(tapered, matrix(0, nrow = nFFT-ndata, ncol = nord))
yk <- mvfft(tapered)
yk <- t(yk)
yk <- cbind( yk[,(nFFT/2+2):nFFT], yk[,1:(nFFT/2+1)] )

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

Fcoef <- t(dw) %*% FrA
FPcoef <- t(H) %*% Fcoef

F1stuff <- FLoop2(nord, mxdeg, H, Fcoef, FPcoef, nFFT)
F1 <- matrix(data= F1stuff$F1, nrow = mxdeg+1, ncol = nFFT)

Fit <- Poly %*% FPcoef
FitC <- Fit + dw %*% ( Fcoef - H %*% FPcoef )

df.ind <- 3
crit.F <- qf(1-1/nFFT, 2, 2*(nord-df.ind)) #should it be 1 - 1/ndata or 1 - 1/nFFT ?
plot(freq, F1[df.ind,], type = 'l')
abline(h = crit.F, lty = 2, col = 2)
#lots of spurious lines. might want to see if I can use Dave's quadratic inverse thing to reduce some of
#them


plot(FMpoly, type = 'l')
lines(Fit[,f.ind], col = 2)
lines(FitC[,f.ind], col = 3)
lines(FrA[,f.ind], col = 4)
