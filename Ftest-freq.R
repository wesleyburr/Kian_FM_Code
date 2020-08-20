source("utility.R")
##########################################################################
#set up frequency modulation and generate data
ndata <- 2000 #N
nppb <- ndata  #use Dave's block function to set this if I do this over blocks
nord <- 10 #This is K = 2NW - 1
mxdeg <- 6 #max degree of polynomial modulation
#CONSIDER ADDING A VARIABLE P=mxdeg+1
js1 <- 69
js2 <- 5318008
seed <- sample(js1:js2, size = 1)
set.seed(seed)

f <- 0.2 #center frequency
W <- (nord + 1)/(2*ndata)
fbw <- 0.8
stm11 <- 2.0/(ndata-1)
Shift <- 0.0
dt <- 1
t <- 0:(ndata-1)        #time
tt <- (t + Shift) * stm11 - 1.0  # this runs from -1 to 1
FMpoly <- W * fbw * (1.0 - 2.0 * tt^2)
Fmod <- f + FMpoly  # quadratic frequency modulation
PhMod <- 2*pi*f*dt*t + 2*pi*cumsum(FMpoly)*dt
#PhMod <- cumsum(Fmod)*dt*2*pi - Fmod*dt*2*pi
cd <- cos(PhMod) + 1i * sin(PhMod) #data without noise

snr <- 5.0									#signal to noise ratio
sigmat2 <- ndata/(2*snr*nord)  # looks like a variance, but of what?   - WHY THIS VARIANCE?
sigmat <- sqrt(sigmat2)
CDN <- cd + rnorm(ndata, sd = sigmat) + rnorm(ndata, sd = sigmat) * 1i

##########################################################################
#slepians, derivatives and special polynomials
DW <- dpss(ndata, nord, ndata * W) #slepian sequences
dw <- DW$v
ev <- DW$eigen

dwp <- dpssp7(dw,nppb,nord,W,ev)

HeSmx <- 0.9
system.time({
  hp <- dpsshp(dw,nppb,nord,mxdeg,HeSmx,norm = TRUE)
})
Poly <- hp$Poly
Ukz <- hp$A

######################################################################################
## NOW HERE I NEED TO TAKE THE COMPLEX DEMODULATE WITH ALL FREQUENCIES. ZERO PAD A LOT
## AND CREATE MATRIX WITH COLUMNS BEING THE COMPLEX DEMODULATE DATA. CARRY OUT THE PROCESS ON ALL 
## COLUMNS OF THE DATA TO GET A MATRIX OF F STATISTICS - DEGREE BY FREQUENCY

nFFT <- 2^ceiling(log2(2*ndata))
freq <- seq(-0.5/dt, 0.5/dt, length.out = nFFT + 1)
freq <- freq[1:nFFT]

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

#pc = c, Poly = G, yk = r = Y - H*c, Ukz = H in '09 paper
pc <- t(Ukz) %*% yk
r <- yk - Ukz %*% pc #this matches the paper, but the original code maybe doesn't
#FitC <- Poly %*% pc + dw %*% r # G*c + V*r. line 11

#P2  <- Arg(FitC)

#phi <- c(1/2,0,-1/2)
#Fr1 <- filter(PhD, phi, method = 'convolution') / (2*pi)
#FrC <- filter(P2, phi, method = 'convolution') / (2*pi)
#FrC[1,] <- rep(0,nFFT)
#FrC[ndata,] <- rep(0,nFFT)




Fcoef <- t(dw) %*% FrA
SSQC <- colSums( Fcoef^2 )
FPcoef <- t(Ukz) %*% Fcoef
#REPLACE THE NEXT BIT FOR USING ALL FREQUENCIES. MATRICES MIGHT NOT WORK HERE

F1stuff <- FLoop2(nord, mxdeg, H, Fcoef, FPcoef, nFFT)
F1 <- matrix(data= F1stuff$F1, nrow = mxdeg+1, ncol = nFFT)

df.ind <- 3
crit.F <- qf(1-1/ndata, 2, 2*(nord-df.ind))
plot(freq, F1[df.ind,], type = 'l')
abline(h = crit.F, lty = 2, col = 2)

# WRITE SUBROUTINE IN FORTRAN. USE FRsd = Fcoef[,i]
output <- Floop(nord, mxdeg, Ukz, Fcoef, FPcoef, dw, SSQC, nFFT, ndata)

  #output <- Fcomp(i,Fcoef[,i],nord,mxdeg,Ukz,FPcoef[,i],dw,SSQC[i],nFFT,ndata)
  
F1 <- matrix(data=output$F1, nrow = mxdeg+1, ncol = nFFT)
SSQRtd <- matrix(data=output$SSQRtd, nrow = mxdeg+1, ncol = nFFT)
  
  #FRsd <- matrix(data=Fcoef[,i], nrow = nord, ncol = mxdeg+1)
  #FRsd <- FRsd - Ukz %*% diag( FPcoef[,i] )
  #ssqr <- colSums(FRsd^2)
  #SSQFK <- SSQC[i] - ssqr
  #F1[,i] <- ( SSQFK/(1:(mxdeg+1)) ) / ( ssqr /(nord - 1:(mxdeg+1)) ) #why do I have two F stats?
  #Rsd <- dw %*% FRsd
  #SSQRtd[,i] <- colSums(Rsd^2)


SSQFtd <- matrix(data=0,nrow = mxdeg+1, ncol = nFFT)
for(L in 1:(mxdeg+1)){
  Fit <- Poly[,1:L,drop=FALSE] %*% FPcoef[1:L,,drop=FALSE]
  SSQFtd[L,] <- colSums( Fit^2 )
}


df1 <- matrix(data = 1:(mxdeg+1), nrow = mxdeg+1, ncol = nFFT, byrow = FALSE)
df2 <- matrix(data = nord - 1:(mxdeg+1), nrow = mxdeg+1, ncol = nFFT, byrow = FALSE)
F2 <- (SSQFtd / df1) / (SSQRtd / df2)

##############################################################

#IN FP TRY REPLACING THE PARTS DEPENDING ON EIGENCOEFFICIENTS OR DATA BY THE FrA AS IN THE
#OTHER TWO F STATISTICS
FP <- matrix(data=0, nrow = mxdeg+1, ncol = nFFT)

for(p in 1:(mxdeg+1)){
  FP[p,] <- (colSums(Mod(Ukz[,1:p,drop=FALSE] %*% pc[1:p,,drop=FALSE])^2) / (2*p) ) /
            (colSums(Mod(yk - Ukz[,1:p,drop=FALSE] %*% pc[1:p, ,drop = FALSE])^2) / (2*(nord - p)) )
}

df.ind <- 3
crit.F <- qf(1-1/ndata, 2*df.ind, 2*(nord-df.ind))

plot(freq,FP[df.ind,], type = 'l')
abline(h = crit.F, lty = 2, col = 2)

plot(freq, F1[df.ind,], type = 'l')
abline(h = crit.F, lty = 2, col = 2)
       
plot(freq, F2[df.ind,], type = 'l')
abline(h = crit.F, lty = 2, col = 2)

spec <- spec.mtm(CDN, k=nord, nw = ndata*W, nFFT = nFFT, dpssIN = dw, Ftest = TRUE, deltat = 1)
plot(spec, Ftest = TRUE)
abline(h = qf(1-1/ndata, 2, 2*(nord-1)), lty = 2, col = 2)
