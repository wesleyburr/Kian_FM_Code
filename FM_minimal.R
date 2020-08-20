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
PhMod <- 2*pi*f*dt*t + cumsum(FMpoly)*2*pi*dt
cd <- cos(PhMod) + 1i * sin(PhMod) #data without noise


snr <- 5.0 											#signal to noise ratio
sigmat2 <- ndata/(2*snr*nord)  # looks like a variance, but of what?   - WHY THIS VARIANCE?
sigmat <- sqrt(sigmat2)
CDN <- cd + rnorm(ndata, sd = sigmat) + rnorm(ndata, sd = sigmat) * 1i
f.demod <- f
CDN <- CDN * exp( -1i*2*pi*f.demod*dt*t )
##########################################################################
#slepians, derivatives and special polynomials
DW <- dpss(ndata, nord, ndata * W) #slepian sequences
dw <- DW$v
ev <- DW$eigen

dwp <- dpssp7(dw,ndata,nord,W,ev)
#remove the part where it computes the eigenvalues and use the ones from dpss

HeSmx <- 0.9
hp <- dpsshp(dw,ndata,nord,mxdeg,HeSmx,norm = TRUE)
Poly <- hp$Poly
H <- hp$A

#ap <- dpssap(dw,mxdeg)
#Poly <- ap$R
#H <- ap$U

#########################################################################
#inversion and fitting - messy still

yk <- t(dw) %*% CDN

Fit1 <- dw %*% yk #standard inverse
Zdot <- dwp %*% yk #time derivative of Fit1.

amp <- Mod(Fit1)
PhD <- Arg(Fit1) #phase of standard inverse
PhD <- unwrap(PhD)
x1 <- Re(Fit1)
y1 <- Im(Fit1)
Xdot <- Re(Zdot)
Ydot <- Im(Zdot)
FrA <- ( x1*Ydot - Xdot*y1 ) / ((2*pi)*amp^2) #derivative of phase? -> instantaneous frequency
#this is for the standard inverse fit
####

#pc = c, Poly = G, yk = r = Y - H*c, Ukz = H in '09 paper
pc <- t(H) %*% yk
r <- yk - H %*% pc #this matches the paper, but the original code maybe doesn't
#for (jd in 0:mxdeg){
#  pc[jd+1] <- sum( r * H[,jd+1]) 
#  r <- r - H[,jd+1]*pc[jd+1]
#}

FitC <- Poly %*% pc + dw %*% r # G*c + V*r. line 11  

P2  <- Arg(FitC)
P2 <- unwrap(P2)

phi <- c(1/2,0,-1/2)
Fr1 <- filter(PhD, phi, method = 'convolution') / (2*pi)
FrC <- filter(P2, phi, method='convolution') / (2*pi)

#FrA[click(FrA,0.1)] <- NA
#Fr1[click(Fr1,0.1)] <- NA
#FrC[click(FrC,0.1)] <- NA

#FrA <- na.interp(FrA) #not sure this helps with anything.
#Fr1 <- na.interp(Fr1)
#FrC <- na.interp(FrC)

Fcoef <- t(dw) %*% FrA
#Fcoef <- Fcoef[,1]
SSQC <- colSums( Fcoef^2 )
FPcoef <- t(H) %*% Fcoef


#Rsd <- matrix(data=0, nrow = ndata, ncol = mxdeg+1)
#ssqr <- numeric(mxdeg+1)
#for(L in 1:(mxdeg+1)){
#  FRsd <- Fcoef - H[,L] * FPcoef[L]
#  ssqr[L] <- sum(FRsd^2)
#  Rsd[,L] <- dw %*% FRsd
#}

#start.time <- Sys.time()
ssqr <- numeric(mxdeg+2)
ssqr[1] <- sum(Fcoef^2)
F1 <- numeric(mxdeg+1)

for(L in 1:(mxdeg+1)){
  projFcoef <- H[,1:L, drop = FALSE] %*% FPcoef[1:L, drop = FALSE]
  ssqr[L+1] <- sum( (Fcoef - projFcoef)^2 )
  F1[L] <- ((ssqr[L] - ssqr[L+1])/ 2) / ( ssqr[L+1] / (2*(nord-L)))
}
#end.time <- Sys.time()
#F1oldruntime <- end.time - start.time
#F1oldruntime


start.time <- Sys.time()
F1stuff <- F1comp(mxdeg,nord,H,FPcoef,Fcoef)
end.time <- Sys.time()
F1runtime <- end.time - start.time
F1runtime
F1 <- F1stuff$F1
res <- matrix(data = F1stuff$res, nrow = nord, ncol = mxdeg+1)

F2stuff <- F2comp(dw, ndata, nord, res, Fcoef, mxdeg)
F2 <- F2stuff$F2

yPk <- t(H) %*% yk
Ukz <- apply(dw, MARGIN = 2, FUN = sum)
FPstuff <- FPcomp(yk=yk,H=H,Ukz=Ukz,nord=nord,mxdeg=mxdeg,yPk=yPk)


alpha <- 0.99
sigF <- qf(alpha, 2, 2*(nord - 1:(mxdeg+1)))
plot(0:mxdeg, F1, type = 'h', lwd = 3, ylim = c(0, max(c(F1,sigF))))
points(0:mxdeg, sigF, pch = 8, col = 2) 
title("New one, projected thing, differenced ssqr")


SSQFK <- numeric(mxdeg+1)
ssqr <- numeric(mxdeg+1)
for(L in 1:(mxdeg+1)){
  projFcoef <- H[,1:L, drop = FALSE] %*% FPcoef[1:L, drop = FALSE]
  SSQFK[L] <- sum( projFcoef^2 )
  ssqr[L] <- sum( (Fcoef - projFcoef)^2 )
}
F1 <- ( SSQFK / 1:(mxdeg+1) ) / (ssqr / (nord - 1:(mxdeg+1)))
sigF <- qf(alpha, 2*(1:mxdeg+1), 2*(nord - 1:(mxdeg+1)))
plot(0:mxdeg, F1, type = 'h', lwd = 3, ylim = c(0, max(c(F1,sigF))))
points(0:mxdeg, sigF, pch = 8, col = 2)
title("new one, projected thing, regular ssq")


#SSQFK <- cumsum(FPcoef^2) #MAKE NOTE OF THIS. CHANGED TO A MORE CONVENTIONAL SSQ


FRsd <- matrix(data=Fcoef, nrow = nord, ncol = mxdeg+1)
FRsd <- FRsd - H %*% diag( FPcoef[,1] )
ssqr <- colSums(FRsd^2)
SSQFK <- SSQC - ssqr 
F1 <- ( SSQFK/(1:(mxdeg+1)) ) / ( ssqr /(nord - 1:(mxdeg+1)) ) 
plot(0:mxdeg, F1, type = 'h', lwd = 3, ylim = c(0, max(c(F1,sigF))))
points(0:mxdeg, sigF, pch = 8, col = 2) 
title("old one, projected thing maybe, weird numerator")


Rsd <- dw %*% FRsd
SSQRtd <- colSums(Rsd^2)

SSQyk <- sum( Mod(yk)^2 )
FPRsd <- matrix(data=yk, nrow = nord, ncol = mxdeg+1)
FPRsd <- FPRsd - H %*% diag( pc[,1] )
ssqr <- colSums(Mod(FPRsd)^2)
SSQFP <- SSQyk - ssqr 
FP <- ( SSQFP/(1:(mxdeg+1)) ) / ( ssqr /(nord - 1:(mxdeg+1)) ) 
#FP <- numeric(mxdeg+1)
#for(p in 1:(mxdeg+1)){
##  ykproj <- H %*% diag(pc[,1])
#  SSQFP <- sum( Mod(ykproj)^2 )
#  res <- yk - ykproj
#  SSQRP <- sum( Mod(res)^2 )
#  FP[p] <- ( SSQFP / (2*p) ) / ( SSQRP / (2*(nord-p)))
#}

plot(0:mxdeg, FP, type = 'h', lwd = 3, ylim = c(0, max(c(FP,sigF))))
points(0:mxdeg, sigF, pch = 8, col = 2)
title("FP in the spirit of Dave's F1")

FP <- numeric(mxdeg+1)
for(p in 1:(mxdeg+1)){
  ykproj <- H[,1:p, drop = FALSE] %*% pc[1:p, drop = FALSE]
  SSQFP <- sum( Mod(ykproj)^2 )
  res <- yk - ykproj
  SSQRP <- sum( Mod(res)^2 )
  FP[p] <- ( SSQFP / (2*p) ) / ( SSQRP / (2*(nord-p)))
}

plot(0:mxdeg, FP, type = 'h', lwd = 3, ylim = c(0, max(c(FP,sigF))))
points(0:mxdeg, sigF, pch = 8, col = 2)
title("my original FP")

FP <- numeric(mxdeg+1)
SSQRP <- numeric(mxdeg+2)
SSQRP[1] <- sum( Mod(yk)^2 )
for(p in 1:(mxdeg+1)){
  ykproj <- H[,1:p, drop = FALSE] %*% pc[1:p, drop = FALSE]
  res <- yk - ykproj
  SSQRP[p+1] <- sum( Mod(res)^2 )
  FP[p] <- (( SSQRP[p] - SSQRP[p+1] ) / 2 ) / ( SSQRP[p+1] / (2*(nord-p)))
}
sigF <- qf(alpha, 2, 2*(nord-1:(mxdeg+1)))
plot(0:mxdeg, FP, type = 'h', lwd = 3, ylim = c(0, max(c(FP,sigF))))
points(0:mxdeg, sigF, pch = 8, col = 2)
title("FP in the style of my new F1")

Fit <- 0 
#Fit2 <- 0
SSQFtd <- numeric(mxdeg+1)
#SSQRtd <- SSQFtd
for(L in 1:(mxdeg+1)){#compute F statistic for each degree
  Fit <- Poly[,1:L, drop = FALSE] %*% FPcoef[1:L, drop = FALSE]
  FitC <- Fit + dw %*% (Fcoef - H[,1:L, drop = FALSE] %*% FPcoef[1:L, drop = FALSE])
  #Fit2 <- Poly[,1:(L+1), drop = FALSE] %*% FPcoef[1:(L+1), drop = FALSE]
  #SSQFtd[L+1] <- sum( ( dw %*% t(dw) %*% Fit )^2 )
  SSQFtd[L] <- sum( Fit^2 )
  
}


sigF <- qf(alpha, 2*(0:mxdeg+1), 2*(nord-1:(mxdeg+1)))
F2 <- ( SSQFtd/(1:(mxdeg+1)) ) / ( SSQRtd/(nord - 1:(mxdeg+1)) )
plot(0:mxdeg, F2, type = 'h', lwd = 3)
points(0:mxdeg, sigF, pch = 8, col = 2)

Fit.reg <- lm(Fit ~ poly(tt, degree = mxdeg))

##################################################################
#plotting and checking
#alpha <- max(0.9, 1 - 1/ndata)
#sigF <- qf(alpha, 2*(1:(mxdeg+1)), 2*(nord - 1:(mxdeg+1)))
#plot(0:mxdeg, F1, type = 'h', lwd = 3)
#points(0:mxdeg, sigF, pch = 8, col = 2)
#plot(0:mxdeg, F2, type = 'h', lwd = 3)
#points(0:mxdeg, sigF, pch = 8, col = 2)
#plot(0:mxdeg, FP, type = 'h', lwd = 3)
#points(0:mxdeg, sigF, pch = 8, col = 2)



#plot(t, FrA, type = 'l')
#lines(t,Fr1, col = 2)
#lines(t,FrC, col = 4)

plot(t, FMpoly, type = 'l')
lines(Fit, col = 2)
lines(FitC, col = 4)

sum( (FMpoly - Fit)^2 )
sum( (FMpoly - FitC)^2 )

plot(t, FrA, type = 'l')
lines(Fr1, col = 2)
lines(FrC, col = 4)
#lines(FrC + f, col = 4)
#PhFit <- cumsum(f+Fit)*dt*2*pi - (f+Fit)*dt*2*pi
#fit.sig <- cos(PhFit) + 1i * sin(PhFit)
#plot(Re(fit.sig), type = 'l', xlim = c(250,500))
#lines(Re(cd), col = 2)


#Z <- CDN * exp( 1i * PhFit)
#plot(Re(Z), type = 'l')
#ft <- spec.mtm(Z, nw = (nord+1)/2, k = nord, dpssIN = dw, Ftest = TRUE, plot = TRUE, deltat = 1)
#plot(ft, Ftest = TRUE)
#abline(h = qf(1 - 1/ndata, 2, 2*(nord-1)), lty = 2, col = 2)
#ft$freq[which.max(ft$mtm$Ftest)]