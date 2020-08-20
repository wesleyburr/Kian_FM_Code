source("utility.R")

ndata <- 2000 #N
nppb <- ndata
nord <- 10 #This is K
mxdeg <- 6 #max degree of polynomial modulation

js1 <- 131
js2 <- 17107
seed <- sample(js1:js2, size = 1)
set.seed(seed)

RtoD <- 180/pi           #radians to degrees
DtoR <- pi/180
FrScl <- 1/360

f <- 0.2 #center frequency. can't seem to find a way to incorporate this

TNW <- nord+1
W <- TNW/(2*ndata)
fbw <- 0.8
stm11 <- 2.0/(ndata-1)
Shift <- 0.0

dt <- 1.0
t <- 0:(ndata-1)        #time
cent <- (1+ndata)/2.0
scl <- 1.0/cent
tc <- (t - cent)*scl    #centered time indices. not used here
tt <- (t + Shift) * stm11 - 1.0  # this runs from -1 to 1
FMpoly <- W * fbw * (1.0 - 2.0 * tt^2)
Fmod <- f + FMpoly # quadratic frequency modulation
PhMod <- cumsum(Fmod)*dt*2*pi
cd <- cos(PhMod) #+ 1i * sin(PhMod)
#cd <- cd * exp( -1i*2*pi*f*t )
Pref <- RtoD * cumsum(FMpoly)*dt*2*pi#Pref is the prespecified phase modultation

DW <- dpss(ndata, nord, ndata * W)
dw <- DW$v
ev <- DW$eigen
#ev <- dw$eigen       #I feel like the eigenvalues are never used for anything
#Ukz0 <- apply(dw, MAR = 2, FUN = sum) 
#Don't think I need this for anything
############################################
#start_time <- Sys.time()
dwp <- dpssp7(dw,nppb,nord,W,ev)
#runtime <- Sys.time() - start_time
#runtime


##########################################################################
norm <- TRUE
HeSmx <- 0.90

hp <- dpsshp(dw,nppb,nord,mxdeg,HeSmx,norm)

Poly <- hp$Poly
Ukz <- hp$A

#these are a check that the columns of Ukz are orthonormal
HtH <- t(Ukz) %*% Ukz
norms <- colSums(Ukz^2)
#######################################################################
snr <- 5.0 											#signal to noise ratio
sigmat2 <- ndata/(2*snr*nord)  # looks like a variance, but of what?   - WHY THIS VARIANCE?
sigmat <- sqrt(sigmat2)
CDN <- cd + rnorm(ndata, sd = sigmat)# + rnorm(ndata, sd = sigmat) * 1i
CDN <- CDN * exp( -1i*2*pi*f*t )
pwrn <- Mod(CDN - cd)^2     #probably don't need this
###################
yk <- t(dw) %*% CDN

Fit1 <- dw %*% yk #standard inverse
Zdot <- dwp %*% yk #time derivative of Fit1.


amp <- Mod(Fit1)
PhD <- RtoD*Arg(Fit1) #phase of standard inverse
x1 <- Re(Fit1)
y1 <- Im(Fit1)
Xdot <- Re(Zdot)
Ydot <- Im(Zdot)
FrA <- ( x1*Ydot - Xdot*y1 ) / ((2*pi)*amp^2) #derivative of phase? -> instantaneous frequency
                                              #this is for the standard inverse fit
Fcoef <- t(dw) %*% FrA
FPcoef <- t(Ukz) %*% Fcoef
Fit <- Poly %*% FPcoef

plot(FMpoly, type = 'l')
lines(FrA, col = 2)
lines(Fit, col = 4)


#I may not even need this curvature stuff
FaC2 <- curvature(FrA)

Fcoef <- t(dw) %*% FrA
Fcoef <- Fcoef[,1]

####

pc <- numeric(mxdeg+1)
for (jd in 0:mxdeg){
  pc[jd+1] <- sum( yk * Ukz[,jd+1]) 
  yk <- yk - Ukz[,jd+1]*pc[jd+1]
}

#pc = c, Poly = G, yk = r = Y - H*c, Ukz = H in '09 paper

FitP <- Poly %*% pc  #line 11
FitC <- FitP + dw %*% yk # G*c + V*r? I think at this point, yk is yk - H*c

A1 <- Mod(FitP)
P1  <- RtoD*Arg(FitP)
A2  <- Mod(FitC)
P2  <- RtoD*Arg(FitC)
x1  <-  Re(FitP)
y1  <- Im(FitP)
x2  <-  Re(FitC)
y2  <- Im(FitC)

#bunch of plotting and cm4pp stuff if I want. otherwise, most of this stuff isn't used


phi <- c(1/2,0,-1/2)
Fr1 <- FrScl * filter(PhD, phi, method = 'convolution')
FrP <- FrScl * filter(P1, phi, method = 'convolution')
FrC <- FrScl * filter(P2, phi, method='convolution')

Zhat <- t(dw) %*% FrC #this is just verification that the eigencoefficients are in fact the same
plot(yk)
points(Zhat, pch = 2, col = 2)
#more plotting stuff here
Frsd <- Fcoef
SSQC <- sum( Fcoef^2 )

FPcoef <- t(Ukz) %*% Fcoef


Frsd <- matrix(data=Frsd, nrow = nord, ncol = mxdeg+1)
FRsd <- Frsd - Ukz %*% diag( FPcoef[,1] )
ssqr <- colSums(FRsd^2)
SSQFK <- SSQC - ssqr
F1 <- ( SSQFK/(1:(mxdeg+1)) ) / ( ssqr /(nord - 1:(mxdeg+1)) ) 
Rsd <- dw %*% FRsd
SSQRtd <- colSums(Rsd^2)
Rcurv2 <- apply(Rsd, MARGIN = 2, FUN = curvature)
#Fit <- Poly %*% FPcoef # G times Px1 vector
#this does the same thing, but with the F test and stuff

Fit <- 0 
SSQFtd <- numeric(mxdeg+1)
Fcurv2 <- numeric(mxdeg+1)
AERR2 <- numeric(mxdeg+1)

for(L in 0:mxdeg){#compute F statistic for each degree
  Fit <- Fit + Poly[,L+1] * FPcoef[L+1]
  SSQFtd[L+1] <- sum( Fit^2 )
  AERR2[L+1] <- sum( (Fit[2:(ndata-1)] - Fmod[2:(ndata-1)])^2 )
  Fcurv2[L+1] <- curvature(Fit)
}
F2 <- ( SSQFtd/(1:(mxdeg+1)) ) / ( SSQRtd/(nord - 1:(mxdeg+1)) )

Fit.reg <- lm(Fit ~ poly(tt, degree = 2))


#more cm4pp and plotting
plot(t, Fmod, type = 'l')
plot(t, Pref, type='l')
plot(t, Re(cd), type='l')
lines(t, Im(cd), type='l', col=2)

plot(t, Re(CDN), type='l')
lines(t, Im(CDN), col=2)

#bunch of plotting here. probably don't need any of it
plot(t,x1, type='l')
lines(t,y1, col=2)
plot(t,Xdot, type='l')
lines(t,Ydot, col=2)
plot(t,amp,type='l')
plot(t,PhD,type='l')
lines(t,Pref,col=2)
plot(t, FrA, type='l')

plot(0:mxdeg, F1, type = 'h', xlab = 'Degree')
plot(0:mxdeg, F2, type = 'h', xlab = 'Degree')

plot(FrA, type= 'l')
lines(Fr1, col=2)
threshhold <- 0.5
clicksFr1 <- which(abs(Fr1) > threshhold)
abline(v=clicksFr1, col = 4) #clicks in phase

plot(FrA, type = 'l')
lines(FrC, col =2)
clicksFrC <- which(abs(FrC) > threshhold)
abline(v=clicksFrC, col = 4)

plot(Fmod, type = 'l')
lines(Fit + f, col =2)



#Z <- CDN * exp( -1i * 2 * pi * Fit) #* 0:(ndata-1))
#plot(Re(Z), type = 'l')
#ft <- spec.mtm(Z, nw = TNW/2, k = nord, dpssIN = dw, Ftest = TRUE, plot = FALSE, deltat = 1)
#plot(ft, Ftest = TRUE)
#abline(h = qf(1 - 1/ndata, 2, 2*(nord-1)), lty = 2, col = 2)
