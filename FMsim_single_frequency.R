source("utility.R")

FBW <- (0:5) / 5
fbw.sim <- length(FBW)
num.sim <- 100


F1array <- array(0, dim = c(mxdeg+1, fbw.sim, num.sim))
F2array <- array(0, dim = c(mxdeg+1, fbw.sim, num.sim))
FParray <- array(0, dim = c(mxdeg+1, fbw.sim, num.sim))



F1countermatrix <- matrix(data = 0, nrow = fbw.sim, ncol = mxdeg+1)
F2countermatrix <- matrix(data = 0, nrow = fbw.sim, ncol = mxdeg+1)
FPcountermatrix <- matrix(data = 0, nrow = fbw.sim, ncol = mxdeg+1)


###PARAMETERS AND STUFF THAT IS THE SAME FOR EVERY ROUND OF THE SIMULATION
##########################################################################
#set up frequency modulation and generate data
ndata <- 20000 #N
nppb <- ndata  #use Dave's block function to set this if I do this over blocks. not currently used.
nord <- 10 #This is K = 2NW - 1
mxdeg <- 6 #max degree of polynomial modulation
#CONSIDER ADDING A VARIABLE P=mxdeg+1
js1 <- 69
js2 <- 5318008 #used for the seed

f <- 0.2 #center frequency
W <- (nord + 1)/(2*ndata) #bandwidth

stm11 <- 2.0/(ndata-1) #stuff used to generate modulating polynomial
Shift <- 0.0
dt <- 1
t <- 0:(ndata-1)        #time
tt <- (t + Shift) * stm11 - 1.0  # this runs from -1 to 1

snr <- 5										#signal to noise ratio
sigmat2 <- ndata/(2*snr*nord)  # looks like a variance, but of what?   - WHY THIS VARIANCE?
sigmat <- sqrt(sigmat2)

##########################################################################
#slepians, derivatives and special polynomials
DW <- dpss(ndata, nord, ndata * W) #slepian sequences
dw <- DW$v
ev <- DW$eigen

dwp <- dpssp7(dw,ndata,nord,W) #so-called time derivatives of the Slepians

HeSmx <- 0.9
hp <- dpsshp(dw,nppb,nord,mxdeg,HeSmx,norm = TRUE) 
Poly <- hp$Poly
Ukz <- hp$A
#the above bit generates the associated polynomials, G (in this case Hermite), and the matrix of 
#inner products, H, between Slepians and Hermite polynomials 

######################################################################################
FP <- numeric(mxdeg+1)
SEED <- NULL
start.time <- Sys.time()
for(i in 1:fbw.sim){
  fbw <- FBW[i]
  FPcounter <- rep(0,mxdeg+1) 
  F1counter <- rep(0,mxdeg+1)
  F2counter <- rep(0,mxdeg+1)
  
  for(j in 1:num.sim){
    seed <- sample(js1:js2, size = 1)
    set.seed(seed)
    SEED <- c(SEED, seed)
    FMpoly <- W * fbw * (1.0 - 2.0 * tt^2)
    Fmod <- f + FMpoly  # quadratic frequency modulation at center frequency f
    PhMod <- cumsum(Fmod)*dt*2*pi - Fmod*dt*2*pi #modulated phase
    cd <- cos(PhMod) + 1i * sin(PhMod) #data without noise
    
    
    CDN <- cd + rnorm(ndata, sd = sigmat) + rnorm(ndata, sd = sigmat) * 1i #data with noise
    CDN <- CDN * exp( -1i*2*pi*f*dt*t )
    
    yk <- t(dw) %*% CDN
    
    #############################################################################
    Fit1 <- dw %*% yk # "standard inverse" Z (Z tilde maybe) from the paper. Equation 3 I think
    Zdot <- dwp %*% yk #time derivative of Fit1.
    
    amp <- Mod(Fit1)
    PhD <- Arg(Fit1) #phase of standard inverse
    x1 <- Re(Fit1)
    y1 <- Im(Fit1)
    Xdot <- Re(Zdot)
    Ydot <- Im(Zdot)
    #FrA is the instantaneous frequency from section 3 of the paper
    FrA <- ( x1*Ydot - Xdot*y1 ) / ((2*pi)*amp^2) 
    
    #pc = c, Poly = G, yk = r = Y - H*c, Ukz = H in '09 paper
    pc <- t(Ukz) %*% yk
    r <- yk - Ukz %*% pc #this matches the paper, but the original code maybe doesn't
    
    Fcoef <- t(dw) %*% FrA
    
    SSQC <- colSums( Fcoef^2 )
    FPcoef <- t(Ukz) %*% Fcoef
    FRsd <- matrix(data=Fcoef, nrow = nord, ncol = mxdeg+1)
    FRsd <- FRsd - Ukz %*% diag( FPcoef[,1] )
    ssqr <- colSums(FRsd^2)
    SSQFK <- SSQC - ssqr
    F1 <- ( SSQFK/(1:(mxdeg+1)) ) / ( ssqr /(nord - 1:(mxdeg+1)) ) #why do I have two F stats?
    Rsd <- dw %*% FRsd
    SSQRtd <- colSums(Rsd^2)
    
    Fit <- 0 
    SSQFtd <- numeric(mxdeg+1)
    
    for(L in 0:mxdeg){#compute F statistic for each degree
      Fit <- Fit + Poly[,L+1] * FPcoef[L+1]
      SSQFtd[L+1] <- sum( Fit^2 )
    }
    
    F2 <- ( SSQFtd/(1:(mxdeg+1)) ) / ( SSQRtd/(nord - 1:(mxdeg+1)) )
    
    for(p in 1:(mxdeg+1)){
      FPn <- sum( Mod( Ukz[,1:p, drop = FALSE] %*% pc[1:p, drop = FALSE] )^2 ) / (2*p) 
      FPd <- sum( Mod( yk - Ukz[,1:p, drop = FALSE] %*% pc[1:p, drop = FALSE] )^2 ) / (2*(nord - p))
      FP[p] <- FPn / FPd
    }
    
    
    alpha <- 0.9
    sigF <- qf(alpha, 2*(1:(mxdeg+1)), 2*(nord - 1:(mxdeg+1)))
    
    for(m in 1:(mxdeg+1)){
      if(F1[m] >= sigF[m]){
        F1counter[m] <- F1counter[m] + 1
      }
      if(F2[m] >= sigF[m]){
        F2counter[m] <- F2counter[m] + 1
      }
      if(FP[m] >= sigF[m]){
        FPcounter[m] <- FPcounter[m] + 1
      }
    }
    
    FParray[,i,j] <- FP
    F1array[,i,j] <- F1
    F2array[,i,j] <- F2
    
  }
  F1countermatrix[i,] <- F1counter
  F2countermatrix[i,] <- F2counter
  FPcountermatrix[i,] <- FPcounter
}
end.time <- Sys.time()

total.time <- end.time - start.time

F1countermatrix
F2countermatrix
FPcountermatrix

