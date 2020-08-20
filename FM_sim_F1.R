source("utility.R")

FBW <- (0:5) / 5
fbw.sim <- length(FBW)
num.sim <- 100

F1counterlist <- NULL
F1DJTcounterlist <- NULL
Fcounterlist <- NULL

###PARAMETERS AND STUFF THAT IS THE SAME FOR EVERY ROUND OF THE SIMULATION
##########################################################################
#set up frequency modulation and generate data
ndata <- 2000 #N
nppb <- ndata  #use Dave's block function to set this if I do this over blocks. not currently used.
nord <- 10 #This is K = 2NW - 1
mxdeg <- 3 #max degree of polynomial modulation
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

alpha <- 1 - 1/ndata
snr <- 20										#signal to noise ratio
sigmat2 <- ndata/(2*snr*nord)  # looks like a variance, but of what?   - WHY THIS VARIANCE?
sigmat <- sqrt(sigmat2)

##########################################################################
#slepians, derivatives and special polynomials
DW <- dpss(ndata, nord, ndata * W) #slepian sequences
dw <- DW$v
ev <- DW$eigen

dwp <- dpssp7(dw,ndata,nord,W,ev) #so-called time derivatives of the Slepians

HeSmx <- 0.9
hp <- dpsshp(dw,nppb,nord,mxdeg,HeSmx,norm = TRUE) 
Poly <- hp$Poly
H <- hp$A
#the above bit generates the associated polynomials, G (in this case Hermite), and the matrix of 
#inner products, H, between Slepians and Hermite polynomials 

######################################################################################

nFFT <- 2^ceiling(log2(2*ndata))
freq <- seq(-0.5, 0.5, length.out = nFFT + 1)
freq <- freq[1:nFFT]

#selection of the band to test how accurate(?) the F statistics are
band.ind <- which(freq >= f-W & freq <= f+W) 
band <- freq[band.ind]

F1array <- array(0, dim = c(mxdeg+1, nFFT, fbw.sim, num.sim))
F1DJTarray <- array(0, dim = c(mxdeg+1, nFFT, fbw.sim, num.sim))
Farray <- array(0, dim = c(nFFT, fbw.sim, num.sim))
######################################################################################

start.time <- Sys.time()
for(i in 1:fbw.sim){
  fbw <- FBW[i] #set the fbw for the next 100 simulations
  #these counters keep track of the number of times the F statistic correctly has a maximum greater
  #than a critical F with level of significance 1 - 1/ndata (or 1/ndata. I can never remember which
  #it is) inside the band set above (f-W,f+W).
  F1counter <- rep(0,mxdeg+1)
  F1DJTcounter <- F1counter
  Fcounter <- 0
  for(j in 1:num.sim){
    
    
    seed <- sample(js1:js2, size = 1)
    set.seed(seed)
    
    
    FMpoly <- W * fbw * (1.0 - 2.0 * tt^2)
    PMpoly <- cumsum(FMpoly) - FMpoly[1]
    Fmod <- f + FMpoly  # quadratic frequency modulation at center frequency f
    PhMod <- 2*pi*f*dt*t + 2*pi*PMpoly #modulated phase
    cd <- cos(PhMod) + 1i * sin(PhMod) #data without noise
    
    #data with noise
    CDN <- cd + complex(real = rnorm(ndata, sd = sigmat), imaginary = rnorm(ndata, sd = sigmat) ) 
    
    
    ## NOW HERE I TAKE THE COMPLEX DEMODULATE WITH ALL FREQUENCIES. ZERO PAD A LOT
    ## AND CREATE MATRIX WITH COLUMNS BEING THE COMPLEX DEMODULATE DATA. CARRY OUT THE PROCESS ON ALL 
    ## COLUMNS OF THE DATA TO GET A MATRIX OF F STATISTICS - DEGREE BY FREQUENCY
    
    
    #taper and zero pad the data
    tapered <- dw * CDN 
    tapered <- rbind(tapered, matrix(0, nrow = nFFT-ndata, ncol = nord))
    #yk is a K x nFFT matrix of eigencoefficients. the frequencies go from -0.5 to 0.5
    yk <- mvfft(tapered)
    yk <- t(yk)
    yk <- cbind( yk[,(nFFT/2+2):nFFT], yk[,1:(nFFT/2+1)])
    
    
    #############################################################################
    Fit1 <- dw %*% yk # "standard inverse" Z (Z tilde maybe) from the paper. Equation 3 I think
    Zdot <- dwp %*% yk #time derivative of Fit1.
    
    amp <- Mod(Fit1)
    x1 <- Re(Fit1)
    y1 <- Im(Fit1)
    Xdot <- Re(Zdot)
    Ydot <- Im(Zdot)
    #FrA is the instantaneous frequency from section 3 of the paper
    FrA <- ( x1*Ydot - Xdot*y1 ) / ((2*pi)*amp^2) 
    
    Fcoef <- t(dw) %*% FrA
    FPcoef <- t(H) %*% Fcoef
    SSQC <- colSums( Fcoef^2 )
    #this is a Fortran subroutine that computes F1
    F1stuff <- FLoop2(nord,mxdeg,H,Fcoef,FPcoef,nFFT)
    F1DJTstuff <- FLoopDJT(nord, mxdeg, H, Fcoef, FPcoef, SSQC, nFFT)
    
    F1 <- matrix(data=F1stuff$F1, nrow = mxdeg+1, ncol = nFFT)
    F1DJT <- matrix(data=F1DJTstuff$F1mat, nrow = mxdeg+1, ncol = nFFT)
    #DON'T NEED TO COMPUTE THE FITS FOR NOW
    #Fit is the estimated modulating polynomial. It is actually computed using all the eigencoefficients
    #at each frequency. So we have a ndata x nFFT matrix
    #for(L in 1:(mxdeg+1)){
    #  Fit <- Poly[,1:L,drop=FALSE] %*% FPcoef[1:L,,drop=FALSE]
    #}
    
    ##############################################################
    
    #should probably compute this first and maybe use adaptively weighted eigencoefficients
    # and also use the usual harmonic F test for the degree 0 one in F1
    spec <- spec.mtm(CDN,  Ftest = TRUE, deltat = 1, plot = FALSE, dpssIN = dw, nFFT = nFFT)
    
    HarmF <- spec$mtm$Ftest #this is the regular Harmonic F statistic
    
    #update the counters across each degree. probably a cleaner way to write the if statements.
    for(m in 1:(mxdeg+1)){
      if(sum(which.max(F1[m,]) == band.ind) == 1 & max(F1[m,]) >= qf(alpha,1,nord-m)){
        F1counter[m] <- F1counter[m] + 1
      }
      if(sum(which.max(F1DJT[m,]) == band.ind) == 1 & max(F1DJT[m,]) >= qf(alpha,2,2*(nord-m))){
        F1DJTcounter[m] <- F1DJTcounter[m] + 1
      }
      
    }
    
    if(sum(which.max(HarmF) == band.ind) == 1 & max(HarmF) >= qf(alpha,2,2*(nord-1))){
      Fcounter <- Fcounter + 1
    }
    
    #this stores all of the F statistics in a vector to be plugged into an array
    F1array[,,i,j] <- F1
    F1DJTarray[,,i,j] <- F1DJT
    Farray[,i,j] <- HarmF
    
  }
  #store the counters from the 100 simulations at each level of FM in a vector
  F1counterlist <- c(F1counterlist, F1counter)
  F1DJTcounterlist <- c(F1DJTcounterlist, F1DJTcounter)
  Fcounterlist <- c(Fcounterlist, Fcounter)
}
end.time <- Sys.time()

total.time <- end.time - start.time

#this puts the counters in a (hopefully) easy to read matrix (level of FM by degree) that shows
#how the performance of these things changes with increased modulation and different degrees
F1countermatrix <- matrix(data = F1counterlist, nrow = fbw.sim, ncol = mxdeg+1, byrow = TRUE)
F1DJTcountermatrix <- matrix(data = F1DJTcounterlist, nrow = fbw.sim, ncol = mxdeg+1, byrow = TRUE)
#the next bit will sum all the F statistics over the band (f-W,f+W) to compare total power (abuse
#of terminology?) in the band. Maybe there is a better way to do it than an apply within for loops
F1band <- F1array[,band.ind,,]
F1DJTband <- F1DJTarray[,band.ind,,]
Fband <- Farray[band.ind,,]

F1bandtotal <- array(0, dim = c(mxdeg+1, fbw.sim, num.sim))
F1DJTbandtotal <- array(0, dim = c(mxdeg+1, fbw.sim, num.sim))
Fbandtotal <- matrix(0, nrow = fbw.sim, ncol = num.sim)

for(i in 1:(mxdeg+1)){
  for(j in 1:fbw.sim){
    F1bandtotal[i,j,] <- apply(F1band[i,,j,], MARGIN = 2, FUN = sum)
    F1DJTbandtotal[i,j,] <- apply(F1DJTband[i,,j,], MARGIN = 2, FUN = sum)
  }
}

for(j in 1:fbw.sim){
  Fbandtotal[j,] <- apply(Fband[,j,], MARGIN = 2, FUN = sum)
}

#the next part computes the 2nd, 50th and 98th quantiles of the simulated F statistics for each
#degree, frequency and level of FM
F1sort <- array(0, dim = c(mxdeg+1, nFFT, fbw.sim, num.sim))
F1DJTsort <- array(0, dim = c(mxdeg+1, nFFT, fbw.sim, num.sim))

probs <- c(0.02, 0.5, 0.98)
F1quant <- array(0, dim = c(mxdeg+1, nFFT, fbw.sim, length(probs)))
F1DJTquant <- array(0, dim = c(mxdeg+1, nFFT, fbw.sim, length(probs)))

for(i in 1:(mxdeg+1)){
  for(j in 1:nFFT){
    for(k in 1:fbw.sim){
      F1sort[i,j,k,] <- sort(F1array[i,j,k,])
      F1quant[i,j,k,] <- quantile(F1sort[i,j,k,], probs = probs)
      
      F1DJTsort[i,j,k,] <- sort(F1DJTarray[i,j,k,])
      F1DJTquant[i,j,k,] <- quantile(F1DJTsort[i,j,k,], probs = probs)
    }
  }
}


#the next bit averages the F statistics over the simulations for each level of FM, degree and
#frequency
F1ave <- array(data=0, dim = c(fbw.sim,mxdeg+1,nFFT))
F1DJTave <- array(data=0, dim = c(fbw.sim, mxdeg+1, nFFT))
for(i in 1:fbw.sim){
  for(j in 1:(mxdeg+1)){
    F1ave[i,j,] <- apply(F1array[j,,i,], MARGIN = 1, FUN = mean)
    F1DJTave[i,j,] <- apply(F1DJTarray[j,,i,], MARGIN = 1, FUN = mean)
  }
}



save.image(file = "SimulationDataF1compare.RData")