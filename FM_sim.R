source("utility.R")

FBW <- (0:5) / 5
fbw.sim <- length(FBW)
num.sim <- 100

FPlist <- NULL
F1list <- NULL
F2list <- NULL
HarmFlist <- NULL


FPcounterlist <- NULL
F1counterlist <- NULL
F2counterlist <- NULL
Fcounterlist <- NULL

###PARAMETERS AND STUFF THAT IS THE SAME FOR EVERY ROUND OF THE SIMULATION
##########################################################################
#set up frequency modulation and generate data
ndata <- 2000 #N
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

dwp <- dpssp7(dw,ndata,nord,W,ev) #so-called time derivatives of the Slepians

HeSmx <- 0.9
hp <- dpsshp(dw,nppb,nord,mxdeg,HeSmx,norm = TRUE) 
Poly <- hp$Poly
Ukz <- hp$A
#the above bit generates the associated polynomials, G (in this case Hermite), and the matrix of 
#inner products, H, between Slepians and Hermite polynomials 

######################################################################################

nFFT <- 2^ceiling(log2(2*ndata))
freq <- seq(-0.5, 0.5, length.out = nFFT + 1)
freq <- freq[1:nFFT]

#selection of the band to test how accurate(?) the F statistics are
band.ind <- which(freq >= f-W & freq <= f+W) 
band <- freq[band.ind]

#degrees of freedom to be used in the F statistics. df1 = numerator, df2 = denominator
df1 <- matrix(data = 1:(mxdeg+1), nrow = mxdeg+1, ncol = nFFT, byrow = FALSE)
df2 <- matrix(data = nord - 1:(mxdeg+1), nrow = mxdeg+1, ncol = nFFT, byrow = FALSE)

#set some 'empty' matrices to be filled later
SSQFtd <- matrix(data=0,nrow = mxdeg+1, ncol = nFFT) #numerator ssq for F2
FP <- matrix(data=0, nrow = mxdeg+1, ncol = nFFT) #this will be the FP from the paper
######################################################################################

start.time <- Sys.time()
for(i in 1:fbw.sim){
  fbw <- FBW[i] #set the fbw for the next 100 simulations
  #these counters keep track of the number of times the F statistic correctly has a maximum greater
  #than a critical F with level of significance 1 - 1/ndata (or 1/ndata. I can never remember which
  #it is) inside the band set above (f-W,f+W).
  FPcounter <- rep(0,mxdeg+1) 
  F1counter <- rep(0,mxdeg+1)
  F2counter <- rep(0,mxdeg+1)
  Fcounter <- 0
  for(j in 1:num.sim){

    
    seed <- sample(js1:js2, size = 1)
    set.seed(seed)
    
   
    FMpoly <- W * fbw * (1.0 - 2.0 * tt^2)
    Fmod <- f + FMpoly  # quadratic frequency modulation at center frequency f
    PhMod <- 2*pi*f*dt*t + cumsum(FMpoly)*2*pi*dt #modulated phase
    cd <- cos(PhMod) + 1i * sin(PhMod) #data without noise
    
    
    CDN <- cd + rnorm(ndata, sd = sigmat) + rnorm(ndata, sd = sigmat) * 1i #data with noise
    
    
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
    
    #this is a pair of Fortran subroutines that compute F1 and SSQRtd (the denominator ssq of F2)
    # DJT seemed to think the 'td' here is for time domain, which makes sense.
    output <- Floop(Fcoef, nord, mxdeg, Ukz, FPcoef, dw, SSQC, nFFT, ndata)
    
    F1 <- matrix(data=output$F1, nrow = mxdeg+1, ncol = nFFT)
    SSQRtd <- matrix(data=output$SSQRtd, nrow = mxdeg+1, ncol = nFFT)
    
    #Fit is the estimated modulating polynomial. It is actually computed using all the eigencoefficients
    #at each frequency. So we have a ndata x nFFT matrix
    for(L in 1:(mxdeg+1)){
      Fit <- Poly[,1:L,drop=FALSE] %*% FPcoef[1:L,,drop=FALSE]
      SSQFtd[L,] <- colSums( Fit^2 )
    }
    
    
    
    F2 <- (SSQFtd / df1) / (SSQRtd / df2)
    
    ##############################################################
    
    #IN FP TRY REPLACING THE PARTS DEPENDING ON EIGENCOEFFICIENTS OR DATA BY THE FrA AS IN THE
    #OTHER TWO F STATISTICS - IGNORE THIS FOR NOW: POSSIBLY DUMB OR SHORTSIGHTED
    
    #Fp from the paper. Doesn't seem to work very well, so it may be coded incorrectly. Maybe
    #it just isn't a great statistic.
    for(p in 1:(mxdeg+1)){
      FP[p,] <- (colSums(Mod(Ukz[,1:p,drop=FALSE] %*% pc[1:p,,drop=FALSE])^2) / (2*p) ) /
        (colSums(Mod(yk - Ukz[,1:p,drop=FALSE] %*% pc[1:p, ,drop = FALSE])^2) / (2*(nord - p)) )
    }
    
    spec <- spec.mtm(CDN,  Ftest = TRUE, deltat = 1, plot = FALSE, dpssIN = DW, nFFT = nFFT)
    
    HarmF <- spec$mtm$Ftest #this is the regular Harmonic F statistic
    
    #update the counters across each degree. probably a cleaner way to write the if statements.
    for(m in 1:(mxdeg+1)){
      if(sum(which.max(FP[m,]) == band.ind) == 1 & max(FP[m,]) >= qf(1-1/ndata,2*m,2*(nord-m))){
        FPcounter[m] <- FPcounter[m] + 1
      }
      if(sum(which.max(F1[m,]) == band.ind) == 1 & max(F1[m,]) >= qf(1-1/ndata,2*m,2*(nord-m))){
        F1counter[m] <- F1counter[m] + 1
      }
      if(sum(which.max(F2[m,]) == band.ind) == 1 & max(F1[m,]) >= qf(1-1/ndata,2*m,2*(nord-m))){
        F2counter[m] <- F2counter[m] + 1
      }
    }
    
    if(sum(which.max(HarmF) == band.ind) == 1 & max(HarmF[m,]) >= qf(1-1/ndata,2,2*(nord-1))){
      Fcounter <- Fcounter + 1
    }
    
    #this stores all of the F statistics in a vector to be plugged into an array (and a dataframe
    #for some crappy plotting I was doing) later.
    FPlist <- c(FPlist, FP)
    F1list <- c(F1list, F1)
    F2list <- c(F2list, F2)
    HarmFlist <- c(HarmFlist, HarmF)
    
  }
  #store the counters from the 100 simulations at each level of FM in a vector
  FPcounterlist <- c(FPcounterlist, FPcounter)
  F1counterlist <- c(F1counterlist, F1counter)
  F2counterlist <- c(F2counterlist, F2counter)
  Fcounterlist <- c(Fcounterlist, Fcounter)
}
end.time <- Sys.time()

total.time <- end.time - start.time

#this puts the counters in a (hopefully) easy to read matrix (level of FM by degree) that shows
#how the performance of these things changes with increased modulation and different degrees
F1countermatrix <- matrix(data = F1counterlist, nrow = fbw.sim, ncol = mxdeg+1, byrow = TRUE)
F2countermatrix <- matrix(data = F2counterlist, nrow = fbw.sim, ncol = mxdeg+1, byrow = TRUE)
FPcountermatrix <- matrix(data = FPcounterlist, nrow = fbw.sim, ncol = mxdeg+1, byrow = TRUE)

#these are just the arrays mentioned above
F1array <- array(F1list, dim = c(mxdeg+1, nFFT, fbw.sim, num.sim))
F2array <- array(F2list, dim = c(mxdeg+1, nFFT, fbw.sim, num.sim))
FParray <- array(FPlist, dim = c(mxdeg+1, nFFT, fbw.sim, num.sim))
Farray <- array(HarmFlist, dim = c(nFFT, fbw.sim, num.sim))

#the next bit will sum all the F statistics over the band (f-W,f+W) to compare total power (abuse
#of terminology?) in the band. Maybe there is a better way to do it than an apply within for loops
F1band <- F1array[,band.ind,,]
F2band <- F2array[,band.ind,,]
FPband <- FParray[,band.ind,,]
Fband <- Farray[band.ind,,]

F1bandtotal <- array(0, dim = c(mxdeg+1, fbw.sim, num.sim))
F2bandtotal <- F1bandtotal
FPbandtotal <- F1bandtotal
Fbandtotal <- matrix(0, nrow = fbw.sim, ncol = num.sim)

for(i in 1:(mxdeg+1)){
  for(j in 1:fbw.sim){
    F1bandtotal[i,j,] <- apply(F1band[i,,j,], MARGIN = 2, FUN = sum)
    F2bandtotal[i,j,] <- apply(F2band[i,,j,], MARGIN = 2, FUN = sum)
    FPbandtotal[i,j,] <- apply(FPband[i,,j,], MARGIN = 2, FUN = sum)
  }
}

for(j in 1:fbw.sim){
  Fbandtotal[j,] <- apply(Fband[,j,], MARGIN = 2, FUN = sum)
}

#the next part computes the 2nd, 50th and 98th quantiles of the simulated F statistics for each
#degree, frequency and level of FM
F1sort <- array(0, dim = c(mxdeg+1, nFFT, fbw.sim, num.sim))
F2sort <- F1sort
FPsort <- F1sort

probs <- c(0.02, 0.5, 0.98)
F1quant <- array(0, dim = c(mxdeg+1, nFFT, fbw.sim, length(probs)))
F2quant <- F1quant
FPquant <- F1quant

start.time2 <- Sys.time()
for(i in 1:(mxdeg+1)){
  for(j in 1:nFFT){
    for(k in 1:fbw.sim){
      F1sort[i,j,k,] <- sort(F1array[i,j,k,])
      F1quant[i,j,k,] <- quantile(F1sort[i,j,k,], probs = probs)
      
      F2sort[i,j,k,] <- sort(F2array[i,j,k,])
      F2quant[i,j,k,] <- quantile(F2sort[i,j,k,], probs = probs)
      
      FPsort[i,j,k,] <- sort(FParray[i,j,k,])
      FPquant[i,j,k,] <- quantile(FPsort[i,j,k,], probs = probs)
    }
  }
}
end.time2 <- Sys.time()

total.time2 <- end.time2 - start.time2


#the next bit averages the F statistics over the simulations for each level of FM, degree and
#frequency
F1ave <- array(data=0, dim = c(fbw.sim,mxdeg+1,nFFT))
for(i in 1:fbw.sim){
  for(j in 1:(mxdeg+1)){
    F1ave[i,j,] <- apply(F1array[j,,i,], MARGIN = 1, FUN = mean)
  }
}

F2ave <- array(data=0, dim = c(fbw.sim,mxdeg+1,nFFT))
for(i in 1:fbw.sim){
  for(j in 1:(mxdeg+1)){
    F2ave[i,j,] <- apply(F2array[j,,i,], MARGIN = 1, FUN = mean)
  }
}

FPave <- array(data=0, dim = c(fbw.sim,mxdeg+1,nFFT))
for(i in 1:fbw.sim){
  for(j in 1:(mxdeg+1)){
    FPave[i,j,] <- apply(FParray[j,,i,], MARGIN = 1, FUN = mean)
  }
}

save.image(file = "SimulationData2.RData")

#the rest of this stuff is for dataframes I was using to make terrible plots. Ignore for now.
frequency <- rep(freq, (mxdeg+1)*fbw.sim*num.sim)
degree <- rep(0:mxdeg, nFFT)
degree <- sort(degree)
degree <- rep(degree, fbw.sim*num.sim)
sim <- rep(1:num.sim, (mxdeg+1)*fbw.sim*nFFT)
sim <- sort(sim)
fbW <- rep(FBW, (mxdeg+1)*nFFT*num.sim)
fbW <- sort(fbW)

F1data <- cbind(fbW, sim, degree, frequency, F1list)
F1data <- as.data.frame(F1data)
colnames(F1data)[5] <- "F1"

F2data <- cbind(fbW, sim, degree, frequency, F2list)
F2data <- as.data.frame(F2data)
colnames(F2data)[5] <- "F2"

FPdata <- cbind(fbW, sim, degree, frequency, FPlist)
FPdata <- as.data.frame(FPdata)
colnames(FPdata)[5] <- "FP"

wireframe(F1 ~ degree * frequency, data = F1data[1:( (mxdeg+1)*nFFT ),])
wireframe(F2 ~ degree * frequency, data = F2data[1:( (mxdeg+1)*nFFT ),])
wireframe(FP ~ degree * frequency, data = FPdata[1:( (mxdeg+1)*nFFT ),])

frequency <- rep(freq, (mxdeg+1)*fbw.sim)
degree <- rep(0:mxdeg, nFFT)
degree <- sort(degree)
degree <- rep(degree, fbw.sim)
fbW <- rep(FBW, (mxdeg+1)*nFFT)
fbW <- sort(fbW)

F1ave.data <- cbind(fbW, degree, frequency, unlist(F1ave))
F1ave.data <- as.data.frame(F1ave.data)
colnames(F1ave.data)[4] <- "F1"

wireframe(F1 ~ degree*frequency, data = F1ave.data[which(F1ave.data$fbW == 0.4),])
