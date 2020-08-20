library(parallel)
library(multitaper)
library(VGAM)
dyn.load("F95lib")

dpssap <- function(V, maxdeg, alpha=0.75) {
  
  # Sanity checks
  stopifnot(is.matrix(V), is.numeric(maxdeg), maxdeg>=0)
  N <- length(V[, 1])
  K <- length(V[1, ])
  P <- maxdeg + 1
  timeArr <- 1:N
  
  R <- matrix(data=0, nrow=N, ncol=P)
  U <- matrix(data=0, nrow=K, ncol=P)
  
  # Setup centered time index
  midTime <- (1+N) / 2
  scl <- 2/(N-1)
  timeArrC <- (timeArr - midTime) * scl
  
  # Start with Gegenbauer polynomials; convergence is faster
  #alpha <- 0.75
  R[, 1] <- 1.0
  if(maxdeg > 0) {
    R[, 2] <- 2 * alpha * timeArrC
    if(maxdeg > 1) {
      for(j in 2:maxdeg) {
        A1 <- 2 * ( (j-1) + alpha ) / j
        A2 <- ( (j-2) + 2 * alpha ) / j
        
        R[, (j+1)] <- A1 * timeArrC * R[, j] - A2 * R[, (j-1)]
      } # end of loop on higher orders
    } # end of maxdeg > 1
  } # end of maxdeg > 0
  
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
  
  Hn <- colSums(R^2)
  return(list(U=U,R=R,Hn=Hn))
}

gegenbauer <- function(maxdeg, N, alpha = 0.75){
  
  P <- maxdeg + 1
  timeArr <- 1:N
  
  R <- matrix(data=0, nrow=N, ncol=P)
  # Setup centered time index
  midTime <- (1+N) / 2
  scl <- 2/(N-1)
  timeArrC <- (timeArr - midTime) * scl
  
  # Start with Gegenbauer polynomials; convergence is faster
  R[, 1] <- 1.0
  if(maxdeg > 0) {
    R[, 2] <- 2 * alpha * timeArrC
    if(maxdeg > 1) {
      for(j in 2:maxdeg) {
        A1 <- 2 * ( (j-1) + alpha ) / j
        A2 <- ( (j-2) + 2 * alpha ) / j
        
        R[, (j+1)] <- A1 * timeArrC * R[, j] - A2 * R[, (j-1)]
      } # end of loop on higher orders
    } # end of maxdeg > 1
  } # end of maxdeg > 0
  return(R=R)
}


topM <- function(x, M=1){
  return( order(x, decreasing=TRUE)[1:M] )
}


reconstruct <- function(time.series, num = 1, NW = 4, K = 7, nFFT = 2^ceiling(log2(10*length(time.series)))){
  
  x <- time.series - mean(time.series) 
  N <- length(x)
  v <- dpss(N,K,NW, returnEigenvalues = FALSE)$v
  
  xbar <- v %*% t(v) %*% x
  temp.x <- x - xbar
  
  tapered <- rbind(v * temp.x[,1], matrix(rep(0, (nFFT - N)*K), ncol = K))
  y <- mvfft(tapered)
  
  spec <- spec.mtm(temp.x, nw=NW, k=K, nFFT = nFFT, dpssIN = v, Ftest = TRUE, 
                   deltat = 1,  plot=FALSE, centre='none')
  
  f.mesh <- seq(0,0.5, by = 2/N)
  L <- length(f.mesh)
  F.max <- NULL
  F.indices <- NULL
  for(l in 1:(L-1)){
    temp.spec <- dropFreqs(spec,f.mesh[l],f.mesh[l+1])
    F.max <- c(F.max, max(temp.spec$mtm$Ftest))
    F.indices <- c(F.indices, which(spec$mtm$Ftest == F.max[l]) )
  }
  peak.freqs <- spec$freq[F.indices]
  topF <- topM(spec$mtm$Ftest[F.indices], num)
  f <- peak.freqs[topF] 
  #need a way to reject frequencies selected that are too close together or are not at peaks,
  #i.e. two frequencies that are on either side of a peak that falls near the end point of 
  #one of the frequency bins
  
  
  h <- nFFT * f + 1
  Yf <- y[h, , drop = FALSE]
  mu <- spec$mtm$cmv[h]
  Z <- v %*% t(Yf)
  Z.tilde <- matrix(data=0, nrow = N, ncol = num)
  
  for(j in 1:num){
    Z.tilde[,j] <- Z[,j] * exp(1i * 2 * pi * f[j] * 0:(N-1) )
  }
  
  Z.hat <- apply(Z.tilde, MARGIN = 1, FUN = sum)
  Z.hat <- Z.hat + (1 + 1i) * ( xbar + mean(time.series) )
  
  #shift the imaginary part up in case it is used for something 
  plot(time.series , type='l')
  lines(Re(Z.hat), col=2)
  #lines(Im(Z.hat),col=4)
  
  data <- list( Z.hat = Z.hat, Z.tilde = Z.tilde, freq = f)
  return(data)
}

block <- function(n, blockSize, overlap){ 
  increment <- ceiling(blockSize * (1-overlap)) 
  sectIdx <- seq(1, n-blockSize+1, by=increment) 
  numSect <- length(sectIdx) 
  list(startIdx = sectIdx, increment = increment, nBlock = numSect, blockSize = blockSize)
}

curvature <- function(x){
  curv <- filter(x,c(1,-2,1), method = 'convolution')
  curv <- sum( curv[2:(ndata-1)]^2 )
  return(curv)
}

dpsshp <- function(dw,nppb,nord,mxdeg,HeScl,norm){
  
  
  #cntr <- (1+nppb)/2
  #scl <- HeScl*2/(nppb-1)
  #xx <- (1:nppb - cntr) * scl
  xx <- HeScl * seq(-1,1, length.out = nppb)
  
  #                   Start with Hermite Polynomials for reasonable match
  Poly <- HePolA(xx, mxdeg)
  #PolyP <- Poly$HeP
  #Poly <- Poly$He
  
  
  if(norm == TRUE){ #normalizes the columns of Poly
    sum <- colSums( Poly^2 )
    scl <- 1/sqrt(sum)
    scl <- diag(scl)
    Poly <- Poly %*% scl
    #PolyP <- PolyP %*% scl
  }
  #         Weighted Inner Products of Polynomials and Slepian Sequences
  A <- matrix(data=0, nrow=nord, ncol=mxdeg+1)
  for(L in 0:mxdeg){
    Kmin <- L %% 2 + 1
    A[seq(Kmin, nord, 2), L+1] <- t(dw[,seq(Kmin, nord, 2)]) %*% Poly[,L+1]
  }
  #                         Degree 0, 1
  m <- min(1,mxdeg) + 1
  sum <- colSums(A[,1:m]^2)
  scl <- 1/sqrt(sum)
  scl <- diag(scl)
  A[,1:m] <- A[,1:m] %*% scl
  Poly[,1:m] <- Poly[,1:m] %*% scl
  
  
  
  
  if(mxdeg >= 2){#        Do Gram-Schmidt
    for(L in 2:mxdeg){
      Kmin <- L %% 2 
      for(j in seq(Kmin, L-2, 2)){
        scl <- sum( A[,L+1] * A[,j+1] ) / sum( A[,j+1]^2 )
        A[,L+1] <- A[,L+1] - scl * A[,j+1]
        Poly[,L+1] <- Poly[,L+1] - scl * Poly[,j+1]
      }
    }                    #Orthonormalize
    
    for(L in 3:(mxdeg+1)){
      scl <- 1/ sqrt( sum( A[,L]^2 ) )
      A[,L] <- A[,L] * scl
      Poly[,L] <- Poly[,L] * scl
    }
  }
  Hn <- colSums(Poly^2)
  output <- list(A = A, Poly = Poly, Hn = Hn)#, PolyP = PolyP)
  return(output)
}

HePolA <- function(x, mxdeg){
  #compute probabilists' Hermite polynomials
  ndata <- length(x)
  He <- matrix(data=0, nrow = ndata, ncol = mxdeg + 1)
  
  
  He[,1] <- 1
  if(mxdeg > 0){
    He[,2] <- x
    if(mxdeg > 1){
      for(k in 1:(mxdeg-1)){
        He[,k+2] <- x*He[,k+1] - k*He[,k]
      }
    }
  }
  return(He)
}
comp.evals <- function(k,efn,nhalf,ndata,x,ev){
  ev8 <- 0
  lmx <- which.max(Mod(efn[1:nhalf,k]))
  j <- Mod(1:ndata - lmx) + 1
  ev8 <- sum( efn[,k]*x[j] ) / efn[lmx,k]
  ev[k] <- ev8
  return(ev[k])
}
looping <- function(k,efn,efnp,ndata,y,nefn,ev){
  output <- .Fortran("loop",ndata=as.integer(ndata),y=as.double(y),efn=as.double(efn),
                     efnp=as.double(efnp),k=as.integer(k),nord=as.integer(nefn))
  efnp <- matrix(output$efnp, nrow = ndata, ncol = nefn)
  efnp[,k] <- efnp[,k] / ev[k]
  return(efnp[,k])
}

dpssp <- function(dw,ndata,nord,W,ev){
  output <- .Fortran("dpssp7", dw = as.double(dw),
                     ndata = as.integer(ndata),
                     nord = as.integer(nord),
                     W = as.double(W),
                     ev = as.double(ev),
                     dwp = double(ndata * nord))
}

#NOTE: CAN REMOVE THE PART ABOUT COMPUTING EIGENVALUES UNLESS THEY BECOME NECESSARY FOR SOMETHING. 
#THEY ARE COMPUTED IN THE DPSS FUNCTION AND CAN EASILY BE TAKEN FROM THERE.
dpssp7 <- function(efn,ndata,nefn,W,ev){
  x <- c(2*W, sin( 2*pi*W*(1:ndata) ) / (pi*(1:ndata)) )
  y <- c(0, ( 2*W*cos( 2*pi*W*(1:ndata) ) - x[2:(ndata+1)])  / (1:ndata) )
  
  
  nhalf <- 1 + ndata/2
  efnp <- matrix(data=0, nrow = ndata, ncol = nefn)
  
  num.cores <- nefn
  cl <- makeCluster(num.cores, type = 'FORK')
  mat <- parLapply(cl,1:nefn, fun = function(k) looping(k,efn,efnp,ndata,y,nefn,ev))
  stopCluster(cl)
  efnp <- matrix(unlist(mat), nrow = ndata, ncol = nefn) 
  
  return(efnp)
}

click <- function(time.series, threshold){
  which(abs(time.series) > threshold)
}

bottomM <- function(x, M=1){
  order(x)[1:M]
}

Fcomp <- function(i,FRsd,nord,mxdeg,Ukz,Fcoef,FPcoef,dw,SSQC,nFFT,ndata){
  # F1 and SSQRtd are mxdeg+1 by nFFT matrices
  output <- .Fortran("Fcomp",i=as.integer(i),FRsd=double(nord),nord=as.integer(nord),
                     mxdeg=as.integer(mxdeg),Ukz=as.double(Ukz),Fcoef=as.double(Fcoef),
                     FPcoef=as.double(FPcoef), dw=as.double(dw),F1=double(mxdeg+1),
                     SSQC=as.double(SSQC),nFFT=as.integer(nFFT),
                     SSQRtd=double(mxdeg+1),ndata=as.integer(ndata))
}

Floop <- function(nord,mxdeg,Ukz,Fcoef,FPcoef,dw,SSQC,nFFT,ndata){
  output <- .Fortran("Floop",FRsd=double(nord*nFFT),nord=as.integer(nord),mxdeg=as.integer(mxdeg),
                     Ukz=as.double(Ukz),Fcoef=as.double(Fcoef),FPcoef=as.double(FPcoef),
                     dw=as.double(dw),F1mat=double((mxdeg+1)*nFFT),SSQC=as.double(SSQC),
                     nFFT=as.integer(nFFT),SSQRtdmat=double((mxdeg+1)*nFFT),ndata=as.integer(ndata))
  output[c("F1mat","SSQRtdmat")]
}

FLoop2 <- function(nord, mxdeg, H, Fcoef, FPcoef, nFFT){
  output <- .Fortran("Floop2", nord = as.integer(nord),
                     mxdeg = as.integer(mxdeg),
                     H = as.double(H),
                     Fcoef = as.double(Fcoef),
                     FPcoef = as.double(FPcoef),
                     F1 = double( (mxdeg+1)*nFFT ),
                     nFFT = as.integer(nFFT),
                     projFcoef = double( nord*(mxdeg+1)*nFFT ),
                     res = double( nord*(mxdeg+1)*nFFT ))
}

unwrap <- function(phase, tol = 6){
  N <- length(phase)
  
  for(n in 1:(N-1)){
    if(phase[n+1] - phase[n] > tol){
      phase[(n+1):N] <- phase[(n+1):N] - 2*pi
    } else if (phase[n] - phase[n+1] > tol){
      phase[(n+1):N] <- phase[(n+1):N] + 2*pi
    }
  }
  return(phase)
}

F1comp <- function(mxdeg, nord, H, FPcoef, Fcoef){
  output <- .Fortran("F1comp",mxdeg=as.integer(mxdeg),
                     nord=as.integer(nord),
                     H=as.double(H),
                     FPcoef=as.double(FPcoef),
                     Fcoef=as.double(Fcoef),
                     projFcoef=double(nord*(mxdeg+1)),
                     F1=double(mxdeg+1),
                     res=double(nord*(mxdeg+1)))
}

F2comp <- function(dw, ndata, nord, res, Fcoef, mxdeg){
  output <- .Fortran("F2comp", dw=as.double(dw),ndata=as.integer(ndata),nord=as.integer(nord),
                     res=as.double(res),F2=double(mxdeg+1),Fcoef=as.double(Fcoef),
                     mxdeg=as.integer(mxdeg),Fres=double(ndata*(mxdeg+2)))
}

FPcomp <- function(yk, H, Ukz, nord, mxdeg, yPk){
  output <- .Fortran("FPcomp", yk=as.complex(yk),H=as.double(H),FP=double(mxdeg+1),Ukz=as.double(Ukz),
                     nord=as.integer(nord),mxdeg=as.integer(mxdeg),
                     res=complex(nord*(mxdeg+1)),yPk=as.complex(yPk))
}

FLoopDJT <- function(nord, mxdeg, Ukz, Fcoef, FPcoef, SSQC, nFFT){
  output <- .Fortran("FLoopDJT", FRsd = double(nord*nFFT),
                     nord = as.integer(nord),
                     mxdeg = as.integer(mxdeg),
                     Ukz = as.double(Ukz),
                     Fcoef = as.double(Fcoef),
                     FPcoef = as.double(FPcoef),
                     F1mat = double((mxdeg+1)*nFFT),
                     SSQC = as.double(SSQC),
                     nFFT = as.integer(nFFT))
}

jackknifeF <- function(mxdeg, nord, H, Fcoef){
  output <- .Fortran("jackknifeF", mxdeg = as.integer(mxdeg),
                     nord = as.integer(nord),
                     H = as.double(H),
                     Fcoef = as.double(Fcoef),
                     jkVar = double(mxdeg),
                     FPave = double(mxdeg))
}

FPLoop <- function(nord, mxdeg, H, Psi, nFFT){
  output <- .Fortran("FPLoop", nord = as.integer(nord),
                     mxdeg = as.integer(mxdeg),
                     H = as.double(H),
                     Psi = as.double(Psi),
                     FP = double( (mxdeg+1)*nFFT ),
                     nFFT = as.integer(nFFT))
}

dpsstp <- function(V, maxdeg) {
  
  # Sanity checks
  stopifnot(is.matrix(V), is.numeric(maxdeg), maxdeg>=0)
  N <- length(V[, 1])
  K <- length(V[1, ])
  P <- maxdeg + 1
  timeArr <- 1:N
  
  R <- matrix(data=0, nrow=N, ncol=P)
  U <- matrix(data=0, nrow=K, ncol=P)
  
  # Setup centered time index
  midTime <- (1+N) / 2
  scl <- 2/(N-1)
  timeArrC <- (timeArr - midTime) * scl
  
  # Chebyshev polynomials of the first kind
  R[, 1] <- 1.0
  if(maxdeg > 0) {
    R[, 2] <- timeArrC
    if(maxdeg > 1) {
      for(j in 2:maxdeg) {
        R[, (j+1)] <- 2*timeArrC * R[, j] - R[, (j-1)]
      } # end of loop on higher orders
    } # end of maxdeg > 1
  } # end of maxdeg > 0
  
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
  
  Hn <- colSums(R^2)
  return(list(U=U,R=R,Hn=Hn))
}

dpsshp2 <- function(V, maxdeg) {
  
  # Sanity checks
  stopifnot(is.matrix(V), is.numeric(maxdeg), maxdeg>=0)
  N <- length(V[, 1])
  K <- length(V[1, ])
  P <- maxdeg + 1
  timeArr <- 1:N
  
  R <- matrix(data=0, nrow=N, ncol=P)
  U <- matrix(data=0, nrow=K, ncol=P)
  
  # Setup centered time index
  midTime <- (1+N) / 2
  scl <- 2/(N-1)
  timeArrC <- (timeArr - midTime) * scl
  
  # Hermite polynomials
  R[, 1] <- 1.0
  if(maxdeg > 0) {
    R[, 2] <- timeArrC
    if(maxdeg > 1) {
      for(j in 2:maxdeg) {
        R[, (j+1)] <- timeArrC * R[, j] - (j-1) * R[, (j-1)]
      } # end of loop on higher orders
    } # end of maxdeg > 1
  } # end of maxdeg > 0
  
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
  
  Hn <- colSums(R^2)
  return(list(U=U,R=R,Hn=Hn))
}

dpssd <- function(n, nw = 4, k = 7){
  
  DW <- dpss(n, k, nw=nw, returnEigenvalues=TRUE)
  dw <- DW$v
  ev <- DW$eigen 
  
  
  W <- nw/n
  dwp <- dpssp7(dw,n,k,W,ev)
  
  D <- 2*dw[1,] - dw[2,]
  
  dd <- matrix(data=0, nrow = n, ncol = k - 2)
  ddp <- dd
  evd <- numeric(k-2)
  for(j in 1:(k-2)){
    L <- j + 2
    g <- D[j]/D[L]
    gsK <- 1/sqrt(1 + g^2)
    gsK2 <- - gsK*g
    
    dd[,j] <- gsK*dw[,j] + gsK2*dw[,L]
    ddp[,j] <- gsK * dwp[,j] + gsK2 * dwp[,L]
    evd[j] <- gsK^2 * ev[j] + gsK2^2 * ev[L]
    
    
  }
  DD <- list(v = dd, eigen = evd, deriv = ddp)
  class(DD) <- "dpss"
  return(DD)
}


ModulatedF <- function(X, nw = 4, k = 7, mxdeg = 1){
  
  #might need to add something in if the sampling rate is not 1
  X <- as.numeric(X)
  N <- length(X)
  w <- nw/N
  nFFT <- 2^ceiling(log2(2*N))
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  DW <- dpss(n=N, k=K, nw=nw) #slepian sequences
  V <- DW$v
  ev <- DW$eigen
  
  Vdot <- dpssp7(V,N,K,w,ev)
  
  ap <- dpssap(V,mxdeg,alpha=1)
  H <- ap$U
  G <- ap$R
  
  ######################################################################################
  
  tapered <- V * X
  tapered <- rbind(tapered, matrix(0, nrow = 2*nFFT-N, ncol = K))
  yk <- mvfft(tapered)
  yk <- t(yk)
  yk <- yk[,1:nFFT]
  
  #############################################################################
  Z <- V %*% yk #standard inverse
  Zdot <- Vdot %*% yk #time derivative of Z
  
  amp <- Mod(Z)
  #ampd <- Mod(Zdot)  #not needed
  #PhD <- Arg(Z) #phase of standard inverse   #also not used here
  U <- Re(Z)
  W <- Im(Z)
  Udot <- Re(Zdot)
  Wdot <- Im(Zdot)
  phi <- ( U*Wdot - Udot*W ) / ((2*pi)*amp^2)
  
  Phi <- t(V) %*% phi
  PhiP <- t(H) %*% Phi
  
  Foutput <- FLoop2(k,mxdeg,H,Phi,PhiP,nFFT)
  ModF <- matrix(data = Foutput$F1, nrow = mxdeg + 1, ncol = nFFT)
  
  output <- list(ModF = ModF, phi = phi, Phi = Phi, PhiP = PhiP, H = H, G = G)
  
}

#' Linear regression using a singular value decomposition
#' 
#' Performs a multivariate linear regression using a singular value decomposition.
#' 
#' @param x A \code{matrix} or \code{vector} containing the predictors.
#' @param y A \code{vector} containing the response.
#' 
#' @return A \code{list} containing the regression coefficients corresponding to the 
#' columns of x, the standard error estimates (*not* the variance) on those coefficients, 
#' and the eigenvalues of the decomposition (singular values squared).
#' 
#' @details Uses the process described in Mandel (1982), "Use of the Singular value 
#' Decomposition in Regression Analysis".
#' 
#' @export
svdRegression <- function(x, y){
  # X = U D V'
  sng <- svd(x)
  
  beta <- sng$v %*% ((Conj(t(sng$u)) %*% y) / sng$d)
  rownames(beta) <- colnames(x)
  colnames(beta) <- NULL
  
  # stdErr <- apply(t(apply(sng$v, 1, "/", sng$d)), 1, sum) * sd.complex(y)
  
  list(coef = t(beta), stdErr = NA, ev = sng$d^2)
}


# designing filters ..
# e-mail to Emily - July 22, 2016

# Generate Slepians (v^{(k)}_{n})
# 1.1) Set W to the low-pass cutoff when generating, use first 2NW-3 (ish)
# 1.2) Set N to the number of filter coefficients you want (an odd number since you want the filter to be symmetric.
# 2.1) Fourier transform the v_{n}^{(k)} (zero-pad out to 16,000-ish)
# 2.2) Use only the real parts of the even tapers and imaginery part of the odd tapers
# 3) Regress the value 1 onto 2.1 from (-W + \Delta f, W + \Delta f) , where \Delta f is about 10% of W.
# 4) Regression coefficients are the weights.
# 5) Add up the V's with their regression coefficients and inverse Fourier transform.
# 6) Take the first N values.  ... should be your filter coefficients.

# nflt - number of points in the filter
# ndec - decimation ratio (2 = 2:1, 5 = 5:1, etc...)
# M - number of frequency bins
filterDesign <- function(nflt, ndec = NULL, M = 2^14, wflt = NULL){
  if ((is.null(ndec) & is.null(wflt)) | (!is.null(ndec) & !is.null(wflt))){
    stop("You must set one and only one of ndec or wflt.")
  }
  
  fudge <- 1.1 # deals with passband of the filter
  if (is.null(wflt)){
    wflt <- 0.5 * fudge / ndec
  } else {
    wflt <- wflt * fudge
  }
  nw <- floor(nflt*wflt)
  k <- 2*nw - 3
  nfreq <- 1 + M/2
  # generate slepians, keep even ordered ones
  slep.tmp <- multitaper::dpss(n = nflt, k = k, nw = nw, returnEigenvalues = FALSE)$v[, seq(1, k, by = 2)]
  
  neh <- (nflt-1)/2
  nc <- neh + 1
  
  # # how *I* think it should be ... 
  slep <- matrix(0, nrow = M, ncol = ncol(slep.tmp))
  
  slep[1, ] <- slep.tmp[nc, ]
  slep[2:(neh+1), ] <- slep[M:(M-neh+1), ] <- slep.tmp[(nc + 1):nflt, ]
  
  taper <- mvfft(slep)[1:nfreq, ]
  taper.real <- Re(taper)
  freq <- seq(0, 0.5, by = 1/M)
  
  fCut1 <- tail(which(freq <= wflt), 1)
  fCut <- trunc(min(fCut1, 0.85 * M / (2*ndec))) ### This is the important piece!
  
  d.qr <- qr(taper.real[1:(fCut), ])
  coefs <- qr.coef(d.qr, rep(1, fCut))
  
  fitted <- qr.fitted(d.qr, rep(1, fCut))
  
  filter1 <- slep.tmp %*% coefs
  filter2 <- filter1 / sum(filter1)
  
  # H <- fft(c(filter2, rep(0, M - nflt)))[1:(M/2+1)]
  # plot(abs(H)[1:1000]^2, type='l', log='y')
  
  filter2
}
