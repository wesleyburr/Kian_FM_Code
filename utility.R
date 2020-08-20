library(parallel)
library(multitaper)
#library(Rcpp)
#library(RcppArmadillo)
#library(RcppEigen)
#sourceCpp("MatMult.cpp", rebuild = TRUE)
dyn.load("lib1.so") #main one
#dyn.load("lib2.so") #just F1comp.f95 and F2comp.f95

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
  ap <- list(U=U,R=R,Hn=Hn)
  class(ap) <- "ass.poly"
  return(ap)
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
  
  num.cores <- min(nefn, detectCores())
  cl <- makeCluster(num.cores, type = 'FORK')
  mat <- parLapply(cl,1:nefn, fun = function(k) looping(k,efn,efnp,ndata,y,nefn,ev))
  stopCluster(cl)
  efnp <- matrix(unlist(mat), nrow = ndata, ncol = nefn) 
  class(efnp) <- "dpssp"
  
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
    ndec <- ceiling(0.5/wflt)
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

findGaps <- function(data, missingSignature=NA){
  nas <- data %in% missingSignature
  N <- length(data)
  
  if (N == 0){
    stop("Data of length 0.")
  }
  
  x <- 1:N
  gaps <- c(0, 0, 0)
  numGaps <- 0
  
  i <- 1
  while (i <= N){
    gapLength <- 0
    
    while(nas[i+gapLength] && (i+gapLength <= length(nas))){
      gapLength <- gapLength + 1
    }
    
    if (gapLength > 0){
      gaps <- rbind(gaps, c(i, i+gapLength-1, gapLength), deparse.level=0)
      i <- i + gapLength
      numGaps <- numGaps + 1
    } else { i <- i + 1 }
  }
  
  if (numGaps == 0){ 
    t(as.matrix(gaps))
  } else if (numGaps == 1) {
    t(as.matrix(gaps[-1,]))
  } else { as.matrix(gaps[-1,])}
}

detrend <- function(x, W, deltat = 1){
  vars <- multitaper::multitaperTrend(as.numeric(x), B = W, deltat = deltat, t.in = 1:length(x))
  
  x - (vars[[1]] + vars[[2]]*vars[[3]])
}

findLocalFMax <- function(obj, cutoff){
  # Check whether this is a spec.mtm() object, or from my own CMV code.
  if (any(class(obj) == "Ftest")){
    Fval <- obj$mtm$Ftest
    k <- obj$mtm$k
  }  else {
    stop("obj needs to be of class 'Ftest'.")
  }
  
  # Thanks to Frank Marshall for catching the missing "-2"
  fMaxInd <- which(Fval > qf(cutoff, 2, 2*k-2))
  maxes <- c()
  
  if (length(fMaxInd) == 0){
    return(maxes)
  }
  
  for (i in 1:length(fMaxInd)){
    if (fMaxInd[i] == 1 || fMaxInd[i] == length(Fval)){
      next
    }
    
    if (Fval[fMaxInd[i]] > Fval[fMaxInd[i]-1] && 
        Fval[fMaxInd[i]] > Fval[fMaxInd[i]+1]){
      maxes <- c(maxes, fMaxInd[i])
    }
  }
  
  maxes
}

filterDesign2 <- function(nflt, ndec = NULL, M = 2^14, wflt = NULL){
  if ((is.null(ndec) & is.null(wflt)) | (!is.null(ndec) & !is.null(wflt))){
    stop("You must set one and only one of ndec or wflt.")
  }
  
  fudge <- 1.1 # deals with passband of the filter
  if (is.null(wflt)){
    wflt <- 0.5 * fudge / ndec
  } else {
    ndec <- ceiling(0.5/wflt)
    wflt <- wflt * fudge
  }
  nw <- floor(nflt*wflt)
  k <- 2*nw - 3
  nfreq <- 1 + M/2
  # generate modified slepians, keep even ordered ones
  slep.tmp <- dpssm(n = nflt, k = k, nw = nw)$v[, seq(1, k-2, by = 2)]
  
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

prewhiten <- function(series, nw = 10, k = 19, order = 1, deltat = 1, rmLineComponents = TRUE){
  #INPUTS
  #series = time series to be whitened
  #nw = time bandwidth product to use to compute the pilot spectrum
  #k = number of tapers to use to compute the pilot spectrum
  #order = order of the AR process to be fit to the residuals. should typically be "small"
  #deltat = sampling period in seconds
  #note that the confidence level is set to 1-1/L. this could be changed.
  
  #RETURNS
  #pw = the prewhitened series
  #sp.pw = the CORRECTED mtm object for the prewhitened series
  #extras = list of some other stuff possibly of interest
  #       tf = AR transfer function
  #       sigF = indices of significant periodic components from the initial Harmonic F test
  #       sp.int = intermediate (flat) spectrum
  
  stopifnot(order >= 1, deltat > 0, is.logical(rmLineComponents), 
            is.numeric(series))
  
  L <- length(series)
  nFFT <- 2^ceiling(log2(10*L))
  
  if(rmLineComponents){
    sp <- multitaper::spec.mtm(series, deltat = deltat, dtUnits = 'second', Ftest = TRUE, 
                               plot = FALSE, nFFT = nFFT)
    conf <- 1-1/L
    sigF <- findLocalFMax(sp, conf)
    
    cmv.tmp <- sp$mtm$cmv
    cmv.tmp[-sigF] <- 0
    cmv <- c(cmv.tmp, Conj(rev(cmv.tmp[-1])[-1]))
    p.recon <- Re(fft(cmv, inverse = TRUE))[1:L]
    res <- series - p.recon
  } else {
    res <- series
  }
  
  sp.pilot.tmp <- multitaper::spec.mtm(res, nw = nw, k = k, deltat = deltat, dtUnits = 'second',
                                       plot = FALSE)
  sp.pilot <- c(sp.pilot.tmp$spec, rev(sp.pilot.tmp$spec[-1])[-1])
  
  acvf <- Re(fft(sp.pilot, inverse = TRUE))/length(sp.pilot)
  gam <- acvf[1:(order+1)]
  Gam <- toeplitz(gam[-order])
  QR <- qr(Gam)
  phi <- qr.coef(QR, y = gam[-1])
  
  
  fit <- filter(res, filter = c(0,phi), sides = 1)
  if(rmLineComponents){
    pw.tmp <- fit - (res + p.recon)
  } else {
    pw.tmp <- fit - res
  }
  
  pw <- pw.tmp[-which(is.na(pw.tmp))]
  
  
  sp.int <- multitaper::spec.mtm(pw, deltat = deltat, dtUnits = 'second', plot = FALSE, 
                                 nFFT = nFFT, Ftest = TRUE)
  tf <- fft(c(-1,phi, rep(0,nFFT-order-1)))[1:sp.int$mtm$nfreqs]
  sp.pw <- sp.int
  sp.pw$spec <- sp.int$spec / Mod(tf)^2
  
  if(rmLineComponents){
    extra <- list(tf = tf, sigF = sigF, sp.int = sp.int, phi = phi)
  } else {
    extra <- list(tf = tf, sp.int = sp.int, phi = phi)
  }
  
  return(list(pw = pw, sp.pw = sp.pw, extra = extra))
}

dpssm <- function(n,nw = 4,k = 7){
  
  DW <- dpss(n=n, k=k, nw=nw, returnEigenvalues=TRUE)
  dw <- DW$v
  ev <- DW$eigen
  D <- 2*dw[1,] - dw[2,]
  
  dm <- matrix(data=0, nrow = n, ncol = k - 2)
  em <- numeric(k-2)
  for(j in 1:(k-2)){
    L <- j + 2
    g <- D[j]/D[L]
    gsK <- 1/sqrt(1 + g^2)
    gsK2 <- - gsK*g
    
    dm[,j] <- gsK*dw[,j] + gsK2*dw[,L]
    em[j] <- gsK^2 * ev[j] + gsK2^2 * ev[L]
  }
  DM <- list(v = dm, eigen = em)
  class(DM) <- "dpss"
  return(DM)
}

ModulatedF2 <- function(spec, mxdeg = 1, adaptive = TRUE, subBand = NULL){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive))
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt[subBand,] * spec$mtm$eigenCoefs[subBand,])
    } else {
      yk <- t(spec$mtm$eigenCoefs[subBand,])
    }
    nfreqs <- ncol(yk)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs)
    } else {
      yk <- t(spec$mtm$eigenCoefs)
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  N <- spec$mtm$nFFT - spec$pad
  k <- spec$mtm$k
  w <- spec$mtm$nw / N
  
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  DW <- spec$mtm$dpss
  V <- DW$v
  ev <- DW$eigen
  rm(DW)
  
  Vdot <- dpssp7(V,N,k,w,ev)
  rm(ev)
  
  ap <- dpssap(V,mxdeg,alpha=1)
  H <- ap$U
  G <- ap$R
  rm(ap)
  #############################################################################
  
  #these two operations are the problem. N x nfreqs matrix requires too much memory
  Z <- V %*% yk #standard inverse
  amp <- Mod(Z)
  U <- Re(Z)
  W <- Im(Z)
  rm(Z)
  
  Zdot <- Vdot %*% yk #time derivative of Z
  Udot <- Re(Zdot)
  Wdot <- Im(Zdot)
  rm(Zdot, yk)
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  phi <- Y / ((2*pi)*amp^2)
  rm(Y,amp)
  
  Phi <- t(V) %*% phi
  rm(V,phi)
  PhiP <- t(H) %*% Phi
  
  Foutput <- FLoop2(k,mxdeg,H,Phi,PhiP,nfreqs)
  ModF <- matrix(data = Foutput$F1, nrow = mxdeg + 1, ncol = nfreqs)
  class(ModF) <- "Ftest"
  
  output <- list(ModF = ModF, Phi = Phi, PhiP = PhiP, H = H, G = G, freq = freq)
  
}

matMult.Fortran <- function(A,B){
  stopifnot(ncol(A) == nrow(B))
  output <- .Fortran("matMult", A = as.double(A),
                     B = as.double(B),
                     m1 = nrow(A),
                     l = ncol(A),
                     n2 = ncol(B),
                     C = double(nrow(A)*ncol(B)))
  
  C <- matrix(data = output$C, nrow = nrow(A), ncol = ncol(B))
}

ModulatedFcmv <- function(spec, mxdeg = 1, subBand = NULL){
  
  #this version returns uses the cmv for inversion instead. just trying it out.
  #this shit doesn't seem to work. NaNs and stuff in phi. not sure why. 0/0
  stopifnot(any(class(spec) == 'Ftest'), mxdeg >= 1)
  
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  
  ###################################################
  # eigencoefficients and parameters
  
  
  N <- spec$mtm$nFFT - spec$pad
  k <- spec$mtm$k
  w <- spec$mtm$nw / N
  if(!is.null(subBand)){
    cmv <- spec$mtm$cmv[subBand]
    freq <- spec$freq[subBand]
    nfreqs <- length(subBand)
  } else {
    cmv <- spec$mtm$cmv[,1]
    freq <- spec$freq
    nfreqs <- spec$mtm$nfreqs
  }
  
  expo <- vandermonde.matrix(freq, N)
  Z <- expo * cmv
  Zdot <- expo * (1i * 2 * pi * freq * cmv)
  rm(expo)
  
  amp <- Mod(Z)
  U <- Re(Z)
  W <- Im(Z)
  rm(Z)
  
  Udot <- Re(Zdot)
  Wdot <- Im(Zdot)
  rm(Zdot)
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  Y <- X1 - X2
  rm(X1,X2)
  
  phi <- Y / ((2*pi)*amp^2)
  rm(Y,amp)
  phi <- t(phi)
  ##########################################################################
  V <- spec$mtm$dpss$v
  
  ap <- dpssap(V,mxdeg,alpha=1)
  H <- ap$U
  #G <- ap$R
  rm(ap)
  #############################################################################
  
  Phi <- crossprod(V, phi)
  rm(V,phi)
  PhiP <- crossprod(H, Phi)
  
  Foutput <- FLoop2(k,mxdeg,H,Phi,PhiP,nfreqs)
  ModF <- matrix(data = Foutput$F1, nrow = mxdeg + 1, ncol = nfreqs)
  class(ModF) <- "Ftest"
  
  #for now don't return phi, Z, G
  output <- list(ModF=ModF, Phi=Phi, PhiP=PhiP, H=H, freq=freq)
  
}

vandermonde.matrix <- function (alpha, n){
  #"borrowed" from matrixcalc package
  if (!is.vector(alpha)) 
    stop("argument alpha is not a vector")
  if (!is.numeric(alpha)) 
    stop("argument n is not a numeric vector")
  m <- length(alpha)
  V <- matrix(0, nrow = m, ncol = n)
  V[, 1] <- rep(1, m)
  j <- 2
  while (j <= n) {
    x <- alpha^(j - 1)
    V[, j] <- x
    j <- j + 1
  }
  return(V)
}

transferFunction <- function(imp, nFFT = 2^ceiling(log2(10*length(imp)))){
  tf <- fft(c(imp, rep(0, nFFT-length(imp))))[1:(nFFT/2)]
  fr <- seq(0, 0.5, length.out = length(tf))
  return(list(tf = tf, fr = fr))
}

prewhitenM <- function(series, nw = 10, k = 19, order = 1, deltat = 1, rmLineComponents = TRUE){
  #this uses the modified Slepians instead
  
  #INPUTS
  #series = time series to be whitened
  #nw = time bandwidth product to use to compute the pilot spectrum
  #k = number of tapers to use to compute the pilot spectrum
  #order = order of the AR process to be fit to the residuals. should typically be "small"
  #deltat = sampling period in seconds
  #note that the confidence level is set to 1-1/L. this could be changed.
  
  #RETURNS
  #pw = the prewhitened series
  #sp.pw = the CORRECTED mtm object for the prewhitened series
  #extras = list of some other stuff possibly of interest
  #       tf = AR transfer function
  #       sigF = indices of significant periodic components from the initial Harmonic F test
  #       sp.int = intermediate (flat) spectrum
  
  stopifnot(order >= 1, deltat > 0, is.logical(rmLineComponents), 
            is.numeric(series))
  
  L <- length(series)
  nFFT <- 2^ceiling(log2(10*L))
  
  if(rmLineComponents){
    DW <- dpss(n = length(series), nw = 4, k = 7)
    sp <- multitaper::spec.mtm(series, deltat = deltat, dtUnits = 'second', Ftest = TRUE, 
                               plot = FALSE, nFFT = nFFT, dpssIN = DW)
    conf <- 1-1/L
    sigF <- findLocalFMax(sp, conf)
    
    cmv.tmp <- sp$mtm$cmv
    cmv.tmp[-sigF] <- 0
    cmv <- c(cmv.tmp, Conj(rev(cmv.tmp[-1])[-1]))
    p.recon <- Re(fft(cmv, inverse = TRUE))[1:L]
    res <- series - p.recon
  } else {
    res <- series
  }
  
  DM <- dpssm(n = length(res), nw = nw + 1, k = k + 2)
  sp.pilot.tmp <- multitaper::spec.mtm(res, nw = nw+1, k = k, deltat = deltat, dtUnits = 'second',
                                       plot = FALSE, dpssIN = DM)
  sp.pilot <- c(sp.pilot.tmp$spec, rev(sp.pilot.tmp$spec[-1])[-1])
  
  acvf <- Re(fft(sp.pilot, inverse = TRUE))/length(sp.pilot)
  gam <- acvf[1:(order+1)]
  Gam <- toeplitz(gam[-order])
  QR <- qr(Gam)
  phi <- qr.coef(QR, y = gam[-1])
  
  
  fit <- filter(res, filter = c(0,phi), sides = 1)
  if(rmLineComponents){
    pw.tmp <- fit - (res + p.recon)
  } else {
    pw.tmp <- fit - res
  }
  
  pw <- pw.tmp[-which(is.na(pw.tmp))]
  
  DM <- dpssm(n = length(pw), nw = 5, k = 9)
  sp.int <- multitaper::spec.mtm(pw, deltat = deltat, dtUnits = 'second', plot = FALSE, nw = 5,
                                 nFFT = nFFT, dpssIN = DM, Ftest = TRUE, returnInternals = TRUE)
  tf <- fft(c(-1,phi, rep(0,nFFT-order-1)))[1:sp.int$mtm$nfreqs]
  sp.pw <- sp.int
  sp.pw$spec <- sp.int$spec / Mod(tf)^2
  
  if(rmLineComponents){
    extra <- list(tf = tf, sigF = sigF, sp.int = sp.int, phi = phi)
  } else {
    extra <- list(tf = tf, sp.int = sp.int, phi = phi)
  }
  
  return(list(pw = pw, sp.pw = sp.pw, extra = extra))
}

ModulatedF4 <- function(spec, mxdeg = 1, adaptive = TRUE, subBand = NULL,
                        derivIN = NULL, apIN = NULL, dpssIN = NULL){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive))
  
  N <- spec$mtm$nFFT - spec$pad
  k <- spec$mtm$k
  if(.is.ap(apIN)){
    stopifnot(ncol(apIN$U) == mxdeg+1)
  }
  
  if(.is.deriv(derivIN)){
    stopifnot(nrow(derivIN) == N)
  }
    
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt[subBand,] * spec$mtm$eigenCoefs[subBand,])
    } else {
      yk <- t(spec$mtm$eigenCoefs[subBand,])
    }
    nfreqs <- ncol(yk)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs)
    } else {
      yk <- t(spec$mtm$eigenCoefs)
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  
  
  
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    V <- dpssIN$v
    ev <- dpssIN$eigen
  } else {
    V <- dpssIN$v
    ev <- dpssIN$eigen
    stopifnot(nrow(V) == N)
  }
  
  
  if(!.is.deriv(derivIN)){
    w <- spec$mtm$nw / N
    derivIN <- dpssp11R(dpssIN,w)
  } 
  
  
  if(!.is.ap(apIN)){
    apIN <- dpssap(V,mxdeg,alpha=1)
    H <- apIN$U
    G <- apIN$R
  } else {
    H <- apIN$U
    G <- apIN$R
  }
  
  #############################################################################
  
  #these two operations are the problem. N x nfreqs matrix requires too much memory
  #also dpssp function was slowing things down. could maybe rewrite so that fewer transposes
  #are used and tcrossprod and crossprod are utilized instead of %*%
  Z <- V %*% yk #standard inverse
  amp <- Mod(Z)
  U <- Re(Z)
  W <- Im(Z)
  rm(Z)
  
  Zdot <- derivIN %*% yk #time derivative of Z
  Udot <- Re(Zdot)
  Wdot <- Im(Zdot)
  rm(Zdot, yk)
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  phi <- Y / ((2*pi)*amp^2)
  rm(Y,amp)
  
  Phi <- crossprod(V,phi)
  rm(phi)
  PhiP <- crossprod(H,Phi)
  
  Foutput <- FLoop2(k,mxdeg,H,Phi,PhiP,nfreqs)
  ModF <- matrix(data = Foutput$F1, nrow = mxdeg + 1, ncol = nfreqs)
  class(ModF) <- "Ftest"
  
  output <- list(ModF = ModF)#will add other stuff to this list at another time.
  
}
  
.is.deriv <- function(obj){
  class(obj) == "dpssp"
}

.is.dpss <- function(obj){
  class(obj) == "dpss"
}

.is.ap <- function(obj){
  class(obj) == "ass.poly"
}

dpssp8R <- function(DW,W){
  efn <- DW$v
  ev <- DW$eigen
  ndata <- nrow(efn)
  K <- ncol(efn)
  tn <- 1:(ndata-1)
  
  b <- c(0, sin( 2*pi*W*tn) / (pi * tn^2) )
  B <- toeplitz(b)
  B[upper.tri(B)] <- -B[upper.tri(B)]
  
  a <- 2 * W * c(0, cos( 2*pi*W*tn ) / tn )
  A <- toeplitz(a)
  A[upper.tri(A)] <- -A[upper.tri(A)]
  
  C <- A - B
  
  efnp <- C %*% efn
  for(k in 1:K){
    efnp[,k] <- efnp[,k] / ev[k]
  }
  
  class(efnp) <- "dpssp"
  return(efnp)
}

skew.toeplitz <- function(x, upper = TRUE){
  X <- toeplitz(x)
  if(upper){
    X[upper.tri(X)] <- -X[upper.tri(X)]
  } else {
    X[lower.tri(X)] <- -X[lower.tri(X)]
  }
  
  return(X)
}


dpssp9R <- function(DW,W){
  efn <- DW$v
  ev <- DW$eigen
  ndata <- nrow(efn)
  K <- ncol(efn)
  tn <- 1:(ndata-1)
  
  b <- c(0, sin( 2*pi*W*tn) / (pi * tn^2) )
  a <- 2 * W * c(0, cos( 2*pi*W*tn ) / tn )
  c <- a - b
  C <- skew.toeplitz(c)
  
  
  efnp <- C %*% efn
  for(k in 1:K){
    efnp[,k] <- efnp[,k] / ev[k]
  }
  
  class(efnp) <- "dpssp"
  return(efnp)
}

dpssp10R <- function(DW,W){
  efn <- DW$v
  ev <- DW$eigen
  ndata <- nrow(efn)
  K <- ncol(efn)
  tn <- 1:(ndata-1)
  
  b <- c(0, sin( 2*pi*W*tn) / (pi * tn^2) )
  a <- 2 * W * c(0, cos( 2*pi*W*tn ) / tn )
  c <- a - b
  C <- skew.toeplitz(c)
  
  #replace this part by skew.toepmult(C, efn)
  efnp <- skew.toepmult(C,efn)
  for(k in 1:K){
    efnp[,k] <- efnp[,k] / ev[k]
  }
  
  class(efnp) <- "dpssp"
  return(efnp)
}

toepmult <- function(A,B){
  stopifnot(is.matrix(A))
  if(is.null(ncol(B))){
    B <- matrix(B, ncol = 1)
  }
  n <- nrow(A)
  m <- ncol(B)
  out <- NULL
  x <- as.matrix(c(A[1,],0,A[1,][n:2]))
  for(i in 1:m){
    p <- c(B[,i],rep(0,n))
    h <- as.vector(fft(p)*fft(x))
    out <- cbind(out, Re(fft(h, inverse = TRUE)[1:n] / length(h)))
  }
  return( out )
}

#keep this thing around for reference. "borrowed" from 
# https://stackoverflow.com/questions/39810967/toeplitz-matrix-vector-multiplication-in-r
#skew.toepmult <- function(A,v){
#  n <- nrow(A)
#  x <- as.matrix(c(A[1,1],-A[1,2:n],0,A[1,n:2]))
#  p <- c(v,rep(0,n))
#  h <- as.vector(fft(p)*fft(x))
#  out <- Re(fft(h, inverse = TRUE)[1:n] / length(h))
#  return( matrix(out,n) )
#}

skew.toepmult <- function(A,B){
  stopifnot(is.matrix(A))
  if(is.null(ncol(B))){
    B <- matrix(B, ncol = 1)
  }
  n <- nrow(A)
  m <- ncol(B)
  out <- NULL
  x <- as.matrix(c(A[1,1],-A[1,2:n],0,A[1,n:2]))
  for(i in 1:m){
    p <- c(B[,i],rep(0,n))
    h <- as.vector(fft(p)*fft(x))
    out <- cbind(out,Re(fft(h, inverse = TRUE)[1:n] / length(h)))
  }
  return( out )
}

dpssp11R <- function(DW,NW){
  stopifnot(.is.dpss(DW))
  efn <- DW$v
  ev <- DW$eigen
  ndata <- nrow(efn)
  K <- ncol(efn)
  W <- NW / ndata
  tn <- 1:(ndata-1)
  
  b <- c(0, sin( 2*pi*W*tn) / (pi * tn^2) )
  a <- 2 * W * c(0, cos( 2*pi*W*tn ) / tn )
  y <- a - b
  
  #this is multiplication by a skew-symmetric toeplitz matrix.
  # the matrix is skew.toeplitz(y, upper = TRUE), meaning that 
  # we actually need to negate different stuff.justification on 9 June 2020 notes. 
  #if the skew-symmetric matrix generated by y is C, then the first row of C is 
  # C[1,] = c(y[1], -y[2:n]). then the relevant line from skew.toepmult is 
  #x <- as.matrix(c(C[1,1], -C[1,2:n], 0, C[1,n:2])), which is equivalent to
  #x <- as.matrix(c(y[1], y[2:n], 0, -y[n:2])), or x <- as.matrix(c(y,0,-y[n:2]))
  efnp <- matrix(data = NA, nrow = ndata, ncol = K)
  x <- as.matrix(c(y,0,-y[ndata:2]))
  for(k in 1:K){
    p <- c(efn[,k],rep(0,ndata))
    h <- as.vector(fft(p)*fft(x))
    efnp[,k] <- Re(fft(h, inverse = TRUE)[1:ndata] / length(h)) / ev[k] 
  }
  
  class(efnp) <- "dpssp"
  return(efnp)
}

findLocalFMaxM <- function(obj, k, cutoff){
  
  M <- nrow(obj)
  stopifnot(k > M)
  MAXES <- list()
  for(m in 1:M){
    Fval <- obj[m,]
    fMaxInd <- which(Fval > qf(cutoff, 1, k-m))
    maxes <- c()
    
    if (length(fMaxInd) == 0){
      next
    }
    
    for (i in 1:length(fMaxInd)){
      if (fMaxInd[i] == 1 || fMaxInd[i] == length(Fval)){
        next
      }
      
      if (Fval[fMaxInd[i]] > Fval[fMaxInd[i]-1] && 
          Fval[fMaxInd[i]] > Fval[fMaxInd[i]+1]){
        maxes <- c(maxes, fMaxInd[i])
      }
    }
    MAXES[[m]] <- maxes
  }
  return(MAXES)
}

ModulatedF5 <- function(spec, mxdeg = 1, adaptive = TRUE, subBand = NULL,
                        derivIN = NULL, apIN = NULL, dpssIN = NULL, dsIN = NULL,
                        downsample = FALSE){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive))
  
  N <- spec$mtm$nFFT - spec$pad
  k <- spec$mtm$k
  nw <- spec$mtm$nw
  deltat <- spec$mtm$deltaT
  if(.is.ap(apIN)){
    stopifnot(ncol(apIN$U) == mxdeg+1)
  }
  
  if(.is.deriv(derivIN)){
    stopifnot(nrow(derivIN) == N)
  }
  
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  # can downsample to 10*K. pass in dpss object with length 10*K
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt[subBand,] * spec$mtm$eigenCoefs[subBand,])
    } else {
      yk <- t(spec$mtm$eigenCoefs[subBand,])
    }
    nfreqs <- ncol(yk)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs)
    } else {
      yk <- t(spec$mtm$eigenCoefs)
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  
  
  
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    #don't actually think the deltat matters. it should cancel out everywhere
    V <- dpssIN$v * sqrt(deltat)
    ev <- dpssIN$eigen
  } else {
    V <- dpssIN$v * sqrt(deltat)
    ev <- dpssIN$eigen
    stopifnot(nrow(V) == N, ncol(V) == K)
  }
    
  if(!.is.deriv(derivIN)){
    derivIN <- dpssp11R(dpssIN,nw)
  } 
  
  
  #############################################################################
  
  #these two operations are the problem. N x nfreqs matrix requires too much memory
  #also dpssp function was slowing things down. could maybe rewrite so that fewer transposes
  #are used and tcrossprod and crossprod are utilized instead of %*%
  
  if(downsample){
    L <- 10*k
    deltat <- floor(N/L)
    ds <- deltat * (1:L)
    V <- V[ds,] 
    derivIN <- derivIN[ds,]
    N <- L
  }
  
  Z <- V %*% yk #standard inverse
  amp <- Mod(Z)
  U <- Re(Z)
  W <- Im(Z)
  rm(Z)
  
  Zdot <- derivIN %*% yk #time derivative of Z
  Udot <- Re(Zdot)
  Wdot <- Im(Zdot)
  rm(Zdot, yk)
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  phi <- Y / ((2*pi)*amp^2)
  #remove the mean from instantaneous frequency? otherwise two noncentral chisq
  phi <- apply(phi, MARGIN = 2, FUN = function(x) x - mean(x))
  rm(Y,amp)
  
  if(downsample){
    stopifnot(nrow(dsIN$v) == N)
    if(!.is.dpss(dsIN)){
      dsIN <- dpss(n = N, nw = nw, k = k)
    }
    
    V <- dsIN$v * sqrt(deltat)
  }
  
  if(!.is.ap(apIN) || nrow(apIN$R) != N){
    apIN <- dpssap(V,mxdeg,alpha=1)
  }
  H <- apIN$U
  G <- apIN$R
  
  Phi <- crossprod(V,phi)
  rm(phi)
  PhiP <- crossprod(H,Phi)
  
  ModF <- modulatedFs(G,H,V,Phi,PhiP,nfreqs,mxdeg+1,k)
  class(ModF) <- "Ftest"
  
  output <- list(ModF = ModF)#will add other stuff to this list at another time.
  
}

modulatedFs <- function(G,H,V,Phi,PhiP,nfreqs,P,K){
  
  F1 <- F2 <- F3 <- matrix(data = NA, nrow = P, ncol = ncol(Phi))
  
  for(p in 1:P){
    #rewrite these using simplified forms.
    projPhi <- H[,1:p, drop = FALSE] %*% PhiP[1:p,, drop = FALSE]
    res <- Phi - projPhi
    Res <- V %*% res
    
    Fit <- G[,1:p, drop = FALSE] %*% PhiP[1:p,, drop = FALSE]
    
    ssq2 <- colSums(Fit^2)
    ssq3 <- colSums(Res^2)
    ssq4 <- colSums(res^2)
    ssq5 <- colSums(projPhi^2)
    
    F1[p,] <- (K/p - 1) * ssq2 / ssq3
    
    F2[p,] <- (K/p - 1) * ssq5 / ssq4
  }
  F3.temp <- FLoop2(K,P-1,H,Phi,PhiP,nfreqs)
  F3 <- matrix(data = F3.temp$F1, nrow = P, ncol = nfreqs)
  
  return(list(F1=F1,F2=F2,F3=F3))
}

bigMatMult <- function(A,B){
  stopifnot(ncol(A) == nrow(B))
  
  N <- nrow(A)
  numBlox <- detectCores()-1
  n <- ceiling(N/numBlox)
  m <- -N %% numBlox
  M <- N + m
  idx <- seq(1, M, by = n)
  A <- rbind(A, matrix(0, nrow = m, ncol = ncol(A)))
  clust <- makeCluster(numBlox)
  temp <- parLapply(cl = clust, X = idx, fun = function(x) blockMult(x,A,B,n))
  stopCluster(clust)
  C.temp <- matrix(data = unlist(temp), nrow = n*numBlox, ncol = ncol(B))
  C <- C.temp[1:N,]
}

blockMult <- function(ind,A,B,len){
  A[ind+len-1,] %*% B
}

ModulatedF6 <- function(spec, mxdeg = 1, adaptive = TRUE, subBand = NULL,
                        derivIN = NULL, apIN = NULL, dpssIN = NULL, dsIN = NULL){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive))
  
  N <- spec$mtm$nFFT - spec$pad
  k <- spec$mtm$k
  nw <- spec$mtm$nw
  deltat <- spec$mtm$deltaT
  if(.is.ap(apIN)){
    stopifnot(ncol(apIN$U) == mxdeg+1)
  }
  
  if(.is.deriv(derivIN)){
    stopifnot(nrow(derivIN) == N)
  }
  
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  # can downsample to 10*K. pass in dpss object with length 10*K
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt[subBand,] * spec$mtm$eigenCoefs[subBand,])
    } else {
      yk <- t(spec$mtm$eigenCoefs[subBand,])
    }
    nfreqs <- ncol(yk)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs)
    } else {
      yk <- t(spec$mtm$eigenCoefs)
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  
  
  
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    #don't actually think the deltat matters. it should cancel out everywhere
    V <- dpssIN$v 
    ev <- dpssIN$eigen
  } else {
    V <- dpssIN$v
    ev <- dpssIN$eigen
    stopifnot(nrow(V) == N, ncol(V) == K)
  }
  
  if(!.is.deriv(derivIN)){
    derivIN <- dpssp11R(dpssIN,nw)
  } 
  
  
  #############################################################################
  
  #N x nfreqs matrix requires too much memory
  #maybe change the amount of downsampling.
  
  if(!is.null(dsIN)){
    if(!.is.dpss(dsIN)){
      L <- 10*k
      dsIN <- dpss(n = L, nw = nw, k = k)
    }
    L <- nrow(dsIN$v)
    deltat <- floor(N/L)
    ds <- deltat * (1:L)
    V <- V[ds,] 
    derivIN <- derivIN[ds,]
    N <- L
  }
  
  Z <- V %*% yk #standard inverse
  amp <- Mod(Z)
  U <- Re(Z)
  W <- Im(Z)
  rm(Z)
  
  Zdot <- derivIN %*% yk #time derivative of Z
  Udot <- Re(Zdot)
  Wdot <- Im(Zdot)
  rm(Zdot, yk)
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  phi <- Y / ((2*pi)*amp^2)
  #remove the mean from instantaneous frequency? otherwise two noncentral chisq
  phi <- apply(phi, MARGIN = 2, FUN = function(x) x - mean(x))
  rm(Y,amp)
  
  if(!is.null(dsIN)){
    V <- dsIN$v 
  }
  
  if(!.is.ap(apIN) || nrow(apIN$R) != N){
    apIN <- dpssap(V,mxdeg,alpha=1)
  }
  H <- apIN$U
  G <- apIN$R
  
  
  Phi <- crossprod(V,phi)
  PhiP <- crossprod(H,Phi)
  
  ModF <- modulatedFs3(G,Phi,PhiP)

  auxiliary <- list(phi=phi, 
                    Phi=Phi,
                    PhiP=PhiP,
                    nfreqs=nfreqs,
                    mxdeg=mxdeg,
                    subBand=subBand,
                    dpssp=derivIN,
                    dsss=dsIN,
                    dpssap=apIN
  )
  class(ModF) <- "Ftest"
  
  output <- list(ModF = ModF, aux = auxiliary)
  
}

modulatedFs2 <- function(G,Phi,PhiP){
  P <- ncol(G)
  K <- nrow(Phi)
  nfreqs <- ncol(Phi)
 
  
  F1 <- F2 <- F3 <- matrix(data = NA, nrow = P, ncol = ncol(Phi))
  
  ssq1 <- colSums(Phi^2)
  for(p in 1:P){
    Fit <- G[,1:p, drop = FALSE] %*% PhiP[1:p,, drop = FALSE]
    
    ssq2 <- colSums(Fit^2)
    ssq3 <- colSums(PhiP[1:p,, drop = FALSE]^2)
    
    F1[p,] <- (K/p - 1) * ssq1 / (ssq1 - ssq3)
    F2[p,] <- (K/p - 1) * ssq2 / (ssq1 - ssq3)
    #if I can get an easy closed form for the proportionality between F1, F2 
    #can generate these faster
   
    
  }
  F3 <- FLoop3(P-1,K,PhiP,Phi,nfreqs)
  Fs <- list(F1=F1,F2=F2,F3=F3)
  
  return(Fs)
}

FLoop3 <- function(mxdeg, nord, FPcoef, Fcoef, nfreqs){
  output <- .Fortran("Floop3", 
                     mxdeg = as.integer(mxdeg),
                     nord = as.integer(nord),
                     FPcoef = as.double(FPcoef),
                     Fcoef = as.double(Fcoef),
                     Fp = double( (mxdeg+1)*nfreqs ),
                     nfreqs = as.integer(nfreqs)
  )
  Fp <- matrix(data = output$Fp, nrow = mxdeg+1, ncol = nfreqs)
}

modF2 <- function(G,Phi,PhiP){
  P <- ncol(G)
  K <- nrow(Phi)
  nfreqs <- ncol(Phi)
  
  F2 <- matrix(data = NA, nrow = P, ncol = ncol(Phi))
  
  
  for(p in 1:P){
    Fit <- G[,1:p, drop = FALSE] %*% PhiP[1:p,, drop = FALSE]
    
    ssq1 <- colSums(Fit^2)
    ssq2 <- colSums(Phi^2)
    ssq3 <- colSums(PhiP[1:p,, drop = FALSE]^2)
    
    F2[p,] <- (K/p - 1) * ssq1 / (ssq2 - ssq3)
    
  }
  
  
  return(F2)
}

ModulatedF2Only <- function(spec, mxdeg = 1, adaptive = TRUE, subBand = NULL,
                            derivIN = NULL, apIN = NULL, dpssIN = NULL, dsIN = NULL){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive))
  
  N <- spec$mtm$nFFT - spec$pad
  k <- spec$mtm$k
  nw <- spec$mtm$nw
  deltat <- spec$mtm$deltaT
  if(.is.ap(apIN)){
    stopifnot(ncol(apIN$U) == mxdeg+1)
  }
  
  if(.is.deriv(derivIN)){
    stopifnot(nrow(derivIN) == N)
  }
  
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  # can downsample to 10*K. pass in dpss object with length 10*K
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt[subBand,] * spec$mtm$eigenCoefs[subBand,])
    } else {
      yk <- t(spec$mtm$eigenCoefs[subBand,])
    }
    nfreqs <- ncol(yk)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs)
    } else {
      yk <- t(spec$mtm$eigenCoefs)
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  
  
  
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    #don't actually think the deltat matters. it should cancel out everywhere
    V <- dpssIN$v #* sqrt(deltat)
    ev <- dpssIN$eigen
  } else {
    V <- dpssIN$v #* sqrt(deltat)
    ev <- dpssIN$eigen
    stopifnot(nrow(V) == N, ncol(V) == K)
  }
  
  if(!.is.deriv(derivIN)){
    derivIN <- dpssp11R(dpssIN,nw)
  } 
  
  
  #############################################################################
  
  #N x nfreqs matrix requires too much memory
  #maybe change the amount of downsampling.
  
  if(!is.null(dsIN)){
    if(!.is.dpss(dsIN)){
      L <- 10*k
      dsIN <- dpss(n = L, nw = nw, k = k)
    }
    L <- nrow(dsIN$v)
    deltat <- floor(N/L)
    ds <- deltat * (1:L)
    V <- V[ds,] 
    derivIN <- derivIN[ds,]
    N <- L
  }
  
  Z <- V %*% yk #standard inverse
  amp <- Mod(Z)
  U <- Re(Z)
  W <- Im(Z)
  rm(Z)
  
  Zdot <- derivIN %*% yk #time derivative of Z
  Udot <- Re(Zdot)
  Wdot <- Im(Zdot)
  rm(Zdot, yk)
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  phi <- Y / ((2*pi)*amp^2)
  #remove the mean from instantaneous frequency? otherwise two noncentral chisq
  phi <- apply(phi, MARGIN = 2, FUN = function(x) x - mean(x))
  rm(Y,amp)
  
  if(!is.null(dsIN)){
    V <- dsIN$v #* sqrt(deltat)
  }
  
  if(!.is.ap(apIN) || nrow(apIN$R) != N){
    apIN <- dpssap(V,mxdeg,alpha=1)
  }
  H <- apIN$U
  G <- apIN$R
  
  
  Phi <- crossprod(V,phi)
  rm(phi)
  PhiP <- crossprod(H,Phi)
  
  ModF <- modF2(G,Phi,PhiP)
  class(ModF) <- "Ftest"
  
  output <- list(ModF = ModF)#will add other stuff to this list at another time.
  
}

ModulatedF3Only <- function(spec, mxdeg = 1, adaptive = TRUE, subBand = NULL,
                            derivIN = NULL, apIN = NULL, dpssIN = NULL, dsIN = NULL){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive))
  
  N <- spec$mtm$nFFT - spec$pad
  k <- spec$mtm$k
  nw <- spec$mtm$nw
  deltat <- spec$mtm$deltaT
  if(.is.ap(apIN)){
    stopifnot(ncol(apIN$U) == mxdeg+1)
  }
  
  if(.is.deriv(derivIN)){
    stopifnot(nrow(derivIN) == N)
  }
  
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  # can downsample to 10*K. pass in dpss object with length 10*K
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt[subBand,] * spec$mtm$eigenCoefs[subBand,])
    } else {
      yk <- t(spec$mtm$eigenCoefs[subBand,])
    }
    nfreqs <- ncol(yk)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs)
    } else {
      yk <- t(spec$mtm$eigenCoefs)
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  
  
  
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    #don't actually think the deltat matters. it should cancel out everywhere
    V <- dpssIN$v #* sqrt(deltat)
    ev <- dpssIN$eigen
  } else {
    V <- dpssIN$v #* sqrt(deltat)
    ev <- dpssIN$eigen
    stopifnot(nrow(V) == N, ncol(V) == K)
  }
  
  if(!.is.deriv(derivIN)){
    derivIN <- dpssp11R(dpssIN,nw)
  } 
  
  
  #############################################################################
  
  #N x nfreqs matrix requires too much memory
  #maybe change the amount of downsampling.
  
  if(!is.null(dsIN)){
    if(!.is.dpss(dsIN)){
      L <- 10*k
      dsIN <- dpss(n = L, nw = nw, k = k)
    }
    L <- nrow(dsIN$v)
    deltat <- floor(N/L)
    ds <- deltat * (1:L)
    V <- V[ds,] 
    derivIN <- derivIN[ds,]
    N <- L
  }
  
  Z <- V %*% yk #standard inverse
  amp <- Mod(Z)
  U <- Re(Z)
  W <- Im(Z)
  rm(Z)
  
  Zdot <- derivIN %*% yk #time derivative of Z
  Udot <- Re(Zdot)
  Wdot <- Im(Zdot)
  rm(Zdot, yk)
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  phi <- Y / ((2*pi)*amp^2)
  #remove the mean from instantaneous frequency? otherwise two noncentral chisq
  phi <- apply(phi, MARGIN = 2, FUN = function(x) x - mean(x))
  rm(Y,amp)
  
  if(!is.null(dsIN)){
    V <- dsIN$v #* sqrt(deltat)
  }
  
  if(!.is.ap(apIN) || nrow(apIN$R) != N){
    apIN <- dpssap(V,mxdeg,alpha=1)
  }
  H <- apIN$U
  G <- apIN$R
  
  
  Phi <- crossprod(V,phi)
  rm(phi)
  PhiP <- crossprod(H,Phi)
  
  ModF <- FLoop3(mxdeg,K,PhiP,Phi,nfreqs)
  class(ModF) <- "Ftest"
  
  output <- list(ModF = ModF)#will add other stuff to this list at another time.
  
}

ModulatedF1Only <- function(spec, mxdeg = 1, adaptive = TRUE, subBand = NULL,
                            derivIN = NULL, apIN = NULL, dpssIN = NULL, dsIN = NULL){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive))
  
  N <- spec$mtm$nFFT - spec$pad
  k <- spec$mtm$k
  nw <- spec$mtm$nw
  deltat <- spec$mtm$deltaT
  if(.is.ap(apIN)){
    stopifnot(ncol(apIN$U) == mxdeg+1)
  }
  
  if(.is.deriv(derivIN)){
    stopifnot(nrow(derivIN) == N)
  }
  
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  # can downsample to 10*K. pass in dpss object with length 10*K
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt[subBand,] * spec$mtm$eigenCoefs[subBand,])
    } else {
      yk <- t(spec$mtm$eigenCoefs[subBand,])
    }
    nfreqs <- ncol(yk)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- t(spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs)
    } else {
      yk <- t(spec$mtm$eigenCoefs)
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  
  
  
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    #don't actually think the deltat matters. it should cancel out everywhere
    V <- dpssIN$v #* sqrt(deltat)
    ev <- dpssIN$eigen
  } else {
    V <- dpssIN$v #* sqrt(deltat)
    ev <- dpssIN$eigen
    stopifnot(nrow(V) == N, ncol(V) == K)
  }
  
  if(!.is.deriv(derivIN)){
    derivIN <- dpssp11R(dpssIN,nw)
  } 
  
  
  #############################################################################
  
  #N x nfreqs matrix requires too much memory
  #maybe change the amount of downsampling.
  
  if(!is.null(dsIN)){
    if(!.is.dpss(dsIN)){
      L <- 10*k
      dsIN <- dpss(n = L, nw = nw, k = k)
    }
    L <- nrow(dsIN$v)
    deltat <- floor(N/L)
    ds <- deltat * (1:L)
    V <- V[ds,] 
    derivIN <- derivIN[ds,]
    N <- L
  }
  
  Z <- V %*% yk #standard inverse
  amp <- Mod(Z)
  U <- Re(Z)
  W <- Im(Z)
  rm(Z)
  
  Zdot <- derivIN %*% yk #time derivative of Z
  Udot <- Re(Zdot)
  Wdot <- Im(Zdot)
  rm(Zdot, yk)
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  phi <- Y / ((2*pi)*amp^2)
  #remove the mean from instantaneous frequency? otherwise two noncentral chisq
  phi <- apply(phi, MARGIN = 2, FUN = function(x) x - mean(x))
  rm(Y,amp)
  
  if(!is.null(dsIN)){
    V <- dsIN$v #* sqrt(deltat)
  }
  
  if(!.is.ap(apIN) || nrow(apIN$R) != N){
    apIN <- dpssap(V,mxdeg,alpha=1)
  }
  H <- apIN$U
  G <- apIN$R
  
  
  Phi <- crossprod(V,phi)
  rm(phi)
  PhiP <- crossprod(H,Phi)
  
  ModF <- modF1(Phi,PhiP)
  class(ModF) <- "Ftest"
  
  output <- list(ModF = ModF)#will add other stuff to this list at another time.
  
}

modF1 <- function(Phi,PhiP){
  P <- nrow(PhiP)
  K <- nrow(Phi)
  nfreqs <- ncol(Phi)
  
  F1 <- matrix(data = NA, nrow = P, ncol = ncol(Phi))
  
  
  for(p in 1:P){
    ssq1 <- colSums(Phi^2)
    ssq2 <- colSums(PhiP[1:p,, drop = FALSE]^2)
    
    F1[p,] <- (K/p - 1) * ssq2 / (ssq1 - ssq2)
    
  }
  
  
  return(F1)
}

modulatedFs3 <- function(G,Phi,PhiP){
  P <- ncol(G)
  K <- nrow(Phi)
  N <- nrow(G)
  nfreqs <- ncol(Phi)
  
  F1 <- F2 <- F3 <- matrix(data = NA, nrow = P, ncol = ncol(Phi))
  ssq1 <- colSums(Phi^2)
  ssq3 <- PhiP[1,]^2
  
  F1[1,] <- (K-1) * ssq3 / (ssq1 - ssq3)
  F2[1,] <- N*G[1,1]^2*F1[1,]
  
  for(p in 2:P){
    Gc <- G[,1:p] %*% PhiP[1:p,]
    ssq2 <- colSums(Gc^2)
    
    ssq3 <- ssq3 + PhiP[p,]^2
    
    
    F1[p,] <- (K/p - 1) * ssq3 / (ssq1 - ssq3)
    
    F2[p,] <- (K/p - 1) * ssq2 / (ssq1 - ssq3)
    
    
  }
  F3 <- FLoop3(P-1,K,PhiP,Phi,nfreqs)
  
  return(list(F1=F1,F2=F2,F3=F3))
}

ModulatedF8 <- function(spec, mxdeg = 1, adaptive = TRUE, subBand = NULL,
                        derivIN = NULL, apIN = NULL, dpssIN = NULL, dsIN = NULL,
                        dsTapers = "dpss", BoxCox = FALSE){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive))
  
  N <- spec$mtm$nFFT - spec$pad
  k <- spec$mtm$k
  nw <- spec$mtm$nw
  deltat <- spec$mtm$deltaT
  if(.is.ap(apIN)){
    stopifnot(ncol(apIN$U) == mxdeg+1)
  }
  
  if(.is.deriv(derivIN)){
    stopifnot(nrow(derivIN) == N)
  }
  
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  # can downsample to 10*K. pass in dpss object with length 10*K
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt[subBand,] * spec$mtm$eigenCoefs[subBand,]
    } else {
      yk <- spec$mtm$eigenCoefs[subBand,]
    }
    nfreqs <- ncol(yk)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs
    } else {
      yk <- spec$mtm$eigenCoefs
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  
  
  
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    #don't actually think the deltat matters. it should cancel out everywhere
    V <- dpssIN$v #* sqrt(deltat)
    ev <- dpssIN$eigen
  } else {
    V <- dpssIN$v 
    ev <- dpssIN$eigen
    stopifnot(nrow(V) == N, ncol(V) == K)
  }
  
  if(!.is.deriv(derivIN)){
    if(dsTapers == "sine"){
      derivIN <- sDeriv(N,K)
    } else {
      derivIN <- dpssp11R(dpssIN,nw)
    }
  } 
  
  
  #############################################################################
  
  #N x nfreqs matrix requires too much memory
  # what length should I downsample to? K corresponds to the actual new nyquist of W,
  # but it just doesn't work. maybe 10K or N/100 or something
  
  if(!is.null(dsIN)){
    if(!.is.dpss(dsIN)){
      L <- 10*k
      if(dsTapers == "sine"){
        dsIN <- sTaper(N = L, K = K)
      } else {
        dsIN <- dpss(n = L, nw = nw, k = k)
      }
    }
    #low <- ceiling(N*0.05)
    #hi <- floor(N-N*0.05)
    L <- nrow(dsIN$v)
    #ds <- seq(low, hi, length.out = L)
    ds <- seq(1,N, length.out = L)
    V <- V[ds,] 
    derivIN <- derivIN[ds,]
    N <- L
  }
  
  Z <- tcrossprod(V,yk) #standard inverse
  amp <- Mod(Z)
  U <- Re(Z)
  W <- Im(Z)
  rm(Z)
  
  
  Udot <- tcrossprod(derivIN, Re(yk))
  Wdot <- tcrossprod(derivIN, Im(yk))
  
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  phi <- Y / ((2*pi)*amp^2)
  #remove the mean from instantaneous frequency? otherwise two noncentral chisq
  phi <- apply(phi, MARGIN = 2, FUN = function(x) x - mean(x))
  if(BoxCox){
    bc <- matrixBoxCox(phi, lambdaIN = rep(0.5, nfreqs))
    phi <- bc$A2
  }
  rm(Y,amp)
  
  if(!is.null(dsIN)){
    V <- dsIN$v
  }
  
  if(!.is.ap(apIN) || nrow(apIN$R) != N){
    apIN <- dpssap(V,mxdeg,alpha=1)
  }
  H <- apIN$U
  G <- apIN$R
  
  
  Phi <- crossprod(V,phi)
  rm(phi)
  PhiP <- crossprod(H,Phi)
  
  
  ModF <- modulatedFs3(G,Phi,PhiP)
  
  class(ModF) <- "Ftest"
  
  return(ModF)
  
}

ModulatedF10 <- function(spec, mxdeg = 1, adaptive = FALSE, subBand = NULL,
                         derivIN = NULL, apIN = NULL, dpssIN = NULL, dsIN = NULL,
                         BoxCox = FALSE, smooth = -1){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive), smooth %% 1 == 0)
  
  N <- spec$mtm$nFFT - spec$pad
  k <- spec$mtm$k
  nw <- spec$mtm$nw
  
  if(.is.ap(apIN)){
    stopifnot(ncol(apIN$U) == mxdeg+1)
  }
  
  if(.is.deriv(derivIN)){
    stopifnot(nrow(derivIN) == N)
  }
  
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  # can downsample to 10*K. pass in dpss object with length 10*K
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt[subBand,] * spec$mtm$eigenCoefs[subBand,]
    } else {
      yk <- spec$mtm$eigenCoefs[subBand,]
    }
    nfreqs <- ncol(yk)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs
    } else {
      yk <- spec$mtm$eigenCoefs
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  
  
  
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    V <- dpssIN$v
    ev <- dpssIN$eigen
  } else {
    V <- dpssIN$v 
    ev <- dpssIN$eigen
    stopifnot(nrow(V) == N, ncol(V) == K)
  }
  
  if(!.is.deriv(derivIN)){
    derivIN <- dpssp11R(dpssIN,nw)
  } 
  
  if(!.is.ap(apIN) || nrow(apIN$R) != N){
    apIN <- dpssap(V,mxdeg,alpha=1)
  }
  H <- apIN$U
  G <- apIN$R
  D <- apIN$D
  
  
  #############################################################################
  
  #N x nfreqs matrix requires too much memory
  # what length should I downsample to? K corresponds to the actual new nyquist of W,
  # but it just doesn't work. maybe 10K or N/100 or something
  
  if(!is.null(dsIN)){
    if(!.is.dpss(dsIN)){
      L <- 10*k
      dsIN <- dpss(n = L, nw = nw, k = k)
    }
    #low <- ceiling(N*0.05)
    #hi <- floor(N-N*0.05)
    L <- nrow(dsIN$v)
    #ds <- seq(low, hi, length.out = L)
    ds <- seq(1,N, length.out = L)
    V <- V[ds,] 
    derivIN <- derivIN[ds,]
    N <- L
  }
  
  if(smooth >= 0){
    Hproj <- diag(x=1, nrow = k) - tcrossprod(H[,1:(smooth+1)])
    Gproj <- tcrossprod(G[,1:(smooth+1)], H[,1:(smooth+1)])
    A <- tcrossprod(V, Hproj) + Gproj
    U <- tcrossprod(A, Re(yk))
    W <- tcrossprod(A, Im(yk))
    amp2 <- U^2 + W^2
    
    
    Dproj <- tcrossprod(D[,1:(smooth+1)], H[,1:(smooth+1)])
    Adot <- tcrossprod(derivIN, Hproj) + Dproj
    
    Udot <- tcrossprod(Adot, Re(yk))
    Wdot <- tcrossprod(Adot, Im(yk))
  } else { 
    U <- tcrossprod(V, Re(yk)) #standard inverse
    W <- tcrossprod(V, Im(yk))
    amp2 <- U^2 + W^2
    
    Udot <- tcrossprod(derivIN, Re(yk))
    Wdot <- tcrossprod(derivIN, Im(yk))
  }
  
  
  
  
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  phi <- Y / ((2*pi)*amp2)
  #remove the mean from instantaneous frequency? otherwise two noncentral chisq
  phi <- apply(phi, MARGIN = 2, FUN = function(x) x - mean(x))
  if(BoxCox){
    bc <- matrixBoxCox(phi, lambdaIN = rep(0.5, nfreqs))
    phi <- bc$A2
  }
  rm(Y,amp2)
  
  if(!is.null(dsIN)){
    V <- dsIN$v
  }
  
  
  
  
  Phi <- crossprod(V,phi)
  PhiP <- crossprod(H,Phi)
  
  ModF <- modulatedFs3(G,Phi,PhiP)
  
  class(ModF) <- "Ftest"
  
  return(ModF)
  
}

ModulatedF12 <- function(spec, mxdeg = 1, adaptive = FALSE, subBand = NULL,
                         derivIN = NULL, apIN = NULL, dpssIN = NULL, dsIN = NULL,
                         BoxCox = FALSE, smooth = -1, scale = FALSE, scalesIN = NULL,
                         smoothIN = NULL){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive), smooth %% 1 == 0)
  
  N <- spec$mtm$nFFT - spec$pad
  k <- spec$mtm$k
  nw <- spec$mtm$nw
  
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  # can downsample to 10*K. pass in dpss object with length 10*K
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt[subBand,,drop=FALSE] * spec$mtm$eigenCoefs[subBand,,drop=FALSE]
    } else {
      yk <- spec$mtm$eigenCoefs[subBand,,drop=FALSE]
    }
    nfreqs <- ncol(yk)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs
    } else {
      yk <- spec$mtm$eigenCoefs
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  
  
  
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    V <- dpssIN$v
    ev <- dpssIN$eigen
  } else {
    V <- dpssIN$v 
    ev <- dpssIN$eigen
    stopifnot(nrow(V) == N, ncol(V) == k)
  }
  
  if(!.is.deriv(derivIN)){
    derivIN <- dpssp11R(dpssIN,nw)
  } 
  
  if(!.is.ap(apIN)){
    apIN <- dpssap(V,mxdeg,alpha=1)
  }
  H <- apIN$U
  G <- apIN$R
  D <- apIN$D
  
  
  #############################################################################
  
  #N x nfreqs matrix requires too much memory
  # what length should I downsample to? K corresponds to the actual new nyquist of W,
  # but it just doesn't work. maybe 10K or N/100 or something
  
  if(!is.null(dsIN)){
    if(!.is.dpss(dsIN)){
      L <- 10*k
      dsIN <- dpss(n = L, nw = nw, k = k)
    }
    V <- dsIN$v
    #L <- nrow(dsIN$v)
    #ds <- seq(1,N, length.out = L)
    #V <- V[ds,] 
    #derivIN <- derivIN[ds,]
    #scalesIN$scl <- scalesIN$scl[ds]
    #N <- L
  }
  
  if(smooth >= 0){
    if(is.null(smoothIN)){
      smoothIN <- assPolySmoother(apIN, V = V, Vdot = derivIN, smooth=smooth)
    } 
    
    A <- smoothIN$A
    Adot <- smoothIN$Adot
    
    U <- tcrossprod(A, Re(yk)) 
    W <- tcrossprod(A, Im(yk)) 
    
    amp2 <- U^2 + W^2 
    
    Udot <- tcrossprod(Adot, Re(yk))
    Wdot <- tcrossprod(Adot, Im(yk)) 
    
    
  } else { 
    U <- tcrossprod(V, Re(yk))  #standard inverse
    W <- tcrossprod(V, Im(yk)) 
    
    amp2 <- U^2 + W^2
    
    Udot <- tcrossprod(derivIN, Re(yk)) 
    Wdot <- tcrossprod(derivIN, Im(yk)) 
    
  }
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  if(scale){
    if(is.null(scalesIN)){
      if(smooth >= 0){
          scl <- inverseScales(A,Adot)$scl
      } else {
        scl <- inverseScales(V,derivIN)$scl
      }
    } else {
      scl <- scalesIN$scl
    }
    phi <- scl * Y / ((2*pi)*amp2)
  } else {
    phi <- Y / (2*pi*amp2)
  }
  
   
  #remove the mean from instantaneous frequency? otherwise two noncentral chisq
  #phi <- apply(phi, MARGIN = 2, FUN = function(x) x - mean(x))
  if(BoxCox){
    bc <- matrixBoxCox(phi, lambdaIN = rep(0.5, nfreqs))
    phi <- bc$A2
  }
  rm(Y)
  
  if(!is.null(dsIN)){
    V <- dsIN$v
  }
  
  Phi <- crossprod(V,phi)
  PhiP <- crossprod(H,Phi)
  
  ModF <- modulatedFs3(G,Phi,PhiP)
  
  class(ModF) <- "Ftest"
  
  
  return(ModF=ModF)
  
}

inverseScales <- function(V, Vdot){
  scl1 <- rowSums(V^2)
  scl2 <- rowSums(Vdot^2)
  scl3 <- rowSums(V*Vdot)^2 / (scl1 * scl2)
  
  scl <- 2*pi*sqrt(scl1) / sqrt( scl2 * (1 - scl3) )
  out <- list(scl = scl, scl1 = scl1, scl2 = scl2, scl3 = scl3)
  return(out)
}

assPolySmoother <- function(AP, V, Vdot, smooth = 0){
  H <- AP$U
  G <- AP$R
  D <- AP$D
  k <- nrow(H)
  
  Hproj <- diag(x=1, nrow = k) - tcrossprod(H[,1:(smooth+1), drop = FALSE])
  Gproj <- tcrossprod(G[,1:(smooth+1), drop = FALSE], H[,1:(smooth+1), drop = FALSE])
  A <- tcrossprod(V, Hproj) + Gproj
  
  Dproj <- tcrossprod(D[,1:(smooth+1), drop = FALSE], H[,1:(smooth+1), drop = FALSE])
  Adot <- tcrossprod(Vdot, Hproj) + Dproj
  
  out = list(A=A, Adot=Adot)
  return(out)
}

dpssapD <- function(V, maxdeg, alpha=0.75, deriv = FALSE) {
  
  # Sanity checks
  stopifnot(is.matrix(V), is.numeric(maxdeg), maxdeg>=0)
  N <- length(V[, 1])
  K <- length(V[1, ])
  P <- maxdeg + 1
  timeArr <- 1:N
  
  if(!deriv){
    R <- matrix(data=0, nrow=N, ncol=P)
    U <- matrix(data=0, nrow=K, ncol=P)
    D <- NULL
    
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
  }  else {
    R <- D <- matrix(data=0, nrow=N, ncol=P)
    U <- matrix(data=0, nrow=K, ncol=P)
    
    
    # Setup centered time index
    midTime <- (1+N) / 2
    tscl <- 2/(N-1)
    timeArrC <- (timeArr - midTime) * tscl
    
    # Start with Gegenbauer polynomials; convergence is faster
    #alpha <- 0.75
    R[, 1] <- 1.0
    if(maxdeg > 0) {
      R[, 2] <- 2 * alpha * timeArrC
      D[, 2] <- tscl 
      if(maxdeg > 1) {
        for(j in 2:maxdeg) {
          A1 <- 2 * ( (j-1) + alpha ) / j
          A2 <- ( (j-2) + 2 * alpha ) / j
          
          R[, (j+1)] <- A1 * timeArrC * R[, j] - A2 * R[, (j-1)]
          D[, (j+1)] <- tscl * j * R[, j]
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
      D[, L] <- D[, L] * scl
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
          D[, L] <- D[, L] - scl * D[, j]
        }
        scl <- 1 / sqrt(sum(U[, L]^2))
        U[, L] <- U[, L] * scl  # orthonormalize
        R[, L] <- R[, L] * scl
        D[, L] <- D[, L] * scl
      }
    }
  }  
  
  Hn <- colSums(R^2)
  ap <- list(U=U,R=R,D=D,Hn=Hn)
  class(ap) <- "ass.poly"
  return(ap)
}

analSig <- function(spec,dpssIN = NULL, derivIN = NULL, scale = FALSE, scaleIN = NULL,
                    subBand = NULL){
  stopifnot(.is.dpss(dpssIN), .is.deriv(derivIN))
  V <- dpssIN$v
  sp <- spec$spec
  
  if(scale){
    if(is.null(scaleIN)){
      scaleIN <- inverseScales(V,derivIN)
    }
    scl <- scaleIN$scl
    scl1 <- scaleIN$scl1
    scl2 <- scaleIN$scl2
    scl3 <- scaleIN$scl3
    
    if(!is.null(subBand)){
      yk <- spec$mtm$eigenCoefs[subBand,,drop = FALSE] / sqrt(0.5*sp[subBand,drop=FALSE])
    } else {
      yk <- spec$mtm$eigenCoefs / sqrt(0.5*sp)
    }
    
    U <- tcrossprod(V, Re(yk)) / sqrt(scl1) #standard inverse
    W <- tcrossprod(V, Im(yk)) / sqrt(scl1)
      
    amp2 <- (U^2 + W^2)
      
    Udot <- tcrossprod(derivIN, Re(yk)) / sqrt(scl2)
    Wdot <- tcrossprod(derivIN, Im(yk)) / sqrt(scl2)
      
    
    X1 <- U*Wdot
    rm(U,Wdot)
    X2 <- Udot*W
    rm(Udot,W)
    
    Y <- (X1 - X2) / sqrt(1 - scl3)
    
    phi <- Y / (2*pi*amp2)
  } else {
    if(!is.null(subBand)){
      yk <- spec$mtm$eigenCoefs[subBand,, drop = FALSE]
    } else {
      yk <- spec$mtm$eigenCoefs
    }
    
    U <- tcrossprod(V, Re(yk))  
    W <- tcrossprod(V, Im(yk)) 
    
    amp2 <- (U^2 + W^2)
    
    Udot <- tcrossprod(derivIN, Re(yk)) 
    Wdot <- tcrossprod(derivIN, Im(yk)) 
    
    
    X1 <- U*Wdot
    rm(U,Wdot)
    X2 <- Udot*W
    rm(Udot,W)
    
    Y <- (X1 - X2) 
    
    phi <- Y / (2*pi*amp2)
  }
  out <- list(num = Y, amp2 = amp2, phi = phi)
  return(out)
}

ModulatedF13 <- function(spec, mxdeg = 1, adaptive = FALSE, subBand = NULL,
                         derivIN = NULL, apIN = NULL, dpssIN = NULL,
                         smooth = -1, scale = FALSE, scalesIN = NULL,
                         smoothIN = NULL){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive), smooth %% 1 == 0)
  
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  # can downsample to 10*K. pass in dpss object with length 10*K
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt[subBand,,drop=FALSE] * spec$mtm$eigenCoefs[subBand,,drop=FALSE]
    } else {
      yk <- spec$mtm$eigenCoefs[subBand,,drop=FALSE]
    }
    nfreqs <- nrow(yk)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs
    } else {
      yk <- spec$mtm$eigenCoefs
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    V <- dpssIN$v
  } else {
    V <- dpssIN$v
  }
  
  if(!.is.deriv(derivIN)){
    derivIN <- dpssp11R(dpssIN,spec$mtm$nw)
  } 
  
  if(!.is.ap(apIN)){
    apIN <- dpssap(V,mxdeg,alpha=1)
  }
  H <- apIN$U
  G <- apIN$R
  D <- apIN$D
  
  
  #############################################################################
  
  if(smooth >= 0){
    if(is.null(smoothIN)){
      smoothIN <- assPolySmoother(apIN, V = V, Vdot = derivIN, smooth=smooth)
    } 
    
    A <- smoothIN$A
    Adot <- smoothIN$Adot
    
    U <- tcrossprod(A, Re(yk)) 
    W <- tcrossprod(A, Im(yk)) 
    
    amp2 <- U^2 + W^2 
    
    Udot <- tcrossprod(Adot, Re(yk))
    Wdot <- tcrossprod(Adot, Im(yk)) 
    
    
  } else { 
    U <- tcrossprod(V, Re(yk))  #standard inverse
    W <- tcrossprod(V, Im(yk)) 
    
    amp2 <- U^2 + W^2
    
    Udot <- tcrossprod(derivIN, Re(yk)) 
    Wdot <- tcrossprod(derivIN, Im(yk)) 
    
  }
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  if(scale){
    if(is.null(scalesIN)){
      if(smooth >= 0){
        scl <- inverseScales(A,Adot)$scl
      } else {
        scl <- inverseScales(V,derivIN)$scl
      }
    } else {
      scl <- scalesIN$scl
    }
    phi <- scl * Y / ((2*pi)*amp2)
  } else {
    phi <- Y / (2*pi*amp2)
  }
  
  rm(Y)
  phi <- apply(phi, MARGIN = 2, FUN = function(x) x - mean(x))
  
  Phi <- crossprod(V,phi)
  PhiP <- crossprod(H,Phi)
  
  ModF <- modulatedFs3(G,Phi,PhiP)
  
  class(ModF) <- "Ftest"
  
  
  return(ModF)
  
}

modulatedF1 <- function(Phi,PhiP){
  P <- nrow(PhiP)
  K <- nrow(Phi)
  
  F1 <- matrix(data = NA, nrow = P, ncol = ncol(Phi))
  ssq1 <- colSums(Phi^2)
  ssq3 <- PhiP[1,]^2
  
  F1[1,] <- (K-1) * ssq3 / (ssq1 - ssq3)
  
  for(p in 2:P){
    ssq3 <- ssq3 + PhiP[p,]^2
    F1[p,] <- (K/p - 1) * ssq3 / (ssq1 - ssq3)
     
  }
  
  return(F1)
}
ModulatedF14 <- function(spec, mxdeg = 1, adaptive = FALSE, subBand = NULL,
                         derivIN = NULL, apIN = NULL, dpssIN = NULL,
                         smooth = -1, scale = FALSE, scalesIN = NULL,
                         smoothIN = NULL){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive), smooth %% 1 == 0)
  
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  # can downsample to 10*K. pass in dpss object with length 10*K
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt[subBand,,drop=FALSE] * spec$mtm$eigenCoefs[subBand,,drop=FALSE]
    } else {
      yk <- spec$mtm$eigenCoefs[subBand,,drop=FALSE]
    }
    nfreqs <- length(subBand)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs
    } else {
      yk <- spec$mtm$eigenCoefs
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    V <- dpssIN$v
  } else {
    V <- dpssIN$v
  }
  
  if(!.is.deriv(derivIN)){
    derivIN <- dpssp11R(dpssIN,spec$mtm$nw)
  } 
  
  if(!.is.ap(apIN)){
    apIN <- dpssap(V,mxdeg,alpha=1)
  }
  H <- apIN$U
  G <- apIN$R
  D <- apIN$D
  
  
  #############################################################################
  
  if(smooth >= 0){
    if(is.null(smoothIN)){
      smoothIN <- assPolySmoother(apIN, V = V, Vdot = derivIN, smooth=smooth)
    } 
    
    A <- smoothIN$A
    Adot <- smoothIN$Adot
    
    U <- tcrossprod(A, Re(yk)) 
    W <- tcrossprod(A, Im(yk)) 
    
    amp2 <- U^2 + W^2 
    
    Udot <- tcrossprod(Adot, Re(yk))
    Wdot <- tcrossprod(Adot, Im(yk)) 
    
    
  } else { 
    U <- tcrossprod(V, Re(yk))  #standard inverse
    W <- tcrossprod(V, Im(yk)) 
    
    amp2 <- U^2 + W^2
    
    Udot <- tcrossprod(derivIN, Re(yk)) 
    Wdot <- tcrossprod(derivIN, Im(yk)) 
    
  }
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  if(scale){
    if(is.null(scalesIN)){
      if(smooth >= 0){
        scl <- inverseScales(A,Adot)$scl
      } else {
        scl <- inverseScales(V,derivIN)$scl
      }
    } else {
      scl <- scalesIN$scl
    }
    phi <- scl * Y / ((2*pi)*amp2)
  } else {
    phi <- Y / (2*pi*amp2)
  }
  
  rm(Y)
  Phi <- crossprod(V,phi)
  PhiP <- crossprod(H,Phi)
  
  K <- spec$mtm$k
  F3 <- FLoop3(mxdeg,K,PhiP,Phi,nfreqs)
  
  phi <- apply(phi, MARGIN = 2, FUN = function(x) x - mean(x))
  Phi <- crossprod(V,phi)
  PhiP <- crossprod(H,Phi)
  
  F1 <- modulatedF1(Phi,PhiP)
  
  
  ModF <- list(F1=F1,F3=F3)
  
  class(ModF) <- "Ftest"
  
  
  return(ModF)
  
}

ModulatedF15 <- function(spec, mxdeg = 1, adaptive = FALSE, subBand = NULL,
                         derivIN = NULL, apIN = NULL, dpssIN = NULL,
                         smooth = -1, scale = FALSE, scalesIN = NULL,
                         smoothIN = NULL){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive), smooth %% 1 == 0)
  
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  # can downsample to 10*K. pass in dpss object with length 10*K
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt[subBand,,drop=FALSE] * spec$mtm$eigenCoefs[subBand,,drop=FALSE]
    } else {
      yk <- spec$mtm$eigenCoefs[subBand,,drop=FALSE]
    }
    nfreqs <- length(subBand)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs
    } else {
      yk <- spec$mtm$eigenCoefs
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    V <- dpssIN$v
  } else {
    V <- dpssIN$v
  }
  
  if(!.is.deriv(derivIN)){
    derivIN <- dpssp11R(dpssIN,spec$mtm$nw)
  } 
  
  if(!.is.ap(apIN)){
    apIN <- dpssap(V,mxdeg,alpha=1)
  }
  H <- apIN$U
  G <- apIN$R
  D <- apIN$D
  
  
  #############################################################################
  
  if(smooth >= 0){
    if(is.null(smoothIN)){
      smoothIN <- assPolySmoother(apIN, V = V, Vdot = derivIN, smooth=smooth)
    } 
    
    A <- smoothIN$A
    Adot <- smoothIN$Adot
    
    U <- tcrossprod(A, Re(yk)) 
    W <- tcrossprod(A, Im(yk)) 
    
    amp2 <- U^2 + W^2 
    
    Udot <- tcrossprod(Adot, Re(yk))
    Wdot <- tcrossprod(Adot, Im(yk)) 
    
    
  } else { 
    U <- tcrossprod(V, Re(yk))  #standard inverse
    W <- tcrossprod(V, Im(yk)) 
    
    amp2 <- U^2 + W^2
    
    Udot <- tcrossprod(derivIN, Re(yk)) 
    Wdot <- tcrossprod(derivIN, Im(yk)) 
    
  }
  
  X1 <- U*Wdot
  rm(U,Wdot)
  X2 <- Udot*W
  rm(Udot,W)
  
  Y <- X1 - X2
  rm(X1,X2)
  
  if(scale){
    if(is.null(scalesIN)){
      if(smooth >= 0){
        scl <- inverseScales(A,Adot)$scl
      } else {
        scl <- inverseScales(V,derivIN)$scl
      }
    } else {
      scl <- scalesIN$scl
    }
    phi <- scl * Y / ((2*pi)*amp2)
  } else {
    phi <- Y / (2*pi*amp2)
  }
  
  rm(Y)
  
  Phi <- crossprod(V,phi)
  PhiP <- crossprod(H,Phi)
  
  K <- spec$mtm$k
  #F3 <- FLoop3(mxdeg,K,PhiP,Phi,nfreqs)
  F3stuff <- F3ncp(Phi,PhiP)
  F3 <- F3stuff$F3
  ncp <- F3stuff$ncp
  
  phi <- apply(phi, MARGIN = 2, FUN = function(x) x - mean(x))
  Phi <- crossprod(V,phi)
  PhiP <- crossprod(H,Phi)
  
  F1 <- modulatedF1(Phi,PhiP)
  
  
  ModF <- list(F1=F1,F3=F3, ncp = ncp)
  
  class(ModF) <- "Ftest"
  
  
  return(ModF)
  
}

F3ncp <- function(Phi, PhiP){
  P <- nrow(PhiP)-1
  K <- nrow(Phi)
  
  ssq1 <- colSums(Phi^2)
  ncp <- ssq1 - PhiP[1,]^2
  
  F3 <- matrix(0, nrow = P, ncol = ncol(Phi))
  ssq2 <- 0
  for(p in 1:P){
    ssq2 <- ssq2 + PhiP[p+1,]^2
    F3[p,] <- (ssq1 - ssq2) / ((K - p)*PhiP[p+1,]^2)
  }
  out <- list(F3 = F3, ncp = ncp)
  
  return(out)
}


FcompsNo0 <- function(Phi, PhiP){
  P <- nrow(PhiP)-1
  K <- nrow(Phi)
  
  ssq1 <- colSums(Phi^2)
  
  F1 <- F3 <- matrix(0, nrow = P, ncol = ncol(Phi))
  ssq2 <- 0
  for(p in 1:P){
    ssq2 <- ssq2 + PhiP[p+1,]^2
    F1[p,] <- (K/p - 1)*ssq2 / (ssq1 - ssq2)
    F3[p,] <- (K - p)*PhiP[p+1,]^2 / (ssq1 - ssq2) 
  }
  out <- list(F1=F1, F3=F3)
  return(out)
}

ModulatedF16 <- function(spec, mxdeg = 1, adaptive = FALSE, subBand = NULL,
                         derivIN = NULL, apIN = NULL, dpssIN = NULL,
                         smooth = -1, scale = FALSE, scalesIN = NULL,
                         smoothIN = NULL){
  
  stopifnot(any(class(spec) == 'mtm'), mxdeg >= 1, is.logical(adaptive), smooth %% 1 == 0)
  
  
  # adaptive: use the adaptively weighted eigencoefficients or don't. 
  # subBand: the INDICES of the frequency range to compute the modulated 
  #          F statistic at. get these from the spec object. could also
  #          just use dropFreqs and pass in that spec object.
  #can now pass in a dpss object, dpssp (derivative) object (recommended), or 
  # an associated polynomial object.
  #the dpssp function has also been sped up by like 100 times
  # can downsample to 10*K. pass in dpss object with length 10*K
  ###################################################
  # eigencoefficients and parameters
  if(!is.null(subBand)){
    freq <- spec$freq[subBand]
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt[subBand,,drop=FALSE] * spec$mtm$eigenCoefs[subBand,,drop=FALSE]
    } else {
      yk <- spec$mtm$eigenCoefs[subBand,,drop=FALSE]
    }
    nfreqs <- length(subBand)
  } else {
    freq <- spec$freq
    if(adaptive){
      yk <- spec$mtm$eigenCoefWt * spec$mtm$eigenCoefs
    } else {
      yk <- spec$mtm$eigenCoefs
    }
    nfreqs <- spec$mtm$nfreqs
  }
  
  ##########################################################################
  #slepians, derivatives and special polynomials
  
  if(!.is.dpss(dpssIN)){
    dpssIN <- spec$mtm$dpss
    V <- dpssIN$v
  } else {
    V <- dpssIN$v
  }
  
  if(!.is.deriv(derivIN)){
    derivIN <- dpssp11R(dpssIN,spec$mtm$nw)
  } 
  
  if(!.is.ap(apIN)){
    apIN <- dpssap(V,mxdeg,alpha=1)
  }
  H <- apIN$U
  G <- apIN$R
  D <- apIN$D
  
  
  #############################################################################
  
  if(smooth >= 0){
    if(is.null(smoothIN)){
      smoothIN <- assPolySmoother(apIN, V = V, Vdot = derivIN, smooth=smooth)
    } 
    
    A <- smoothIN$A
    Adot <- smoothIN$Adot
    
    U <- tcrossprod(A, Re(yk)) 
    W <- tcrossprod(A, Im(yk)) 
    Udot <- tcrossprod(Adot, Re(yk))
    Wdot <- tcrossprod(Adot, Im(yk)) 
    num <- U*Wdot - Udot*W
    rm(Udot, Wdot)
    
    amp2 <- U^2 + W^2 
    rm(U,W)
    
  } else { 
    U <- tcrossprod(V, Re(yk))  #standard inverse
    W <- tcrossprod(V, Im(yk)) 
    
    Udot <- tcrossprod(derivIN, Re(yk)) 
    Wdot <- tcrossprod(derivIN, Im(yk)) 
    
    num <- U*Wdot - Udot*W
    rm(Udot, Wdot)
    amp2 <- U^2 + W^2
    rm(U,W)
  }
  
  if(scale){
    if(is.null(scalesIN)){
      if(smooth >= 0){
        scl <- inverseScales(A,Adot)$scl
      } else {
        scl <- inverseScales(V,derivIN)$scl
      }
    } else {
      scl <- scalesIN$scl
    }
    phi <- scl * num / ((2*pi)*amp2)
  } else {
    phi <- num / (2*pi*amp2)
  }
  
  rm(num, amp2)
  
  Phi <- crossprod(V,phi)
  PhiP <- crossprod(H,Phi)
  
  Fs <- modulatedFs5(Phi,PhiP)
  F1 <- Fs$F1
  F3 <- Fs$F3
  
  mFs <- FcompsNo0(Phi,PhiP)
  mF1 <- mFs$F1
  mF3 <- mFs$F3
  
  
  ModF <- list(F1=F1,F3=F3,mF1=mF1,mF3=mF3)
  
  class(ModF) <- "Ftest"
  
  
  return(ModF)
  
}

modulatedFs4 <- function(Phi,PhiP){
  P <- nrow(PhiP)
  K <- nrow(Phi)
  nfreqs <- ncol(Phi)
  
  F1 <- F3 <- matrix(data = NA, nrow = P, ncol = ncol(Phi))
  ssq1 <- colSums(Phi^2)
  ssq2 <- 0
  
  for(p in 1:P){
    ssq2 <- ssq2 + PhiP[p,]^2
    F1[p,] <- (K/p - 1) * ssq2 / (ssq1 - ssq2)
  }
  F3 <- FLoop3(P-1,K,PhiP,Phi,nfreqs)
  
  return(list(F1=F1,F3=F3))
}

modulatedFs5 <- function(Phi,PhiP){
  P <- nrow(PhiP)
  K <- nrow(Phi)
  
  F1 <- F3 <- matrix(data = NA, nrow = P, ncol = ncol(Phi))
  ssq1 <- colSums(Phi^2)
  ssq2 <- 0
  
  for(p in 1:P){
    ssq2 <- ssq2 + PhiP[p,]^2
    F1[p,] <- (K/p - 1) * ssq2 / (ssq1 - ssq2)
    F3[p,] <- (K - p) * PhiP[p,]^2 / (ssq1 - ssq2)
  }
  
  return(list(F1=F1,F3=F3))
}

findLocalFMaxM1 <- function(obj, k, cutoff){
  
  M <- nrow(obj)
  stopifnot(k > M)
  MAXES <- list()
  for(m in 1:M){
    Fval <- obj[m,]
    fMaxInd <- which(Fval > qf(cutoff, m, k-m))
    maxes <- c()
    
    if (length(fMaxInd) == 0){
      next
    }
    
    for (i in 1:length(fMaxInd)){
      if (fMaxInd[i] == 1 || fMaxInd[i] == length(Fval)){
        next
      }
      
      if (Fval[fMaxInd[i]] > Fval[fMaxInd[i]-1] && 
          Fval[fMaxInd[i]] > Fval[fMaxInd[i]+1]){
        maxes <- c(maxes, fMaxInd[i])
      }
    }
    MAXES[[m]] <- maxes
  }
  return(MAXES)
}

findLocalFMaxM2 <- function(obj, k, cutoff, m, type = 1){
  
  
  stopifnot(k > m, (type == 1 | type == 3))
  
  Fval <- obj[m,]
  if(type == 1){
    fMaxInd <- which(Fval > qf(cutoff, m, k-m))
  } else {
    fMaxInd <- which(Fval > qf(cutoff, 1, k-m))
  }
  maxes <- c()
    
  if (length(fMaxInd) > 0){
    for (i in 1:length(fMaxInd)){
      if (fMaxInd[i] == 1 || fMaxInd[i] == length(Fval)){
        next
      }
      
      if (Fval[fMaxInd[i]] > Fval[fMaxInd[i]-1] && 
          Fval[fMaxInd[i]] > Fval[fMaxInd[i]+1]){
        maxes <- c(maxes, fMaxInd[i])
      }
    }
  }

  return(maxes)
}

ModulatedF17 <- function(yk, derivIN = NULL, apIN = NULL, dpssIN = NULL){
  
  #speed version of this function. everything must be passed in.
  K <- ncol(yk)
  V <- dpssIN$v
  H <- apIN$U
  P <- ncol(H)
  nfreqs <- nrow(yk)
  
  #############################################################################

  phi <- IFcompute(yk,V,derivIN)
  
  
  Phi <- crossprod(V,phi)
  rm(phi)
  PhiP <- crossprod(H,Phi)
  
  F3 <- FLoop3(P-1,K,PhiP,Phi,nfreqs)
  
  mFs <- modulatedFs5(Phi,PhiP[2:P,,drop = FALSE])
  mF1 <- mFs$F1
  mF3 <- mFs$F3
  
  
  ModF <- list(F3=F3,mF1=mF1,mF3=mF3)
  
  return(ModF)
  
}

HarmonicF <- function(yk, dw){
  
  k <- ncol(yk)
  Ukz <- colSums(dw)
  ssqUkz <- sum(Ukz^2)
  cmv <- (yk %*% Ukz) / ssqUkz
  ssqave <- ssqUkz * Mod(cmv)^2
  Ukz <- as.matrix(Ukz)
  
  ssqres <- apply( Mod(yk - (cmv %*% t(Ukz)))^2, MARGIN = 1, FUN = sum)
  HF <- (k-1) * ssqave / ssqres
  class(HF) <- "Ftest"
  
  out <- list(HF = HF, cmv = cmv)
  
  return(out)
}


ModulatedF18 <- function(yk, derivIN = NULL, apIN = NULL, dpssIN = NULL){
  

  #speed version of this function. everything must be passed in.
  #just computes modified F3. note that H matrix is missing first column,
  #the zero degree column
  K <- ncol(yk)
  V <- dpssIN$v
  H <- apIN$U[,-1, drop = FALSE]
  P <- ncol(H)
  nfreqs <- nrow(yk)
  
  #############################################################################
  
  phi <- IFcompute(yk, V, derivIN)
  
  Phi <- crossprod(V,phi)
  rm(phi)
  PhiP <- crossprod(H,Phi)
  
  mF3 <- FLoop4(P,K,PhiP,Phi,nfreqs)
  
  return(mF3)
  
}

FLoop4 <- function(mxdeg, nord, FPcoef, Fcoef, nfreqs){
  output <- .Fortran("FLoop4", 
                     mxdeg = as.integer(mxdeg),
                     nord = as.integer(nord),
                     FPcoef = as.double(FPcoef),
                     Fcoef = as.double(Fcoef),
                     Fp = double( mxdeg*nfreqs ),
                     nfreqs = as.integer(nfreqs)
  )
  Fp <- matrix(data = output$Fp, nrow = mxdeg, ncol = nfreqs)
}

filterDesign3 <- function(nflt, ndec = NULL, M = 2^14, wflt = NULL){
  if ((is.null(ndec) & is.null(wflt)) | (!is.null(ndec) & !is.null(wflt))){
    stop("You must set one and only one of ndec or wflt.")
  }
  
  fudge <- 1.1 # deals with passband of the filter
  if (is.null(wflt)){
    wflt <- 0.5 * fudge / ndec
  } else {
    ndec <- ceiling(0.5/wflt)
    wflt <- wflt * fudge
  }
  nw <- floor(nflt*wflt)
  k <- 2*nw - 1
  nfreq <- 1 + M/2
  # generate modified slepians, keep even ordered ones
  slep.tmp <- dpssm(n = nflt, k = k, nw = nw)$v[, seq(1, k-2, by = 2)]
  
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
  
  #H <- fft(c(filter2, rep(0, M - nflt)))[1:(M/2+1)]
  #plot(abs(H)[1:1000]^2, type='l', log='y')
  
  filter2
}

IFcompute <- function(yk, V, Vdot){
  U <- tcrossprod(V, Re(yk))
  W <- tcrossprod(V, Im(yk))
  Udot <- tcrossprod(Vdot, Re(yk))
  Wdot <- tcrossprod(Vdot, Im(yk))
  
  num <- U*Wdot - Udot*W
  rm(Udot, Wdot)
  amp2 <- U^2 + W^2
  rm(U,W)
  
  phi <- num / (2 * pi * amp2)
  
  return(phi)
}

jkFreq <- function(Phi, H, freq){
  K <- nrow(H)
  P <- ncol(H)
  
  F3 <- matrix(0, nrow = P, ncol = length(freq))
  fkp <- matrix(0, nrow = K, ncol = P)
  for(k in 1:K){
    ssq1 <- colSums( Phi[-k,]^2 )
    PhiP <- crossprod(H[-k,], Phi[-k,])
    ssq2 <- 0
    for(p in 1:P){
      ssq2 <- ssq2 + PhiP[p,]^2
      F3[p,] <- (K-1-p) * PhiP[p,]^2 / (ssq1 - ssq2)
      fkp[k,p] <- freq[which.max(F3[p,])]
    }
  }
  jkMean <- apply(fkp, MARGIN = 2, FUN = function(x) mean(x))
  jkVar <- numeric(P)
  for(p in 1:P){
    jkVar[p] <- (1-1/K) * sum( (fkp[,p] - jkMean[p])^2 )
  }
  
  jk <- list(jkMean = jkMean, jkVar = jkVar)
  return(jk)
}

ModulatedF19 <- function(yk, mxdeg = 1, derivIN = NULL, apIN = NULL, dpssIN = NULL,
                         deltat = 1, freq){
  
  stopifnot(mxdeg >= 1)
  #speed version of this function. everything must be passed in.
  #just computes modified F3. note that H matrix is missing first column,
  #the zero degree column.
  # also does jackknife estimates of mean and variance of frequency
  K <- ncol(yk)
  V <- dpssIN$v
  H <- apIN$U[,2:(mxdeg+1), drop = FALSE]
  nfreqs <- nrow(yk)
  #############################################################################
  
  phi <- IFcompute(yk, V, derivIN)
  
  Phi <- crossprod(V,phi)
  rm(phi)
  PhiP <- crossprod(H,Phi)
  
  mF3 <- FLoop4(mxdeg,K,PhiP,Phi,nfreqs)
  jk <- jkFreq(Phi,H,freq)
  
  ModF <- list(mF3 = mF3, jk = jk)
  
  return(ModF)
  
}


ModulatedF20 <- function(yk, derivIN = NULL, apIN = NULL, dpssIN = NULL){
  
  #speed version of this function. everything must be passed in.
  #computes all 4 test statistics
  K <- ncol(yk)
  V <- dpssIN$v
  H <- apIN$U
  P <- ncol(H)
  #############################################################################
  
  phi <- IFcompute(yk, V, derivIN)
  
  Phi <- crossprod(V,phi)
  rm(phi)
  PhiP <- crossprod(H,Phi)
  
  Fs <- modulatedFs5(Phi, PhiP)
  mFs <- modulatedFs5(Phi, PhiP[2:P,,drop = FALSE])
  
  ModF <- list(F1 = Fs$F1, F3 = Fs$F3, mF1 = mFs$F1, mF3 = mFs$F3)
  
  return(ModF)
  
}

findLocalFMax2 <- function(obj, cutoff, k){
  # Check whether this is a spec.mtm() object, or from my own CMV code.
  if (any(class(obj) == "Ftest")){
    Fval <- obj
  }  else {
    stop("obj needs to be of class 'Ftest'.")
  }
  
  
  fMaxInd <- which(Fval > qf(cutoff, 2, 2*k-2))
  maxes <- c()
  
  if (length(fMaxInd) == 0){
    return(maxes)
  }
  
  for (i in 1:length(fMaxInd)){
    if (fMaxInd[i] == 1 || fMaxInd[i] == length(Fval)){
      next
    }
    
    if (Fval[fMaxInd[i]] > Fval[fMaxInd[i]-1] && 
        Fval[fMaxInd[i]] > Fval[fMaxInd[i]+1]){
      maxes <- c(maxes, fMaxInd[i])
    }
  }
  
  maxes
}
