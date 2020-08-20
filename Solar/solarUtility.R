library(multitaper)
load("GOLF_PM2.RData")
library(tsinterp)
library(forecast)
dyn.load("solarf95.so")

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

.is.deriv <- function(obj){
  class(obj) == "dpssp"
}

.is.dpss <- function(obj){
  class(obj) == "dpss"
}

.is.ap <- function(obj){
  class(obj) == "ass.poly"
}

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