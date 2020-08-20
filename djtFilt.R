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

decim <- function(dat, file, step){
  filter <- scan(file, skip=2)
  nFilt <- length(filter)
  nDat <- length(dat)
  neh <- (nFilt-1)/2
  st <- neh+1
  lt <- nDat - neh
  filt <- rep(NA,nDat)
  
  for(j in seq(st,lt,step)) {
    filt[j] <- sum(dat[(j-neh):(j+neh)] * filter)
  }
  
  filt
}
