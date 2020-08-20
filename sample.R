
FPcomp <- function(yk, H, Ukz, nord, mxdeg, yPk){
  output <- .Fortran("FPcomp", yk=as.complex(yk),H=as.double(H),Ukz=as.double(Ukz),
                     nord=as.integer(nord),mxdeg=as.integer(mxdeg),FP=double(mxdeg+1),
                     res=complex(nord*(mxdeg+1)),yPk=as.complex(yPk))
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

##########################################################################
#set up frequency modulation and generate data
ndata <- 2000 #N
nppb <- ndata  #use Dave's block function to set this if I do this over blocks
nord <- 15 #This is K = 2NW - 1
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
PhMod <- cumsum(Fmod)*dt*2*pi - Fmod*dt*2*pi
cd <- cos(PhMod) #+ 1i * sin(PhMod) #data without noise


snr <- 20.0 											#signal to noise ratio
sigmat2 <- ndata/(2*snr*nord)  # looks like a variance, but of what?   - WHY THIS VARIANCE?
sigmat <- sqrt(sigmat2)
CDN <- cd + rnorm(ndata, sd = sigmat) #+ rnorm(ndata, sd = sigmat) * 1i
CDN <- CDN * exp( -1i*2*pi*f*dt*t )
##########################################################################
#slepians, derivatives and special polynomials
dw <- dpss(ndata, nord, ndata * W) #slepian sequences
dw <- dw$v

HeSmx <- 0.9
hp <- dpsshp(dw,nppb,nord,mxdeg,HeSmx,norm = TRUE)
Poly <- hp$Poly
H <- hp$A


yk <- t(dw) %*% CDN

yPk <- t(H) %*% yk
Ukz <- apply(dw, MARGIN = 2, FUN = sum)
FPstuff <- FPcomp(yk=yk,H=H,Ukz=Ukz,nord=nord,mxdeg=mxdeg,yPk=yPk)