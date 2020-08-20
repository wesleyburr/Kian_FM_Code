#This simulation is to get a better idea of the distribution of the instantaneous frequency. It's
# hard analytically, so hopefully this will be illuminating.



rm(list = ls())
source("utility.R")
##########################################################################
#set up frequency modulation and generate data
N <- 500 #N
K <- 10 #This is K = 2NW - 1
nFFT <- 2^ceiling(log2(2*N))
A <- 1
f0 <- 0.3 #carrier frequency
B <- 0.005
rho <- 0.8
stm11 <- 2.0/(N-1)
dt <- 1
t <- 0:(N-1)        #time
tt <- t * stm11 - 1.0  # this runs from -1 to 1
phi.prime <- B * rho * (1.0 - 2.0 * tt^2)
Fmod <- f0 + phi.prime  # quadratic frequency modulation
phi <- cumsum(phi.prime)*2*pi*dt
Phi <- 2*pi*f0*dt*t + phi

sig <- A*cos(Phi)

num.sim <- 1000
IFarray <- array(0, dim = c(N,nFFT/2+1,num.sim))
AMParray <- IFarray
NUMarray <- IFarray

low.seed <- 69
high.seed <- 5318008

w <- 0.01
sigma <- 1
snr <- A^2 / (4*w*sigma^2)


seeds <- NULL
dw <- dpss(N,K,N*w)
V <- dw$v
lambda <- dw$eigen
Vdot <- dpssp7(V,N,K,w,lambda)

start.time <- Sys.time()
for(i in 1:num.sim){
  seed <- sample(low.seed:high.seed, size = 1)
  seeds <- c(seeds, seed)
  set.seed(seed)
  
  noise <- rnorm(N, sd = sigma)
  x <- sig + noise
  
  sp.X <- spec.mtm(x, k = K, nw = N*w, deltat = 1, plot = FALSE, dpssIN = V, nFFT = nFFT, returnInternals = TRUE)
  
  Y.unweighted <- t(sp.X$mtm$eigenCoefs)
  dk <- t(sqrt(sp.X$mtm$eigenCoefWt))
  Y <- dk * Y.unweighted
  
  Z <- V %*% Y
  U <- Re(Z)
  W <- Im(Z)
  amp2 <- Mod(Z)^2
  
  
  Zdot <- Vdot %*% Y
  Udot <- Re(Zdot)
  Wdot <- Im(Zdot)
  num <- Wdot*U - Udot*W
  
  psi <- num / (2*pi*amp2)
  IFarray[,,i] <- psi
  AMParray[,,i] <- amp2
  NUMarray[,,i] <- num
  
}
sim.time <- Sys.time() - start.time

sortedIFarray <- array(0, dim = c(N, nFFT/2+1, num.sim))
sortedAMParray <- sortedIFarray
sortedNUMarray <- sortedIFarray

sort.start <- Sys.time()
for(n in 1:N){
  for(f in 1:(nFFT/2+1)){
    sortedIFarray[n,f,] <- sort(IFarray[n,f,])
    sortedAMParray[n,f,] <- sort(AMParray[n,f,])
    sortedNUMarray[n,f,] <- sort(NUMarray[n,f,])
  }
}
sort.time <- Sys.time() - sort.start

sigma.hat <- sd(sortedIFarray[200,69,])
qqplot(sortedIFarray[200,69,], qnorm(seq(0,1, length.out = 1000), mean = 0, sd = sigma.hat), xlab = '',
       ylab = '')
#weird tails. probably not gaussian.

vtk2 <- rowSums( V^2 )/2
vdtk2 <- rowSums( Vdot^2 )/2

plot(sortedAMParray[200,69,], type = 'l', ylab = '', main = 'amplitude of |Z(t)|^2')
plot(sortedNUMarray[200,69,], type = 'l', ylab = '', main = 'numerator of instantaneous frequency')
plot(sortedIFarray[200,69,], type = 'l', ylab = '', main = 'instantaneous frequency')

qqplot(sortedAMParray[200,69,]/vtk2[200], qchisq(seq(0,1, length.out = 1000), 2),
       xlab = 'quantiles of |Z(t)|^2', ylab = expression(paste('quantiles of ', chi[2]^2)))
abline(a = 0, b = 1, col = 2)

freq <- sp.X$freq
f.ind <- which.min(abs(freq-f0))

plot(density(sortedAMParray[200,69,] / vtk2[200]), type = 'l')
plot(density(sortedNUMarray[200,69,] / sqrt(vtk2[200]*vdtk2[200])), type = 'l')
plot(density(sortedIFarray[200,69,]), type = 'l')



prod.dens <- function(x,sig1,sig2){
  besselK(abs(x)/(sig1*sig2),0)/(pi*sig1*sig2)
}
support <- seq(-6,6,length.out = 1000)
plot(support, prod.dens(support,1,1), type = 'l')
