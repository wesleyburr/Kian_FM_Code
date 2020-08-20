N <- 1024
x <- rnorm(N)
nw <- 4
k <- 2*nw - 1
v <- multitaper::dpss(N,k,nw,returnEigenvalues = FALSE)$v

#y <- matrix(rep(0, N*k), nrow = N, ncol = k)

#for( i in 1:k){
#  y[,i] <- fft( v[,i] * x )
#}

y <- mvfft( v*x )

#this recovers the mean, which for standard normal noise is 0
y0 <- y[1,]

z <- v %*% y0

z <- Re(z)
plot(x, type = "l")
lines(z, col = "red")

#this isn't even coming close to giving me the correct thing. Something must be wrong. Nope

f <- 6/1024
x1 <- sin(2*pi*f*( 0:(N-1) )) + x

#get Fourier frequency associated with f -> 6/1024, 7th one
#y <- matrix(rep(0, N*k), nrow = N, ncol = k)

#for( i in 1:k){
#  y[,i] <- fft( v[,i] * x1 )
#}
y <- mvfft(v * x1)
y1 <- y[7,]

z1 <- v %*% y1

for(t in 1:N){
  z1[t] <- z1[t] * exp(2*pi*1i*f*t) / Mod( z1[t] )
}
plot(x1, type='l')
lines(Re(z1),col='red',lwd=2)
lines(Im(z1),col='blue',lwd=2)

#getting something closer to correct. amplitude is off though and signal is still complex
#dividing by the modulus at each point gives a great reconstruction of the signal, though I don't
#know why


#now do this with zero-padding present
set.seed(1224)
N <- 1000
f <- 6/N


sig <- 100*sin( 2*pi*f* 0:(N-1) )
noiz <- rnorm(N, sd = sigmat)

x <- sig + noiz

NW <- 6
K <- 2*NW - 2

v <- multitaper::dpss(N,K,NW,returnEigenvalues = FALSE)$v

#need to zero-pad a lot to get close to the correct frequency and get a decent fit
pad.length <- 6^ceiling( log(2*N,6) )

tapered <- matrix(rep(0,pad.length*K), nrow=pad.length,ncol=K)
for(k in 1:K){
  tapered[,k] <- c( v[,k]*x, rep(0,pad.length-N) )
}



y <- mvfft(tapered)

h <- floor(pad.length*f + 1) #find index corresponding to relevant Fourier frequency
# pad.length*f +1 is actually closer to 48 than 47, and (48-1)/7776 is closer to f=0.006 than
# (47-1)/7776, but using h=47 gives a better estimate of the signal.

yh <- y[h,]

Z <- v %*% yh

#estimated complex mean value for signal at specified frequency and incorporated into
#the inverse. Had a couple ideas but the one that seems to work best is 2|mu|. No idea why.
# I guess for a pure sinusoid with amplitude A, the complex mean value is mu = A/2i
U <- apply(v, MARGIN=2, FUN=sum)
mu <- ( t(U)%*%U )^-1 %*% t(U) %*% yh

for(t in 0:(N-1)){
  Z[t+1] <- 2*Mod(mu)*Z[t+1] * exp(1i*2*pi*(h-1)/pad.length*t) / Mod(Z[t+1])
}

plot(x, type='l')
lines(Re(Z), col='red')
lines(Im(Z), col='blue')


sum( (Re(Z)-x)^2 )
sum( (Re(Z)-x)^2 ) / N

sum( ( Re(Z) - (x-noiz) )^2 )
sum( ( Re(Z) - (x-noiz) )^2 ) / N



####################################
#want to try this stuff on real data

data <- read.csv("SN_y_tot_V2.0.csv", header=FALSE, sep=";")
time <- data[,1]
sunspot <- data[,2]
sunspot <- sunspot - mean(sunspot)

N <- length(sunspot)
NW <- 5
K <- 2*NW - 1
nFFT <- 2^ceiling(log2(10*N))

v <- multitaper::dpss(N,K,NW, returnEigenvalues = FALSE)$v

spec <- multitaper::spec.mtm(ts(sunspot), nw=NW, k=K, nFFT = nFFT, dpssIN = v, Ftest = TRUE, 
                             deltat = 1, dtUnits = 'year')
plot(spec, Ftest = TRUE)
f <- spec$freq[ which( spec$mtm$Ftest == max(spec$mtm$Ftest) ) ]

h <- nFFT*f + 1

tapered <- rbind(v * sunspot, matrix(rep(0, (nFFT - N)*K), ncol = K))

y <- mvfft(tapered)

Yh <- y[h,]

Z <- v %*% Yh

U <- apply(v, MARGIN=2, FUN=sum)
mu <- ( t(U)%*%U )^-1 %*% t(U) %*% Yh

Z.tilde <- numeric(N)
for(t in 0:(N-1)){
  Z.tilde[t+1] <- 2*Mod(mu)*Z[t+1] * exp(1i*2*pi*f*t) / Mod(Z[t+1])
}

plot(time, sunspot, type='l')
lines(time, Re(Z.tilde), col=2)

plot(Re(Z.tilde), type='l', xlim = c(0,50))


zoom <- multitaper::dropFreqs(spec, 0.05,0.15)
f2 <- zoom$freq[ which( zoom$mtm$Ftest == max(zoom$mtm$Ftest) ) ]

h2 <- nFFT*f2 + 1

Yh2 <- y[h2,]

Z2 <- v %*% Yh2

mu2 <- ( t(U)%*%U )^-1 %*% t(U) %*% Yh2
Z.tilde2 <- numeric(N)
for(t in 0:(N-1)){
  Z.tilde2[t+1] <- 2*Mod(mu2)*Z2[t+1] * exp(1i*2*pi*f2*t) / Mod(Z2[t+1])
}

plot(time, sunspot, type = 'l')
lines(time, Re(Z.tilde2), col=2)

Z.hat <- Z.tilde + Z.tilde2
plot(time,sunspot, type='l')
lines(time, Re(Z.hat), col=2)
#lines(time, Im(Z.hat), col=4)

plot(spec, Ftest = TRUE)
abline(h = qf(0.99,2,2*(K-1)), lty=2,col=2)

zoom3 <- multitaper::dropFreqs(spec, 0,0.05)
f3 <- zoom3$freq[ which( zoom3$mtm$Ftest == max(zoom3$mtm$Ftest) ) ]

h3 <- nFFT*f3 + 1

Yh3 <- y[h3,]

Z3 <- v %*% Yh3

mu3 <- ( t(U)%*%U )^-1 %*% t(U) %*% Yh3
Z.tilde3 <- numeric(N)
for(t in 0:(N-1)){
  Z.tilde3[t+1] <- 2*Mod(mu3)*Z3[t+1] * exp(1i*2*pi*f3*t) / Mod(Z3[t+1])
}

zoom4 <- multitaper::dropFreqs(spec, 0.25,0.3)
f4 <- zoom4$freq[ which( zoom4$mtm$Ftest == max(zoom4$mtm$Ftest) ) ]

h4 <- nFFT*f4 + 1

Yh4 <- y[h4,]

Z4 <- v %*% Yh4

mu4 <- ( t(U)%*%U )^-1 %*% t(U) %*% Yh4
Z.tilde4 <- numeric(N)
for(t in 0:(N-1)){
  Z.tilde4[t+1] <- 2*Mod(mu4)*Z4[t+1] * exp(1i*2*pi*f4*t) / Mod(Z4[t+1])
}

Z.hat <- Z.tilde + Z.tilde2 + Z.tilde3 + Z.tilde4
plot(time,sunspot,type='l')
lines(time,Re(Z.hat), col=2)

#########################################################
#reconstruct out of most significant frequencies
#########################################################
data <- read.csv("SN_y_tot_V2.0.csv", header=FALSE, sep=";")
time <- data[,1]
sunspot <- data[,2]


topM <- function(x, M=1){
  return( order(x, decreasing=TRUE)[1:M] )
}

m <- 13

reconstruct <- function(time.series, num = 1, NW = 4, K = 7){
  
  
  x <- time.series - mean(time.series) 
  N <- length(x)
  nFFT <- 2^ceiling(log2(10*N))
  v <- multitaper::dpss(N,K,NW, returnEigenvalues = FALSE)$v
  
  
  xbar <- v %*% t(v) %*% x
  temp.x <- x - xbar
  
  tapered <- rbind(v * temp.x[,1], matrix(rep(0, (nFFT - N)*K), ncol = K))
  y <- mvfft(tapered)
  
  spec <- multitaper::spec.mtm(temp.x, nw=NW, k=K, nFFT = nFFT, dpssIN = v, Ftest = TRUE, 
                               deltat = 1,  plot=FALSE, centre='none')
  
  f.mesh <- seq(0,0.5, by = 2/N)
  L <- length(f.mesh)
  F.max <- NULL
  F.indices <- NULL
  for(l in 1:(L-1)){
    temp.spec <- multitaper::dropFreqs(spec,f.mesh[l],f.mesh[l+1])
    F.max <- c(F.max, max(temp.spec$mtm$Ftest))
    F.indices <- c(F.indices, which(spec$mtm$Ftest == F.max[l]) )
  }
  peak.freqs <- spec$freq[F.indices]
  topF <- topM(spec$mtm$Ftest[F.indices], num)
  f <- peak.freqs[topF]

  
  h <- nFFT * f + 1
  Yf <- y[h,]
  U <- apply(v, MARGIN=2, FUN=sum)
  if(num == 1){
    Z <- v %*% Yf
    mu <- ( t(U)%*%U )^-1 %*% t(U) %*% Yf
  } else {
    Z <- v %*% t(Yf)
    mu <- ( t(U)%*%U )^-1 %*% t(U) %*% t(Yf)
  }
  
  
  

  Z.tilde <- matrix(rep(0,num*N), nrow=N)
  for(j in 1:num){
    for(t in 0:(N-1)){
      #Z.tilde[t+1,j] <-2*Mod(mu[j])* cos(2*pi*f[j]*t + Arg(mu[j]))
      Z.tilde[t+1,j] <-Z[t+1,j] * exp(1i*2*pi*f[j]*t) 
      
    }
  }

  Z.hat <- apply(Z.tilde, MARGIN = 1, FUN = sum)
  Z.hat <- Z.hat  + xbar + mean(time.series) + 1i * ( xbar + mean(time.series) )
  #shift the imaginary part up in case it is used for something 
  plot(time.series , type='l')
  lines(Re(Z.hat), col=2)
  #lines(Im(Z.hat),col=4)
  
  data <- list( Z.hat = Z.hat, Z.tilde = Z.tilde, freq = f)
  return(data)
}
system.time({
  Recon <- reconstruct(sunspot,m,NW = 8, K=15)
})
MSE <- mean( (Re(Recon$Z.hat) - sunspot)^2 )
x <- sunspot - Re( Recon$Z.hat )
plot(x)
qqnorm(x/sd(x))
abline(a=0,b=1)
multitaper::spec.mtm(ts(x),nw=8,k=15)
