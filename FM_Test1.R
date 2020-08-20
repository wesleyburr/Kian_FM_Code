
# simple sinusoid with no FM
set.seed(1250)
t <- 1:10000
N <- length(t)
nw <- 5
w <- nw / N
k <- 2*nw - 1
f <- 0.3
amp <- 1
phi <- 0
alpha <- 1/N

sig <- amp*sin(2*pi*f*t + phi)
noise <- rnorm(N)
z <- sig + noise
z <- ts(z)

s <- multitaper::spec.mtm(z, k=k, nw=nw, Ftest = TRUE)
plot(s, Ftest = TRUE)
abline(h=qf(1-alpha,2,2*(k-1)), lty=2,col=2)




#single line component
#simple FM -W -> W -> -W by line segments
set.seed(11)
t <- 1:10000
N <- length(t)
nw <- 12
w <- nw/N
k <- 2*nw - 1
amp <- 1
phase <- 0
f <- 0.3


f.mod <- function(t) {
  if(0 <= t && t <= N/2){
    return(-w*(1-2*t/N) + 2*t*w/N)
  } else if (N/2 < t && t <= N){
    return( w*2*t/N + 3*w*(1 - 2*t/N) )
  }
  
}

sig <- amp*sin(2*pi*(f + f.mod(t))*t + phase)
noise <- rnorm(N)

z <- sig + noise
z <- ts(z)
plot(z, xlim=c(0,100))

padding <- 6^ceiling(log(10*N,6))
taper <- multitaper::dpss(N, k=k-1, nw=nw, returnEigenvalues =TRUE)

s <- multitaper::spec.mtm(z,nw=nw,k=k-1, Ftest = TRUE, dpssIN = taper$v, nFFT=padding)
plot(s, Ftest = TRUE)

s2 <- multitaper::dropFreqs(s, 0.28, 0.32)
plot(s2)
plot(s2, Ftest=TRUE)

max.f <- s$freq[which(s$mtm$Ftest == max(s$mtm$Ftest))]

#Dave Riegert gave me this block function
block <- function(n, blockSize, overlap){ 
  increment <- ceiling(blockSize * (1-overlap)) 
  sectIdx <- seq(1, n-blockSize+1, by=increment) 
  numSect <- length(sectIdx) 
  list(startIdx = sectIdx, increment = increment, nBlock = numSect, blockSize = blockSize)
}

blok <- block(N, N/3, 0.4)

nBlock <- blok$nBlock
startIdx <- blok$startIdx
blockSize <- blok$blockSize

#system.time({
spec.list <- list()
spec.list2 <- list()
for( i in 1:nBlock ){
  spec.list[[i]] <- multitaper::spec.mtm(ts( z[ startIdx[i]:( startIdx[i] + blockSize - 1) ]), 
                                         k=k-1,nw=nw,Ftest = TRUE, dpssIN = taper$v,nFFT=padding)
  plot(spec.list[[i]], Ftest = TRUE)
  abline(h=qf(0.99,2,2*( (k-1) - 1)), col=2, lty=2)
    
  spec.list2[[i]] <- multitaper::dropFreqs(spec.list[[i]], f - 0.02, f + 0.02)
  plot(spec.list2[[i]])
  plot(spec.list2[[i]], Ftest=TRUE)
  abline(h=qf(0.99,2,2*( (k-1) - 1)), col=2, lty=2)
  
}
#})
  


#quadratic FM with a single line component. Again, -W to W to -W

#time and bandwidth parameters
N <- 10000
t <- 0:(N-1)
nw <- 7
k <- 2*nw - 1
w <- nw / N

#signal parameters
A <- 1
P <- 0
f <- 0.3

#FM parameters
b <- 8*w/(N-1)

quadratic.fm <- function(t) -b*t^2 / (N-1) + b*t - w

#sig <- A * sin( 2*pi*( f + quadratic.fm(t) ) + P)
sig <- A * sin( 2*pi*( f + quadratic.fm(t) )*t + P)

set.seed(420)
noise <- rnorm(N)

x <- sig + noise
x <- ts(x)

plot(t,quadratic.fm(t), type='l')
plot(x, xlim=c(0,200))

taper <- multitaper::dpss(N,nw,k,returnEigenvalues = TRUE)
padding <- 6^ceiling(log(10*N,6))


s.deg2 <- multitaper::spec.mtm(x,nw=nw,k=k,dpssIN=taper$v,Ftest=TRUE,nFFT=padding)

Fs <- c(qf(0.99,2,2*(k-1)),qf(0.999,2,2*(k-1)),qf(0.9999,2,2*(k-1)) )
plot(s.deg2, Ftest=TRUE)
abline(h=Fs,lty=2,col=2)
text(x=c(0.3,0.3,0.3),y=c(7,11,16),labels=c("0.99","0.999","0.9999"))

s.deg2.zoom <- multitaper::dropFreqs(s.deg2,f - 0.02,f + 0.02)
plot(s.deg2.zoom)
plot(s.deg2.zoom, Ftest = TRUE)

#check the quadratic over different blocks as before. also figure out signal to noise ratio and 
#how to generate the signal

spec.list3 <- list()
spec.list4 <- list()
for( i in 1:nBlock ){
  spec.list3[[i]] <- multitaper::spec.mtm(ts( x[ startIdx[i]:( startIdx[i] + blockSize - 1) ]), 
                                         k=k-1,nw=nw,Ftest = TRUE, dpssIN = taper$v,nFFT=padding)
  plot(spec.list3[[i]], Ftest = TRUE)
  abline(h=Fs, col=2, lty=2)
  
  spec.list4[[i]] <- multitaper::dropFreqs(spec.list3[[i]], f - 0.02, f + 0.02)
  plot(spec.list4[[i]])
  plot(spec.list4[[i]], Ftest=TRUE)
  abline(h=Fs, col=2, lty=2)
  
}




#more than one line component. center frequencies separated by more than W
#time and bandwidth parameters
N <- 10000
t <- 0:N
nw <- 7
k <- 2*nw - 1
w <- nw / N

#signal parameters
A1 <- 3
A2 <- 4
P <- 0
f1 <- 0.3
f2 <- 0.15

#FM parameters
b <- 8*w/N

quadratic.fm <- function(t) -b*t^2 / N + b*t - w

sig1 <- A1 * sin( 2*pi*( f1 + quadratic.fm(t) )*t + P)
sig2 <- A2 * sin( 2*pi*( f2 + quadratic.fm(t) )*t + P)

set.seed(428)
noiz <- rnorm(N+1)

z <- ts( sig1 + sig2 + noiz)

taper <- multitaper::dpss(N,nw,k,returnEigenvalues = TRUE)
padding <- 6^ceiling(log(10*N,6))

strum <- multitaper::spec.mtm(z,nw=nw,k=k,dpssIN=taper$v,Ftest=TRUE,nFFT=padding)
plot(strum, Ftest=TRUE)

#what if the center frequencies are less than W apart?
f2 <- f1 - w/2
sig2 <- A2 * sin( 2*pi*( f2 + quadratic.fm(t) )*t + P)

set.seed(433)
noiz <- rnorm(N+1)

z <- ts( sig1 + sig2 + noiz)

strum <- multitaper::spec.mtm(z,nw=nw,k=k,dpssIN=taper$v,Ftest=TRUE,nFFT=padding)
plot(strum, Ftest=TRUE)

spec.list3 <- list()
spec.list4 <- list()
for( i in 1:nBlock ){
  spec.list3[[i]] <- multitaper::spec.mtm( z[ startIdx[i]:( startIdx[i] + blockSize - 1) ], 
                                          k=k-1,nw=nw,Ftest = TRUE, dpssIN = taper$v,nFFT=padding)
  plot(spec.list3[[i]], Ftest = TRUE)
  abline(h=c( qf(0.99,2,2*(k-1) ) ,qf(0.999,2,2*(k-1) )), col=2, lty=2)
  
  spec.list4[[i]] <- multitaper::dropFreqs(spec.list3[[i]], f - 0.02, f + 0.02)
  plot(spec.list4[[i]])
  plot(spec.list4[[i]], Ftest=TRUE)
  abline(h=c( qf(0.99,2,2*(k-1) ) ,qf(0.999,2,2*(k-1) )), col=2, lty=2)
  
}





#now to see how much FM it takes to ruin your day
N <- 10000
t <- 0:(N-1)
nw <- 7
k <- 2*nw - 1
w <- nw / N

#signal parameters
A <- 1
P <- 0
f <- 0.3

#FM parameters. prop = proportion of half bandwidth involved in FM. FM goes from 
#f-pW to f+pW. b is just a coefficient determined by creating a suitable quadratic
prop <- 0.69 
pW <- prop*w
b <- 8*pW/(N-1)


quadratic.fm <- function(t) -b*t^2 / (N-1) + b*t - pW

#sig <- A * sin( 2*pi*( f + quadratic.fm(t) ) + P)
sig <- A * sin( 2*pi*( f + quadratic.fm(t) )*t + P)

plot(sig[1:200],type='l')
lines(1:200,sin(2*pi*f*(1:200)),col=2)


set.seed(420)
noise <- rnorm(N)

x <- sig + noise
x <- ts(x)

plot(t,quadratic.fm(t), type='l')
plot(x, xlim=c(0,200))

taper <- multitaper::dpss(N,nw,k,returnEigenvalues = TRUE)
padding <- 6^ceiling(log(10*N,6))

s.deg2 <- multitaper::spec.mtm(x,nw=nw,k=k,dpssIN=taper$v,Ftest=TRUE,nFFT=padding)

Fs <- c(qf(0.99,2,2*(k-1)),qf(0.999,2,2*(k-1)),qf(0.9999,2,2*(k-1)) )
plot(s.deg2, Ftest=TRUE)
abline(h=Fs,lty=2,col=2)
text(x=c(0.3,0.3,0.3),y=c(7,11,16),labels=c("0.99","0.999","0.9999"))

s.deg2.zoom <- multitaper::dropFreqs(s.deg2,f - 0.02,f + 0.02)
plot(s.deg2.zoom)
plot(s.deg2.zoom, Ftest = TRUE)
abline(h=Fs,lty=2,col=2)


block <- function(n, blockSize, overlap){ 
  increment <- ceiling(blockSize * (1-overlap)) 
  sectIdx <- seq(1, n-blockSize+1, by=increment) 
  numSect <- length(sectIdx) 
  list(startIdx = sectIdx, increment = increment, nBlock = numSect, blockSize = blockSize)
}

blok <- block(N, N/4, 0.5)

nBlock <- blok$nBlock
startIdx <- blok$startIdx
blockSize <- blok$blockSize



spec.list3 <- list()
spec.list4 <- list()
for( i in 1:nBlock ){
  spec.list3[[i]] <- multitaper::spec.mtm(ts( x[ startIdx[i]:( startIdx[i] + blockSize - 1) ]), 
                                          k=k-1,nw=nw,Ftest = TRUE, dpssIN = taper$v,nFFT=padding)
  plot(spec.list3[[i]], Ftest = TRUE)
  abline(h=Fs, col=2, lty=2)
  
  spec.list4[[i]] <- multitaper::dropFreqs(spec.list3[[i]], f - 0.02, f + 0.02)
  plot(spec.list4[[i]])
  plot(spec.list4[[i]], Ftest=TRUE)
  abline(h=Fs, col=2, lty=2)
  
}
# early blocks pick up the center frequency, but later blocks completely miss it. 
# lower p, around 0.2 it starts being able to pick up the center frequency in later blocks
# and even almost in the original F test with no segmenting. The data is extremely contrived
# though and I am specifically looking for a predetermined frequency. 




N <- 10000
t <- 0:(N-1)
nw <- 10
k <- 2*nw - 1
w <- nw / N

#signal parameters
A <- 1
P <- 0
f <- 0.3

#FM parameters. prop = proportion of half bandwidth involved in FM. FM goes from 
#f-pW to f+pW. b is just a coefficient determined by creating a suitable quadratic
prop <- 0.69 
pW <- prop*w
b <- 8*pW/(N-1)

quadratic.fm <- function(t) -b*t^2 / (N-1) + b*t - pW
sig <- A * cos(2*pi*( f + quadratic.fm(t) ) * t + P ) + A * sin(2*pi*( f + quadratic.fm(t) ) * t + P )
noiz <- rnorm(N) + 1i * rnorm(N)

X <- sig + noiz
v <- dpss(N,k,nw)
v <- v$v

nFFT <- 2^ceiling(log2(2*N))

tapered.X <- rbind(v*X, matrix(data=0,nrow = nFFT-N, ncol=k))
tt <- 0:(nFFT-1)
Y <- mvfft(tapered.X)

d <- 3
phi <- -pW * cos( 2*pi*d/(2*(N-1))*t )
Phi <- f + phi
plot(t,Phi, type = 'l')
points(c(0, (N-1)/3, 2*(N-1)/3, N-1), c(f-pW,f+pW,f-pW,f+pW))
#problem: don't have local mins and maxs at interior f+-pW points
#solution? use derivative constraints at interior points instead. this doesn't work

set.seed(111)
ndata <- 200
f0 <- 0.3
sig <- cos((f+Fmod)*0:(ndata-1)) + 1i*sin((f0+Fmod)*0:(ndata-1)) 
noiz <- rnorm(ndata, sd = sigmat) + 1i * rnorm(ndata, sd = sigmat)
CDN <- sig + noiz
