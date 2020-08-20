#Moving stuff from the fortran code to R. I can't remember what everything is
#or how it all works

#Set parameters. Need to remember what all these mean
ndmax <- 200
ndata <- ndmax #size of simulated data?
ndstak <- 15600
nDrop <- 40
nppb <- ndmax
nord <- 10 #I think this is K. Actually this is definitely K
mxdeg <- 6 #max degree of polynomial modulation?
nHesmn <- 30
nHesmx <- 900
mxfit <- 5
nfbnds <- 5
mxcoef <- 1+mxfit
#CLEANED UP TO HERE
#real-valued variables
xd <- numeric(ndmax) #not sure
t <- numeric(ndmax) #time?
x <- numeric(ndmax) # real part of signal?
#dw <- matrix(rep(0,nppb*nord),nrow=nppb,ncol=nord) #data window. Should this be a matrix?
#theta <- numeric(nord) #eigenvalues?
dwp <- matrix(rep(0,nppb*nord),nrow=nppb,ncol=nord) # Slepians, but adjusted for polynomial passing
#ev <- numeric(nord) # not sure of these. eigenvalues probably
evp<- numeric(nord)
evdx <- numeric(nord)
Fmod <- numeric(ndmax) #something with frequency modulation
PhMod <- numeric(ndmax) #something with phase modulation
y <- numeric(ndmax) #imaginary part of signal?
amp <- numeric(ndmax) #amplitude?
PhD <- numeric(ndmax) #phase something?
Ukz(nord,0:mxdeg) #what type of object is this? matrix, but the column indices go from 0 to mxdeg
#Ukz is probably the U_k at 0. But then what is Ukz0? 0th order?
base <- numeric(ndmax)
Hn(0:mxdeg) #???
Ukz0 <- numeric(nord) #???
Poly(nppb,0:mxdeg) #maybe the polynomial matrix - **THIS IS G IN THE PAPER**
SSQUkz(0:mxdeg) #sum of squares of columns of Ukz?
yd <- numeric(ndata)
x1 <- numeric(ndata) #real and imaginary parts of something
y1 <- numeric(ndata)
x2 <- numeric(ndata)
y2 <- numeric(ndata)
tc <- numeric(ndmax) #centered time?
A1 <- numeric(ndata)
P1 <- numeric(ndata) #amplitude and phase of something?
A2 <- numeric(ndata)
P2 <- numeric(ndata)
Pref <- numeric(ndata)
Fr1 <- numeric(ndata)
FrP <- numeric(ndata)
FrC <- numeric(ndata)
FrA <- numeric(ndata) #frequency as derivative of phase line 15 in '09 paper
Fbounds <- numeric(nfbnds)
Xdot <- numeric(ndmax)
Ydot <- numeric(ndmax)
HePol(ndmax,0:mxdeg) #some more vectors with different indices
HeCo(0:mxdeg)
HeCHat(0:mxdeg)
Xmat(ndmax,0:mxdeg) #and a matrix
Fit <- numeric(ndmax)
Rsd <- numeric(ndmax)
Fcoef <- numeric(nord)
FPcoef(0:mxdeg)
FRsd <- numeric(nord)
Cdpss2 <- numeric(nord)
Csp2(0:mxdeg)
SSQFrsd(0:mxdeg)

#Now complex-valued variables

cd <- complex(ndata) #complex demodulate? complex data?
yout <- complex(ndata)
yk <- complex(nord) #eigencoefficients
pc(0:mxdeg) #polynomial coefficients?
CDN <- complex(ndata)
cdotr <- complex(1) #do I even need to bother with this? Likely not
zz <- complex(1)
Fit1 <- complex(ndata)
FitP <- complex(ndata) #something with the polynomials
FitC <- complex(ndata) # cumulative version of above thing?
pca(0:mxdeg)
Cum <- complex(mxdeg) #sum something ( ͡° ͜ʖ ͡°)
Zdot <- complex(ndata) #inverted data I think. Line 3 in paper

#Don't think I need these character variables
#character pname*16, Gfile*16, ttl*96, tlbl*64, stl*96 

#also don't think I need this
#double precision dstak(ndstak)
#common/cstak/dstak

#probably don't need this
#data  pname/"CspcdmFMH3.f"/, Gfile/"G-CspcdmFMH3g.ps"/
#  1, dt/1.0/, tlbl/"Sample Number"/
#  2, ttl/"Quadratic FM test"/
#  4,  js1/131/, js2/17107/
#  3, f1mn/-1.499/, f1mx/1.499/, norm/1/

dt <- 1.0
js1 <- 131
js2 <- 17107
f1mn <- -1.499 #these are just +- nyquist for plotting?
f1mx <- 1.499
norm <- 1

#couple of functions. one computes |z|^2 , other computes the phase
# cabs2 <- function(zz) Re(zz)^2 + Im(zz)^2 #  this is just Mod(zz)^2
phse <- function(zz) { RtoD*atan2(Im(zz), Re(zz)) }

#need to reorder some of this stuff. this function is defined before RtoD is assigned a value

#call istkin(ndstak,4)
#Data file Input
#call prid(pname, "Test projection filter")
#call Gbegin(pname, Gfile)
tpi <- 2*pi
RtoD <- 360/tpi           #radians to degrees
DtoR <- tpi/360           #degrees to radians
#call ranset(js1,js2)  # sets the seed for the RNG
set.seed(sample(seq(js1, js2, 1)))

nCalLo <- nDrop + 1
nCompare <- ndata - 2*nDrop
#CX      TNW = float(nord+2)
TNW <- nord+1 #time bandwidth product
W <- TNW/(2*ndata) #double check this
#write(*,*) " 2 N W", TNW," nord",nord," W",W
fbw <- 0.8 #fraction of bandwidth for FM apparently
#write(*,*) "Fraction of bandwidth for FM", fbw
#                       Frequency bounds for plots. Probably don't need these
#Fbounds(1) = -W
#Fbounds(2) = -fbw
#Fbounds(3) = 0.0
#Fbounds(4) = fbw
#Fbounds(5) = W

stm11 <- 2.0/(ndata-1)
# Pcum <- 0.0
Shift <- 0.0
#write(*,*) " Shift",Shift                        #DELETE ALL THIS
#write(stl,"('K',i3,' mxdeg',i2,' W',1pe11.5,' shift',0pf4.0
#      1 ,'  ')") nord,mxdeg,W,shift
#call blank1(stl,stl)
#call Gbplin(stl)
cent <- (1+ndata)/2.0
scl <- 1.0/cent

# Rewrite, using vectorization
t <- 0:(ndata-1)
tc <- (t - cent)*scl
tt <- (t + Shift) * stm11 - 1.0  # this runs from -1 to 1
Fmod <- W * fbw * (1.0 - 2.0 * tt^2)  # quadratic frequency modulation
for(j in 1:ndata) {
  Pcum <- sum(Fmod[1:j] * dt * tpi)
  PhMod[j] <- Pcum
  cd[j] <- cos(Pcum) + 1i * sin(Pcum)
}
Pref <- RtoD * PhMod

#for (n in 1:ndata){
#  t[n] <- n-1
#  tc[n] <- (t[n]-cent)*scl
#  tt <- ( t[n] + Shift)*stm11 - 1.0
#  Fmod[n] <- W * fbw *( 1.0 - 2.0*(tt^2) )
#  Pcum <- Pcum + Fmod[n]* dt * tpi
#  PhMod[n] <- Pcum
#  Pref[n] <- RtoD*PhMod[n]
#  x[n] <- cos(Pcum)
#  y[n] <- sin(Pcum)
#  cd[n] <-  x[n] + 1i*y[n] 
#}

#plotting stuff. replace?
#call sdplot(t,Fmod,ndata,ttl,tlbl,"Frequency", 0.5)
#call Ghline(Fbounds,nfbnds,2, 1.0,"[12 12] 4","red")
#call sdplot(t,Pref,ndata,ttl,tlbl,"Phase",0.5)
#call sdplos(t,x,ndata,ttl,tlbl,"Re{CD}",0.8
           # 1, t(1),t(ndata), -1.0, 1.0)
#call Gline(t,y,ndata, 1,2, 2.0,"[4 4] 0","red")
 
# this matters
#call dpss(dw,ndmax,nppb,nord,W,theta)      does this compute the dpss and put them in dw?
            #U_k (0), rms ripple            what is theta doing? I forget what nppb is
dw <- multitaper::dpss(ndata, nord, ndata * W)
# grab the evals, dw$eigen
ev <- dw$eigen
dw <- dw$v


Ukz0 <- apply(dw, MAR = 2, FUN = sum) 
ssq <- sum(Ukz0^2)

#ssq <- 0.0
#for(k in 1:nord){
#  Ukz0[k] <- ssum(nppb,dw[1,k],1) #what is ssum?
#  ssq <- ssq + Ukz0[k]^2
#}

#write(*,*) "Ukz0",Ukz0," ssq Ukz0", ssq
#call dpssev(dw,ndmax,nppb,nord,W,ev,evp)   # yes, computes eigenvalues - probably doesn't matter (ev set above)
#write(*,*) "ev", ev," evp", evp
#call dpssp(dw,ndmax,nppb,nord,W,evdx,dwp)  # MATTERS - this is another function that computes the dpss-polynomials

#don't know if I need this loop. It's just writing eigenvalues and stuff
#do 48 k = 1, nord
#48 write(*,*) "k,ev,evx,evp,theta",k,ev(k),evdx(k), evp(k),theta(k)

#no clue what this stuff is doing
nprint <- 2
Trmx <- 0.0 #this is the only place that Trmx appears
HeSmx <- 0.90
# ------------------------
#call dpsshp(dw,nppb,nppb,nord,Poly,nppb,nppb,mxdeg    # MATTERS - NEED THIS
#            HeSmx,norm,Ukz,Hn,nprint)
# call cm4pp(Fmod,ndata,ave,var,"Modulating Frequency")  # not important

######################################################
snr <- 5.0 											#signal to noise ratio
sigmat2 <- ndata/(2*snr*nord)  # looks like a variance, but of what?   - WHY THIS VARIANCE?
sigmat <- sqrt(sigmat2)

#likely don't need this write
#write(*,*) "snr", snr," sigma each component", sigmat

#loop looks like it's making noise, accumulating the squared noise, and adding it to the signal
#for(n in 1:ndata){
#  x1[n] <- rnorm(n)*sigmat      #this is probably generating random normal. or using an 
#  y1[n] <- rnorm(n+17)*sigmat #??? already generated random normal vector? indices are weird. seed?
#  pwrn <- pwrn + x1[n]^2 + y1[n]^2 # PoWeR in Noise
#  x1[n] <- x[n] + x1[n]  # x is data, x1 is current noies
#  y1[n] <- y[n] + y1[n]  # y is data, y1 is current noise
#  CDN[n] <- x1[n] + 1i*y1[n]  # CDN is the combination, additively
#}
CDN <- cd + rnorm(ndata, sd = sigmat) + rnorm(ndata, sd = sigmat) * 1i
pwrn <- Mod(CDN - cd)^2

#more writing and plotting (?) 
#write(*,*) "Pwrn",pwrn
#call sdplot(t,x1,ndata,ttl,tlbl,"Re + noise", 0.9)
#call sdplot(t,y1,ndata,ttl,tlbl,"Im + noise", 0.9)
#ssqyk <- 0.0        #sum of squares of eigencoefficients

# ------------------------------------------------------
# START OF THE INTERESTING STUFF
yk <- t(CDN) %*% dw
yksq <- Mod(yk)^2
ssqyk <- sum(yksq)
#for(k in 1:nord){
  #yk[k] <- cdotr(nppb,CDN,1,dw[1,k],1) #cdotr is complex(c) times real(r)
#  yk[k] <- sum(CDN * dw[,k]) #yk at frequency zero?
#  yksq <- Mod(yk[k])^2
#  ssqyk <- ssqyk + yksq
  #write(*,"(i3,1p8e14.6)") k,yk(k),yksq,ssqyk,theta(k),ev(k), evp(k)
#}
#do I need these?
#call setc(ndata, (0.0, 0.0), Fit1)
#call setc(ndata, (0.0, 0.0), Zdot)
#call setc(ndata, (0.0, 0.0), FitP)
#call setc(ndata, (0.0, 0.0), FitC)


Fit1 <- dw %*% yk
Zdot <- dwp %*% yk

amp <- Mod(Fit1)
PhD <- phse(Fit1)
x1 <- Re(Fit1)
y1 <- Im(Fit1)
Xdot <- Re(Zdot)
Ydot <- Im(Zdot)
FrA <- ( x1*Ydot - Xdot*y1 ) / (tpi*amp^2)

#for(n in 1:ndata){
#  for(k in 1:nord){
#    Fit1[n] <- Fit1[n] + yk[k]*dw[n,k]   #or is this the standard inverse?
#    Zdot[n] <- Zdot[n] + yk[k]*dwp[n,k] #inverse of data. But is dwp even defined yet?
#  }
#  amp[n] <- abs(Fit1[n]) 
#  PhD[n] <- phse(Fit1[n])
#  x1[n] <-  Re(Fit1[n])
#  y1[n] <- Im(Fit1[n])
#  Xdot[n] <-  Re(Zdot[n])
#  Ydot[n] <- Im(Zdot[n])
#  FrA[n] <- ( x1[n]*Ydot[n] - Xdot[n]*y1[n])/(tpi*(amp[n]^2)) # frequency line 15
#}

system.time({
curv.filter <- c(-2,1,rep(0,ndata-3),1)
FaC2 <- Re( fft( fft(curv.filter) * fft(FrA), inverse = TRUE) / ndata )
FaC2 <- sum( FaC2[2:(ndata-1)]^2 )
}
)



#system.time({
#FaC2 <- 0.0
#for(n in 2:(ndata-1) ){
#  FaC2 <- FaC2 +  (FrA[n-1]-2.0*FrA[n]+FrA[n+1])^2
#}
#})
#write(*,*) "Curvature of FrA", FaC2

#ssqFc <- 0.0 #sum of squares of whatever Fcoef is
Fcoef <- t(FrA) %*% dw
ssqFc <- sum( Fcoef^2 )
#for(k in 1:nord){
#  Fcoef[k] <- sdot(ndata, FrA,1, dw[1,k],1)  #sdot?
#  ssqFc <- ssqFc + Fcoef[k]^2
# write(*,*) "k, Fcoef, ssqFc", k,Fcoef(k), ssqFc
#}

#likely don't need all this plotting stuff. Maybe replace it?
#call sdplos(t,x1,ndata,"Fit with dpss expansion"
#            1, tlbl,"Fit1",0.8,t(1),t(ndata),f1mn,f1mx)
#call Gline(t,y1,ndata,1,2, 2.,"[4 4] 0","red")
#call a2mnmx(Xdot,Ydot,ndata,XYDmn,XYDmx)
#call sdplos(t,Xdot,ndata,"Analytic Derivative"
#            1,tlbl,"Xdot,Ydot", 0.9,t(1),t(ndata),XYDmn,XYDmx)
#call Gline(t,Ydot,ndata, 1,2, 2.,"[5 5] 0","red")
#call sdplot(t,amp,ndata,"Amp, dpss expansion"
#            1,tlbl,"Dpss AMP", 0.8)
#call adstoa(amp,yd,ndata, -1.0)
#call cm4pp(yd,ndata,ave, var,"DPSS amplitude - 1.0")
#call cm4pp(yd(nCalLo),nCompare,ave,var,"DPSS Amp -1., Center")
#call sphsed(PhD,ndata)
#call sdplot(t,PhD,ndata,"Phse, dpss expansion"
#            1,tlbl,"Dpss Phase", 0.8)
#call Gline(t,Pref,ndata, 1,2, 2.,"[4 4] 0","blue")
yd <- PhD - Pref
#for(n in 1:ndata){
#  yd[n] <- PhD[n]-Pref[n]  #not sure. something with phase
#}

#possibly don't need
#call cm4pp(yd,ndata,ave, var,"DPSS Phse - Ref")
#call cm4pp(yd(nCalLo),nCompare,ave,var,"DPSS Phse - Ref, Cntr")
#write(*,*) " ssqyk", ssqyk

#this part might have to stay as is since it updates the eigencoefficients each step
ssqpc <- 0.0
for (jd in 0:mxdeg){ #mxdeg = max degree? pc must be the polynomial part
  pc[jd+1] <- sum( yk * Ukz[,jd+1]) #cdotr again
  ssqpc <- ssqpc + Mod(pc[jd+1])^2
  #write(*,*) "pc",jd, pc[jd]," ssqpc",ssqpc
  ssqR <- 0.0
  for(k in 1:nord){
    yk[k] <- yk[k] - Ukz[k,jd+1]*pc[jd+1] #eigencoefficients - U_k and some polynomials
    ssqR <- ssqR + Mod(yk[k])^2
  }
  #write(*,*) "jd",jd," ssq yk rsd",ssqR

}
#don't need this
#write(*,*) "yrsd ", (yk(n),n=1,nord)

#for(n in 1:ndata){
#  for(ip in 0:mxdeg){
#    FitP[n] <- FitP[n] + Poly[n,ip+1]*pc[ip+1]
#  }
#  FitC[n] <- FitP[n]
#  for(k in 1:nord){
#    FitC[n] <- FitC[n] +  dw[n,k]*yk[k]
#  }
  #A1[n] <- abs(FitP[n])  #cabs again?
  #P1[n] <- phse(FitP[n])
  #A2[n] <- abs(FitC[n])
  #P2[n] <- phse(FitC[n])
  #x1[n] <-  Re(FitP[n])
  #y1[n] <- Im(FitP[n])
  #x2[n] <-  Re(FitC[n])
  #y2[n] <- Im(FitC[n])
#}
FitP <- FitP + Poly %*% pc
FitC <- FitP
FitC <- FitC + dw %*% yk

A1 <- abs(FitP)
P1  <- phse(FitP)
A2  <- abs(FitC)
P2  <- phse(FitC)
x1  <-  Re(FitP)
y1  <- Im(FitP)
x2  <-  Re(FitC)
y2  <- Im(FitC)
#these next few (7) lines are things he commented out
#CX      call a2mnmx(x1,y1,ndata, f1mn,f1mx)
#CX      call sdplos(t,x1,ndata,"Fit with Poly expansion"
#                    CX     1, tlbl,"Fit1",0.8,t(1),t(ndata),f1mn,f1mx)
#CX      call Gline(t,y1,ndata,1,2, 2.,"[4 4] 0","red")
#CX      call sdplot(t,A1,ndata,"Poly Amplitude"
#                    CX     1,tlbl,"Poly Amp", 0.8)
#CX      call a2mnmx(x2,y2,ndata, f1mn,f1mx)
#


#more plotting stuff. possibly replace
#call sdplos(t,x2,ndata,"Fit with composite expansion"
#            1, tlbl,"FitC",0.8,t(1),t(ndata),f1mn,f1mx)
#call Gline(t,y2,ndata,1,2, 2.,"[4 4] 0","red")
#call sdplot(t,A2,ndata,"Comp Amplitude"
#            1,tlbl,"Comp Amp", 0.8)
#call adstoa(A2,yd,ndata, -1.0)
#call cm4pp(yd,ndata,ave, var,"Comp amplitude - 1.0")
#call cm4pp(yd(nCalLo),nCompare,ave,var,"Comp Amp -1., Center")
#call sphsed(P1,ndata)
#call sphsed(P2,ndata)
#CX      call sdplot(t,P1,ndata,"Phase with Poly expansion"
#                    CX     1,tlbl,"P phse", 0.8)
#CX      call Gline(t,Pref,ndata, 1,2, 2.,"[4 4] 0","blue")
#call sdplot(t,P2,ndata,"Phase with Comp expansion"
#            1,tlbl,"C phse", 0.8)
#call Gline(t,Pref,ndata, 1,2, 2.,"[4 4] 0","blue")
yd <- P2 - Pref
#for(n in 1:ndata){
#  yd[n] = P2[n]-Pref[n]      #more phase stuff. 
#}
# have at least two occurences of stuff being assigned to yd. I think it's just assigned
# and then plotted and then assigned something new and then plotted again
#call cm4pp(yd,ndata,ave, var,"COMP Phse - Ref")
#call cm4pp(yd(nCalLo),nCompare,ave,var,"COMP Phse - Ref, Cntr")


FrScl <- 1.0/360.0
nplt <- ndata-1

system.time({
phi <- c(-1,rep(0,ndata-2),1)
Fr1 <- Re( FrScl*fft( fft(phi) * fft(PhD), inverse=TRUE) / ndata )
Fr1 <- Fr1[1:nplt]

FrP <- Re( FrScl*fft( fft(phi) * fft(P1), inverse=TRUE) / ndata )
FrP <- FrP[1:nplt]

FrC <- Re( FrScl*fft( fft(phi) * fft(P2), inverse=TRUE) / ndata )
FrC <- FrC[1:nplt]
})

#system.time({
#for(n in 1:nplt){
#  Fr1[n] = FrScl*(PhD[n+1]-PhD[n])
#  FrP[n] = FrScl*(P1[n+1]-P1[n])
#  FrC[n] = FrScl*(P2[n+1]-P2[n])
#}
#})
#more plotting
#call sdplot(t,Fr1,nplt,"Frequency, DPSS",tlbl,"Fr1",0.8)
#call Ghline(Fbounds,nfbnds,2, 1.0,"[12 12] 4","red")
#call Gline(t,Fmod,ndata,1,2, 2.,"[8 8] 0","blue")
#call sdplot(t,FrA,ndata,"Analytic Deriv",tlbl,"dpss FrA", 0.9)
#call Ghline(Fbounds,nfbnds,2, 1.0,"[12 12] 4","red")
#call Gline(t,Fmod,ndata,1,2, 2.,"[8 8] 0","blue")
#c                      Project Polynomials


#what do these do?
call setr(ndata, 0.0, Fit)  # set real
call movefr(nord, Fcoef, Frsd)  # move from one array to the other (copy)
#so is Frsd = Fcoef?
Frsd <- Fcoef
SSQC <- sum( Fcoef^2 )
#SSQC <- sdot(nord, Fcoef,1, Fcoef,1) #sdot is dot product. 1 is an increment
#write(*,*) "nord,SSQC", nord,SSQC

#nested loops. the end of the thing looks like plotting
#I'm currently butchering the following bits
FPcoef <- t(Ukz) %*% Fcoef

F1 <- numeric(mxdeg+1)
F2 <- F1
SSQFK <- numeric(mxdeg+1)
ssqr <- numeric(mxdeg+1)
SSQFtd <- numeric(mxdeg+1)
SSQRtd <- numeric(mxdeg+1)
Rsd <- matrix(rep(0,ndmax*(mxdeg+1)), nrow=ndmax)

Fcurv2 <- numeric(mxdeg+1)
Rcurv2 <- numeric(mxdeg+1)
AERR2 <- numeric(mxdeg+1)

FRsd <- matrix(rep(0,nord*(mxdeg+1)), nrow=nord)
for(L in 0:mxdeg){
  FRsd[,L+1] <- Frsd - Ukz[,L+1] * FPcoef[L+1]  
  ssqr <- sum( FRsd[,L+1]^2 )
  SSQFK <- SSQC - ssqr
  F1[L+1] <- (SSQFK/(L+1)) / (ssqr/(nord-L-1))
  
  #there may be something weird going on with this next part. 
  #it accumulates over L so that SSQFtd will be the sum of squares of
  Fit <- Fit + Poly[,L+1] * FPcoef[L+1]
  SSQFtd[L+1] <- sum( Fit^2 )
  
  
  Rsd[,L+1] <- dw %*% FRsd[,L+1]
  SSQRtd[L+1] <- sum( Rsd[,L+1]^2 )
  
  
  AERR2[L+1] <- sum( (Fit[2:(ndata-1)] - Fmod[2:(ndata-1)])^2 )
  Fcurv <- Re( fft( fft(curv.filter) * fft(Fit), inverse = TRUE) / ndata )
  Fcurv2[L+1] <- sum( Fcurv[2:(ndata-1)]^2 )
  Rcurv <- Re( fft( fft(curv.filter) * fft(Rsd[,L+1]), inverse = TRUE) / ndata )
  Rcurv2[L+1] <- sum( Rcurv[2:(ndata-1)]^2 )
  
  
  F2[L+1] <- (SSQFtd[L+1]/(L+1)) / (SSQRtd[L+1]/(nord-1-L))
}




#for(L in 0:mxdeg){
  
  
  #FPcoef[L] <- sdot(nord, Fcoef,1, Ukz[1,L],1) #sdot again
  #write(*,*) "L, FPcoef(L)", L, FPcoef(L)
  #ssqr <- 0.0
  #SSQFK <- 0.0
  #for(k in 1:nord){
  #  FRsd[k] <- Frsd[k] - FPcoef[L]*Ukz[k,L]
  #  ssqr <- ssqr + FRsd[k]^2
  #}
  #SSQFK <- SSQC - ssqr
  #F1 <- (SSQFK/(L+1)) / (ssqr/(nord-L-1))
  #write(*,*) "ssqr,SSQFK",ssqr,SSQFK," F ", F1
  #SSQFtd <- 0.0
  #for(n in 1:ndata){
  #  Fit[n] <- Fit[n] + Poly[n,L]* FPcoef[L]
  #  SSQFtd <- SSQFtd + Fit[n]^2
  #}
  #SSQRtd <- 0.0
  #for(n in 1:ndata){
  #  Rsd[n] <- 0.0
  #  for(k in 1:nord){
  #    Rsd[n] <- Rsd[n] + dw[n,k]*FRsd[k]
  #    SSQRtd <- SSQRtd + Rsd[n]^2
  #  }
  #}
  #Fcurv2 <- 0.0
  #Rcurv2 <- 0.0
  #AERR2 <- 0.0
  #for(n in 2:(ndata-1)){
  #  AERR2 <- AERR2 + (Fit[n]- Fmod[n])^2
  #  Fcurv2 <- Fcurv2 + (Fit[n-1]-2.0*Fit[n]+Fit[n+1])^2
  #  Rcurv2 <- Rcurv2 + (Rsd[n-1]-2.0*Rsd[n]+Rsd[n+1])^2
  #}
  #write(*,*) "SSQRtd",SSQRtd," SSQFtd",SSQFtd," Rcurv2",Rcurv2
  #1," Fcurv2",Fcurv2," AERR2",AERR2
  
  #if statement and line that follows are plotting stuff?
#  if(L > 1){ 
#    call sdplot(t,Fit,ndata,"Fit",tlbl,"Fit", 0.9)
#    call Gline(t,Fmod,ndata,1,2, 2.,"[8 8] 0","blue")
#  }
#  call sdplot(t,Rsd,ndata,"Rsd",tlbl,"Rsd", 0.9)
  #F2 <- (SSQFtd/(L+1)) / (SSQRtd/(nord-1-L))  #is this an F statistic?
  #write(*,*) "L,C, ssq Rsd", L,FPcoef(L), ssqr," F-TD", F2
#}
#for(n in 1:(ndata-1) ){
#  yd[n] = Fr1[n]-Fmod[n]
#}
yd <- Fr1 - Fmod
#no idea what this next stuff is. cm4pp shows up a lot
#call cm4pp(yd,ndata-1,ave, var,"DPSS Freq - Ref")
#call cm4pp(yd(nCalLo),nCompare,ave,var,"DPSS Freq - Ref, Cntr")
#CX      call sdplot(t,FrP,nplt,"Frequency, Poly",tlbl,"FrP",0.8)
#CX      call Gline(t,Fmod,ndata,1,2, 2.,"[8 8] 0","blue")
#call sdplot(t,FrC,nplt,"Frequency, Comp",tlbl,"FrC",0.8)
#call Gline(t,Fmod,ndata,1,2, 2.,"[8 8] 0","blue")


#for(n in 1:(ndata-1)){
#  yd[n] = FrC[n] - Fmod[n]
#}
yd <- FrC - Fmod

#more cm4pp. no clue. likely remove all this stuff.
#call cm4pp(yd,ndata-1,ave, var,"COMP  Freq - Ref")
#call cm4pp(yd(nCalLo),nCompare,ave,var,"COMP Freq - Ref, Cntr")
#9999 call Gclose(pname)
#write(*,*) " maximum stack size used =",istkst(3)
#stop "normal"
#end#
